import tensorflow as tf
import tensorflow_datasets as tfds


def normalize(x, y):
  return tf.cast(x, tf.float32) / 255., y


(data_train, data_test), data_info = tfds.load(
    'mnist',
    split=['train', 'test'],
    shuffle_files=True,
    as_supervised=True,
    with_info=True,
)

data_train = data_train.map(normalize)
data_train = data_train.cache()
data_train = data_train.shuffle(data_info.splits['train'].num_examples)
data_train = data_train.batch(128)

data_test = data_test.map(normalize)
data_test = data_test.batch(128)
data_test = data_test.cache()


model = tf.keras.models.Sequential([
  tf.keras.layers.Flatten(input_shape=(28, 28)),
  tf.keras.layers.Dense(100, activation='relu'),
  tf.keras.layers.Dense(50, activation='relu'),
  tf.keras.layers.Dense(10)
])

model.compile(
    optimizer=tf.keras.optimizers.Adam(0.005),
    loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    metrics=[tf.keras.metrics.SparseCategoricalAccuracy()],
)

model.fit(
    data_train,
    epochs=6,
    validation_data=data_test,
)
