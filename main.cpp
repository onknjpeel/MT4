#define NOMINMAX

#include <Novice.h>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <utility>
#include <iostream>

const char kWindowTitle[] = "MT4_01_01_Basic";

struct Vector3 {
	float x, y, z;
};

struct Matrix4x4 {
	float m[4][4];
};



float Dot(const Vector3& v1, const Vector3& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }

float Length(const Vector3& v) { return std::sqrt(Dot(v, v)); }

Vector3 Normalize(const Vector3& v) {
	float length = Length(v);
	assert(length != 0.0f);
	return { v.x / length, v.y / length, v.z / length };
}

Matrix4x4 MakeRotateAxisAngle(const Vector3& axis, float angle) {

	float cos = std::cos(angle);
	float sin = std::sin(angle);

	Matrix4x4 rotateMatrix = {};
	rotateMatrix.m[0][0] = (axis.x * axis.x) * (1 - cos) + cos;
	rotateMatrix.m[0][1] = (axis.x * axis.y) * (1 - cos) + (axis.z * sin);
	rotateMatrix.m[0][2] = (axis.x * axis.z) * (1 - cos) - (axis.y * sin);
	rotateMatrix.m[0][3] = 0;

	rotateMatrix.m[1][0] = (axis.x * axis.y) * (1 - cos) - (axis.z * sin);
	rotateMatrix.m[1][1] = (axis.y * axis.y) * (1 - cos) + cos;
	rotateMatrix.m[1][2] = (axis.y * axis.z) * (1 - cos) + (axis.x * sin);
	rotateMatrix.m[1][3] = 0;

	rotateMatrix.m[2][0] = (axis.x * axis.z) * (1 - cos) + (axis.y * sin);
	rotateMatrix.m[2][1] = (axis.y * axis.z) * (1 - cos) - (axis.x * sin);
	rotateMatrix.m[2][2] = (axis.z * axis.z) * (1 - cos) + cos;
	rotateMatrix.m[2][3] = 0;

	rotateMatrix.m[3][0] = 0;
	rotateMatrix.m[3][1] = 0;
	rotateMatrix.m[3][2] = 0;
	rotateMatrix.m[3][3] = 1;

	return rotateMatrix;
}

Matrix4x4 Inverse(const Matrix4x4& m) {
	float determinant =
		+m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3]
		+ m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1]
		+ m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]

		- m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1]
		- m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3]
		- m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]

		- m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3]
		- m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1]
		- m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]

		+ m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1]
		+ m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3]
		+ m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]

		+ m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3]
		+ m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1]
		+ m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]

		- m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1]
		- m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3]
		- m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]

		- m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0]
		- m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0]
		- m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]

		+ m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0]
		+ m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0]
		+ m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0];

	Matrix4x4 result = {};
	float recpDeterminant = 1.0f / determinant;
	result.m[0][0] = (m.m[1][1] * m.m[2][2] * m.m[3][3] + m.m[1][2] * m.m[2][3] * m.m[3][1] +
		m.m[1][3] * m.m[2][1] * m.m[3][2] - m.m[1][3] * m.m[2][2] * m.m[3][1] -
		m.m[1][2] * m.m[2][1] * m.m[3][3] - m.m[1][1] * m.m[2][3] * m.m[3][2]) * recpDeterminant;
	result.m[0][1] = (-m.m[0][1] * m.m[2][2] * m.m[3][3] - m.m[0][2] * m.m[2][3] * m.m[3][1] -
		m.m[0][3] * m.m[2][1] * m.m[3][2] + m.m[0][3] * m.m[2][2] * m.m[3][1] +
		m.m[0][2] * m.m[2][1] * m.m[3][3] + m.m[0][1] * m.m[2][3] * m.m[3][2]) * recpDeterminant;
	result.m[0][2] = (m.m[0][1] * m.m[1][2] * m.m[3][3] + m.m[0][2] * m.m[1][3] * m.m[3][1] +
		m.m[0][3] * m.m[1][1] * m.m[3][2] - m.m[0][3] * m.m[1][2] * m.m[3][1] -
		m.m[0][2] * m.m[1][1] * m.m[3][3] - m.m[0][1] * m.m[1][3] * m.m[3][2]) * recpDeterminant;
	result.m[0][3] = (-m.m[0][1] * m.m[1][2] * m.m[2][3] - m.m[0][2] * m.m[1][3] * m.m[2][1] -
		m.m[0][3] * m.m[1][1] * m.m[2][2] + m.m[0][3] * m.m[1][2] * m.m[2][1] +
		m.m[0][2] * m.m[1][1] * m.m[2][3] + m.m[0][1] * m.m[1][3] * m.m[2][2]) * recpDeterminant;

	result.m[1][0] = (-m.m[1][0] * m.m[2][2] * m.m[3][3] - m.m[1][2] * m.m[2][3] * m.m[3][0] -
		m.m[1][3] * m.m[2][0] * m.m[3][2] + m.m[1][3] * m.m[2][2] * m.m[3][0] +
		m.m[1][2] * m.m[2][0] * m.m[3][3] + m.m[1][0] * m.m[2][3] * m.m[3][2]) * recpDeterminant;
	result.m[1][1] = (m.m[0][0] * m.m[2][2] * m.m[3][3] + m.m[0][2] * m.m[2][3] * m.m[3][0] +
		m.m[0][3] * m.m[2][0] * m.m[3][2] - m.m[0][3] * m.m[2][2] * m.m[3][0] -
		m.m[0][2] * m.m[2][0] * m.m[3][3] - m.m[0][0] * m.m[2][3] * m.m[3][2]) * recpDeterminant;
	result.m[1][2] = (-m.m[0][0] * m.m[1][2] * m.m[3][3] - m.m[0][2] * m.m[1][3] * m.m[3][0] -
		m.m[0][3] * m.m[1][0] * m.m[3][2] + m.m[0][3] * m.m[1][2] * m.m[3][0] +
		m.m[0][2] * m.m[1][0] * m.m[3][3] + m.m[0][0] * m.m[1][3] * m.m[3][2]) * recpDeterminant;
	result.m[1][3] = (m.m[0][0] * m.m[1][2] * m.m[2][3] + m.m[0][2] * m.m[1][3] * m.m[2][0] +
		m.m[0][3] * m.m[1][0] * m.m[2][2] - m.m[0][3] * m.m[1][2] * m.m[2][0] -
		m.m[0][2] * m.m[1][0] * m.m[2][3] - m.m[0][0] * m.m[1][3] * m.m[2][2]) * recpDeterminant;

	result.m[2][0] = (m.m[1][0] * m.m[2][1] * m.m[3][3] + m.m[1][1] * m.m[2][3] * m.m[3][0] +
		m.m[1][3] * m.m[2][0] * m.m[3][1] - m.m[1][3] * m.m[2][1] * m.m[3][0] -
		m.m[1][1] * m.m[2][0] * m.m[3][3] - m.m[1][0] * m.m[2][3] * m.m[3][1]) * recpDeterminant;
	result.m[2][1] = (-m.m[0][0] * m.m[2][1] * m.m[3][3] - m.m[0][1] * m.m[2][3] * m.m[3][0] -
		m.m[0][3] * m.m[2][0] * m.m[3][1] + m.m[0][3] * m.m[2][1] * m.m[3][0] +
		m.m[0][1] * m.m[2][0] * m.m[3][3] + m.m[0][0] * m.m[2][3] * m.m[3][1]) * recpDeterminant;
	result.m[2][2] = (m.m[0][0] * m.m[1][1] * m.m[3][3] + m.m[0][1] * m.m[1][3] * m.m[3][0] +
		m.m[0][3] * m.m[1][0] * m.m[3][1] - m.m[0][3] * m.m[1][1] * m.m[3][0] -
		m.m[0][1] * m.m[1][0] * m.m[3][3] - m.m[0][0] * m.m[1][3] * m.m[3][1]) * recpDeterminant;
	result.m[2][3] = (-m.m[0][0] * m.m[1][1] * m.m[2][3] - m.m[0][1] * m.m[1][3] * m.m[2][0] -
		m.m[0][3] * m.m[1][0] * m.m[2][1] + m.m[0][3] * m.m[1][1] * m.m[2][0] +
		m.m[0][1] * m.m[1][0] * m.m[2][3] + m.m[0][0] * m.m[1][3] * m.m[2][1]) * recpDeterminant;

	result.m[3][0] = (-m.m[1][0] * m.m[2][1] * m.m[3][2] - m.m[1][1] * m.m[2][2] * m.m[3][0] -
		m.m[1][2] * m.m[2][0] * m.m[3][1] + m.m[1][2] * m.m[2][1] * m.m[3][0] +
		m.m[1][1] * m.m[2][0] * m.m[3][2] + m.m[1][0] * m.m[2][2] * m.m[3][1]) * recpDeterminant;
	result.m[3][1] = (m.m[0][0] * m.m[2][1] * m.m[3][2] + m.m[0][1] * m.m[2][2] * m.m[3][0] +
		m.m[0][2] * m.m[2][0] * m.m[3][1] - m.m[0][2] * m.m[2][1] * m.m[3][0] -
		m.m[0][1] * m.m[2][0] * m.m[3][2] - m.m[0][0] * m.m[2][2] * m.m[3][1]) * recpDeterminant;
	result.m[3][2] = (-m.m[0][0] * m.m[1][1] * m.m[3][2] - m.m[0][1] * m.m[1][2] * m.m[3][0] -
		m.m[0][2] * m.m[1][0] * m.m[3][1] + m.m[0][2] * m.m[1][1] * m.m[3][0] +
		m.m[0][1] * m.m[1][0] * m.m[3][2] + m.m[0][0] * m.m[1][2] * m.m[3][1]) * recpDeterminant;
	result.m[3][3] = (m.m[0][0] * m.m[1][1] * m.m[2][2] + m.m[0][1] * m.m[1][2] * m.m[2][0] +
		m.m[0][2] * m.m[1][0] * m.m[2][1] - m.m[0][2] * m.m[1][1] * m.m[2][0] -
		m.m[0][1] * m.m[1][0] * m.m[2][2] - m.m[0][0] * m.m[1][2] * m.m[2][1]) * recpDeterminant;

	return result;
}

Vector3 Transform(Vector3 vector, Matrix4x4 matrix) {
	Vector3 result;
	result = {
		vector.x * matrix.m[0][0] + vector.y * matrix.m[1][0] + vector.z * matrix.m[2][0] + matrix.m[3][0],
		vector.x * matrix.m[0][1] + vector.y * matrix.m[1][1] + vector.z * matrix.m[2][1] + matrix.m[3][1],
		vector.x * matrix.m[0][2] + vector.y * matrix.m[1][2] + vector.z * matrix.m[2][2] + matrix.m[3][2]
	};
	float w = vector.x * matrix.m[0][3] + vector.y * matrix.m[1][3] + vector.z * matrix.m[2][3] + matrix.m[3][3];

	assert(w != 0.0f);
	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
}

static const int kRowHeight = 20;
static const int kColumnWidth = 60;

void VectorScreenPrintf(float width, float height, Vector3 v, const char* label) {
	Novice::ScreenPrintf(int(width), int(height), "%.02f", v.x);
	Novice::ScreenPrintf(int(width + kColumnWidth), int(height), "%.02f", v.y);
	Novice::ScreenPrintf(int(width + kColumnWidth * 2), int(height), "%.02f", v.z);
	Novice::ScreenPrintf(int(width + kColumnWidth * 3), int(height), "%s", label);
}

void MatrixScreenPrintf(int x, int y, const Matrix4x4& matrix, const char* label) {
	Novice::ScreenPrintf(x, y, "%s", label);
	for (int row = 0; row < 4; ++row) {
		for (int column = 0; column < 4; ++column) {
			Novice::ScreenPrintf(
				x + column * kColumnWidth, y + (row + 1) * kRowHeight, "%6.03f",
				matrix.m[row][column]);
		}
	}
}

#pragma region Quaternion
/*四元数*/
struct Quaternion {
	float x;
	float y;
	float z;
	float w;
};

/*四元数同士の積*/
Quaternion Multiply(const Quaternion& lhs, const Quaternion& rhs)
{
	Quaternion result = {
		((lhs.w * rhs.x) + (lhs.x * rhs.w) + (lhs.y * rhs.z) - (lhs.z * rhs.y)),
		((lhs.w * rhs.y) - (lhs.x * rhs.z) + (lhs.y * rhs.w) + (lhs.z * rhs.x)),
		((lhs.w * rhs.z) + (lhs.x * rhs.y) - (lhs.y * rhs.x) + (lhs.z * rhs.w)),
		(lhs.w * rhs.w) - (lhs.x * rhs.x) - (lhs.y * rhs.y) - (lhs.z * rhs.z)
	};
	return result;
}

/*単位四元数*/
Quaternion IdentityQuaternion()
{
	Quaternion result = { 0,0,0,1 };
	return result;
}

/*共役四元数*/
Quaternion Conjugate(const Quaternion& quaternion)
{
	Quaternion result = {
		quaternion.x * -1,
		quaternion.y * -1,
		quaternion.z * -1,
		quaternion.w
	};
	return result;
}

/*四元数のノルム*/
float Norm(const Quaternion& quaternion)
{
	float result = std::sqrtf((quaternion.x * quaternion.x) + (quaternion.y * quaternion.y) + (quaternion.z * quaternion.z) + (quaternion.w * quaternion.w));
	return result;
}

/*四元数の正規化*/
Quaternion Normalize(const Quaternion& quaternion)
{
	float length = Norm(quaternion);
	assert(length != 0.0f);

	return { quaternion.x / length,quaternion.y / length,quaternion.z / length ,quaternion.w / length };
}

/*逆四元数*/
Quaternion Inverse(const Quaternion& quaternion)
{
	float normSquared = quaternion.x * quaternion.x +
		quaternion.y * quaternion.y +
		quaternion.z * quaternion.z +
		quaternion.w * quaternion.w;

	Quaternion conjugated = Conjugate(quaternion);
	return Quaternion{
		conjugated.x / normSquared,
		conjugated.y / normSquared,
		conjugated.z / normSquared,
		conjugated.w / normSquared
	};
}

/*四元数を表示*/
void QuaternionScreenPrintf(int x, int y, Quaternion quaternion, const char* label) {
	Novice::ScreenPrintf(x, y, "%.02f", quaternion.x);
	Novice::ScreenPrintf(x + kColumnWidth, y, "%.02f", quaternion.y);
	Novice::ScreenPrintf(x + kColumnWidth * 2, y, "%.02f", quaternion.z);
	Novice::ScreenPrintf(x + kColumnWidth * 3, y, "%.02f", quaternion.w);
	Novice::ScreenPrintf(x + kColumnWidth * 4, y, label);
}

//任意軸回転を表すQuaternionの生成
Quaternion MakeRotateAxisAngleQuaternion(const Vector3& axis, float angle) {
	Quaternion result = {
		axis.x * std::sin(angle / 2),axis.y * std::sin(angle / 2),axis.z * std::sin(angle / 2),std::cos(angle / 2)
	};
	return result;
}

//ベクトルをQuaternionで回転させた結果のベクトルを求める
Vector3 RotateVector(const Vector3& vector, const Quaternion& q) {
	// ベクトルをクォータニオンに変換 (スカラー部分が 0)
	Quaternion vecQuat(vector.x, vector.y, vector.z, 0);

	// 回転クォータニオンとその共役を使って回転
	Quaternion qConjugate = Conjugate(q);

	// 回転を計算
	Quaternion rotatedQuat = Multiply(Multiply(q, vecQuat), qConjugate);

	// 結果をベクトルとして返す
	return Vector3(rotatedQuat.x, rotatedQuat.y, rotatedQuat.z);
}

//Quaternionから回転行列を求める
Matrix4x4 MakeRotateMatrix(const Quaternion q) {
	Matrix4x4 result = {
		(q.w * q.w) + (q.x * q.x) - (q.y * q.y) - (q.z * q.z),2 * ((q.x * q.y) + (q.w * q.z)),2 * ((q.x * q.z) - (q.w * q.y)),0,
		2 * ((q.x * q.y) - (q.w * q.z)),(q.w * q.w) - (q.x * q.x) + (q.y * q.y) - (q.z * q.z),2 * ((q.y * q.z) + (q.w * q.x)),0,
		2 * ((q.x * q.z) + (q.w * q.y)),2 * ((q.y * q.z) - (q.w * q.x)),(q.w * q.w) - (q.x * q.x) - (q.y * q.y) + (q.z * q.z),0,
		0,0,0,1
	};
	return result;
}

Quaternion operator-(Quaternion q) {
	Quaternion result;
	result = { -q.x,-q.y ,-q.z ,-q.w };
	return result;
}

Quaternion operator*(const Quaternion& q, float scalar) {
	return { q.x * scalar, q.y * scalar, q.z * scalar, q.w * scalar };
}

Quaternion operator*(float scalar, const Quaternion& q) {
	return q * scalar;  // 既存のオーバーロードを利用
}

Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
	return { q1.x + q2.x, q1.y + q2.y, q1.z + q2.z, q1.w + q2.w };
}

float Dot(const Quaternion& q0, const Quaternion& q1) {
	return q0.x * q1.x + q0.y * q1.y + q0.z * q1.z + q0.w * q1.w;
}

Quaternion Slerp(const Quaternion& q0, const Quaternion& q1, float t) {
	
	Quaternion p0 = q0;
	Quaternion p1 = q1;

	float dot = Dot(q0, q1);
	if (dot < 0) {
		p0 = -p0;
		dot = -dot;
	}

	float theta = std::acos(dot);

	float scale0 = std::sin((1 - t) * theta) / std::sin(theta);
	float scale1 = std::sin(t * theta) / std::sin(theta);

	return scale0 * p0 + scale1 * p1;
}


#pragma endregion

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	const int kWindowWidth = 1280;
	const int kWindowHeight = 720;
	Novice::Initialize(kWindowTitle, kWindowWidth, kWindowHeight);

	Quaternion rotation0 = MakeRotateAxisAngleQuaternion({ 0.71f,0.71f,0.0f }, 0.3f);
	Quaternion rotation1 = MakeRotateAxisAngleQuaternion({ 0.71f,0.0f,0.71f }, 3.141592f);

	Quaternion interpolate0 = Slerp(rotation0, rotation1, 0.0f);
	Quaternion interpolate1 = Slerp(rotation0, rotation1, 0.3f);
	Quaternion interpolate2 = Slerp(rotation0, rotation1, 0.5f);
	Quaternion interpolate3 = Slerp(rotation0, rotation1, 0.7f);
	Quaternion interpolate4 = Slerp(rotation0, rotation1, 1.0f);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		QuaternionScreenPrintf(0, kRowHeight * 0, interpolate0, ": interpolate0, Slerp(q0, p1, 0.0f)");
		QuaternionScreenPrintf(0, kRowHeight * 1, interpolate1, ": interpolate1, Slerp(q0, p1, 0.3f)");
		QuaternionScreenPrintf(0, kRowHeight * 2, interpolate2, ": interpolate2, Slerp(q0, p1, 0.5f)");
		QuaternionScreenPrintf(0, kRowHeight * 3, interpolate3, ": interpolate3, Slerp(q0, p1, 0.7f)");
		QuaternionScreenPrintf(0, kRowHeight * 4, interpolate4, ": interpolate4, Slerp(q0, p1, 1.0f)");

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}