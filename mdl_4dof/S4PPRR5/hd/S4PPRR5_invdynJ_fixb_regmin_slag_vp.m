% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x12]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:53
% EndTime: 2019-12-31 16:19:54
% DurationCPUTime: 0.20s
% Computational Cost: add. (107->47), mult. (256->80), div. (0->0), fcn. (181->6), ass. (0->34)
t11 = sin(pkin(6));
t12 = cos(pkin(6));
t24 = -g(1) * t11 + g(2) * t12 + qJDD(2);
t10 = qJDD(1) - g(3);
t14 = sin(qJ(3));
t16 = cos(qJ(3));
t17 = qJD(4) ^ 2;
t41 = -(2 * qJDD(3) * pkin(3)) + pkin(5) * t17 + t10 * t14 - t24 * t16;
t32 = t14 * qJDD(3);
t18 = qJD(3) ^ 2;
t38 = t17 + t18;
t40 = t38 * t16 + t32;
t15 = cos(qJ(4));
t13 = sin(qJ(4));
t8 = t13 ^ 2;
t39 = -t15 ^ 2 + t8;
t37 = qJD(3) * pkin(3);
t6 = -qJD(1) * t14 + qJD(2) * t16;
t36 = t6 * qJD(3);
t33 = qJDD(4) * t13;
t31 = t15 * qJDD(3);
t30 = t16 * qJDD(3);
t29 = qJD(3) * qJD(4);
t28 = t13 * t29;
t27 = g(1) * t12 + g(2) * t11;
t25 = -qJD(2) * qJD(3) - t10;
t23 = -qJDD(4) * t16 + 0.2e1 * t14 * t29;
t22 = qJD(1) * qJD(3) - t24;
t3 = -t6 - t37;
t21 = -pkin(5) * qJDD(4) + (t3 + t6 - t37) * qJD(4);
t19 = -qJDD(3) * pkin(5) - t3 * qJD(3) - t10 * t16 - t24 * t14 - t36;
t5 = -t16 * t18 - t32;
t4 = t14 * t18 - t30;
t1 = [t10, t10, 0, t5, t4, 0, 0, 0, 0, 0, t23 * t13 - t40 * t15, t40 * t13 + t23 * t15; 0, t24, 0, -t4, t5, 0, 0, 0, 0, 0, (-0.2e1 * t28 + t31) * t16 + (-t38 * t15 - t33) * t14, (-qJDD(4) * t14 - 0.2e1 * t16 * t29) * t15 + (t38 * t14 - t30) * t13; 0, 0, qJDD(3), t25 * t14 - t22 * t16 + (t16 * qJD(1) + t14 * qJD(2)) * qJD(3), t22 * t14 + t25 * t16 + t36, qJDD(3) * t8 + 0.2e1 * t15 * t28, 0.2e1 * t13 * t31 - 0.2e1 * t39 * t29, t15 * t17 + t33, qJDD(4) * t15 - t13 * t17, 0, t21 * t13 - t41 * t15, t41 * t13 + t21 * t15; 0, 0, 0, 0, 0, -t13 * t18 * t15, t39 * t18, t13 * qJDD(3), t31, qJDD(4), t19 * t13 - t27 * t15, t27 * t13 + t19 * t15;];
tau_reg = t1;
