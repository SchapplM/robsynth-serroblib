% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPR4
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% tau_reg [4x14]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:01
% EndTime: 2019-12-31 16:22:02
% DurationCPUTime: 0.12s
% Computational Cost: add. (103->35), mult. (167->49), div. (0->0), fcn. (83->4), ass. (0->33)
t16 = qJD(2) ^ 2;
t8 = pkin(6) + qJ(2);
t7 = cos(t8);
t5 = g(2) * t7;
t36 = -t16 * qJ(3) + t5;
t6 = sin(t8);
t35 = g(1) * t6;
t24 = qJDD(3) - t35;
t14 = -pkin(2) - pkin(5);
t13 = cos(qJ(4));
t10 = t13 ^ 2;
t12 = sin(qJ(4));
t34 = -t12 ^ 2 + t10;
t33 = t12 * t13;
t15 = qJD(4) ^ 2;
t32 = -t15 - t16;
t30 = pkin(2) * qJDD(2);
t11 = qJDD(1) - g(3);
t29 = qJDD(4) * t12;
t28 = qJDD(4) * t13;
t27 = t13 * qJDD(2);
t26 = qJ(3) * qJDD(2);
t25 = qJD(2) * qJD(4);
t23 = g(1) * t7 + g(2) * t6;
t22 = -t24 + t30;
t21 = 0.2e1 * qJ(3) * t25 + qJDD(4) * t14;
t20 = (2 * qJD(2) * qJD(3)) - t23;
t19 = -t14 * qJDD(2) - t24 - t36;
t18 = t20 + 0.2e1 * t26;
t17 = -t14 * t15 + t18;
t2 = -t15 * t12 + t28;
t1 = -t15 * t13 - t29;
t3 = [t11, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, t1, -t2; 0, qJDD(2), -t5 + t35, t23, t24 + t5 - 0.2e1 * t30, t18, (t22 - t5) * pkin(2) + (t20 + t26) * qJ(3), t10 * qJDD(2) - 0.2e1 * t25 * t33, -0.2e1 * t12 * t27 - 0.2e1 * t34 * t25, t2, t1, 0, t17 * t12 + t21 * t13, -t21 * t12 + t17 * t13; 0, 0, 0, 0, qJDD(2), -t16, -t22 + t36, 0, 0, 0, 0, 0, t32 * t12 + t28, t32 * t13 - t29; 0, 0, 0, 0, 0, 0, 0, t16 * t33, t34 * t16, t27, -t12 * qJDD(2), qJDD(4), -t11 * t12 - t19 * t13, -t11 * t13 + t19 * t12;];
tau_reg = t3;
