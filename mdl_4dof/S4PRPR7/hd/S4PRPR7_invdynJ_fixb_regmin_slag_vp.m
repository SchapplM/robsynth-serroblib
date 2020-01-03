% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRPR7
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:56
% EndTime: 2019-12-31 16:25:57
% DurationCPUTime: 0.22s
% Computational Cost: add. (138->69), mult. (293->88), div. (0->0), fcn. (190->6), ass. (0->48)
t16 = sin(qJ(2));
t18 = cos(qJ(2));
t13 = sin(pkin(6));
t14 = cos(pkin(6));
t34 = g(1) * t14 + g(2) * t13;
t60 = -g(3) * t18 + t34 * t16;
t19 = -pkin(2) - pkin(5);
t40 = t18 * qJDD(1);
t35 = qJDD(3) - t40;
t39 = qJD(1) * qJD(2);
t8 = t16 * t39;
t30 = t35 + t8;
t45 = qJD(2) * qJ(3);
t48 = t16 * qJD(1);
t7 = t45 + t48;
t50 = t7 * qJD(2);
t59 = -t19 * qJDD(2) - t30 + t50 + t60;
t57 = (t45 + t7 - t48) * qJD(4) + qJDD(4) * t19;
t46 = qJDD(1) - g(3);
t56 = -t46 * t16 + t34 * t18;
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t53 = t15 * t17;
t12 = t17 ^ 2;
t52 = t15 ^ 2 - t12;
t20 = qJD(4) ^ 2;
t21 = qJD(2) ^ 2;
t51 = t20 + t21;
t49 = qJDD(2) * pkin(2);
t47 = t18 * qJD(1);
t44 = qJDD(2) * t16;
t43 = qJDD(4) * t15;
t41 = t17 * qJDD(2);
t38 = qJD(2) * qJD(3);
t37 = qJD(2) * qJD(4);
t36 = qJDD(2) * qJ(3);
t33 = g(1) * t13 - g(2) * t14;
t32 = (-qJD(2) * pkin(2) + qJD(3) - t47) * t16 + t7 * t18;
t25 = t51 * t18 + t44;
t24 = -qJDD(4) * t18 + 0.2e1 * t16 * t37;
t23 = t35 - t60;
t2 = t36 + t16 * qJDD(1) + (qJD(3) + t47) * qJD(2);
t22 = -g(3) * t16 + t36 + t38 - t19 * t20 + t2 + (-t34 - t39) * t18;
t9 = qJDD(4) * t17;
t5 = -t18 * qJDD(2) + t21 * t16;
t4 = t21 * t18 + t44;
t3 = t30 - t49;
t1 = [t46, 0, -t5, -t4, t5, t4, t32 * qJD(2) + t2 * t16 - t3 * t18 - g(3), 0, 0, 0, 0, 0, t25 * t15 + t24 * t17, -t24 * t15 + t25 * t17; 0, qJDD(2), t40 + t60, t56, t23 - 0.2e1 * t49, 0.2e1 * t36 + 0.2e1 * t38 - t56, t2 * qJ(3) + t7 * qJD(3) - t3 * pkin(2) - g(3) * (t18 * pkin(2) + t16 * qJ(3)) - t32 * qJD(1) + t34 * (pkin(2) * t16 - qJ(3) * t18), t12 * qJDD(2) - 0.2e1 * t37 * t53, -0.2e1 * t15 * t41 + 0.2e1 * t52 * t37, -t20 * t15 + t9, -t20 * t17 - t43, 0, t22 * t15 + t57 * t17, -t57 * t15 + t22 * t17; 0, 0, 0, 0, qJDD(2), -t21, t23 - t49 + t8 - t50, 0, 0, 0, 0, 0, -t51 * t15 + t9, -t51 * t17 - t43; 0, 0, 0, 0, 0, 0, 0, t21 * t53, -t52 * t21, t41, -t15 * qJDD(2), qJDD(4), t33 * t15 - t59 * t17, t59 * t15 + t33 * t17;];
tau_reg = t1;
