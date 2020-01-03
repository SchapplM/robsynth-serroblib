% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RPRP3
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tau_reg [4x13]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:49
% EndTime: 2019-12-31 16:42:50
% DurationCPUTime: 0.30s
% Computational Cost: add. (285->82), mult. (576->120), div. (0->0), fcn. (309->8), ass. (0->57)
t28 = sin(pkin(6));
t18 = t28 * pkin(1) + pkin(5);
t56 = qJ(4) + t18;
t33 = cos(qJ(3));
t20 = t33 * pkin(3) + pkin(2);
t29 = cos(pkin(6));
t59 = t29 * pkin(1);
t43 = -t20 - t59;
t31 = sin(qJ(3));
t51 = qJD(1) * qJD(3);
t50 = t31 * t51;
t62 = pkin(3) * t50 + t43 * qJDD(1) + qJDD(4);
t48 = t56 * qJD(1);
t4 = t33 * qJD(2) - t48 * t31;
t57 = qJD(3) * pkin(3);
t3 = t4 + t57;
t61 = t3 - t4;
t60 = g(3) * t33;
t26 = t31 ^ 2;
t27 = t33 ^ 2;
t58 = t26 - t27;
t19 = -pkin(2) - t59;
t16 = qJD(1) * t19;
t55 = qJDD(2) - g(3);
t54 = t31 * qJDD(1);
t53 = t33 * qJDD(1);
t49 = qJD(3) * t56;
t25 = qJ(1) + pkin(6);
t21 = sin(t25);
t22 = cos(t25);
t47 = g(1) * t22 + g(2) * t21;
t46 = -g(1) * t21 + g(2) * t22;
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t45 = g(1) * t32 - g(2) * t34;
t5 = t31 * qJD(2) + t48 * t33;
t44 = t3 * t31 - t33 * t5;
t42 = t48 * qJD(3);
t13 = t18 * qJDD(1);
t41 = -t16 * qJD(1) - t13 + t47;
t40 = 0.2e1 * t16 * qJD(3) - qJDD(3) * t18;
t39 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(3) * qJD(2) + t13;
t35 = qJD(3) ^ 2;
t38 = -0.2e1 * qJDD(1) * t19 - t18 * t35 - t46;
t36 = qJD(1) ^ 2;
t30 = -qJ(4) - pkin(5);
t23 = t33 * qJDD(2);
t12 = qJDD(3) * t33 - t35 * t31;
t11 = qJDD(3) * t31 + t35 * t33;
t10 = t56 * t33;
t9 = t56 * t31;
t8 = t43 * qJD(1) + qJD(4);
t7 = -t31 * qJD(4) - t33 * t49;
t6 = t33 * qJD(4) - t31 * t49;
t2 = (qJDD(2) - t42) * t31 + t39 * t33;
t1 = qJDD(3) * pkin(3) - t39 * t31 - t33 * t42 + t23;
t14 = [qJDD(1), t45, g(1) * t34 + g(2) * t32, (t45 + (t28 ^ 2 + t29 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t26 * qJDD(1) + 0.2e1 * t33 * t50, 0.2e1 * t31 * t53 - 0.2e1 * t58 * t51, t11, t12, 0, t40 * t31 + t38 * t33, -t38 * t31 + t40 * t33, (-qJD(3) * t3 + qJDD(1) * t10 + t2 + (qJD(3) * t9 + t6) * qJD(1)) * t33 + (-qJD(3) * t5 + qJDD(1) * t9 - t1 + (-qJD(3) * t10 - t7) * qJD(1)) * t31 - t47, t2 * t10 + t5 * t6 - t1 * t9 + t3 * t7 + t8 * t31 * t57 - g(1) * (-t32 * pkin(1) - t21 * t20 - t22 * t30) - g(2) * (t34 * pkin(1) + t22 * t20 - t21 * t30) + t62 * t43; 0, 0, 0, t55, 0, 0, 0, 0, 0, t12, -t11, 0, -t44 * qJD(3) + t1 * t33 + t2 * t31 - g(3); 0, 0, 0, 0, -t31 * t36 * t33, t58 * t36, t54, t53, qJDD(3), t41 * t31 + t23 - t60, -t55 * t31 + t41 * t33, -pkin(3) * t54 + (-t57 + t61) * t33 * qJD(1), t61 * t5 + (-t60 + t1 + (-qJD(1) * t8 + t47) * t31) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t26 - t27) * t36, t44 * qJD(1) + t46 + t62;];
tau_reg = t14;
