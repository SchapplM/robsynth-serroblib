% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4PRRP3
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tau_reg [4x15]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:05
% EndTime: 2021-01-14 22:27:07
% DurationCPUTime: 0.37s
% Computational Cost: add. (296->111), mult. (603->140), div. (0->0), fcn. (318->4), ass. (0->69)
t76 = qJDD(1) - g(3);
t55 = qJD(1) * qJD(3);
t62 = qJ(4) + pkin(5);
t38 = -t62 * qJDD(2) - t55;
t75 = -qJD(2) * qJD(4) + t38;
t74 = -2 * pkin(2);
t34 = cos(qJ(3));
t27 = t34 * qJD(1);
t33 = sin(qJ(3));
t51 = qJD(2) * t62;
t5 = -t33 * t51 + t27;
t60 = qJD(3) * pkin(3);
t3 = t5 + t60;
t73 = t3 - t5;
t30 = t33 ^ 2;
t72 = pkin(3) * t30;
t71 = pkin(3) * t34;
t29 = pkin(6) + qJ(2);
t23 = sin(t29);
t70 = g(1) * t23;
t24 = cos(t29);
t69 = g(1) * t24;
t68 = g(2) * t23;
t67 = g(2) * t24;
t66 = g(3) * t34;
t65 = t23 * t34;
t64 = t24 * t33;
t36 = qJD(2) ^ 2;
t63 = t34 * t36;
t31 = t34 ^ 2;
t61 = t30 - t31;
t59 = qJD(2) * t33;
t58 = qJDD(3) * pkin(3);
t22 = pkin(2) + t71;
t57 = qJDD(2) * t22;
t25 = t33 * qJDD(2);
t56 = t34 * qJDD(2);
t54 = qJD(2) * qJD(3);
t13 = t62 * t34;
t52 = t33 * t54;
t50 = qJD(3) * t62;
t4 = pkin(3) * t52 + qJDD(4) - t57;
t49 = t4 - t57;
t48 = 0.2e1 * t34 * t54;
t47 = t34 * t50;
t46 = t68 + t69;
t45 = -t67 + t70;
t7 = t33 * qJD(1) + qJD(2) * t13;
t44 = t3 * t33 - t34 * t7;
t35 = qJD(3) ^ 2;
t43 = pkin(5) * t35 + qJDD(2) * t74;
t42 = pkin(2) * t36 - pkin(5) * qJDD(2);
t26 = t34 * qJDD(1);
t41 = g(1) * t64 + t33 * t68 + t26 - t66;
t40 = -pkin(5) * qJDD(3) + t54 * t74;
t20 = pkin(5) * t52;
t39 = g(2) * t65 - t76 * t33 + t34 * t69 + t20;
t9 = -qJD(2) * t22 + qJD(4);
t37 = (-qJD(4) - t9) * qJD(2) + t38;
t17 = g(1) * t65;
t16 = g(2) * t64;
t12 = t62 * t33;
t11 = qJDD(3) * t34 - t33 * t35;
t10 = qJDD(3) * t33 + t34 * t35;
t8 = -t33 * qJD(4) - t47;
t6 = t34 * qJD(4) - t33 * t50;
t2 = -t20 + (-qJ(4) * t54 + qJDD(1)) * t33 - t75 * t34;
t1 = -qJD(2) * t47 + t75 * t33 + t26 + t58;
t14 = [t76, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, t11, -t10, 0, -t44 * qJD(3) + t1 * t34 + t2 * t33 - g(3); 0, qJDD(2), t45, t46, qJDD(2) * t30 + t33 * t48, 0.2e1 * t33 * t56 - 0.2e1 * t61 * t54, t10, t11, 0, t17 + t40 * t33 + (-t43 - t67) * t34, t16 + t40 * t34 + (t43 - t70) * t33, -qJDD(3) * t12 + t17 + (-t49 - t67) * t34 + (t8 + (t9 + (-t22 - t71) * qJD(2)) * t33) * qJD(3), -qJDD(3) * t13 + t16 + (t49 - t70) * t33 + (t34 * t9 - t6 + (-t22 * t34 + t72) * qJD(2)) * qJD(3), (-qJD(3) * t3 + qJDD(2) * t13 + t2 + (qJD(3) * t12 + t6) * qJD(2)) * t34 + (-qJD(3) * t7 + qJDD(2) * t12 - t1 + (-qJD(3) * t13 - t8) * qJD(2)) * t33 - t46, t2 * t13 + t7 * t6 - t1 * t12 + t3 * t8 - t4 * t22 + t9 * t33 * t60 - g(1) * (-t22 * t23 + t24 * t62) - g(2) * (t22 * t24 + t23 * t62); 0, 0, 0, 0, -t33 * t63, t61 * t36, t25, t56, qJDD(3), t42 * t33 + t41, (-pkin(5) * t59 + t27) * qJD(3) + (t42 - t55) * t34 + t39, 0.2e1 * t58 + (-t34 * t51 + t7) * qJD(3) + (pkin(3) * t63 + t37) * t33 + t41, -t36 * t72 + (qJ(4) * t59 + t5) * qJD(3) + t37 * t34 + t39, -pkin(3) * t25 + (-t60 + t73) * t34 * qJD(2), t73 * t7 + (-t66 + t1 + (-qJD(2) * t9 + t46) * t33) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t52 - t56, t25 + t48, (-t30 - t31) * t36, t44 * qJD(2) + t4 - t45;];
tau_reg = t14;
