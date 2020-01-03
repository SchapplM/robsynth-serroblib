% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:47
% EndTime: 2019-12-31 18:02:50
% DurationCPUTime: 0.90s
% Computational Cost: add. (637->160), mult. (1342->257), div. (0->0), fcn. (796->6), ass. (0->96)
t46 = cos(qJ(4));
t79 = t46 * qJD(1);
t29 = qJD(5) + t79;
t111 = qJD(5) - t29;
t40 = sin(pkin(8));
t48 = qJD(4) ^ 2;
t49 = qJD(1) ^ 2;
t110 = t40 * (t48 + t49);
t44 = sin(qJ(4));
t38 = t44 ^ 2;
t43 = sin(qJ(5));
t109 = (qJD(1) * t38 - t29 * t46) * t43;
t47 = -pkin(1) - pkin(2);
t41 = cos(pkin(8));
t77 = qJD(1) * qJD(2);
t65 = t41 * t77;
t28 = qJD(1) * t47 + qJD(2);
t78 = qJD(1) * qJ(2);
t18 = t40 * t28 + t41 * t78;
t14 = -qJD(1) * pkin(6) + t18;
t9 = qJD(3) * t44 + t14 * t46;
t4 = t9 * qJD(4) + t44 * t65;
t108 = t4 * t43;
t45 = cos(qJ(5));
t107 = t4 * t45;
t76 = qJD(1) * qJD(4);
t66 = t46 * t76;
t89 = qJD(1) * t44;
t69 = t43 * t89;
t75 = qJD(4) * qJD(5);
t11 = qJD(5) * t69 + (-t66 + t75) * t45;
t106 = t11 * t43;
t105 = t11 * t46;
t83 = qJD(5) * t45;
t70 = t44 * t83;
t81 = t43 * qJD(4);
t12 = -t43 * t75 + (t46 * t81 + t70) * qJD(1);
t104 = t12 * t46;
t80 = t45 * qJD(4);
t20 = t69 + t80;
t103 = t20 * t29;
t21 = t45 * t89 - t81;
t102 = t21 * t29;
t101 = t43 * t29;
t100 = t43 * t46;
t99 = t44 * t11;
t98 = t44 * t20;
t97 = t44 * t21;
t96 = t45 * t29;
t95 = t45 * t46;
t94 = t48 * t44;
t93 = t48 * t46;
t92 = t41 * qJ(2) + t40 * t47;
t91 = -t46 ^ 2 + t38;
t88 = qJD(2) * t41;
t23 = -pkin(6) + t92;
t87 = qJD(4) * t23;
t86 = qJD(4) * t44;
t85 = qJD(4) * t46;
t84 = qJD(5) * t43;
t82 = qJD(5) * t46;
t74 = t29 * t95;
t73 = t29 * t84;
t72 = t44 * t84;
t71 = t29 * t83;
t68 = 0.2e1 * t77;
t6 = qJD(4) * pkin(7) + t9;
t67 = t23 * t29 + t6;
t64 = t44 * t76;
t17 = t41 * t28 - t40 * t78;
t63 = -t40 * qJ(2) + t41 * t47;
t62 = t29 * t70;
t61 = t40 * t68;
t60 = 0.2e1 * t64;
t22 = pkin(3) - t63;
t13 = qJD(1) * pkin(3) - t17;
t59 = pkin(4) * t46 + pkin(7) * t44;
t58 = -pkin(4) * t44 + pkin(7) * t46;
t7 = qJD(1) * t59 + t13;
t2 = t43 * t7 + t45 * t6;
t57 = t43 * t6 - t45 * t7;
t56 = t17 * t40 - t18 * t41;
t8 = qJD(3) * t46 - t14 * t44;
t54 = qJD(1) * (t13 - t88);
t53 = t45 * t38 * t76 + t29 * t72;
t52 = -t23 * t48 + t61;
t51 = qJD(4) * (-qJD(1) * t22 - t13 - t88);
t3 = t8 * qJD(4) + t46 * t65;
t5 = -qJD(4) * pkin(4) - t8;
t50 = -qJD(4) * t5 - qJD(5) * t7 - t29 * t88 - t3;
t19 = t40 * qJD(2) + qJD(4) * t58;
t24 = t58 * qJD(1);
t16 = t22 + t59;
t15 = t19 * qJD(1);
t10 = t45 * t15;
t1 = [0, 0, 0, 0, t68, qJ(2) * t68, t61, 0.2e1 * t65, ((-t40 * t63 + t41 * t92) * qJD(1) - t56) * qJD(2), t46 * t60, -0.2e1 * t91 * t76, -t93, t94, 0, t44 * t51 + t46 * t52, -t44 * t52 + t46 * t51, -t45 * t99 + (t46 * t80 - t72) * t21, (-t20 * t45 - t21 * t43) * t85 + (t106 - t12 * t45 + (t20 * t43 - t21 * t45) * qJD(5)) * t44, t105 + (-t74 + t97) * qJD(4) + t53, t62 + t104 + (-t98 - t109) * qJD(4), (-t29 - t79) * t86, (-t16 * t84 + t45 * t19) * t29 + (-t20 * t87 + t43 * t50 - t67 * t83 + t10) * t46 + (-t20 * t88 - t5 * t83 - t23 * t12 - t108 + (t23 * t101 - (-t100 * t23 + t16 * t45) * qJD(1) + t57) * qJD(4)) * t44, -(t16 * t83 + t43 * t19) * t29 + (-t21 * t87 + (qJD(5) * t67 - t15) * t43 + t50 * t45) * t46 + (-t21 * t88 + t5 * t84 + t23 * t11 - t107 + (t23 * t96 + (t16 * t43 + t23 * t95) * qJD(1) + t2) * qJD(4)) * t44; 0, 0, 0, 0, -t49, -t49 * qJ(2), -t40 * t49, -t41 * t49, t56 * qJD(1), 0, 0, 0, 0, 0, -t110 * t46 + t41 * t60, t110 * t44 + 0.2e1 * t41 * t66, 0, 0, 0, 0, 0, t41 * t73 + ((t44 * t81 - t45 * t82) * t29 - t20 * t85 - t44 * t12) * t40 + (-(-t100 * t41 + t40 * t45) * t29 + (-(-t100 * t40 - t41 * t45) * qJD(4) + t41 * t20) * t44) * qJD(1), t41 * t71 + (-(-t43 * t82 - t44 * t80) * t29 - t21 * t85 + t99) * t40 + ((t40 * t43 + t41 * t95) * t29 + ((t40 * t95 - t41 * t43) * qJD(4) + t41 * t21) * t44) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t93, 0, 0, 0, 0, 0, -t62 + t104 + (-t98 + t109) * qJD(4), -t105 + (-t74 - t97) * qJD(4) + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t49 * t46, t91 * t49, 0, 0, 0, t44 * t54, t46 * t54, -t21 * t96 + t106, (t11 + t103) * t45 + (t12 + t102) * t43, t71 + (t74 + (-t21 - t81) * t44) * qJD(1), -t73 + (-t29 * t100 + (t20 - t80) * t44) * qJD(1), t29 * t89, pkin(4) * t12 - t107 - (t24 * t45 - t43 * t8) * t29 + t9 * t20 + (-pkin(7) * t96 + t43 * t5) * qJD(5) + (-t57 * t44 + (pkin(7) * t86 + t46 * t5) * t43) * qJD(1), -pkin(4) * t11 + t108 + (t24 * t43 + t45 * t8) * t29 + t9 * t21 + (pkin(7) * t101 + t45 * t5) * qJD(5) + (t5 * t95 + (pkin(7) * t80 - t2) * t44) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t20, -t20 ^ 2 + t21 ^ 2, t11 - t103, t12 - t102, -t64, -t111 * t2 + t5 * t21 - t43 * t3 + t10, t111 * t57 - t43 * t15 - t5 * t20 - t45 * t3;];
tauc_reg = t1;
