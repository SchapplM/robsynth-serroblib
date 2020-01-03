% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR11_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:54
% EndTime: 2019-12-31 18:05:57
% DurationCPUTime: 0.98s
% Computational Cost: add. (1113->185), mult. (2208->275), div. (0->0), fcn. (1142->4), ass. (0->111)
t32 = qJD(1) * qJ(2) + qJD(3);
t28 = -pkin(6) * qJD(1) + t32;
t42 = sin(qJ(4));
t79 = qJD(1) * qJD(2);
t44 = cos(qJ(4));
t89 = qJD(4) * t44;
t15 = t28 * t89 + t42 * t79;
t62 = pkin(4) * t44 + pkin(7) * t42;
t20 = t62 * qJD(4) + qJD(3);
t17 = t20 * qJD(1);
t41 = sin(qJ(5));
t43 = cos(qJ(5));
t40 = pkin(1) + qJ(3);
t25 = t42 * pkin(4) - t44 * pkin(7) + t40;
t16 = t25 * qJD(1) - qJD(2);
t98 = t42 * t28;
t18 = qJD(4) * pkin(7) + t98;
t5 = t43 * t16 - t41 * t18;
t1 = qJD(5) * t5 + t43 * t15 + t41 * t17;
t82 = t42 * qJD(1);
t30 = qJD(5) + t82;
t121 = -t5 * t30 + t1;
t78 = qJD(1) * qJD(4);
t69 = t42 * t78;
t91 = qJD(1) * t44;
t71 = t43 * t91;
t83 = t41 * qJD(4);
t23 = t71 + t83;
t88 = qJD(5) * t23;
t10 = -t41 * t69 + t88;
t104 = t23 * t30;
t120 = t104 - t10;
t81 = t43 * qJD(4);
t21 = t41 * t91 - t81;
t53 = t21 * t30;
t84 = qJD(5) * t44;
t73 = t41 * t84;
t50 = t42 * t81 + t73;
t9 = t50 * qJD(1) - qJD(5) * t81;
t119 = -t9 + t53;
t6 = t41 * t16 + t43 * t18;
t2 = -qJD(5) * t6 - t41 * t15 + t43 * t17;
t118 = -t6 * t30 - t2;
t117 = qJD(1) * t40;
t37 = t42 ^ 2;
t38 = t44 ^ 2;
t93 = t37 + t38;
t116 = t93 * qJD(2);
t113 = t9 * t41;
t112 = t9 * t42;
t111 = t10 * t42;
t110 = t10 * t43;
t90 = qJD(4) * t42;
t14 = t28 * t90 - t44 * t79;
t109 = t14 * t41;
t108 = t14 * t43;
t94 = t44 * t28;
t19 = -qJD(4) * pkin(4) - t94;
t107 = t19 * t41;
t106 = t19 * t43;
t105 = t23 * t21;
t103 = t23 * t44;
t102 = t30 * t41;
t101 = t30 * t43;
t46 = qJD(1) ^ 2;
t100 = t37 * t46;
t99 = t41 * t42;
t97 = t42 * t43;
t96 = t43 * t44;
t95 = t44 * t10;
t45 = qJD(4) ^ 2;
t92 = -t45 - t46;
t87 = qJD(5) * t41;
t86 = qJD(5) * t42;
t85 = qJD(5) * t43;
t29 = -qJD(2) + t117;
t80 = qJD(2) - t29;
t77 = t30 * t99;
t76 = t30 * t97;
t75 = t44 * t46 * t42;
t74 = t30 * t87;
t72 = t43 * t84;
t35 = 0.2e1 * t79;
t70 = 0.2e1 * qJD(3) * qJD(1);
t68 = t44 * t78;
t67 = t29 + t117;
t66 = t30 + t82;
t65 = t80 * qJD(1);
t64 = qJD(1) + t86;
t63 = t42 * t68;
t39 = -pkin(6) + qJ(2);
t61 = -t39 * t86 + t20;
t60 = t30 * t85 + t41 * t68;
t59 = t41 * t6 + t43 * t5;
t58 = t41 * t5 - t43 * t6;
t56 = qJD(2) + t67;
t55 = qJD(1) * t38 - t30 * t42;
t52 = -t39 * t45 + t70;
t51 = -pkin(7) * t89 + t19 * t42;
t49 = qJD(2) * t42 + qJD(5) * t25 + t39 * t89;
t48 = t58 * qJD(5) - t1 * t41 - t2 * t43;
t47 = -t59 * qJD(5) + t1 * t43 - t2 * t41;
t34 = t38 * t46;
t24 = t62 * qJD(1);
t12 = t41 * t25 + t39 * t97;
t11 = t43 * t25 - t39 * t99;
t8 = t41 * t24 + t43 * t94;
t7 = t43 * t24 - t41 * t94;
t4 = -t49 * t41 + t61 * t43;
t3 = t61 * t41 + t49 * t43;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, qJ(2) * t35, 0, 0, 0, 0, 0, 0, 0, t35, t70, t32 * qJD(2) + t29 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t40) * qJD(1), -0.2e1 * t63, 0.2e1 * (t37 - t38) * t78, -t45 * t42, 0.2e1 * t63, -t45 * t44, 0, t52 * t42 + t56 * t89, t52 * t44 - t56 * t90, -t93 * t35, t67 * qJD(3) + (qJD(1) * t39 + t28) * t116, -t50 * t23 - t9 * t96, (t21 * t43 + t23 * t41) * t90 + (-t110 + t113 + (t21 * t41 - t23 * t43) * qJD(5)) * t44, -t30 * t73 - t112 + (t55 * t43 + t103) * qJD(4), t41 * t95 + (-t42 * t83 + t72) * t21, -t30 * t72 - t111 + (-t21 * t44 - t55 * t41) * qJD(4), t66 * t89, t4 * t30 + (t2 + (t21 * t39 - t107) * qJD(4)) * t42 + (t19 * t85 - qJD(2) * t21 - t10 * t39 + t109 + (qJD(1) * t11 + t5) * qJD(4)) * t44, -t3 * t30 + (-t1 + (t23 * t39 - t106) * qJD(4)) * t42 + (-t19 * t87 - qJD(2) * t23 + t108 + t39 * t9 + (-qJD(1) * t12 - t6) * qJD(4)) * t44, -t12 * t10 + t11 * t9 - t3 * t21 - t4 * t23 + t48 * t44 + t59 * t90, -t14 * t44 * t39 + t1 * t12 + t2 * t11 + t6 * t3 + t5 * t4 + (-qJD(2) * t44 + t39 * t90) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t46 * qJ(2), 0, 0, 0, 0, 0, 0, 0, -t46, 0, (-qJD(3) - t32) * qJD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t68, 0.2e1 * t69, t34 + t100, (-t93 * t28 - qJD(3)) * qJD(1), 0, 0, 0, 0, 0, 0, t74 + (t77 + (t21 - t81) * t44) * qJD(1), (t76 + t103) * qJD(1) + t60, t119 * t43 - t120 * t41, (t19 * t44 + t58 * t42) * qJD(1) + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t65, 0, 0, 0, 0, 0, 0, t92 * t42, t92 * t44, 0, (-t29 + t116) * qJD(1), 0, 0, 0, 0, 0, 0, -t95 - t64 * t101 + (-t66 * t44 * t41 + t21 * t42) * qJD(4), t44 * t9 + t64 * t102 + (-t30 * t96 + (t23 - t71) * t42) * qJD(4), (-t21 * t89 + t64 * t23 - t111) * t43 + (t64 * t21 + t23 * t89 - t112) * t41, -t59 * qJD(1) + (-t58 * qJD(4) - t14) * t44 + (qJD(4) * t19 + t47) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t34 - t100, 0, -t75, 0, 0, t44 * t65, -t80 * t82, 0, 0, t23 * t101 - t113, (-t9 - t53) * t43 + (-t10 - t104) * t41, (t76 - t103) * qJD(1) + t60, t41 * t53 - t110, -t74 + (-t77 + (t21 + t81) * t44) * qJD(1), -t30 * t91, -t21 * t98 - pkin(4) * t10 - t108 - t7 * t30 + (-pkin(7) * t101 + t107) * qJD(5) + (t51 * t41 - t44 * t5) * qJD(1), -t23 * t98 + pkin(4) * t9 + t109 + t8 * t30 + (pkin(7) * t102 + t106) * qJD(5) + (t51 * t43 + t44 * t6) * qJD(1), t8 * t21 + t7 * t23 + ((-t10 + t88) * pkin(7) + t121) * t43 + ((qJD(5) * t21 - t9) * pkin(7) + t118) * t41, -t14 * pkin(4) + t47 * pkin(7) - t19 * t98 - t5 * t7 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, -t21 ^ 2 + t23 ^ 2, t119, -t105, t120, t68, -t19 * t23 - t118, t19 * t21 - t121, 0, 0;];
tauc_reg = t13;
