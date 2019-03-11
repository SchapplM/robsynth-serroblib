% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:39:51
% EndTime: 2019-03-09 00:39:56
% DurationCPUTime: 1.53s
% Computational Cost: add. (2860->194), mult. (6966->333), div. (0->0), fcn. (7109->12), ass. (0->129)
t93 = cos(qJ(5));
t136 = qJD(5) * t93;
t89 = sin(qJ(4));
t146 = t89 * t93;
t88 = sin(qJ(5));
t94 = cos(qJ(4));
t156 = ((t88 * t94 + t146) * qJD(4) + t89 * t136) * pkin(3);
t92 = cos(qJ(6));
t84 = t92 ^ 2;
t87 = sin(qJ(6));
t141 = t87 ^ 2 - t84;
t117 = t141 * qJD(6);
t154 = pkin(8) + pkin(9);
t90 = sin(qJ(3));
t68 = t154 * t90;
t95 = cos(qJ(3));
t69 = t154 * t95;
t110 = -t68 * t94 - t69 * t89;
t120 = qJD(3) * t154;
t63 = t90 * t120;
t64 = t95 * t120;
t33 = -qJD(4) * t110 + t94 * t63 + t89 * t64;
t62 = t89 * t95 + t90 * t94;
t37 = -pkin(10) * t62 + t110;
t109 = t68 * t89 - t69 * t94;
t61 = t89 * t90 - t94 * t95;
t38 = -pkin(10) * t61 - t109;
t30 = t37 * t88 + t38 * t93;
t133 = t95 * qJD(3);
t134 = t90 * qJD(3);
t138 = qJD(4) * t94;
t139 = qJD(4) * t89;
t116 = t133 * t89 + t134 * t94 + t138 * t90 + t139 * t95;
t97 = -pkin(10) * t116 - t33;
t34 = qJD(4) * t109 + t89 * t63 - t94 * t64;
t49 = (-qJD(3) - qJD(4)) * t61;
t98 = t49 * pkin(10) - t34;
t11 = qJD(5) * t30 + t88 * t97 + t93 * t98;
t29 = -t37 * t93 + t38 * t88;
t82 = qJD(6) * t92;
t153 = t11 * t87 + t29 * t82;
t47 = t61 * t93 + t62 * t88;
t21 = -qJD(5) * t47 - t116 * t88 + t93 * t49;
t48 = -t61 * t88 + t62 * t93;
t152 = t48 * t21;
t151 = t48 * t87;
t150 = t48 * t92;
t85 = sin(pkin(6));
t91 = sin(qJ(2));
t149 = t85 * t91;
t96 = cos(qJ(2));
t148 = t85 * t96;
t22 = qJD(5) * t48 + t116 * t93 + t88 * t49;
t147 = t87 * t22;
t145 = t92 * t21;
t144 = t92 * t22;
t137 = qJD(5) * t88;
t79 = pkin(3) * t94 + pkin(4);
t44 = t137 * t79 + t156;
t132 = t88 * t89 * pkin(3);
t55 = -t79 * t93 - pkin(5) + t132;
t143 = t44 * t87 + t55 * t82;
t126 = pkin(4) * t137;
t78 = -pkin(4) * t93 - pkin(5);
t142 = t126 * t87 + t78 * t82;
t140 = qJD(2) * t91;
t135 = qJD(6) * t87;
t131 = -0.2e1 * pkin(2) * qJD(3);
t130 = pkin(5) * t135;
t129 = pkin(5) * t82;
t81 = pkin(3) * t134;
t128 = pkin(3) * t139;
t127 = pkin(3) * t138;
t125 = pkin(4) * t136;
t123 = t85 * t140;
t122 = qJD(2) * t148;
t121 = t87 * t82;
t80 = -pkin(3) * t95 - pkin(2);
t119 = -0.4e1 * t87 * t150;
t52 = t55 * t135;
t118 = -t44 * t92 + t52;
t54 = pkin(4) * t61 + t80;
t28 = pkin(5) * t47 - pkin(11) * t48 + t54;
t115 = t28 * t92 - t30 * t87;
t114 = t28 * t87 + t30 * t92;
t86 = cos(pkin(6));
t57 = -t149 * t90 + t86 * t95;
t58 = t149 * t95 + t86 * t90;
t39 = t57 * t94 - t58 * t89;
t40 = t57 * t89 + t58 * t94;
t32 = t39 * t88 + t40 * t93;
t56 = pkin(3) * t146 + t79 * t88 + pkin(11);
t113 = t47 * t56 - t48 * t55;
t77 = pkin(4) * t88 + pkin(11);
t112 = t47 * t77 - t48 * t78;
t65 = t78 * t135;
t107 = -t126 * t92 + t65;
t106 = t148 * t92 + t32 * t87;
t105 = t148 * t87 - t32 * t92;
t104 = t21 * t87 + t48 * t82;
t103 = t135 * t48 - t145;
t102 = t135 * t47 - t144;
t42 = pkin(4) * t116 + t81;
t43 = -t79 * t136 - t93 * t127 + (qJD(4) + qJD(5)) * t132;
t101 = t21 * t55 - t22 * t56 + t43 * t47 + t44 * t48;
t50 = qJD(3) * t57 + t122 * t95;
t51 = -qJD(3) * t58 - t122 * t90;
t100 = qJD(4) * t40 + t89 * t50 - t94 * t51;
t99 = t21 * t78 - t22 * t77 + (-t47 * t93 + t48 * t88) * qJD(5) * pkin(4);
t72 = 0.2e1 * t121;
t60 = -0.2e1 * t117;
t46 = t48 ^ 2;
t31 = -t39 * t93 + t40 * t88;
t26 = t29 * t135;
t20 = -qJD(4) * t39 - t94 * t50 - t89 * t51;
t16 = t47 * t82 + t147;
t14 = -t117 * t48 + t145 * t87;
t13 = pkin(5) * t22 - pkin(11) * t21 + t42;
t12 = qJD(6) * t119 - t141 * t21;
t10 = -t136 * t37 + t137 * t38 + t88 * t98 - t93 * t97;
t8 = qJD(5) * t32 + t100 * t93 - t88 * t20;
t7 = t100 * t88 - t136 * t39 + t137 * t40 + t20 * t93;
t6 = t135 * t31 - t8 * t92;
t5 = t31 * t82 + t8 * t87;
t4 = qJD(6) * t105 + t123 * t92 + t87 * t7;
t3 = qJD(6) * t106 - t123 * t87 + t92 * t7;
t2 = -qJD(6) * t114 + t87 * t10 + t92 * t13;
t1 = -qJD(6) * t115 + t92 * t10 - t87 * t13;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t123, -t122, 0, 0, 0, 0, 0 (-t134 * t96 - t140 * t95) * t85 (-t133 * t96 + t140 * t90) * t85, 0, 0, 0, 0, 0 (-t116 * t96 + t140 * t61) * t85 (t140 * t62 - t49 * t96) * t85, 0, 0, 0, 0, 0 (t140 * t47 - t22 * t96) * t85 (t140 * t48 - t21 * t96) * t85, 0, 0, 0, 0, 0, t104 * t31 - t106 * t22 + t151 * t8 + t4 * t47, -t103 * t31 + t105 * t22 + t150 * t8 + t3 * t47; 0, 0, 0, 0, 0.2e1 * t90 * t133, 0.2e1 * (-t90 ^ 2 + t95 ^ 2) * qJD(3), 0, 0, 0, t90 * t131, t95 * t131, 0.2e1 * t62 * t49, -0.2e1 * t116 * t62 - 0.2e1 * t49 * t61, 0, 0, 0, 0.2e1 * t116 * t80 + 0.2e1 * t61 * t81, 0.2e1 * t49 * t80 + 0.2e1 * t62 * t81, 0.2e1 * t152, -0.2e1 * t21 * t47 - 0.2e1 * t22 * t48, 0, 0, 0, 0.2e1 * t22 * t54 + 0.2e1 * t42 * t47, 0.2e1 * t21 * t54 + 0.2e1 * t42 * t48, -0.2e1 * t121 * t46 + 0.2e1 * t152 * t84, 0.2e1 * t117 * t46 + t119 * t21, -0.2e1 * t103 * t47 + 0.2e1 * t144 * t48, -0.2e1 * t104 * t47 - 0.2e1 * t147 * t48, 0.2e1 * t47 * t22, 0.2e1 * t104 * t29 + 0.2e1 * t11 * t151 + 0.2e1 * t115 * t22 + 0.2e1 * t2 * t47, 0.2e1 * t1 * t47 - 0.2e1 * t103 * t29 + 0.2e1 * t11 * t150 - 0.2e1 * t114 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t50, 0, 0, 0, 0, 0, -t100, t20, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, t133, -t134, 0, -pkin(8) * t133, pkin(8) * t134, 0, 0, t49, -t116, 0, t34, t33, 0, 0, t21, -t22, 0, -t11, t10, t14, t12, t16, -t102, 0, t26 + (-qJD(6) * t113 - t11) * t92 + t101 * t87, t101 * t92 + t113 * t135 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t128, -0.2e1 * t127, 0, 0, 0, 0, 0, -0.2e1 * t44, 0.2e1 * t43, t72, t60, 0, 0, 0, 0.2e1 * t118, 0.2e1 * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t20, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t116, 0, t34, t33, 0, 0, t21, -t22, 0, -t11, t10, t14, t12, t16, -t102, 0, t26 + (-qJD(6) * t112 - t11) * t92 + t99 * t87, t112 * t135 + t92 * t99 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t127, 0, 0, 0, 0, 0 (-pkin(4) - t79) * t137 - t156, t43 - t125, t72, t60, 0, 0, 0, t52 + t65 + (-t44 - t126) * t92, t142 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t126, -0.2e1 * t125, t72, t60, 0, 0, 0, 0.2e1 * t107, 0.2e1 * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t22, 0, -t11, t10, t14, t12, t16, -t102, 0, t26 + (-pkin(5) * t21 - pkin(11) * t22) * t87 + (-t11 + (-pkin(5) * t48 - pkin(11) * t47) * qJD(6)) * t92, pkin(5) * t103 + pkin(11) * t102 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43, t72, t60, 0, 0, 0, t118 - t130, -t129 + t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t125, t72, t60, 0, 0, 0, t107 - t130, -t129 + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t60, 0, 0, 0, -0.2e1 * t130, -0.2e1 * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, -t104, t22, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t135, 0, t43 * t87 - t56 * t82, t135 * t56 + t43 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t135, 0, -t125 * t87 - t77 * t82, -t125 * t92 + t135 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t135, 0, -pkin(11) * t82, pkin(11) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
