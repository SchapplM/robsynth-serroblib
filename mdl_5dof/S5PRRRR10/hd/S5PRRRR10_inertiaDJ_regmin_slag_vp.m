% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:06
% EndTime: 2019-12-05 17:26:12
% DurationCPUTime: 1.56s
% Computational Cost: add. (1237->212), mult. (3963->426), div. (0->0), fcn. (3886->12), ass. (0->127)
t69 = sin(qJ(4));
t147 = -0.4e1 * t69;
t64 = sin(pkin(6));
t74 = cos(qJ(3));
t139 = t64 * t74;
t116 = pkin(8) * t139;
t70 = sin(qJ(3));
t144 = pkin(2) * t70;
t66 = cos(pkin(6));
t42 = t116 + (pkin(9) + t144) * t66;
t43 = (-pkin(3) * t74 - pkin(9) * t70 - pkin(2)) * t64;
t73 = cos(qJ(4));
t146 = t73 * t42 + t69 * t43;
t72 = cos(qJ(5));
t62 = t72 ^ 2;
t68 = sin(qJ(5));
t131 = t68 ^ 2 - t62;
t103 = t131 * qJD(5);
t44 = (pkin(3) * t70 - pkin(9) * t74) * t64 * qJD(3);
t128 = qJD(3) * t70;
t109 = t64 * t128;
t127 = qJD(3) * t74;
t45 = -t66 * pkin(2) * t127 + pkin(8) * t109;
t14 = -qJD(4) * t146 + t73 * t44 + t69 * t45;
t145 = 0.2e1 * t64;
t143 = pkin(9) * t64;
t142 = t68 * pkin(9);
t108 = t64 * t127;
t140 = t64 * t70;
t48 = t69 * t140 - t73 * t66;
t30 = -qJD(4) * t48 + t73 * t108;
t49 = t73 * t140 + t69 * t66;
t32 = t72 * t139 + t68 * t49;
t16 = -qJD(5) * t32 + t68 * t109 + t72 * t30;
t141 = t16 * t68;
t138 = t69 * t72;
t71 = sin(qJ(2));
t137 = t70 * t71;
t75 = cos(qJ(2));
t136 = t70 * t75;
t135 = t71 * t74;
t134 = t72 * t73;
t133 = t74 * t75;
t61 = t69 ^ 2;
t130 = -t73 ^ 2 + t61;
t65 = sin(pkin(5));
t129 = qJD(2) * t65;
t126 = qJD(4) * t68;
t125 = qJD(4) * t69;
t124 = qJD(4) * t72;
t123 = qJD(4) * t73;
t122 = qJD(4) * t74;
t121 = qJD(5) * t68;
t120 = qJD(5) * t72;
t119 = qJD(5) * t73;
t118 = t73 * t142;
t117 = pkin(9) * t134;
t115 = -0.2e1 * pkin(3) * qJD(4);
t114 = -0.2e1 * pkin(4) * qJD(5);
t113 = t68 * t139;
t59 = t64 ^ 2;
t112 = t59 * t127;
t111 = t68 * t119;
t110 = t72 * t119;
t107 = t71 * t129;
t106 = t68 * t120;
t105 = t69 * t123;
t104 = t72 * t123;
t102 = t130 * qJD(4);
t101 = 0.2e1 * t105;
t100 = t59 * t107;
t99 = t64 * t107;
t98 = t70 * t112;
t97 = t68 * t104;
t96 = -t73 * pkin(4) - t69 * pkin(10);
t95 = pkin(4) * t69 - pkin(10) * t73;
t67 = cos(pkin(5));
t89 = t66 * t136 + t135;
t29 = t67 * t140 + t65 * t89;
t47 = -t65 * t75 * t64 + t67 * t66;
t22 = t29 * t73 + t47 * t69;
t90 = t66 * t133 - t137;
t28 = -t67 * t139 - t65 * t90;
t12 = t22 * t72 + t28 * t68;
t11 = -t22 * t68 + t28 * t72;
t41 = pkin(8) * t140 + (-pkin(2) * t74 - pkin(3)) * t66;
t23 = t48 * pkin(4) - t49 * pkin(10) + t41;
t25 = -pkin(10) * t139 + t146;
t8 = t68 * t23 + t72 * t25;
t21 = t29 * t69 - t47 * t73;
t33 = t72 * t49 - t113;
t94 = -t32 * t72 - t33 * t68;
t92 = -t69 * t42 + t73 * t43;
t55 = -pkin(3) + t96;
t40 = t68 * t55 + t117;
t20 = t67 * t108 + (t90 * qJD(3) + (-t66 * t137 + t133) * qJD(2)) * t65;
t5 = t22 * qJD(4) + t20 * t69 - t73 * t99;
t88 = t21 * t120 + t5 * t68;
t87 = t21 * t121 - t5 * t72;
t86 = t95 * t68;
t10 = -pkin(4) * t109 - t14;
t24 = pkin(4) * t139 - t92;
t85 = t10 * t68 + t24 * t120;
t84 = -t10 * t72 + t24 * t121;
t31 = qJD(4) * t49 + t69 * t108;
t83 = t48 * t120 + t68 * t31;
t82 = t48 * t121 - t72 * t31;
t13 = -t43 * t123 + t42 * t125 - t69 * t44 + t73 * t45;
t81 = t69 * t122 + t73 * t128;
t80 = -t73 * t122 + t69 * t128;
t79 = -t69 * t121 + t104;
t78 = t69 * t124 + t111;
t46 = (t66 * t144 + t116) * qJD(3);
t77 = pkin(10) * t109 - t13;
t76 = t31 * pkin(4) - t30 * pkin(10) + t46;
t39 = t72 * t55 - t118;
t27 = -t40 * qJD(5) + (t69 * t142 + t72 * t95) * qJD(4);
t26 = pkin(9) * t78 - qJD(4) * t86 - t55 * t120;
t19 = t67 * t109 + (t89 * qJD(3) + (t66 * t135 + t136) * qJD(2)) * t65;
t17 = -qJD(5) * t113 - t72 * t109 + t49 * t120 + t68 * t30;
t7 = t72 * t23 - t68 * t25;
t6 = -qJD(4) * t21 + t20 * t73 + t69 * t99;
t4 = t11 * qJD(5) + t19 * t68 + t6 * t72;
t3 = -t12 * qJD(5) + t19 * t72 - t6 * t68;
t2 = -t8 * qJD(5) - t68 * t77 + t72 * t76;
t1 = -t23 * t120 + t25 * t121 - t68 * t76 - t72 * t77;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t107, -t75 * t129, 0, 0, 0, 0, 0, -t74 * t100 + t47 * t109 - t19 * t66, t70 * t100 + t47 * t108 - t20 * t66, 0, 0, 0, 0, 0, t19 * t48 + t28 * t31 + (-t21 * t128 + t5 * t74) * t64, t19 * t49 + t28 * t30 + (-t22 * t128 + t6 * t74) * t64, 0, 0, 0, 0, 0, t11 * t31 + t21 * t17 + t3 * t48 + t5 * t32, -t12 * t31 + t21 * t16 + t5 * t33 - t4 * t48; 0, 0, 0, 0, 0.2e1 * t98, 0.2e1 * (-t70 ^ 2 + t74 ^ 2) * t59 * qJD(3), 0.2e1 * t66 * t108, -0.2e1 * t66 * t109, 0, -0.2e1 * t59 * pkin(2) * t128 - 0.2e1 * t46 * t66, -0.2e1 * pkin(2) * t112 + 0.2e1 * t45 * t66, 0.2e1 * t49 * t30, -0.2e1 * t30 * t48 - 0.2e1 * t49 * t31, (t49 * t128 - t30 * t74) * t145, (-t48 * t128 + t31 * t74) * t145, -0.2e1 * t98, 0.2e1 * t41 * t31 + 0.2e1 * t46 * t48 + 0.2e1 * (t92 * t128 - t14 * t74) * t64, 0.2e1 * t41 * t30 + 0.2e1 * t46 * t49 + 0.2e1 * (-t128 * t146 - t13 * t74) * t64, 0.2e1 * t33 * t16, -0.2e1 * t16 * t32 - 0.2e1 * t33 * t17, 0.2e1 * t16 * t48 + 0.2e1 * t33 * t31, -0.2e1 * t17 * t48 - 0.2e1 * t32 * t31, 0.2e1 * t48 * t31, 0.2e1 * t10 * t32 + 0.2e1 * t24 * t17 + 0.2e1 * t2 * t48 + 0.2e1 * t7 * t31, 0.2e1 * t1 * t48 + 0.2e1 * t10 * t33 + 0.2e1 * t24 * t16 - 0.2e1 * t8 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, 0, 0, 0, 0, t28 * t125 - t19 * t73, t28 * t123 + t19 * t69, 0, 0, 0, 0, 0, (t21 * t126 - t3) * t73 + (qJD(4) * t11 + t88) * t69, (t21 * t124 + t4) * t73 + (-qJD(4) * t12 - t87) * t69; 0, 0, 0, 0, 0, 0, t108, -t109, 0, -t46, t45, t49 * t123 + t30 * t69, t30 * t73 - t69 * t31 + (-t48 * t73 - t49 * t69) * qJD(4), t80 * t64, t81 * t64, 0, -pkin(3) * t31 + t41 * t125 - t80 * t143 - t46 * t73, -pkin(3) * t30 + t41 * t123 - t81 * t143 + t46 * t69, t16 * t138 + t33 * t79, t94 * t123 + (-t141 - t17 * t72 + (t32 * t68 - t33 * t72) * qJD(5)) * t69, (t48 * t124 - t16) * t73 + (qJD(4) * t33 - t82) * t69, (-t48 * t126 + t17) * t73 + (-qJD(4) * t32 - t83) * t69, t48 * t125 - t31 * t73, t27 * t48 + t39 * t31 + (-t2 + (pkin(9) * t32 + t24 * t68) * qJD(4)) * t73 + (pkin(9) * t17 + qJD(4) * t7 + t85) * t69, t26 * t48 - t40 * t31 + (-t1 + (pkin(9) * t33 + t24 * t72) * qJD(4)) * t73 + (pkin(9) * t16 - qJD(4) * t8 - t84) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, -0.2e1 * t102, 0, 0, 0, t69 * t115, t73 * t115, 0.2e1 * t62 * t105 - 0.2e1 * t61 * t106, 0.2e1 * t61 * t103 + t97 * t147, 0.2e1 * t69 * t111 + 0.2e1 * t130 * t124, -0.2e1 * t68 * t102 + 0.2e1 * t69 * t110, -0.2e1 * t105, 0.2e1 * t39 * t125 - 0.2e1 * t27 * t73 + 0.2e1 * (t68 * t101 + t61 * t120) * pkin(9), -0.2e1 * t40 * t125 - 0.2e1 * t26 * t73 + 0.2e1 * (t72 * t101 - t61 * t121) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, t87, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31, t109, t14, t13, t33 * t120 + t141, t94 * qJD(5) + t16 * t72 - t68 * t17, t83, -t82, 0, -pkin(4) * t17 - pkin(10) * t83 + t84, -pkin(4) * t16 + pkin(10) * t82 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, -t125, 0, -pkin(9) * t123, pkin(9) * t125, -t69 * t103 + t97, t106 * t147 - t131 * t123, t68 * t125 - t110, t78, 0, (pkin(10) * t134 + (-t72 * pkin(4) + t142) * t69) * qJD(5) + (t96 * t68 - t117) * qJD(4), (pkin(9) * t138 + t86) * qJD(5) + (t96 * t72 + t118) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t106, -0.2e1 * t103, 0, 0, 0, t68 * t114, t72 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, t31, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t69 * t120 - t68 * t123, t125, t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t121, 0, -pkin(10) * t120, pkin(10) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
