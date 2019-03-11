% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:41:20
% EndTime: 2019-03-09 11:41:23
% DurationCPUTime: 1.37s
% Computational Cost: add. (3814->211), mult. (8199->357), div. (0->0), fcn. (8290->8), ass. (0->126)
t153 = cos(qJ(4));
t144 = -qJ(3) - pkin(7);
t98 = sin(qJ(2));
t80 = t144 * t98;
t100 = cos(qJ(2));
t82 = t144 * t100;
t94 = sin(pkin(10));
t95 = cos(pkin(10));
t47 = t95 * t80 + t94 * t82;
t73 = t94 * t100 + t95 * t98;
t42 = -t73 * pkin(8) + t47;
t48 = t94 * t80 - t95 * t82;
t72 = t95 * t100 - t94 * t98;
t43 = t72 * pkin(8) + t48;
t97 = sin(qJ(4));
t29 = t153 * t43 + t97 * t42;
t99 = cos(qJ(5));
t26 = t99 * t29;
t106 = t153 * t72 - t97 * t73;
t46 = t153 * t73 + t97 * t72;
t124 = -t100 * pkin(2) - pkin(1);
t57 = -t72 * pkin(3) + t124;
t27 = -pkin(4) * t106 - t46 * pkin(9) + t57;
t96 = sin(qJ(5));
t141 = t96 * t27 + t26;
t125 = t95 * pkin(2) + pkin(3);
t155 = t94 * pkin(2);
t157 = -t153 * t125 + t97 * t155;
t90 = qJD(5) * t99;
t127 = t46 * t90;
t105 = t73 * qJD(2);
t133 = qJD(2) * t100;
t136 = qJD(2) * t98;
t71 = t95 * t133 - t94 * t136;
t30 = t106 * qJD(4) - t97 * t105 + t153 * t71;
t109 = t96 * t30 + t127;
t92 = t96 ^ 2;
t93 = t99 ^ 2;
t139 = t92 - t93;
t116 = t139 * qJD(5);
t111 = -qJ(6) * t30 - qJD(6) * t46;
t118 = qJD(2) * t144;
t66 = t100 * qJD(3) + t98 * t118;
t67 = -t98 * qJD(3) + t100 * t118;
t41 = t95 * t66 + t94 * t67;
t101 = -pkin(8) * t105 + t41;
t40 = -t94 * t66 + t95 * t67;
t104 = -t71 * pkin(8) + t40;
t119 = qJD(4) * t153;
t135 = qJD(4) * t97;
t11 = -t153 * t101 - t97 * t104 - t42 * t119 + t43 * t135;
t31 = t46 * qJD(4) + t153 * t105 + t97 * t71;
t88 = pkin(2) * t136;
t52 = pkin(3) * t105 + t88;
t16 = t31 * pkin(4) - t30 * pkin(9) + t52;
t122 = t96 * t11 + t99 * t16;
t138 = qJ(6) * t46;
t1 = t31 * pkin(5) + t111 * t99 + (-t26 + (-t27 + t138) * t96) * qJD(5) + t122;
t131 = -t99 * t11 + t96 * t16 + t27 * t90;
t3 = -qJ(6) * t127 + (-qJD(5) * t29 + t111) * t96 + t131;
t121 = t99 * t27 - t96 * t29;
t91 = t99 * qJ(6);
t7 = -pkin(5) * t106 - t46 * t91 + t121;
t8 = -t96 * t138 + t141;
t156 = -t1 * t99 - t3 * t96 + (t7 * t96 - t8 * t99) * qJD(5);
t154 = t99 * pkin(5);
t152 = t46 * t30;
t151 = t46 * t96;
t150 = t46 * t99;
t148 = t96 * t31;
t147 = t99 * t30;
t146 = t99 * t31;
t61 = t157 * qJD(4);
t145 = t99 * t61;
t143 = -qJ(6) - pkin(9);
t12 = t97 * t101 - t153 * t104 + t43 * t119 + t42 * t135;
t28 = -t153 * t42 + t97 * t43;
t142 = t12 * t96 + t28 * t90;
t103 = t97 * t125 + t153 * t155;
t62 = t103 * qJD(4);
t64 = -pkin(4) + t157;
t140 = t62 * t96 + t64 * t90;
t65 = pkin(9) + t103;
t137 = -qJ(6) - t65;
t134 = qJD(5) * t96;
t132 = -0.2e1 * pkin(1) * qJD(2);
t130 = pkin(4) * t134;
t129 = pkin(4) * t90;
t87 = pkin(5) * t134;
t128 = pkin(5) * t90;
t126 = t96 * t90;
t123 = -0.4e1 * t96 * t150;
t120 = t64 * t134 - t62 * t99;
t117 = qJD(5) * t143;
t115 = qJD(5) * t137;
t112 = -t106 * t65 - t46 * t64;
t108 = t46 * t134 - t147;
t20 = -t106 * t90 + t148;
t107 = -t106 * t134 - t146;
t102 = -t106 * t61 + t30 * t64 - t31 * t65 + t46 * t62;
t89 = t99 * qJD(6);
t86 = -pkin(4) - t154;
t83 = 0.2e1 * t126;
t81 = t99 * pkin(9) + t91;
t79 = t143 * t96;
t77 = -0.2e1 * t116;
t69 = -t96 * qJD(6) + t99 * t117;
t68 = t96 * t117 + t89;
t60 = t68 * t99;
t56 = t64 - t154;
t51 = t87 + t62;
t50 = t99 * t65 + t91;
t49 = t137 * t96;
t44 = t46 ^ 2;
t37 = (-qJD(6) + t61) * t96 + t99 * t115;
t36 = t96 * t115 - t145 + t89;
t35 = t36 * t99;
t22 = t28 * t134;
t18 = pkin(5) * t151 + t28;
t17 = -t46 * t116 + t96 * t147;
t13 = qJD(5) * t123 - t139 * t30;
t6 = t109 * pkin(5) + t12;
t5 = -t141 * qJD(5) + t122;
t4 = t29 * t134 - t131;
t2 = t3 * t99;
t9 = [0, 0, 0, 0.2e1 * t98 * t133, 0.2e1 * (t100 ^ 2 - t98 ^ 2) * qJD(2), 0, 0, 0, t98 * t132, t100 * t132, -0.2e1 * t48 * t105 - 0.2e1 * t40 * t73 + 0.2e1 * t41 * t72 - 0.2e1 * t47 * t71, 0.2e1 * t124 * t88 + 0.2e1 * t47 * t40 + 0.2e1 * t48 * t41, 0.2e1 * t152, 0.2e1 * t106 * t30 - 0.2e1 * t46 * t31, 0, 0, 0, -0.2e1 * t106 * t52 + 0.2e1 * t57 * t31, 0.2e1 * t57 * t30 + 0.2e1 * t52 * t46, -0.2e1 * t126 * t44 + 0.2e1 * t93 * t152, 0.2e1 * t44 * t116 + t123 * t30, 0.2e1 * t106 * t108 + 0.2e1 * t146 * t46, 0.2e1 * t106 * t109 - 0.2e1 * t46 * t148, -0.2e1 * t106 * t31, -0.2e1 * t106 * t5 + 0.2e1 * t109 * t28 + 0.2e1 * t12 * t151 + 0.2e1 * t121 * t31, -0.2e1 * t106 * t4 - 0.2e1 * t108 * t28 + 0.2e1 * t12 * t150 - 0.2e1 * t141 * t31, 0.2e1 * (-t7 * t99 - t8 * t96) * t30 + 0.2e1 * t156 * t46, 0.2e1 * t7 * t1 + 0.2e1 * t18 * t6 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, t133, -t136, 0, -pkin(7) * t133, pkin(7) * t136 (-t94 * t105 - t95 * t71) * pkin(2) (t40 * t95 + t41 * t94) * pkin(2), 0, 0, t30, -t31, 0, -t12, t11, t17, t13, t20, -t107, 0, t22 + (-qJD(5) * t112 - t12) * t99 + t102 * t96, t102 * t99 + t112 * t134 + t142, t2 + (-t30 * t49 - t37 * t46 + (-t46 * t50 - t7) * qJD(5)) * t99 + (-t30 * t50 - t36 * t46 - t1 + (t46 * t49 - t8) * qJD(5)) * t96, t1 * t49 + t18 * t51 + t3 * t50 + t8 * t36 + t7 * t37 + t6 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t62, 0.2e1 * t61, t83, t77, 0, 0, 0, 0.2e1 * t120, 0.2e1 * t140, -0.2e1 * t37 * t96 + 0.2e1 * t35 + 0.2e1 * (-t49 * t99 - t50 * t96) * qJD(5), 0.2e1 * t50 * t36 + 0.2e1 * t49 * t37 + 0.2e1 * t56 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, 0, 0, 0, 0, 0, t31, t30, 0, 0, 0, 0, 0, -t107, -t20 (-t92 - t93) * t30, -t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t96 + t37 * t99 + (-t49 * t96 + t50 * t99) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t31, 0, -t12, t11, t17, t13, t20, -t107, 0, t22 + (-pkin(4) * t30 - pkin(9) * t31) * t96 + (-t12 + (-pkin(4) * t46 + pkin(9) * t106) * qJD(5)) * t99, pkin(4) * t108 + pkin(9) * t107 + t142, t2 + (-t30 * t79 - t46 * t69 + (-t46 * t81 - t7) * qJD(5)) * t99 + (-t30 * t81 - t46 * t68 - t1 + (t46 * t79 - t8) * qJD(5)) * t96, t1 * t79 + t18 * t87 + t3 * t81 + t6 * t86 + t8 * t68 + t7 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t61, t83, t77, 0, 0, 0, t120 - t130, -t129 + t140, t35 + t60 + (-t37 - t69) * t96 + ((-t49 - t79) * t99 + (-t50 - t81) * t96) * qJD(5), t36 * t81 + t37 * t79 + t49 * t69 + t50 * t68 + t51 * t86 + t56 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t68 + t99 * t69 + (-t79 * t96 + t81 * t99) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t77, 0, 0, 0, -0.2e1 * t130, -0.2e1 * t129, -0.2e1 * t69 * t96 + 0.2e1 * t60 + 0.2e1 * (-t79 * t99 - t81 * t96) * qJD(5), 0.2e1 * t81 * t68 + 0.2e1 * t79 * t69 + 0.2e1 * t86 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t109, t31, t5, t4, t108 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t134, 0, t96 * t61 - t65 * t90, t134 * t65 + t145, -t128, t37 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t90, 0, -t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t134, 0, -pkin(9) * t90, pkin(9) * t134, -t128, t69 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
