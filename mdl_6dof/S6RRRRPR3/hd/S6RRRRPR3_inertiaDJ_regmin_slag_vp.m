% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:04:00
% EndTime: 2019-03-09 22:04:04
% DurationCPUTime: 1.65s
% Computational Cost: add. (3592->217), mult. (7997->343), div. (0->0), fcn. (7611->8), ass. (0->127)
t80 = sin(qJ(6));
t78 = t80 ^ 2;
t83 = cos(qJ(6));
t138 = -t83 ^ 2 + t78;
t122 = t138 * qJD(6);
t153 = qJD(2) + qJD(3);
t147 = cos(qJ(3));
t124 = t147 * qJD(3);
t121 = pkin(2) * t124;
t145 = sin(qJ(3));
t123 = t145 * qJD(3);
t146 = cos(qJ(4));
t125 = qJD(4) * t146;
t81 = sin(qJ(4));
t150 = pkin(2) * t81;
t72 = t147 * pkin(2) + pkin(3);
t39 = (t145 * qJD(4) + t123) * t150 - t146 * t121 - t72 * t125;
t86 = 2 * qJD(5);
t152 = pkin(4) + pkin(10);
t151 = pkin(7) + pkin(8);
t82 = sin(qJ(2));
t84 = cos(qJ(2));
t55 = t145 * t84 + t147 * t82;
t103 = t145 * t82 - t147 * t84;
t97 = t103 * qJD(3);
t89 = -t103 * qJD(2) - t97;
t90 = t153 * t55;
t99 = t81 * t103;
t21 = -qJD(4) * t99 + t55 * t125 + t146 * t90 + t81 * t89;
t137 = qJD(4) * t81;
t62 = t151 * t82;
t63 = t151 * t84;
t105 = t145 * t63 + t147 * t62;
t35 = -t55 * pkin(9) - t105;
t104 = t145 * t62 - t147 * t63;
t36 = -t103 * pkin(9) - t104;
t126 = qJD(2) * t151;
t56 = t82 * t126;
t57 = t84 * t126;
t107 = t145 * t57 + t147 * t56;
t87 = -t90 * pkin(9) - t63 * t123 - t62 * t124 - t107;
t106 = t145 * t56 - t147 * t57;
t88 = -t89 * pkin(9) + t62 * t123 - t63 * t124 + t106;
t9 = -t35 * t125 + t36 * t137 - t146 * t87 - t81 * t88;
t6 = -t21 * pkin(5) - t9;
t4 = t6 * t80;
t5 = t6 * t83;
t149 = t82 * pkin(2);
t134 = qJD(6) * t83;
t27 = -t146 * t36 - t81 * t35;
t95 = t146 * t103;
t45 = t81 * t55 + t95;
t19 = -t45 * pkin(5) - t27;
t148 = t19 * t134 + t4;
t144 = t80 * t21;
t143 = t83 * t21;
t38 = -qJD(5) + t39;
t118 = t145 * t146;
t50 = pkin(2) * t118 + t81 * t72 + qJ(5);
t142 = t50 * t134 - t38 * t80;
t119 = pkin(3) * t125;
t66 = t119 + qJD(5);
t69 = t81 * pkin(3) + qJ(5);
t141 = t69 * t134 + t66 * t80;
t130 = qJ(5) * qJD(6);
t139 = qJD(5) * t80 + t83 * t130;
t136 = qJD(6) * t19;
t135 = qJD(6) * t80;
t133 = qJD(6) * t152;
t132 = t82 * qJD(2);
t131 = t84 * qJD(2);
t129 = -0.2e1 * pkin(1) * qJD(2);
t128 = t80 * t143;
t74 = pkin(3) * t137;
t127 = t80 * t134;
t73 = -t84 * pkin(2) - pkin(1);
t120 = pkin(2) * t123;
t71 = -t146 * pkin(3) - pkin(4);
t26 = -t146 * t35 + t81 * t36;
t46 = t146 * t55 - t99;
t48 = t103 * pkin(3) + t73;
t92 = -t46 * qJ(5) + t48;
t17 = t152 * t45 + t92;
t18 = t46 * pkin(5) + t26;
t117 = t83 * t17 + t80 * t18;
t116 = t80 * t17 - t83 * t18;
t115 = -qJ(5) * t21 - qJD(5) * t45;
t51 = t145 * t150 - t146 * t72 - pkin(4);
t20 = qJD(4) * t95 + t55 * t137 - t146 * t89 + t81 * t90;
t114 = t46 * t134 - t80 * t20;
t13 = -t46 * t135 - t83 * t20;
t113 = t45 * t134 + t144;
t112 = t45 * t135 - t143;
t49 = -pkin(10) + t51;
t111 = qJD(6) * (t45 * t50 - t46 * t49);
t68 = -pkin(10) + t71;
t110 = qJD(6) * (t45 * t69 - t46 * t68);
t109 = qJD(6) * (qJ(5) * t45 + t152 * t46);
t40 = t72 * t137 + (qJD(4) * t118 + (t147 * t81 + t118) * qJD(3)) * pkin(2);
t108 = -t50 * t21 + t38 * t45 + t40 * t46;
t102 = t152 * t20 + t115;
t100 = t73 * t55;
t98 = -t69 * t21 - t66 * t45 + t46 * t74;
t96 = -t20 * t49 + t108;
t94 = -t20 * t68 + t98;
t93 = -t39 + t86;
t37 = pkin(2) * t132 + t90 * pkin(3);
t8 = t21 * pkin(4) + t20 * qJ(5) - t46 * qJD(5) + t37;
t10 = t36 * t125 + t35 * t137 - t146 * t88 + t81 * t87;
t77 = qJD(5) * t83;
t65 = -0.2e1 * t127;
t59 = t66 * t83;
t54 = 0.2e1 * t122;
t44 = t45 ^ 2;
t34 = t38 * t83;
t32 = t74 + t40;
t29 = t104 * qJD(3) + t106;
t28 = t105 * qJD(3) + t107;
t25 = t45 * pkin(4) + t92;
t15 = -0.2e1 * t46 * t20;
t12 = -t45 * t122 + t128;
t11 = -0.4e1 * t45 * t127 - t138 * t21;
t7 = -t20 * pkin(5) + t10;
t3 = t21 * pkin(10) + t8;
t2 = -t117 * qJD(6) - t80 * t3 + t83 * t7;
t1 = t116 * qJD(6) - t83 * t3 - t80 * t7;
t14 = [0, 0, 0, 0.2e1 * t82 * t131, 0.2e1 * (-t82 ^ 2 + t84 ^ 2) * qJD(2), 0, 0, 0, t82 * t129, t84 * t129, 0.2e1 * t55 * t89, 0.2e1 * t153 * (t103 ^ 2 - t55 ^ 2) 0, 0, 0, 0.2e1 * qJD(3) * t100 + 0.2e1 * (t103 * t149 + t100) * qJD(2), -0.2e1 * t73 * t97 + 0.2e1 * (-t73 * t103 + t55 * t149) * qJD(2), t15, 0.2e1 * t20 * t45 - 0.2e1 * t46 * t21, 0, 0, 0, 0.2e1 * t48 * t21 + 0.2e1 * t37 * t45, -0.2e1 * t48 * t20 + 0.2e1 * t37 * t46, 0.2e1 * t10 * t46 - 0.2e1 * t26 * t20 + 0.2e1 * t27 * t21 + 0.2e1 * t9 * t45, -0.2e1 * t25 * t21 - 0.2e1 * t8 * t45, 0.2e1 * t25 * t20 - 0.2e1 * t8 * t46, 0.2e1 * t26 * t10 + 0.2e1 * t25 * t8 + 0.2e1 * t27 * t9, 0.2e1 * t78 * t45 * t21 + 0.2e1 * t44 * t127, -0.2e1 * t44 * t122 + 0.4e1 * t45 * t128, 0.2e1 * t114 * t45 + 0.2e1 * t46 * t144, 0.2e1 * t13 * t45 + 0.2e1 * t46 * t143, t15, 0.2e1 * t112 * t19 + 0.2e1 * t116 * t20 + 0.2e1 * t2 * t46 - 0.2e1 * t45 * t5, 0.2e1 * t1 * t46 + 0.2e1 * t113 * t19 + 0.2e1 * t117 * t20 + 0.2e1 * t45 * t4; 0, 0, 0, 0, 0, t131, -t132, 0, -pkin(7) * t131, pkin(7) * t132, 0, 0, t89, -t90, 0, t29, t28, 0, 0, -t20, -t21, 0, -t10, t9, -t51 * t20 + t108, t10, -t9, t10 * t51 + t26 * t40 + t27 * t38 - t9 * t50, t12, t11, t13, -t114, 0, t80 * t111 + t83 * t96 + t148, t5 + t83 * t111 + (-t96 - t136) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t120, -0.2e1 * t121, 0, 0, 0, 0, 0, -0.2e1 * t40, 0.2e1 * t39, 0, 0.2e1 * t40, -0.2e1 * t38, -0.2e1 * t50 * t38 + 0.2e1 * t51 * t40, t65, t54, 0, 0, 0, 0.2e1 * t142, -0.2e1 * t50 * t135 - 0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t90, 0, t29, t28, 0, 0, -t20, -t21, 0, -t10, t9, -t71 * t20 + t98, t10, -t9, t10 * t71 + t26 * t74 - t27 * t66 - t9 * t69, t12, t11, t13, -t114, 0, t110 * t80 + t83 * t94 + t148, t5 + t83 * t110 + (-t94 - t136) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t121, 0, 0, 0, 0, 0, -t32, -t119 + t39, 0, t32, t119 + t93, -t38 * t69 + t40 * t71 + t50 * t66 + t51 * t74, t65, t54, 0, 0, 0, t141 + t142, -t34 + t59 + (-t50 - t69) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t74, -0.2e1 * t119, 0, 0.2e1 * t74, 0.2e1 * t66, 0.2e1 * t69 * t66 + 0.2e1 * t71 * t74, t65, t54, 0, 0, 0, 0.2e1 * t141, -0.2e1 * t69 * t135 + 0.2e1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, 0, -t10, t9, pkin(4) * t20 + t115, t10, -t9, -t10 * pkin(4) - t9 * qJ(5) - t27 * qJD(5), t12, t11, t13, -t114, 0, t102 * t83 + t109 * t80 + t148, t5 + t83 * t109 + (-t102 - t136) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t39, 0, t40, t93, -t40 * pkin(4) - t38 * qJ(5) + t50 * qJD(5), t65, t54, 0, 0, 0, t139 + t142, -t34 + t77 + (-qJ(5) - t50) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t119, 0, t74, t86 + t119, -pkin(4) * t74 + t66 * qJ(5) + t69 * qJD(5), t65, t54, 0, 0, 0, t139 + t141, t59 + t77 + (-qJ(5) - t69) * t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, qJ(5) * t86, t65, t54, 0, 0, 0, 0.2e1 * t139, -0.2e1 * t80 * t130 + 0.2e1 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, t10, 0, 0, 0, 0, 0, t13, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, -t112, -t20, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t134, 0, -t49 * t135 + t83 * t40, -t49 * t134 - t80 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t134, 0, -t68 * t135 + t83 * t74, -t68 * t134 - t80 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t134, 0, t80 * t133, t83 * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, -t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t14;
