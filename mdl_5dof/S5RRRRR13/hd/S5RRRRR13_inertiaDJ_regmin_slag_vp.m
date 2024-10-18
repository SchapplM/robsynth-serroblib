% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR13_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:16
% EndTime: 2024-09-27 17:32:17
% DurationCPUTime: 0.86s
% Computational Cost: add. (1577->176), mult. (4454->250), div. (0->0), fcn. (3596->10), ass. (0->139)
t91 = sin(pkin(5));
t93 = sin(qJ(5));
t94 = sin(qJ(4));
t97 = cos(qJ(5));
t98 = cos(qJ(4));
t67 = (t93 * t94 - t97 * t98) * t91;
t68 = (t93 * t98 + t94 * t97) * t91;
t175 = qJD(4) + qJD(5);
t92 = cos(pkin(5));
t174 = -0.2e1 * t92;
t173 = 0.2e1 * t92;
t172 = pkin(4) * t98;
t171 = pkin(10) * t91;
t89 = t92 * pkin(4);
t163 = t92 * t94;
t95 = sin(qJ(3));
t96 = sin(qJ(2));
t161 = t95 * t96;
t100 = cos(qJ(2));
t87 = pkin(1) * t100 + pkin(2);
t99 = cos(qJ(3));
t72 = -pkin(1) * t161 + t87 * t99 + pkin(3);
t142 = t72 * t163;
t160 = t96 * t99;
t88 = t91 * pkin(9);
t65 = pkin(1) * t160 + t87 * t95 + t88;
t109 = -t65 * t98 - t142;
t164 = t91 * t98;
t85 = pkin(10) * t164;
t31 = -t109 + t85;
t170 = t31 * t97;
t152 = qJD(3) * t99;
t103 = (-t96 * t152 + (-t100 * t95 - t160) * qJD(2)) * pkin(1);
t153 = qJD(3) * t95;
t51 = -t153 * t87 + t103;
t169 = t51 * t94;
t86 = pkin(2) * t99 + pkin(3);
t141 = t86 * t163;
t80 = pkin(2) * t95 + t88;
t108 = -t80 * t98 - t141;
t52 = -t108 + t85;
t168 = t52 * t97;
t145 = pkin(3) * t163;
t107 = -pkin(9) * t164 - t145;
t63 = -t107 + t85;
t167 = t63 * t97;
t90 = t91 ^ 2;
t166 = t90 * t94;
t165 = t90 * t98;
t162 = t92 * t98;
t104 = (t96 * t153 + (-t100 * t99 + t161) * qJD(2)) * pkin(1);
t50 = -t152 * t87 + t104;
t121 = t51 * t162 + t94 * t50;
t13 = qJD(4) * t109 + t121;
t159 = t13 * t92 + t51 * t165;
t151 = qJD(4) * t94;
t126 = t91 * t151;
t81 = pkin(4) * t126;
t35 = -t51 * t91 + t81;
t43 = t175 * t67;
t57 = (-t72 - t172) * t91;
t158 = t35 * t68 - t43 * t57;
t135 = pkin(2) * t153;
t69 = t135 * t91 + t81;
t73 = (-t86 - t172) * t91;
t157 = -t43 * t73 + t69 * t68;
t75 = (-pkin(3) - t172) * t91;
t156 = -t75 * t43 + t68 * t81;
t150 = qJD(4) * t98;
t127 = t92 * t150;
t134 = pkin(2) * t152;
t155 = -t86 * t127 - t98 * t134;
t154 = pkin(1) * qJD(2);
t149 = qJD(5) * t93;
t148 = qJD(5) * t97;
t125 = -t65 - t171;
t112 = t125 * t94;
t138 = t72 * t127 + t51 * t163 - t98 * t50;
t10 = qJD(4) * t112 + t138;
t102 = (t125 * t98 - t142) * qJD(4) + t121;
t24 = t162 * t72 + t112 + t89;
t147 = -t10 * t97 - t102 * t93 - t148 * t24;
t146 = t97 * t89;
t123 = -t93 * t10 + t102 * t97;
t3 = (-t24 * t93 - t170) * qJD(5) + t123;
t44 = t175 * t68;
t144 = t3 * t92 + t35 * t67 + t44 * t57;
t106 = (-t162 * t95 - t94 * t99) * qJD(3) * pkin(2);
t124 = -t80 - t171;
t101 = t106 + (t124 * t98 - t141) * qJD(4);
t115 = t92 * t135;
t30 = (qJD(4) * t124 - t115) * t94 - t155;
t122 = t101 * t97 - t93 * t30;
t47 = t124 * t94 + t162 * t86 + t89;
t6 = (-t47 * t93 - t168) * qJD(5) + t122;
t143 = t44 * t73 + t6 * t92 + t69 * t67;
t131 = t91 * (-pkin(9) - pkin(10));
t105 = (t131 * t98 - t145) * qJD(4);
t113 = t94 * t131;
t82 = pkin(3) * t127;
t61 = qJD(4) * t113 + t82;
t120 = t97 * t105 - t93 * t61;
t60 = pkin(3) * t162 + t113 + t89;
t17 = (-t60 * t93 - t167) * qJD(5) + t120;
t140 = t17 * t92 + t75 * t44 + t67 * t81;
t139 = -t101 * t93 - t47 * t148 - t30 * t97;
t137 = -t93 * t105 - t60 * t148 - t97 * t61;
t136 = t96 * t154;
t133 = pkin(4) * t149;
t132 = pkin(4) * t148;
t130 = t100 * t154;
t129 = t90 * t151;
t128 = t90 * t150;
t84 = t91 * t150;
t119 = qJD(3) * (-pkin(2) - t87);
t118 = qJD(4) * (-pkin(3) - t72);
t117 = qJD(4) * (-pkin(3) - t86);
t116 = qJD(4) * (-t72 - t86);
t114 = t98 * t135;
t79 = 0.2e1 * t94 * t128;
t78 = t135 * t166;
t77 = t84 * t173;
t76 = t126 * t174;
t71 = t107 * qJD(4);
t70 = pkin(9) * t126 - t82;
t66 = 0.2e1 * (-t94 ^ 2 + t98 ^ 2) * t90 * qJD(4);
t64 = t71 * t92;
t39 = t44 * t174;
t38 = t43 * t173;
t37 = qJD(4) * t108 + t106;
t36 = (qJD(4) * t80 + t115) * t94 + t155;
t34 = t37 * t92;
t23 = -0.2e1 * t68 * t43;
t16 = t149 * t63 + t137;
t14 = 0.2e1 * t43 * t67 - 0.2e1 * t44 * t68;
t12 = t151 * t65 - t138;
t5 = t149 * t52 + t139;
t2 = t149 * t31 + t147;
t1 = [0, 0, 0, 0, -0.2e1 * t136, -0.2e1 * t130, 0, 0.2e1 * t51, 0.2e1 * t50, t79, t66, t77, t76, 0, -0.2e1 * t129 * t72 + 0.2e1 * t159, 0.2e1 * t12 * t92 + 0.2e1 * (-t150 * t72 - t169) * t90, t23, t14, -t38, t39, 0, 0.2e1 * t144, 0.2e1 * t2 * t92 + 0.2e1 * t158; 0, 0, 0, 0, -t136, -t130, 0, t119 * t95 + t103, t119 * t99 + t104, t79, t66, t77, t76, 0, t34 + (t116 * t94 - t114) * t90 + t159, t78 + (t12 + t36) * t92 + (t116 * t98 - t169) * t90, t23, t14, -t38, t39, 0, t143 + t144, (t2 + t5) * t92 + t157 + t158; 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t135, -0.2e1 * t134, t79, t66, t77, t76, 0, 0.2e1 * t34 + 0.2e1 * (-t151 * t86 - t114) * t90, -0.2e1 * t128 * t86 + 0.2e1 * t36 * t92 + 0.2e1 * t78, t23, t14, -t38, t39, 0, 0.2e1 * t143, 0.2e1 * t5 * t92 + 0.2e1 * t157; 0, 0, 0, 0, 0, 0, 0, t51, t50, t79, t66, t77, t76, 0, t118 * t166 + t159 + t64, (t12 + t70) * t92 + (t118 * t98 - t169) * t90, t23, t14, -t38, t39, 0, t140 + t144, (t16 + t2) * t92 + t156 + t158; 0, 0, 0, 0, 0, 0, 0, -t135, -t134, t79, t66, t77, t76, 0, t34 + t64 + (t117 * t94 - t114) * t90, t78 + (t36 + t70) * t92 + t117 * t165, t23, t14, -t38, t39, 0, t140 + t143, (t16 + t5) * t92 + t156 + t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t66, t77, t76, 0, -0.2e1 * pkin(3) * t129 + 0.2e1 * t64, -0.2e1 * pkin(3) * t128 + 0.2e1 * t70 * t92, t23, t14, -t38, t39, 0, 0.2e1 * t140, 0.2e1 * t16 * t92 + 0.2e1 * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t126, 0, t13, t12, 0, 0, -t43, -t44, 0, (-t170 + (-t24 - t89) * t93) * qJD(5) + t123, (t31 * t93 - t146) * qJD(5) + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t126, 0, t37, t36, 0, 0, -t43, -t44, 0, (-t168 + (-t47 - t89) * t93) * qJD(5) + t122, (t52 * t93 - t146) * qJD(5) + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t126, 0, t71, t70, 0, 0, -t43, -t44, 0, (-t167 + (-t60 - t89) * t93) * qJD(5) + t120, (t63 * t93 - t146) * qJD(5) + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t133, -0.2e1 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t44, 0, t17, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
