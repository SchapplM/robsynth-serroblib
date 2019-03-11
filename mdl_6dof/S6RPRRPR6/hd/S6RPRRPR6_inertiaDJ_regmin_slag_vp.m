% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRRPR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:43
% EndTime: 2019-03-09 05:17:48
% DurationCPUTime: 1.67s
% Computational Cost: add. (3760->214), mult. (8721->391), div. (0->0), fcn. (8876->10), ass. (0->129)
t124 = sin(qJ(4));
t126 = cos(qJ(4));
t155 = qJD(4) * t126;
t120 = sin(pkin(10));
t127 = cos(qJ(3));
t122 = cos(pkin(10));
t125 = sin(qJ(3));
t160 = t125 * t122;
t98 = t127 * t120 + t160;
t150 = t98 * t155;
t96 = t125 * t120 - t127 * t122;
t88 = t96 * qJD(3);
t178 = -t124 * t88 + t150;
t123 = sin(qJ(6));
t175 = cos(qJ(6));
t119 = sin(pkin(11));
t121 = cos(pkin(11));
t177 = -t119 * t124 + t121 * t126;
t97 = t119 * t126 + t121 * t124;
t60 = t123 * t97 - t175 * t177;
t118 = t126 ^ 2;
t159 = t124 ^ 2 - t118;
t145 = t159 * qJD(4);
t112 = -t122 * pkin(2) - pkin(1);
t176 = 0.2e1 * t112;
t174 = t119 * pkin(4);
t173 = t98 * t88;
t172 = pkin(7) + qJ(2);
t171 = -qJ(5) - pkin(8);
t138 = qJ(5) * t88 - qJD(5) * t98;
t102 = t172 * t120;
t103 = t172 * t122;
t157 = qJD(3) * t127;
t158 = qJD(2) * t127;
t42 = -t122 * t158 + t102 * t157 + (qJD(2) * t120 + qJD(3) * t103) * t125;
t89 = t98 * qJD(3);
t58 = t89 * pkin(3) + t88 * pkin(8);
t148 = t124 * t42 + t126 * t58;
t169 = qJ(5) * t98;
t59 = t96 * pkin(3) - t98 * pkin(8) + t112;
t161 = t125 * t102;
t65 = t127 * t103 - t161;
t62 = t126 * t65;
t14 = t89 * pkin(4) + t138 * t126 + (-t62 + (-t59 + t169) * t124) * qJD(4) + t148;
t153 = t124 * t58 - t126 * t42 + t59 * t155;
t16 = -qJ(5) * t150 + (-qJD(4) * t65 + t138) * t124 + t153;
t7 = t119 * t14 + t121 * t16;
t56 = t126 * t59;
t27 = t96 * pkin(4) - t124 * t65 - t126 * t169 + t56;
t166 = t124 * t98;
t170 = t124 * t59 + t62;
t30 = -qJ(5) * t166 + t170;
t20 = t119 * t27 + t121 * t30;
t147 = qJD(4) * t171;
t84 = t126 * qJD(5) + t124 * t147;
t85 = -t124 * qJD(5) + t126 * t147;
t47 = t119 * t85 + t121 * t84;
t104 = t171 * t124;
t105 = t171 * t126;
t67 = t119 * t104 - t121 * t105;
t165 = qJD(4) * t98;
t162 = t124 * t126;
t156 = qJD(4) * t124;
t154 = -0.2e1 * pkin(3) * qJD(4);
t114 = pkin(4) * t156;
t151 = t98 * t156;
t113 = -t126 * pkin(4) - pkin(3);
t149 = t124 * t155;
t6 = -t119 * t16 + t121 * t14;
t19 = -t119 * t30 + t121 * t27;
t46 = -t119 * t84 + t121 * t85;
t146 = -0.4e1 * t98 * t162;
t64 = t127 * t102 + t125 * t103;
t66 = t121 * t104 + t119 * t105;
t43 = qJD(2) * t160 - qJD(3) * t161 + t103 * t157 + t120 * t158;
t144 = 0.2e1 * (t120 ^ 2 + t122 ^ 2) * qJD(2);
t143 = pkin(3) * t88 - pkin(8) * t89;
t142 = pkin(3) * t98 + pkin(8) * t96;
t44 = pkin(4) * t166 + t64;
t130 = t97 * qJD(4);
t87 = -t119 * t156 + t121 * t155;
t37 = qJD(6) * t60 + t123 * t130 - t175 * t87;
t61 = t123 * t177 + t175 * t97;
t141 = t37 * t96 - t61 * t89;
t140 = t43 * t98 - t64 * t88;
t139 = -t88 * t96 + t98 * t89;
t36 = t178 * pkin(4) + t43;
t52 = t177 * t98;
t10 = t96 * pkin(5) - t52 * pkin(9) + t19;
t51 = t97 * t98;
t11 = -t51 * pkin(9) + t20;
t136 = -t175 * t10 + t123 * t11;
t135 = t123 * t10 + t175 * t11;
t48 = -t97 * pkin(9) + t66;
t49 = pkin(9) * t177 + t67;
t134 = t123 * t49 - t175 * t48;
t133 = t123 * t48 + t175 * t49;
t132 = -t123 * t52 - t175 * t51;
t34 = -t123 * t51 + t175 * t52;
t131 = t124 * t89 + t96 * t155;
t111 = t121 * pkin(4) + pkin(5);
t129 = -t175 * t111 + t123 * t174;
t128 = t123 * t111 + t175 * t174;
t32 = t97 * t165 + t177 * t88;
t4 = t89 * pkin(5) + t32 * pkin(9) + t6;
t31 = -t165 * t177 + t97 * t88;
t5 = t31 * pkin(9) + t7;
t2 = -t135 * qJD(6) - t123 * t5 + t175 * t4;
t1 = t136 * qJD(6) - t123 * t4 - t175 * t5;
t93 = t98 ^ 2;
t76 = t128 * qJD(6);
t75 = t129 * qJD(6);
t70 = -pkin(5) * t177 + t113;
t68 = pkin(5) * t130 + t114;
t63 = 0.2e1 * t96 * t89;
t57 = t126 * t89 - t96 * t156;
t40 = -pkin(9) * t130 + t47;
t39 = -t87 * pkin(9) + t46;
t38 = t61 * qJD(6) + t123 * t87 + t175 * t130;
t35 = t51 * pkin(5) + t44;
t24 = -t38 * t96 - t60 * t89;
t23 = -t170 * qJD(4) + t148;
t22 = t65 * t156 - t153;
t21 = -t31 * pkin(5) + t36;
t18 = -t133 * qJD(6) - t123 * t40 + t175 * t39;
t17 = t134 * qJD(6) - t123 * t39 - t175 * t40;
t9 = t34 * qJD(6) - t123 * t32 - t175 * t31;
t8 = t132 * qJD(6) + t123 * t31 - t175 * t32;
t3 = [0, 0, 0, 0, 0, t144, qJ(2) * t144, -0.2e1 * t173, -0.2e1 * t139, 0, 0, 0, t89 * t176, -t88 * t176, -0.2e1 * t118 * t173 - 0.2e1 * t93 * t149, 0.2e1 * t145 * t93 - t88 * t146, 0.2e1 * t139 * t126 - 0.2e1 * t96 * t151, -0.2e1 * t139 * t124 - 0.2e1 * t96 * t150, t63, 0.2e1 * t64 * t150 + 0.2e1 * t23 * t96 + 0.2e1 * t56 * t89 + 0.2e1 * (-t65 * t89 + t140) * t124, 0.2e1 * t140 * t126 - 0.2e1 * t64 * t151 - 0.2e1 * t170 * t89 + 0.2e1 * t22 * t96, 0.2e1 * t19 * t32 + 0.2e1 * t20 * t31 - 0.2e1 * t7 * t51 - 0.2e1 * t6 * t52, 0.2e1 * t19 * t6 + 0.2e1 * t20 * t7 + 0.2e1 * t44 * t36, 0.2e1 * t34 * t8, 0.2e1 * t132 * t8 - 0.2e1 * t34 * t9, 0.2e1 * t34 * t89 + 0.2e1 * t8 * t96, 0.2e1 * t132 * t89 - 0.2e1 * t9 * t96, t63, -0.2e1 * t132 * t21 - 0.2e1 * t136 * t89 + 0.2e1 * t2 * t96 + 0.2e1 * t35 * t9, 0.2e1 * t1 * t96 - 0.2e1 * t135 * t89 + 0.2e1 * t21 * t34 + 0.2e1 * t35 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, -t88, 0, 0, 0, 0, 0, t57, -t131, t52 * t130 + t177 * t32 + t97 * t31 - t87 * t51, -t19 * t130 + t177 * t6 + t20 * t87 + t7 * t97, 0, 0, 0, 0, 0, t24, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t130 * t177 + 0.2e1 * t97 * t87, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t89, 0, -t43, t42, -t98 * t145 - t88 * t162, qJD(4) * t146 + t159 * t88, t131, t57, 0, -t43 * t126 + t143 * t124 + (t124 * t64 - t142 * t126) * qJD(4), t43 * t124 + t143 * t126 + (t142 * t124 + t126 * t64) * qJD(4), -t20 * t130 + t177 * t7 - t19 * t87 + t67 * t31 + t66 * t32 - t46 * t52 - t47 * t51 - t6 * t97, t36 * t113 + t44 * t114 + t19 * t46 + t20 * t47 + t6 * t66 + t7 * t67, -t34 * t37 + t8 * t61, -t132 * t37 - t34 * t38 - t8 * t60 - t61 * t9, -t141, t24, 0, -t132 * t68 - t134 * t89 + t18 * t96 + t21 * t60 + t35 * t38 + t70 * t9, -t133 * t89 + t17 * t96 + t21 * t61 + t68 * t34 - t35 * t37 + t70 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66 * t130 + t177 * t46 + t97 * t47 + t87 * t67, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t149, -0.2e1 * t145, 0, 0, 0, t124 * t154, t126 * t154, -0.2e1 * t67 * t130 + 0.2e1 * t177 * t47 - 0.2e1 * t46 * t97 - 0.2e1 * t66 * t87, 0.2e1 * t113 * t114 + 0.2e1 * t66 * t46 + 0.2e1 * t67 * t47, -0.2e1 * t61 * t37, 0.2e1 * t37 * t60 - 0.2e1 * t61 * t38, 0, 0, 0, 0.2e1 * t70 * t38 + 0.2e1 * t68 * t60, -0.2e1 * t70 * t37 + 0.2e1 * t68 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t88 - t151, -t178, t89, t23, t22 (t119 * t31 + t121 * t32) * pkin(4) (t119 * t7 + t121 * t6) * pkin(4), 0, 0, t8, -t9, t89, -t129 * t89 - t76 * t96 + t2, -t128 * t89 + t75 * t96 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t155, 0 (t87 * t119 - t121 * t130) * pkin(4), 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t156, 0, -pkin(8) * t155, pkin(8) * t156 (-t119 * t130 - t121 * t87) * pkin(4) (t119 * t47 + t121 * t46) * pkin(4), 0, 0, -t37, -t38, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t76, 0.2e1 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, 0, 0, 0, 0, t38, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, t89, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
