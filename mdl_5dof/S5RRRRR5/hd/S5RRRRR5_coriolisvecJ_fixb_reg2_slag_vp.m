% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:30
% EndTime: 2019-12-05 18:58:35
% DurationCPUTime: 1.47s
% Computational Cost: add. (3956->243), mult. (6562->328), div. (0->0), fcn. (3934->8), ass. (0->174)
t120 = sin(qJ(4));
t117 = t120 ^ 2;
t124 = cos(qJ(4));
t118 = t124 ^ 2;
t180 = t117 + t118;
t125 = cos(qJ(3));
t126 = cos(qJ(2));
t195 = pkin(1) * qJD(1);
t173 = t126 * t195;
t150 = qJD(2) * t173;
t121 = sin(qJ(3));
t122 = sin(qJ(2));
t174 = t122 * t195;
t100 = t121 * t174;
t196 = (qJD(2) + qJD(3)) * t100;
t116 = qJD(1) + qJD(2);
t92 = pkin(2) * t116 + t173;
t46 = t125 * (qJD(3) * t92 + t150) - t196;
t220 = t180 * t46;
t123 = cos(qJ(5));
t176 = qJD(4) * t124;
t183 = t123 * t124;
t219 = -qJD(5) * t183 - t123 * t176;
t112 = qJD(3) + t116;
t119 = sin(qJ(5));
t86 = t119 * t124 + t120 * t123;
t75 = t86 * t112;
t178 = qJD(3) * t125;
t171 = pkin(2) * t178;
t80 = t125 * t173 - t100;
t218 = t171 - t80;
t179 = qJD(3) * t121;
t184 = t122 * t125;
t136 = t121 * t126 + t184;
t79 = t136 * t195;
t148 = pkin(2) * t179 - t79;
t71 = t121 * t92 + t125 * t174;
t63 = pkin(8) * t112 + t71;
t160 = pkin(9) * t112 + t63;
t141 = qJD(4) * t160;
t15 = -t120 * t141 + t124 * t46;
t48 = t160 * t120;
t45 = qJD(4) * pkin(4) - t48;
t217 = (qJD(5) * t45 + t15) * t123;
t115 = qJD(4) + qJD(5);
t106 = pkin(2) * t121 + pkin(8);
t127 = qJD(4) ^ 2;
t216 = t106 * t127 + t112 * t148;
t215 = qJD(2) * t136 + t122 * t178;
t214 = -pkin(9) - pkin(8);
t108 = pkin(1) * t126 + pkin(2);
t189 = pkin(1) * t184 + t121 * t108;
t78 = pkin(8) + t189;
t213 = -pkin(9) - t78;
t212 = pkin(2) * t125;
t211 = pkin(3) * t112;
t210 = pkin(4) * t124;
t109 = -pkin(3) - t210;
t70 = t125 * t92 - t100;
t51 = t109 * t112 - t70;
t209 = t51 * t75;
t186 = t119 * t120;
t168 = t112 * t186;
t73 = -t112 * t183 + t168;
t208 = t75 * t73;
t207 = -pkin(9) - t106;
t81 = t207 * t120;
t114 = t124 * pkin(9);
t82 = t106 * t124 + t114;
t54 = -t119 * t82 + t123 * t81;
t155 = qJD(4) * t207;
t68 = t120 * t155 + t124 * t171;
t69 = -t120 * t171 + t124 * t155;
t85 = -t183 + t186;
t206 = qJD(5) * t54 + t119 * t69 + t123 * t68 + t85 * t80;
t55 = t119 * t81 + t123 * t82;
t205 = -qJD(5) * t55 - t119 * t68 + t123 * t69 + t86 * t80;
t177 = qJD(4) * t120;
t164 = t112 * t177;
t130 = t215 * pkin(1);
t169 = t92 * t179;
t47 = qJD(1) * t130 + t169;
t36 = pkin(4) * t164 + t47;
t61 = t115 * t86;
t204 = t36 * t85 + t51 * t61;
t139 = t115 * t186;
t60 = t139 + t219;
t203 = t36 * t86 - t51 * t60;
t97 = t214 * t120;
t98 = pkin(8) * t124 + t114;
t67 = t119 * t97 + t123 * t98;
t166 = qJD(4) * t214;
t89 = t120 * t166;
t90 = t124 * t166;
t202 = qJD(5) * t67 + t119 * t89 - t123 * t90 - t86 * t70;
t66 = -t119 * t98 + t123 * t97;
t201 = -qJD(5) * t66 - t119 * t90 - t123 * t89 - t85 * t70;
t62 = -t70 - t211;
t199 = t47 * t120 + t62 * t176;
t110 = pkin(4) * t177;
t198 = t110 + t148;
t197 = t219 * t112;
t185 = t121 * t122;
t52 = t108 * t178 + (-t122 * t179 + (t125 * t126 - t185) * qJD(2)) * pkin(1);
t194 = t112 * t52;
t53 = t108 * t179 + t130;
t193 = t112 * t53;
t192 = t112 * t71;
t49 = t160 * t124;
t191 = t119 * t49;
t190 = t123 * t49;
t187 = t112 * t120;
t182 = t127 * t120;
t181 = t117 - t118;
t175 = pkin(4) * t187;
t13 = t123 * t45 - t191;
t14 = t119 * t45 + t190;
t16 = -t120 * t46 - t124 * t141;
t156 = -qJD(5) * t191 + t119 * t16;
t3 = t156 + t217;
t157 = -t119 * t15 + t123 * t16;
t4 = -qJD(5) * t14 + t157;
t170 = t13 * t60 - t14 * t61 - t3 * t85 - t4 * t86;
t111 = t112 ^ 2;
t167 = t120 * t111 * t124;
t161 = -pkin(4) * t115 - t45;
t159 = qJD(4) * t213;
t158 = -t112 * t62 - t46;
t154 = t180 * t70;
t151 = -pkin(1) * t185 + t108 * t125;
t149 = t124 * t164;
t147 = -t71 + t110;
t77 = -pkin(3) - t151;
t146 = (-qJD(2) + t116) * t195;
t145 = pkin(1) * qJD(2) * (-qJD(1) - t116);
t143 = pkin(8) * t127 - t192;
t142 = qJD(4) * (t70 - t211);
t140 = (-pkin(2) * t112 - t92) * qJD(3);
t138 = t127 * t78 + t193;
t64 = t213 * t120;
t65 = t124 * t78 + t114;
t34 = -t119 * t65 + t123 * t64;
t35 = t119 * t64 + t123 * t65;
t137 = qJD(4) * (t112 * t77 - t52);
t135 = t51 * t73 - t156;
t107 = -pkin(3) - t212;
t133 = qJD(4) * (t107 * t112 - t218);
t132 = t218 * t180;
t129 = t215 * t195;
t128 = -t129 - t169;
t113 = t127 * t124;
t94 = t109 - t212;
t84 = -0.2e1 * t149;
t83 = 0.2e1 * t149;
t76 = t77 - t210;
t72 = -0.2e1 * t181 * t112 * qJD(4);
t59 = t61 * t115;
t58 = t60 * t115;
t56 = t62 * t177;
t50 = t110 + t53;
t38 = t61 * t112;
t37 = t112 * t139 + t197;
t29 = -t120 * t52 + t124 * t159;
t28 = t120 * t159 + t124 * t52;
t27 = -t73 ^ 2 + t75 ^ 2;
t21 = -t197 + (-t168 + t73) * t115;
t20 = -t123 * t48 - t191;
t19 = t119 * t48 - t190;
t11 = t38 * t85 + t61 * t73;
t10 = -t37 * t86 - t60 * t75;
t7 = -qJD(5) * t35 - t119 * t28 + t123 * t29;
t6 = qJD(5) * t34 + t119 * t29 + t123 * t28;
t5 = t37 * t85 - t38 * t86 + t60 * t73 - t61 * t75;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 * t145, t126 * t145, 0, 0, 0, 0, 0, 0, 0, 0, t128 - t193, -t46 - t194, 0, -t47 * t151 + t46 * t189 + t71 * t52 - t70 * t53, t83, t72, t113, t84, -t182, 0, t56 + t120 * t137 + (-t138 - t47) * t124, t120 * t138 + t124 * t137 + t199, t180 * t194 + t220, t47 * t77 + t53 * t62 + t180 * (t46 * t78 + t52 * t63), t10, t5, -t58, t11, -t59, 0, t115 * t7 + t38 * t76 + t50 * t73 + t204, -t115 * t6 - t37 * t76 + t50 * t75 + t203, t34 * t37 - t35 * t38 - t6 * t73 - t7 * t75 + t170, t13 * t7 + t14 * t6 + t3 * t35 + t34 * t4 + t36 * t76 + t50 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 * t146, t126 * t146, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t79 + t121 * t140 - t129, t112 * t80 + (t140 - t150) * t125 + t196, 0, t70 * t79 - t71 * t80 + (t121 * t46 - t125 * t47 + (-t121 * t70 + t125 * t71) * qJD(3)) * pkin(2), t83, t72, t113, t84, -t182, 0, t56 + t120 * t133 + (-t216 - t47) * t124, t216 * t120 + t124 * t133 + t199, t112 * t132 + t220, t106 * t220 + t107 * t47 + t132 * t63 + t148 * t62, t10, t5, -t58, t11, -t59, 0, t205 * t115 + t198 * t73 + t38 * t94 + t204, -t206 * t115 + t198 * t75 - t37 * t94 + t203, -t205 * t75 - t206 * t73 + t37 * t54 - t38 * t55 + t170, t205 * t13 + t206 * t14 + t198 * t51 + t3 * t55 + t36 * t94 + t4 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128 + t192, t112 * t70 - t46, 0, 0, t83, t72, t113, t84, -t182, 0, t56 + t120 * t142 + (-t143 - t47) * t124, t120 * t143 + t124 * t142 + t199, -t112 * t154 + t220, -pkin(3) * t47 + pkin(8) * t220 - t63 * t154 - t62 * t71, t10, t5, -t58, t11, -t59, 0, t109 * t38 - t202 * t115 + t147 * t73 + t204, -t109 * t37 + t201 * t115 + t147 * t75 + t203, t201 * t73 + t202 * t75 + t37 * t66 - t38 * t67 + t170, t109 * t36 - t202 * t13 - t201 * t14 + t147 * t51 + t3 * t67 + t4 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t181 * t111, 0, t167, 0, 0, t158 * t120, t158 * t124, 0, 0, t208, t27, t21, -t208, 0, 0, -t73 * t175 - t115 * t19 - t209 + (t161 * t119 - t190) * qJD(5) + t157, -t75 * t175 + t115 * t20 + (t161 * qJD(5) - t15) * t123 + t135, (t14 + t19) * t75 + (-t13 + t20) * t73 + (-t119 * t38 + t123 * t37 + (t119 * t75 - t123 * t73) * qJD(5)) * pkin(4), -t13 * t19 - t14 * t20 + (-t51 * t187 + t119 * t3 + t123 * t4 + (-t119 * t13 + t123 * t14) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t27, t21, -t208, 0, 0, t115 * t14 - t209 + t4, t115 * t13 + t135 - t217, 0, 0;];
tauc_reg = t1;
