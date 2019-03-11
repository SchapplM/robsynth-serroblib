% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPPR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:27:53
% EndTime: 2019-03-09 16:28:03
% DurationCPUTime: 3.26s
% Computational Cost: add. (3982->351), mult. (10803->649), div. (0->0), fcn. (10085->10), ass. (0->160)
t207 = pkin(4) + pkin(9);
t146 = cos(qJ(3));
t143 = sin(qJ(3));
t193 = t143 * qJ(4);
t204 = pkin(3) + qJ(5);
t211 = t146 * t204 + t193;
t137 = sin(pkin(11));
t139 = cos(pkin(11));
t142 = sin(qJ(6));
t145 = cos(qJ(6));
t161 = t137 * t142 - t139 * t145;
t113 = (t137 ^ 2 + t139 ^ 2) * qJD(5);
t185 = qJD(4) * t146;
t210 = qJD(3) * t211 + qJD(5) * t143 - t185;
t138 = sin(pkin(6));
t209 = 0.2e1 * t138;
t208 = 0.2e1 * qJD(4);
t144 = sin(qJ(2));
t206 = pkin(1) * t144;
t205 = pkin(9) * t138;
t203 = -pkin(10) - t204;
t147 = cos(qJ(2));
t189 = qJD(2) * t144;
t131 = qJD(3) * t146;
t187 = qJD(3) * t143;
t140 = cos(pkin(6));
t197 = t138 * t147;
t180 = pkin(8) * t197;
t86 = t180 + (pkin(9) + t206) * t140;
t87 = (-pkin(2) * t147 - pkin(9) * t144 - pkin(1)) * t138;
t92 = (pkin(2) * t144 - pkin(9) * t147) * t138 * qJD(2);
t123 = t138 * t189;
t188 = qJD(2) * t147;
t93 = -pkin(1) * t140 * t188 + pkin(8) * t123;
t31 = -t131 * t86 + t143 * t93 + t146 * t92 - t187 * t87;
t175 = t138 * t188;
t198 = t138 * t144;
t178 = t143 * t198;
t73 = -qJD(3) * t178 + (qJD(3) * t140 + t175) * t146;
t22 = t73 * pkin(4) + (qJD(5) * t147 - t189 * t204) * t138 - t31;
t100 = -t140 * t146 + t178;
t101 = t140 * t143 + t146 * t198;
t94 = (t140 * t206 + t180) * qJD(2);
t148 = -t73 * qJ(4) - t101 * qJD(4) + t94;
t72 = qJD(3) * t101 + t143 * t175;
t23 = t100 * qJD(5) + t204 * t72 + t148;
t8 = t137 * t22 + t139 * t23;
t172 = -t143 * t86 + t146 * t87;
t47 = pkin(3) * t197 - t172;
t36 = pkin(4) * t101 + qJ(5) * t197 + t47;
t85 = pkin(8) * t198 + (-pkin(1) * t147 - pkin(2)) * t140;
t154 = -t101 * qJ(4) + t85;
t39 = t100 * t204 + t154;
t18 = t137 * t36 + t139 * t39;
t202 = t143 * t87 + t146 * t86;
t201 = t73 * t143;
t130 = pkin(9) * t131;
t110 = pkin(4) * t131 + t130;
t171 = pkin(3) * t187 - t143 * qJD(4);
t74 = -t146 * qJD(5) + (-qJ(4) * t146 + qJ(5) * t143) * qJD(3) + t171;
t53 = t110 * t137 + t139 * t74;
t104 = -pkin(2) - t211;
t117 = t207 * t143;
t66 = t104 * t139 + t117 * t137;
t200 = t137 * t143;
t199 = t137 * t146;
t196 = t139 * t143;
t195 = t139 * t146;
t118 = t207 * t146;
t190 = qJ(4) * qJD(4);
t186 = qJD(3) * t147;
t184 = qJD(6) * t142;
t183 = qJD(6) * t145;
t182 = qJD(6) * t146;
t181 = -0.2e1 * pkin(2) * qJD(3);
t179 = pkin(9) * t187;
t135 = t138 ^ 2;
t177 = t135 * t188;
t176 = t137 * t187;
t174 = qJD(4) * t197;
t173 = t139 * t187;
t7 = -t137 * t23 + t139 * t22;
t17 = -t137 * t39 + t139 * t36;
t52 = t110 * t139 - t137 * t74;
t170 = pkin(3) * t123;
t56 = t123 * t137 - t139 * t72;
t169 = pkin(10) * t56 - t8;
t168 = t144 * t177;
t167 = t8 * t137 + t7 * t139;
t166 = -pkin(3) * t146 - t193;
t105 = t137 * t145 + t139 * t142;
t99 = -t137 * t184 + t139 * t183;
t165 = -t101 * t99 - t105 * t73;
t71 = -t100 * t137 + t139 * t197;
t12 = pkin(5) * t101 + pkin(10) * t71 + t17;
t70 = t100 * t139 + t137 * t197;
t13 = pkin(10) * t70 + t18;
t4 = t12 * t142 + t13 * t145;
t26 = t137 * t53 + t139 * t52;
t108 = t139 * t117;
t58 = t143 * pkin(5) + t108 + (pkin(10) * t146 - t104) * t137;
t61 = -pkin(10) * t195 + t66;
t33 = t142 * t58 + t145 * t61;
t164 = t142 * t71 + t145 * t70;
t42 = t142 * t70 - t145 * t71;
t162 = qJD(5) * t101 + t204 * t73;
t111 = t203 * t137;
t112 = t203 * t139;
t69 = t111 * t145 + t112 * t142;
t46 = qJ(4) * t197 - t202;
t159 = -t105 * t131 - t143 * t99;
t30 = -t131 * t87 - t143 * t92 + t146 * t93 + t187 * t86;
t158 = pkin(10) * t173 + t53;
t157 = t143 * t186 + t146 * t189;
t156 = t143 * t189 - t146 * t186;
t40 = -pkin(4) * t100 - t46;
t57 = t123 * t139 + t137 * t72;
t155 = pkin(5) * t73 - pkin(10) * t57 + t7;
t153 = t157 * t138;
t152 = t156 * t138;
t151 = qJD(3) * t166 + t185;
t119 = qJ(4) * t123;
t25 = -t119 + t30 + t174;
t150 = (pkin(5) * t146 - pkin(10) * t200) * qJD(3) + t52;
t28 = -t31 - t170;
t149 = t28 * t143 - t25 * t146 + (t143 * t46 + t146 * t47) * qJD(3);
t24 = -pkin(4) * t72 - t25;
t127 = pkin(5) * t137 + qJ(4);
t122 = 0.2e1 * t143 * t131;
t114 = -pkin(2) + t166;
t109 = t207 * t187;
t98 = t105 * qJD(6);
t97 = pkin(5) * t195 + t118;
t95 = -qJ(4) * t131 + t171;
t90 = t105 * t146;
t89 = t161 * t146;
t84 = (-pkin(5) * t139 - t207) * t187;
t68 = -t111 * t142 + t112 * t145;
t65 = -t104 * t137 + t108;
t62 = -t131 * t161 - t143 * t98;
t60 = t105 * t182 - t142 * t176 + t145 * t173;
t59 = t105 * t187 + t161 * t182;
t55 = 0.2e1 * t101 * t73;
t50 = qJD(5) * t161 - qJD(6) * t69;
t49 = qJD(5) * t105 + t111 * t184 - t112 * t183;
t48 = t101 * t131 + t201;
t45 = t100 * pkin(3) + t154;
t37 = -t101 * t98 - t161 * t73;
t32 = -t142 * t61 + t145 * t58;
t29 = t72 * pkin(3) + t148;
t27 = -pkin(5) * t70 + t40;
t16 = qJD(6) * t42 + t142 * t57 + t145 * t56;
t15 = qJD(6) * t164 - t142 * t56 + t145 * t57;
t14 = pkin(5) * t56 + t24;
t10 = -qJD(6) * t33 - t142 * t158 + t145 * t150;
t9 = -t142 * t150 - t145 * t158 - t183 * t58 + t184 * t61;
t3 = t12 * t145 - t13 * t142;
t2 = -qJD(6) * t4 + t142 * t169 + t145 * t155;
t1 = -t12 * t183 + t13 * t184 - t142 * t155 + t145 * t169;
t5 = [0, 0, 0, 0.2e1 * t168, 0.2e1 * (-t144 ^ 2 + t147 ^ 2) * t135 * qJD(2), 0.2e1 * t140 * t175, -0.2e1 * t140 * t123, 0, -0.2e1 * pkin(1) * t135 * t189 - 0.2e1 * t140 * t94, -0.2e1 * pkin(1) * t177 + 0.2e1 * t140 * t93, t55, -0.2e1 * t100 * t73 - 0.2e1 * t101 * t72 (t101 * t189 - t147 * t73) * t209 (-t100 * t189 + t147 * t72) * t209, -0.2e1 * t168, 0.2e1 * t94 * t100 + 0.2e1 * t85 * t72 + 0.2e1 * (-t31 * t147 + t172 * t189) * t138, 0.2e1 * t94 * t101 + 0.2e1 * t85 * t73 + 0.2e1 * (-t30 * t147 - t189 * t202) * t138, 0.2e1 * t100 * t25 + 0.2e1 * t101 * t28 + 0.2e1 * t46 * t72 + 0.2e1 * t47 * t73, -0.2e1 * t29 * t100 - 0.2e1 * t45 * t72 + 0.2e1 * (-t147 * t28 + t189 * t47) * t138, -0.2e1 * t29 * t101 - 0.2e1 * t45 * t73 + 0.2e1 * (t147 * t25 - t189 * t46) * t138, 0.2e1 * t25 * t46 + 0.2e1 * t28 * t47 + 0.2e1 * t29 * t45, 0.2e1 * t101 * t7 + 0.2e1 * t17 * t73 - 0.2e1 * t24 * t70 + 0.2e1 * t40 * t56, -0.2e1 * t101 * t8 - 0.2e1 * t18 * t73 - 0.2e1 * t24 * t71 + 0.2e1 * t40 * t57, -0.2e1 * t17 * t57 - 0.2e1 * t18 * t56 + 0.2e1 * t7 * t71 + 0.2e1 * t70 * t8, 0.2e1 * t17 * t7 + 0.2e1 * t18 * t8 + 0.2e1 * t24 * t40, 0.2e1 * t42 * t15, 0.2e1 * t15 * t164 - 0.2e1 * t16 * t42, 0.2e1 * t101 * t15 + 0.2e1 * t42 * t73, -0.2e1 * t101 * t16 + 0.2e1 * t164 * t73, t55, 0.2e1 * t101 * t2 - 0.2e1 * t14 * t164 + 0.2e1 * t16 * t27 + 0.2e1 * t3 * t73, 0.2e1 * t1 * t101 + 0.2e1 * t14 * t42 + 0.2e1 * t15 * t27 - 0.2e1 * t4 * t73; 0, 0, 0, 0, 0, t175, -t123, 0, -t94, t93, t48, -t143 * t72 + t73 * t146 + (-t100 * t146 - t101 * t143) * qJD(3), t152, t153, 0, -pkin(2) * t72 - t94 * t146 - t156 * t205 + t187 * t85, -pkin(2) * t73 + t131 * t85 + t94 * t143 - t157 * t205 (t201 - t146 * t72 + (t100 * t143 + t101 * t146) * qJD(3)) * pkin(9) + t149, pkin(9) * t152 - t95 * t100 - t114 * t72 + t29 * t146 - t187 * t45, pkin(9) * t153 - t95 * t101 - t114 * t73 - t131 * t45 - t29 * t143, pkin(9) * t149 + t29 * t114 + t45 * t95, t24 * t195 + t52 * t101 + t109 * t70 + t118 * t56 + t7 * t143 + t65 * t73 + (t146 * t17 - t196 * t40) * qJD(3), -t24 * t199 - t53 * t101 + t109 * t71 + t118 * t57 - t8 * t143 - t66 * t73 + (-t146 * t18 + t200 * t40) * qJD(3), t52 * t71 + t53 * t70 - t66 * t56 - t65 * t57 + (t137 * t7 - t139 * t8) * t146 + (-t137 * t17 + t139 * t18) * t187, -t109 * t40 + t118 * t24 + t17 * t52 + t18 * t53 + t65 * t7 + t66 * t8, -t15 * t90 + t42 * t59, t15 * t89 + t16 * t90 + t164 * t59 + t42 * t60, t101 * t59 + t131 * t42 + t143 * t15 - t73 * t90, t101 * t60 + t131 * t164 - t143 * t16 + t73 * t89, t48, t10 * t101 + t131 * t3 - t14 * t89 + t143 * t2 + t16 * t97 - t164 * t84 - t27 * t60 + t32 * t73, t1 * t143 + t101 * t9 - t131 * t4 - t14 * t90 + t15 * t97 + t27 * t59 - t33 * t73 + t42 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0.2e1 * (-t143 ^ 2 + t146 ^ 2) * qJD(3), 0, 0, 0, t143 * t181, t146 * t181, 0, -0.2e1 * t114 * t187 + 0.2e1 * t146 * t95, -0.2e1 * t114 * t131 - 0.2e1 * t143 * t95, 0.2e1 * t114 * t95, -0.2e1 * t109 * t195 + 0.2e1 * t52 * t143 + 0.2e1 * (-t118 * t196 + t146 * t65) * qJD(3), 0.2e1 * t109 * t199 - 0.2e1 * t53 * t143 + 0.2e1 * (t118 * t200 - t146 * t66) * qJD(3), 0.2e1 * (t137 * t52 - t139 * t53) * t146 + 0.2e1 * (-t137 * t65 + t139 * t66) * t187, -0.2e1 * t109 * t118 + 0.2e1 * t52 * t65 + 0.2e1 * t53 * t66, -0.2e1 * t90 * t59, 0.2e1 * t59 * t89 - 0.2e1 * t60 * t90, -0.2e1 * t131 * t90 + 0.2e1 * t143 * t59, 0.2e1 * t131 * t89 + 0.2e1 * t143 * t60, t122, 0.2e1 * t10 * t143 + 0.2e1 * t131 * t32 - 0.2e1 * t60 * t97 - 0.2e1 * t84 * t89, -0.2e1 * t131 * t33 + 0.2e1 * t143 * t9 + 0.2e1 * t59 * t97 - 0.2e1 * t84 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t72, t123, t31, t30, -pkin(3) * t73 - qJ(4) * t72 - qJD(4) * t100, -t31 - 0.2e1 * t170, 0.2e1 * t119 - t30 - 0.2e1 * t174, -pkin(3) * t28 - qJ(4) * t25 - qJD(4) * t46, qJ(4) * t56 - qJD(4) * t70 + t24 * t137 - t139 * t162, qJ(4) * t57 - qJD(4) * t71 + t137 * t162 + t24 * t139 (-qJD(5) * t71 + t204 * t57 - t7) * t139 + (-qJD(5) * t70 + t204 * t56 - t8) * t137, t24 * qJ(4) + t40 * qJD(4) - t167 * t204 + (-t137 * t18 - t139 * t17) * qJD(5), -t15 * t161 - t42 * t98, -t105 * t15 + t16 * t161 - t164 * t98 - t42 * t99, t37, t165, 0, -qJD(4) * t164 + t101 * t50 + t105 * t14 + t127 * t16 + t27 * t99 + t68 * t73, qJD(4) * t42 + t101 * t49 + t127 * t15 - t14 * t161 - t27 * t98 - t69 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, -t187, 0, -t130, t179, t151, t130, -t179, t151 * pkin(9), -t109 * t137 - t139 * t210, -t109 * t139 + t137 * t210, -t26, -t109 * qJ(4) + t118 * qJD(4) - t26 * t204 + (-t137 * t66 - t139 * t65) * qJD(5), -t161 * t59 + t90 * t98, -t105 * t59 - t161 * t60 - t89 * t98 + t90 * t99, t62, t159, 0, -qJD(4) * t89 + t105 * t84 - t127 * t60 + t131 * t68 + t143 * t50 + t97 * t99, -qJD(4) * t90 + t127 * t59 - t131 * t69 + t143 * t49 - t161 * t84 - t97 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, 0.2e1 * t190, t137 * t208, t139 * t208, 0.2e1 * t113, 0.2e1 * t113 * t204 + 0.2e1 * t190, 0.2e1 * t161 * t98, 0.2e1 * t105 * t98 + 0.2e1 * t161 * t99, 0, 0, 0, 0.2e1 * qJD(4) * t105 + 0.2e1 * t127 * t99, -0.2e1 * qJD(4) * t161 - 0.2e1 * t127 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t123, 0, t28, t73 * t139, -t137 * t73, -t137 * t56 - t139 * t57, t167, 0, 0, 0, 0, 0, t37, t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, 0, t130, t139 * t131, -t137 * t131, 0, t26, 0, 0, 0, 0, 0, t62, t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t57, 0, t24, 0, 0, 0, 0, 0, t16, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t176, 0, -t109, 0, 0, 0, 0, 0, -t60, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), 0, 0, 0, 0, 0, t99, -t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t16, t73, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t60, t131, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t99, 0, t50, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
