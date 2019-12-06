% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:36:13
% EndTime: 2019-12-05 18:36:20
% DurationCPUTime: 2.40s
% Computational Cost: add. (2471->297), mult. (3785->419), div. (0->0), fcn. (2535->14), ass. (0->198)
t127 = qJDD(1) + qJDD(2);
t130 = qJD(1) + qJD(2);
t135 = sin(pkin(9));
t141 = cos(qJ(5));
t142 = cos(qJ(4));
t201 = qJD(4) + qJD(5);
t175 = t201 * t142;
t137 = sin(qJ(5));
t217 = t137 * t142;
t138 = sin(qJ(4));
t231 = t127 * t138;
t205 = qJD(4) * t138;
t185 = t135 * t205;
t218 = t137 * t138;
t190 = t135 * t218;
t254 = -qJD(5) * t190 - t137 * t185;
t15 = t254 * t130 + (t127 * t217 + (t130 * t175 + t231) * t141) * t135;
t128 = t135 ^ 2;
t136 = cos(pkin(9));
t222 = t136 * t130;
t99 = -qJD(4) + t222;
t215 = qJD(4) + t99;
t139 = sin(qJ(2));
t242 = pkin(1) * qJD(1);
t195 = t139 * t242;
t89 = t130 * qJ(3) + t195;
t255 = (t128 * t130 + t136 * t215) * t89;
t143 = cos(qJ(2));
t209 = qJD(1) * t143;
t167 = -pkin(1) * t209 + qJD(3);
t220 = t136 * t142;
t92 = -t136 * pkin(3) - t135 * pkin(7) - pkin(2);
t236 = qJ(3) * t220 + t138 * t92;
t134 = qJ(1) + qJ(2);
t123 = sin(t134);
t125 = cos(t134);
t166 = g(2) * t125 + g(3) * t123;
t246 = t143 * pkin(1);
t212 = -qJD(2) * t195 + qJDD(1) * t246;
t183 = qJDD(3) - t212;
t247 = t127 * pkin(2);
t62 = t183 - t247;
t252 = t166 - t62;
t204 = qJD(4) * t142;
t219 = t136 * t143;
t251 = qJD(3) * t220 + t92 * t204 - (t138 * t139 + t142 * t219) * t242;
t227 = t130 * t135;
t248 = g(1) * t135;
t118 = g(2) * t123;
t224 = t135 * t142;
t191 = t130 * t224;
t51 = t130 * t190 - t141 * t191;
t161 = t141 * t138 + t217;
t52 = t161 * t227;
t245 = t51 * t52;
t115 = t139 * pkin(1) + qJ(3);
t77 = t92 - t246;
t243 = t115 * t220 + t138 * t77;
t202 = qJDD(1) * t139;
t207 = qJD(2) * t143;
t54 = t127 * qJ(3) + t130 * qJD(3) + (qJD(1) * t207 + t202) * pkin(1);
t48 = t128 * t54;
t91 = -qJD(5) + t99;
t241 = t136 * t91;
t194 = t89 * t220;
t47 = t92 * t130 + t167;
t158 = -t138 * t47 - t194;
t225 = t135 * t138;
t199 = pkin(8) * t225;
t26 = -t130 * t199 - t158;
t240 = t141 * t26;
t223 = t136 * t127;
t96 = -qJDD(4) + t223;
t90 = -qJDD(5) + t96;
t239 = t90 * t136;
t238 = t96 * t136;
t235 = qJ(3) * t138;
t234 = t115 * t138;
t233 = t123 * t136;
t232 = t125 * t136;
t230 = t127 * t142;
t126 = t130 ^ 2;
t229 = t128 * t126;
t226 = t130 * t138;
t221 = t136 * t138;
t216 = t138 * t142;
t214 = g(2) * t232 + g(3) * t233;
t213 = g(3) * t125 - t118;
t129 = t136 ^ 2;
t211 = t128 + t129;
t132 = t142 ^ 2;
t210 = t138 ^ 2 - t132;
t208 = qJD(2) * t139;
t206 = qJD(3) * t138;
t203 = qJD(5) * t137;
t39 = t92 * t127 + t183;
t200 = t138 * t39 + t47 * t204 + t54 * t220;
t198 = pkin(8) * t224;
t109 = pkin(1) * t207 + qJD(3);
t196 = pkin(1) * t208;
t197 = t109 * t220 + t138 * t196 + t77 * t204;
t189 = t127 * t224;
t188 = qJ(3) * t204;
t187 = t130 * t208;
t186 = t130 * t204;
t184 = t136 * t205;
t42 = t142 * t47;
t25 = -pkin(8) * t191 - t89 * t221 + t42;
t17 = -t99 * pkin(4) + t25;
t151 = -t89 * t184 + t200;
t155 = t186 + t231;
t9 = -t155 * t135 * pkin(8) + t151;
t182 = qJD(5) * t17 + t9;
t181 = t211 * t54;
t180 = t96 + t223;
t179 = t109 * t211;
t178 = t130 * t215;
t177 = t211 * t127;
t176 = t129 * t54 - t213 + t48;
t174 = t130 * t195;
t173 = -t212 - t166;
t36 = t142 * t39;
t170 = -t54 * t221 + t36;
t81 = t142 * t92;
t40 = -t198 + t81 + (-pkin(4) - t235) * t136;
t169 = qJD(5) * t40 + (-qJ(3) * t221 - t198) * qJD(4) + t251;
t104 = pkin(8) * t185;
t46 = -t199 + t236;
t63 = (-t138 * t219 + t139 * t142) * t242;
t168 = t236 * qJD(4) + qJD(5) * t46 + t136 * t206 - t104 + t63;
t165 = -t137 * t17 - t240;
t74 = t142 * t77;
t30 = -t198 + t74 + (-pkin(4) - t234) * t136;
t37 = -t199 + t243;
t164 = -t137 * t37 + t141 * t30;
t163 = t137 * t30 + t141 * t37;
t162 = qJD(4) * (t99 + t222);
t160 = -t141 * t142 + t218;
t105 = t135 * pkin(4) * t204;
t159 = -t167 * t135 - t105;
t66 = t123 * t221 + t125 * t142;
t68 = -t123 * t142 + t125 * t221;
t157 = -g(2) * t68 - g(3) * t66 + t151 * t136 + t142 * t48;
t67 = t123 * t220 - t125 * t138;
t69 = -t123 * t138 - t125 * t220;
t156 = t128 * t89 * t204 - g(2) * t69 + g(3) * t67 + t138 * t48;
t154 = t174 + t247;
t22 = t26 * t203;
t29 = (t155 * pkin(4) + t54) * t135;
t146 = t201 * t161;
t31 = t146 * t135;
t50 = (pkin(4) * t226 + t89) * t135;
t133 = qJ(4) + qJ(5);
t122 = sin(t133);
t124 = cos(t133);
t57 = t122 * t233 + t125 * t124;
t59 = t122 * t232 - t123 * t124;
t72 = t160 * t135;
t8 = -pkin(8) * t189 - t96 * pkin(4) + (-t194 + (pkin(8) * t227 - t47) * t138) * qJD(4) + t170;
t153 = -g(2) * t59 - g(3) * t57 + (t137 * t8 + t182 * t141 - t22) * t136 - t29 * t72 - t50 * t31;
t120 = -pkin(2) - t246;
t152 = -pkin(1) * t187 - t120 * t127;
t150 = -t99 ^ 2 - t229;
t2 = t165 * qJD(5) - t137 * t9 + t141 * t8;
t32 = t141 * t135 * t175 + t254;
t58 = -t125 * t122 + t124 * t233;
t60 = -t123 * t122 - t124 * t232;
t71 = t161 * t135;
t149 = -g(2) * t60 + g(3) * t58 - t2 * t136 + t29 * t71 + t50 * t32;
t148 = t167 * t211;
t147 = t124 * t248 + t22 + (t26 * t91 - t8) * t137 - g(2) * t58 - g(3) * t60 + t50 * t52;
t145 = -g(2) * t57 + g(3) * t59 + t122 * t248 + t50 * t51 + t2;
t14 = t141 * t189 + (-t127 * t218 - t146 * t130) * t135;
t144 = cos(qJ(1));
t140 = sin(qJ(1));
t114 = t125 * qJ(3);
t112 = pkin(4) * t225;
t107 = t142 * t196;
t85 = t135 * qJ(3) + t112;
t84 = -t130 * pkin(2) + t167;
t75 = t135 * t115 + t112;
t61 = t135 * t109 + t105;
t56 = t62 * t135;
t55 = (t127 * t132 - 0.2e1 * t138 * t186) * t128;
t38 = 0.2e1 * (t210 * t130 * qJD(4) - t127 * t216) * t128;
t28 = (t180 * t138 + t142 * t162) * t135;
t27 = (t138 * t162 - t180 * t142) * t135;
t24 = -t243 * qJD(4) - t109 * t221 + t104 + t107;
t23 = (-t115 * t221 - t198) * qJD(4) + t197;
t16 = t51 ^ 2 - t52 ^ 2;
t13 = t51 * t91 - t15;
t12 = -t52 * t91 + t14;
t11 = t158 * qJD(4) + t170;
t6 = -t14 * t72 + t51 * t31;
t5 = t15 * t136 + t32 * t91 + t71 * t90;
t4 = -t14 * t136 + t31 * t91 + t72 * t90;
t3 = -t14 * t71 + t72 * t15 + t31 * t52 + t51 * t32;
t1 = [qJDD(1), g(2) * t144 + g(3) * t140, -g(2) * t140 + g(3) * t144, t127, (t127 * t143 - t187) * pkin(1) - t173, ((-qJDD(1) - t127) * t139 + (-qJD(1) - t130) * t207) * pkin(1) + t213, (t152 - t62) * t136 + t214, t56 + (-t152 - t166) * t135, t115 * t177 + t130 * t179 + t176, t62 * t120 + t84 * t196 - g(2) * (-t144 * pkin(1) - t125 * pkin(2) - t123 * qJ(3)) - g(3) * (-t140 * pkin(1) - t123 * pkin(2) + t114) + t115 * t181 + t89 * t179, t55, t38, t27, t28, t238, -(-t205 * t77 + t107) * t99 - t74 * t96 + (-(-t109 * t138 - t115 * t204) * t99 + t96 * t234 - t11) * t136 + (t109 * t226 + t115 * t155) * t128 + t156, (-t115 * t184 + t197) * t99 + t243 * t96 + ((t109 * t130 + t115 * t127) * t142 + (-t115 * t130 - t89) * t205) * t128 + t157, t6, t3, t4, t5, t239, -(-qJD(5) * t163 - t137 * t23 + t141 * t24) * t91 - t164 * t90 + t61 * t52 + t75 * t15 + t149, (qJD(5) * t164 + t137 * t24 + t141 * t23) * t91 + t163 * t90 - t61 * t51 + t75 * t14 + t153; 0, 0, 0, t127, -t173 + t174, (-t202 + (-qJD(2) + t130) * t209) * pkin(1) + t213, (t154 - t62) * t136 + t214, t56 + (-t154 - t166) * t135, qJ(3) * t177 + t148 * t130 + t176, -t84 * t195 - g(3) * t114 + t252 * pkin(2) + (t181 + t118) * qJ(3) + t148 * t89, t55, t38, t27, t28, t238, -t81 * t96 + (t205 * t92 + t63) * t99 + (-(-t188 - t206) * t99 + t96 * t235 - t11) * t136 + (qJ(3) * t231 + (t138 * t167 + t188) * t130) * t128 + t156, t236 * t96 + (-qJ(3) * t184 + t251) * t99 + (qJ(3) * t230 - t89 * t205 + (-qJ(3) * t205 + t142 * t167) * t130) * t128 + t157, t6, t3, t4, t5, t239, -(-t137 * t46 + t141 * t40) * t90 + t85 * t15 + (t137 * t169 + t141 * t168) * t91 - t159 * t52 + t149, (t137 * t40 + t141 * t46) * t90 + t85 * t14 + (-t137 * t168 + t141 * t169) * t91 + t159 * t51 + t153; 0, 0, 0, 0, 0, 0, -t223, t135 * t127, -t211 * t126, -t211 * t89 * t130 - t252, 0, 0, 0, 0, 0, t138 * t150 - t142 * t96, t138 * t96 + t142 * t150, 0, 0, 0, 0, 0, t146 * t91 + t160 * t90 + (-t135 * t52 - t161 * t241) * t130, t51 * t227 + t161 * t90 + (t241 * t130 - t201 * t91) * t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216 * t229, -t210 * t229, (-t138 * t178 + t230) * t135, (-t142 * t178 - t231) * t135, -t96, -g(2) * t66 + g(3) * t68 + t36 - t142 * t255 + (-t136 * t54 - t215 * t47 + t248) * t138, g(1) * t224 - g(2) * t67 - g(3) * t69 + t138 * t255 - t42 * t99 - t200, -t245, t16, t12, t13, -t90, (-t137 * t25 - t240) * t91 + (-t141 * t90 - t191 * t52 + t203 * t91) * pkin(4) + t145, (-t25 * t91 - t182) * t141 + (qJD(5) * t141 * t91 + t137 * t90 + t191 * t51) * pkin(4) + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t245, t16, t12, t13, -t90, t165 * t91 + t145, (-t9 + (-qJD(5) - t91) * t17) * t141 + t147;];
tau_reg = t1;
