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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:05:57
% EndTime: 2020-01-03 12:06:03
% DurationCPUTime: 2.43s
% Computational Cost: add. (2471->297), mult. (3785->418), div. (0->0), fcn. (2535->14), ass. (0->198)
t127 = qJDD(1) + qJDD(2);
t130 = qJD(1) + qJD(2);
t135 = sin(pkin(9));
t141 = cos(qJ(5));
t142 = cos(qJ(4));
t202 = qJD(4) + qJD(5);
t175 = t202 * t142;
t137 = sin(qJ(5));
t218 = t137 * t142;
t138 = sin(qJ(4));
t232 = t127 * t138;
t206 = qJD(4) * t138;
t185 = t135 * t206;
t219 = t137 * t138;
t190 = t135 * t219;
t256 = -qJD(5) * t190 - t137 * t185;
t15 = t256 * t130 + (t127 * t218 + (t130 * t175 + t232) * t141) * t135;
t134 = qJ(1) + qJ(2);
t123 = sin(t134);
t125 = cos(t134);
t167 = g(2) * t125 + g(3) * t123;
t139 = sin(qJ(2));
t243 = pkin(1) * qJD(1);
t196 = t139 * t243;
t143 = cos(qJ(2));
t247 = t143 * pkin(1);
t213 = -qJD(2) * t196 + qJDD(1) * t247;
t183 = qJDD(3) - t213;
t248 = t127 * pkin(2);
t62 = t183 - t248;
t259 = t62 + t167;
t128 = t135 ^ 2;
t136 = cos(pkin(9));
t223 = t136 * t130;
t99 = -qJD(4) + t223;
t216 = qJD(4) + t99;
t89 = t130 * qJ(3) + t196;
t257 = (t128 * t130 + t216 * t136) * t89;
t210 = qJD(1) * t143;
t168 = -pkin(1) * t210 + qJD(3);
t221 = t136 * t142;
t92 = -t136 * pkin(3) - t135 * pkin(7) - pkin(2);
t237 = qJ(3) * t221 + t138 * t92;
t205 = qJD(4) * t142;
t220 = t136 * t143;
t254 = qJD(3) * t221 + t92 * t205 - (t138 * t139 + t142 * t220) * t243;
t228 = t130 * t135;
t251 = g(1) * t135;
t118 = g(3) * t125;
t225 = t135 * t142;
t191 = t130 * t225;
t51 = t130 * t190 - t141 * t191;
t162 = t141 * t138 + t218;
t52 = t162 * t228;
t246 = t51 * t52;
t115 = t139 * pkin(1) + qJ(3);
t77 = t92 - t247;
t244 = t115 * t221 + t138 * t77;
t203 = qJDD(1) * t139;
t208 = qJD(2) * t143;
t54 = t127 * qJ(3) + t130 * qJD(3) + (qJD(1) * t208 + t203) * pkin(1);
t48 = t128 * t54;
t91 = -qJD(5) + t99;
t242 = t136 * t91;
t195 = t89 * t221;
t47 = t130 * t92 + t168;
t159 = -t138 * t47 - t195;
t226 = t135 * t138;
t200 = pkin(8) * t226;
t26 = -t130 * t200 - t159;
t241 = t141 * t26;
t224 = t136 * t127;
t96 = -qJDD(4) + t224;
t90 = -qJDD(5) + t96;
t240 = t90 * t136;
t239 = t96 * t136;
t236 = qJ(3) * t138;
t235 = t115 * t138;
t234 = t123 * t136;
t233 = t125 * t136;
t231 = t127 * t142;
t126 = t130 ^ 2;
t230 = t128 * t126;
t227 = t130 * t138;
t222 = t136 * t138;
t217 = t138 * t142;
t215 = t125 * pkin(2) + t123 * qJ(3);
t214 = g(2) * t123 - t118;
t129 = t136 ^ 2;
t212 = t128 + t129;
t132 = t142 ^ 2;
t211 = t138 ^ 2 - t132;
t209 = qJD(2) * t139;
t207 = qJD(3) * t138;
t204 = qJD(5) * t137;
t39 = t127 * t92 + t183;
t201 = t138 * t39 + t47 * t205 + t54 * t221;
t199 = pkin(8) * t225;
t109 = pkin(1) * t208 + qJD(3);
t197 = pkin(1) * t209;
t198 = t109 * t221 + t138 * t197 + t77 * t205;
t193 = t259 * t135;
t189 = t127 * t225;
t188 = qJ(3) * t205;
t187 = t130 * t209;
t186 = t130 * t205;
t184 = t136 * t206;
t42 = t142 * t47;
t25 = -pkin(8) * t191 - t89 * t222 + t42;
t17 = -t99 * pkin(4) + t25;
t152 = -t89 * t184 + t201;
t156 = t186 + t232;
t9 = -pkin(8) * t135 * t156 + t152;
t182 = qJD(5) * t17 + t9;
t181 = t212 * t54;
t180 = t96 + t224;
t179 = t109 * t212;
t178 = t130 * t216;
t177 = t212 * t127;
t176 = t129 * t54 - t214 + t48;
t174 = t130 * t196;
t36 = t142 * t39;
t171 = -t54 * t222 + t36;
t81 = t142 * t92;
t40 = -t199 + t81 + (-pkin(4) - t236) * t136;
t170 = qJD(5) * t40 + (-qJ(3) * t222 - t199) * qJD(4) + t254;
t104 = pkin(8) * t185;
t46 = -t200 + t237;
t63 = (-t138 * t220 + t139 * t142) * t243;
t169 = t237 * qJD(4) + qJD(5) * t46 + t136 * t207 - t104 + t63;
t166 = -t137 * t17 - t241;
t74 = t142 * t77;
t30 = -t199 + t74 + (-pkin(4) - t235) * t136;
t37 = -t200 + t244;
t165 = -t137 * t37 + t141 * t30;
t164 = t137 * t30 + t141 * t37;
t163 = qJD(4) * (t99 + t223);
t161 = -t141 * t142 + t219;
t105 = t135 * pkin(4) * t205;
t160 = -t168 * t135 - t105;
t66 = -t123 * t222 - t125 * t142;
t68 = -t123 * t142 + t125 * t222;
t158 = g(2) * t68 - g(3) * t66 + t152 * t136 + t142 * t48;
t67 = t123 * t221 - t125 * t138;
t69 = t123 * t138 + t125 * t221;
t157 = t128 * t89 * t205 - g(2) * t69 - g(3) * t67 + t138 * t48;
t155 = t167 - t213;
t22 = t26 * t204;
t29 = (pkin(4) * t156 + t54) * t135;
t146 = t202 * t162;
t31 = t146 * t135;
t50 = (pkin(4) * t227 + t89) * t135;
t133 = qJ(4) + qJ(5);
t122 = sin(t133);
t124 = cos(t133);
t57 = -t122 * t234 - t125 * t124;
t59 = t122 * t233 - t123 * t124;
t72 = t161 * t135;
t8 = -pkin(8) * t189 - t96 * pkin(4) + (-t195 + (pkin(8) * t228 - t47) * t138) * qJD(4) + t171;
t154 = g(2) * t59 - g(3) * t57 + (t137 * t8 + t182 * t141 - t22) * t136 - t29 * t72 - t50 * t31;
t120 = -pkin(2) - t247;
t153 = pkin(1) * t187 + t120 * t127;
t151 = -t99 ^ 2 - t230;
t2 = qJD(5) * t166 - t137 * t9 + t141 * t8;
t32 = t135 * t141 * t175 + t256;
t58 = -t125 * t122 + t124 * t234;
t60 = t123 * t122 + t124 * t233;
t71 = t162 * t135;
t150 = -g(2) * t60 - g(3) * t58 - t2 * t136 + t29 * t71 + t50 * t32;
t149 = -t167 + t174;
t148 = t168 * t212;
t147 = t124 * t251 + t22 + g(2) * t58 + (t26 * t91 - t8) * t137 - g(3) * t60 + t50 * t52;
t145 = -g(2) * t57 - g(3) * t59 + t122 * t251 + t50 * t51 + t2;
t14 = t141 * t189 + (-t127 * t219 - t130 * t146) * t135;
t144 = cos(qJ(1));
t140 = sin(qJ(1));
t116 = t123 * pkin(2);
t112 = pkin(4) * t226;
t107 = t142 * t197;
t85 = t135 * qJ(3) + t112;
t84 = -t130 * pkin(2) + t168;
t75 = t135 * t115 + t112;
t61 = t135 * t109 + t105;
t55 = (t127 * t132 - 0.2e1 * t138 * t186) * t128;
t38 = 0.2e1 * (t211 * t130 * qJD(4) - t127 * t217) * t128;
t28 = (t180 * t138 + t142 * t163) * t135;
t27 = (t138 * t163 - t180 * t142) * t135;
t24 = -t244 * qJD(4) - t109 * t222 + t104 + t107;
t23 = (-t115 * t222 - t199) * qJD(4) + t198;
t16 = t51 ^ 2 - t52 ^ 2;
t13 = t51 * t91 - t15;
t12 = -t52 * t91 + t14;
t11 = qJD(4) * t159 + t171;
t6 = -t14 * t72 + t51 * t31;
t5 = t15 * t136 + t32 * t91 + t71 * t90;
t4 = -t14 * t136 + t31 * t91 + t72 * t90;
t3 = -t14 * t71 + t72 * t15 + t31 * t52 + t51 * t32;
t1 = [qJDD(1), -g(2) * t144 - g(3) * t140, g(2) * t140 - g(3) * t144, t127, (t127 * t143 - t187) * pkin(1) - t155, ((-qJDD(1) - t127) * t139 + (-qJD(1) - t130) * t208) * pkin(1) + t214, (-t153 - t259) * t136, t135 * t153 + t193, t115 * t177 + t130 * t179 + t176, t62 * t120 + t84 * t197 - g(2) * (t144 * pkin(1) + t215) - g(3) * (t140 * pkin(1) - t125 * qJ(3) + t116) + t115 * t181 + t89 * t179, t55, t38, t27, t28, t239, -(-t77 * t206 + t107) * t99 - t74 * t96 + (-(-t109 * t138 - t115 * t205) * t99 + t96 * t235 - t11) * t136 + (t109 * t227 + t115 * t156) * t128 + t157, (-t115 * t184 + t198) * t99 + t244 * t96 + ((t109 * t130 + t115 * t127) * t142 + (-t115 * t130 - t89) * t206) * t128 + t158, t6, t3, t4, t5, t240, -(-qJD(5) * t164 - t137 * t23 + t141 * t24) * t91 - t165 * t90 + t61 * t52 + t75 * t15 + t150, (qJD(5) * t165 + t137 * t24 + t141 * t23) * t91 + t164 * t90 - t61 * t51 + t75 * t14 + t154; 0, 0, 0, t127, t149 + t213, (-t203 + (-qJD(2) + t130) * t210) * pkin(1) + t214, (t149 - t62 + t248) * t136, (-t174 - t248) * t135 + t193, qJ(3) * t177 + t148 * t130 + t176, -t62 * pkin(2) - t84 * t196 - g(2) * t215 - g(3) * t116 + (t181 + t118) * qJ(3) + t148 * t89, t55, t38, t27, t28, t239, -t81 * t96 + (t92 * t206 + t63) * t99 + (-(-t188 - t207) * t99 + t96 * t236 - t11) * t136 + (qJ(3) * t232 + (t138 * t168 + t188) * t130) * t128 + t157, t237 * t96 + (-qJ(3) * t184 + t254) * t99 + (qJ(3) * t231 - t89 * t206 + (-qJ(3) * t206 + t168 * t142) * t130) * t128 + t158, t6, t3, t4, t5, t240, -(-t137 * t46 + t141 * t40) * t90 + t85 * t15 + (t137 * t170 + t141 * t169) * t91 - t160 * t52 + t150, (t137 * t40 + t141 * t46) * t90 + t85 * t14 + (-t137 * t169 + t141 * t170) * t91 + t160 * t51 + t154; 0, 0, 0, 0, 0, 0, -t224, t135 * t127, -t212 * t126, -t212 * t89 * t130 + qJDD(3) + t155 - t248, 0, 0, 0, 0, 0, t138 * t151 - t142 * t96, t138 * t96 + t151 * t142, 0, 0, 0, 0, 0, t146 * t91 + t161 * t90 + (-t135 * t52 - t162 * t242) * t130, t51 * t228 + t162 * t90 + (t242 * t130 - t202 * t91) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217 * t230, -t211 * t230, (-t138 * t178 + t231) * t135, (-t142 * t178 - t232) * t135, -t96, -g(2) * t66 - g(3) * t68 + t36 - t142 * t257 + (-t136 * t54 - t216 * t47 + t251) * t138, g(1) * t225 + g(2) * t67 - g(3) * t69 + t138 * t257 - t42 * t99 - t201, -t246, t16, t12, t13, -t90, (-t137 * t25 - t241) * t91 + (-t141 * t90 - t191 * t52 + t204 * t91) * pkin(4) + t145, (-t25 * t91 - t182) * t141 + (qJD(5) * t141 * t91 + t137 * t90 + t191 * t51) * pkin(4) + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, t16, t12, t13, -t90, t166 * t91 + t145, (-t9 + (-qJD(5) - t91) * t17) * t141 + t147;];
tau_reg = t1;
