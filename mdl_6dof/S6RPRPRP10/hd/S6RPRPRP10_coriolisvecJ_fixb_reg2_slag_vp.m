% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP10_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:32
% EndTime: 2019-03-09 03:32:40
% DurationCPUTime: 2.80s
% Computational Cost: add. (3129->359), mult. (6141->440), div. (0->0), fcn. (3209->4), ass. (0->199)
t127 = sin(qJ(5));
t129 = cos(qJ(5));
t130 = cos(qJ(3));
t202 = qJD(1) * t130;
t112 = qJD(5) + t202;
t199 = qJD(3) * t129;
t128 = sin(qJ(3));
t203 = qJD(1) * t128;
t82 = t127 * t203 + t199;
t150 = t112 * t82;
t180 = t129 * t203;
t201 = qJD(3) * t127;
t80 = -t180 + t201;
t151 = t112 * t80;
t193 = qJD(1) * qJD(3);
t178 = t130 * t193;
t195 = qJD(5) * t127;
t44 = qJD(3) * t195 - qJD(5) * t180 - t127 * t178;
t216 = qJD(5) * t82;
t45 = -t129 * t178 + t216;
t271 = (t45 + t150) * t129 - (t44 + t151) * t127;
t149 = (t112 + t202) * t128;
t170 = qJD(5) * t130 + qJD(1);
t215 = t112 * t129;
t229 = t128 * t44;
t270 = (-t127 * t149 + t130 * t82) * qJD(3) + t170 * t215 - t229;
t194 = qJD(5) * t129;
t153 = (qJD(3) * pkin(8) - qJD(4)) * t130;
t122 = qJD(1) * qJD(2);
t177 = t128 * t193;
t186 = pkin(3) * t178 + qJ(4) * t177 + t122;
t32 = qJD(1) * t153 + t186;
t131 = -pkin(3) - pkin(8);
t197 = qJD(3) * t131;
t132 = -pkin(1) - pkin(7);
t105 = t132 * qJD(1) + qJD(2);
t91 = t130 * t105;
t253 = qJD(4) - t91;
t208 = pkin(4) * t202 + t253;
t41 = t208 + t197;
t219 = qJ(4) * t130;
t163 = pkin(8) * t128 - t219;
t206 = pkin(3) * t203 + qJD(1) * qJ(2);
t50 = qJD(1) * t163 + t206;
t200 = qJD(3) * t128;
t84 = t105 * t200;
t55 = -pkin(4) * t177 + t84;
t144 = -t127 * t55 - t129 * t32 - t41 * t194 + t195 * t50;
t16 = -t127 * t50 + t129 * t41;
t267 = -t112 * t16 - t144;
t17 = t127 * t41 + t129 * t50;
t176 = t127 * t32 - t129 * t55 + t50 * t194 + t41 * t195;
t266 = t112 * t17 - t176;
t225 = t129 * t44;
t12 = -t127 * t150 - t225;
t198 = qJD(3) * t130;
t230 = t127 * t80;
t38 = t45 * t127;
t264 = t128 * ((t127 * t82 + t129 * t80) * qJD(5) + t38 + t225) + (-t129 * t82 + t230) * t198;
t125 = t128 ^ 2;
t214 = t112 * t130;
t148 = qJD(1) * t125 - t214;
t181 = t112 * t195;
t263 = qJD(3) * (-t128 * t80 + t129 * t148) + t128 * t181 + t130 * t45;
t123 = qJD(3) * qJ(4);
t90 = t128 * t105;
t73 = -t90 - t123;
t223 = t130 * t73;
t121 = qJD(3) * qJD(4);
t68 = -t105 * t198 - t121;
t235 = qJD(3) * pkin(3);
t173 = -qJD(4) + t235;
t70 = -t173 - t91;
t262 = qJD(3) * (t128 * (-t70 + t91) + t223) + t128 * t68;
t14 = qJ(6) * t112 + t17;
t110 = pkin(5) * t177;
t2 = t110 + t176;
t261 = -t112 * t14 + t2;
t185 = t129 * t202;
t88 = t112 * t194;
t259 = t112 * t185 + t88;
t257 = t45 - t150;
t254 = t128 * pkin(3) + qJ(2);
t75 = t163 + t254;
t239 = pkin(4) - t132;
t95 = t239 * t130;
t236 = t127 * t95 + t129 * t75;
t179 = pkin(3) * t198 + qJ(4) * t200 + qJD(2);
t46 = t153 + t179;
t78 = t239 * t200;
t10 = -qJD(5) * t236 - t127 * t46 - t129 * t78;
t159 = t127 * t16 - t129 * t17;
t248 = -qJD(5) * t159 - t127 * t144 - t129 * t176;
t169 = qJ(6) * t177;
t1 = qJD(6) * t112 - t144 - t169;
t207 = qJD(6) - t16;
t13 = -pkin(5) * t112 + t207;
t161 = t127 * t13 + t129 * t14;
t247 = qJD(5) * t161 + t1 * t127 - t129 * t2;
t246 = t82 ^ 2;
t245 = t112 ^ 2;
t244 = 0.2e1 * t122;
t47 = -pkin(4) * t178 - t68;
t5 = pkin(5) * t45 + qJ(6) * t44 - qJD(6) * t82 + t47;
t243 = t127 * t5;
t242 = t129 * t5;
t64 = -pkin(4) * t203 + t90;
t52 = t123 + t64;
t19 = pkin(5) * t80 - qJ(6) * t82 + t52;
t241 = t19 * t82;
t240 = t82 * t80;
t165 = pkin(5) * t129 + qJ(6) * t127;
t147 = -pkin(4) - t165;
t142 = t147 * t130;
t238 = qJD(1) * t142 - qJD(5) * t165 + qJD(6) * t129 - t253;
t213 = t127 * t131;
t66 = t80 * t194;
t237 = -t131 * t66 - t45 * t213;
t86 = pkin(3) * t202 + qJ(4) * t203;
t63 = pkin(8) * t202 + t86;
t24 = t127 * t64 + t129 * t63;
t231 = t127 * t47;
t226 = t128 * t82;
t224 = t129 * t47;
t221 = t131 * t44;
t220 = t131 * t82;
t218 = qJD(3) * t19;
t217 = qJD(3) * t52;
t212 = t129 * t130;
t133 = qJD(3) ^ 2;
t211 = t133 * t128;
t210 = t133 * t130;
t134 = qJD(1) ^ 2;
t209 = t134 * qJ(2);
t126 = t130 ^ 2;
t205 = t125 - t126;
t204 = t133 + t134;
t196 = qJD(4) * t130;
t192 = t80 * t185 + t38 + t66;
t191 = 0.2e1 * qJD(1);
t190 = t80 ^ 2 - t246;
t67 = t80 * t198;
t189 = t127 * t214;
t188 = t112 * t213;
t187 = t131 * t215;
t184 = t128 * t199;
t183 = t128 * t197;
t182 = t129 * t197;
t71 = -qJ(4) * t202 + t206;
t93 = -t219 + t254;
t175 = qJD(1) * t93 + t71;
t101 = t127 * t177;
t172 = -qJD(3) * t82 + t101;
t168 = t130 * t177;
t164 = pkin(5) * t127 - qJ(6) * t129;
t162 = t127 * t14 - t129 * t13;
t160 = t127 * t17 + t129 * t16;
t23 = -t127 * t63 + t129 * t64;
t30 = -t127 * t75 + t129 * t95;
t43 = -qJD(1) * t196 + t186;
t61 = t179 - t196;
t145 = -qJD(1) * t61 + t132 * t133 - t43;
t9 = -t127 * t78 + t129 * t46 + t95 * t194 - t195 * t75;
t143 = t203 * t80 - t101 + t259;
t140 = -t12 - t192;
t138 = -t129 * t67 + (-t129 * t45 + t195 * t80) * t128;
t137 = t128 * t45 + t67 + t129 * t168 + (t170 * t127 + t184) * t112;
t136 = -t181 - qJD(3) * t80 + (-t184 - t189) * qJD(1);
t135 = -t200 * t230 + (-qJD(1) * t82 + (t45 - t216) * t130) * t127 + (-t82 * t200 + qJD(1) * t80 + (qJD(5) * t80 - t44) * t130) * t129;
t119 = qJ(2) * t244;
t115 = t128 * t132;
t109 = t132 * t198;
t108 = t130 * t134 * t128;
t100 = -0.2e1 * t168;
t99 = 0.2e1 * t168;
t98 = t204 * t130;
t97 = t204 * t128;
t96 = t205 * t134;
t94 = -pkin(4) * t128 + t115;
t92 = qJ(4) + t164;
t87 = t112 * t203;
t79 = -pkin(4) * t198 + t109;
t72 = 0.2e1 * t205 * t193;
t56 = t71 * t202;
t53 = qJD(3) * t149;
t42 = t128 * t147 + t115;
t34 = pkin(5) * t82 + qJ(6) * t80;
t26 = -pkin(5) * t130 - t30;
t25 = qJ(6) * t130 + t236;
t22 = -t44 + t151;
t21 = pkin(5) * t203 - t23;
t20 = -qJ(6) * t203 + t24;
t18 = -t181 + (-t189 + (t82 - t199) * t128) * qJD(1);
t15 = t109 + (qJD(5) * t164 - qJD(6) * t127) * t128 + qJD(3) * t142;
t11 = t194 * t226 + (t198 * t82 - t229) * t127;
t8 = t128 * t88 - t130 * t44 + (-t148 * t127 - t226) * qJD(3);
t7 = pkin(5) * t200 - t10;
t6 = -qJ(6) * t200 + qJD(6) * t130 + t9;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t119, t100, t72, -t211, t99, -t210, 0, -t132 * t211 + (qJ(2) * t198 + qJD(2) * t128) * t191, -t132 * t210 + (-qJ(2) * t200 + qJD(2) * t130) * t191, 0, t119, 0, t211, t210, t100, t72, t99, t262, t128 * t145 - t175 * t198, t130 * t145 + t175 * t200, -t132 * t262 + t43 * t93 + t61 * t71, t11, -t264, t8, t138, -t263, -t53, t10 * t112 + t45 * t94 + t79 * t80 + (-t199 * t52 - t176) * t130 + (t52 * t195 - t224 + (-qJD(1) * t30 - t16) * qJD(3)) * t128, -t112 * t9 - t44 * t94 + t79 * t82 + (t201 * t52 + t144) * t130 + (t52 * t194 + t231 + (qJD(1) * t236 + t17) * qJD(3)) * t128, -t10 * t82 + t30 * t44 - t236 * t45 - t80 * t9 - t159 * t198 + (-qJD(5) * t160 + t127 * t176 - t129 * t144) * t128, t10 * t16 - t144 * t236 + t17 * t9 - t176 * t30 + t47 * t94 + t52 * t79, t11, t8, t264, -t53, t263, t138, -t112 * t7 + t15 * t80 + t42 * t45 + (-t19 * t199 - t2) * t130 + (t19 * t195 - t242 + (qJD(1) * t26 + t13) * qJD(3)) * t128, -t25 * t45 - t26 * t44 - t6 * t80 + t7 * t82 + t161 * t198 + (-qJD(5) * t162 + t1 * t129 + t127 * t2) * t128, t112 * t6 - t15 * t82 + t42 * t44 + (-t19 * t201 + t1) * t130 + (-t19 * t194 - t243 + (-qJD(1) * t25 - t14) * qJD(3)) * t128, t1 * t25 + t13 * t7 + t14 * t6 + t15 * t19 + t2 * t26 + t42 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t209, 0, 0, 0, 0, 0, 0, -t97, -t98, 0, -t209, 0, 0, 0, 0, 0, 0, 0, t97, t98, -qJD(1) * t71 - t262, 0, 0, 0, 0, 0, 0, t137, t270, t135, t159 * qJD(1) + (qJD(3) * t160 + t47) * t128 + (t217 - t248) * t130, 0, 0, 0, 0, 0, 0, t137, t135, -t270, -t161 * qJD(1) + (qJD(3) * t162 + t5) * t128 + (t218 - t247) * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, -t96, 0, -t108, 0, 0, -t130 * t209, t128 * t209, 0, 0, 0, 0, 0, t108, -t96, -t108 ((-t73 - t123) * t130 + (t173 + t70) * t128) * qJD(1), t203 * t86 + t56, 0.2e1 * t121 + (-t128 * t71 + t130 * t86) * qJD(1), -qJ(4) * t68 - qJD(4) * t73 - t71 * t86 + (t223 + (-t70 - t235) * t128) * t105, t12, -t271, t18, t192, -t143, t87, qJ(4) * t45 - t112 * t23 + t231 + t208 * t80 + (t129 * t52 - t188) * qJD(5) + (t52 * t212 + (t16 - t182) * t128) * qJD(1), -qJ(4) * t44 + t112 * t24 + t224 + t208 * t82 + (-t127 * t52 - t187) * qJD(5) + (-t128 * t17 + (-t130 * t52 + t183) * t127) * qJD(1), t23 * t82 + t24 * t80 + (t221 - t266) * t129 + (t16 * t202 + t144 + (t16 + t220) * qJD(5)) * t127 + t237, qJ(4) * t47 + t131 * t248 - t16 * t23 - t17 * t24 + t208 * t52, t12, t18, t271, t87, t143, t129 * t151 + t38, t112 * t21 + t243 + t45 * t92 - t238 * t80 + (t129 * t19 - t188) * qJD(5) + (t19 * t212 + (-t13 - t182) * t128) * qJD(1), t20 * t80 - t21 * t82 + (t221 + t261) * t129 + (-t13 * t202 - t1 + (-t13 + t220) * qJD(5)) * t127 + t237, -t112 * t20 - t242 + t44 * t92 + t238 * t82 + (t127 * t19 + t187) * qJD(5) + (t128 * t14 + (t130 * t19 - t183) * t127) * qJD(1), -t13 * t21 + t131 * t247 - t14 * t20 - t238 * t19 + t5 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t126 * t134 - t133, qJD(3) * t73 + t56 + t84, 0, 0, 0, 0, 0, 0, t136, -t245 * t129 + t172, t140, t267 * t127 + t129 * t266 - t217, 0, 0, 0, 0, 0, 0, t136, t140, -t172 + t259, -t218 - t261 * t129 + (t112 * t13 + t1) * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, -t190, t22, -t240, -t257, -t177, -t52 * t82 + t266, t52 * t80 - t267, 0, 0, t240, t22, t190, -t177, t257, -t240, -t34 * t80 - 0.2e1 * t110 - t241 + t266, pkin(5) * t44 - qJ(6) * t45 + (t14 - t17) * t82 + (t13 - t207) * t80, -0.2e1 * t169 - t19 * t80 + t34 * t82 + (0.2e1 * qJD(6) - t16) * t112 - t144, -pkin(5) * t2 + qJ(6) * t1 - t13 * t17 + t14 * t207 - t19 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177 + t240, t22, -t245 - t246, t241 + t261;];
tauc_reg  = t3;
