% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP4_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:12:18
% EndTime: 2019-03-08 20:12:27
% DurationCPUTime: 3.25s
% Computational Cost: add. (5664->385), mult. (14656->491), div. (0->0), fcn. (11623->10), ass. (0->193)
t140 = cos(pkin(11));
t250 = cos(qJ(4));
t189 = t250 * t140;
t138 = sin(pkin(11));
t143 = sin(qJ(4));
t208 = t143 * t138;
t160 = t189 - t208;
t139 = sin(pkin(6));
t146 = cos(qJ(2));
t211 = t139 * t146;
t151 = t160 * t211;
t243 = pkin(8) + qJ(3);
t123 = t243 * t138;
t124 = t243 * t140;
t262 = -t250 * t123 - t143 * t124;
t238 = -qJD(1) * t151 + t160 * qJD(3) + t262 * qJD(4);
t184 = qJD(4) * t208;
t110 = -qJD(4) * t189 + t184;
t116 = t250 * t138 + t143 * t140;
t111 = t116 * qJD(4);
t144 = sin(qJ(2));
t205 = qJD(1) * t139;
t187 = t144 * t205;
t277 = pkin(4) * t111 + pkin(9) * t110 - t187;
t145 = cos(qJ(5));
t126 = qJD(2) * t184;
t129 = qJD(2) * t189;
t157 = -qJD(4) * t129 + t126;
t197 = t145 * qJD(4);
t142 = sin(qJ(5));
t199 = qJD(5) * t142;
t202 = qJD(2) * t116;
t61 = -qJD(5) * t197 + t145 * t157 + t199 * t202;
t90 = qJD(4) * t142 + t145 * t202;
t82 = t90 * t199;
t268 = -t145 * t61 - t82;
t106 = qJD(2) * t208 - t129;
t101 = qJD(5) + t106;
t198 = qJD(5) * t145;
t203 = qJD(1) * t146;
t186 = t139 * t203;
t114 = (qJD(3) + t186) * qJD(2);
t265 = t160 * t114;
t118 = qJD(2) * qJ(3) + t187;
t141 = cos(pkin(6));
t204 = qJD(1) * t141;
t128 = t140 * t204;
t235 = pkin(8) * qJD(2);
t80 = t128 + (-t118 - t235) * t138;
t93 = t140 * t118 + t138 * t204;
t81 = t140 * t235 + t93;
t269 = -t143 * t81 + t250 * t80;
t29 = qJD(4) * t269 + t265;
t48 = t143 * t80 + t250 * t81;
t43 = qJD(4) * pkin(9) + t48;
t176 = qJD(3) - t186;
t132 = -pkin(3) * t140 - pkin(2);
t201 = qJD(2) * t132;
t100 = t176 + t201;
t52 = pkin(4) * t106 - pkin(9) * t202 + t100;
t200 = qJD(2) * t144;
t185 = t139 * t200;
t125 = qJD(1) * t185;
t99 = qJD(2) * t111;
t60 = t99 * pkin(4) + pkin(9) * t157 + t125;
t183 = t142 * t29 - t145 * t60 + t43 * t198 + t52 * t199;
t19 = t142 * t52 + t145 * t43;
t276 = t101 * t19 - t183;
t275 = t202 * qJD(4);
t88 = t142 * t202 - t197;
t220 = t145 * t88;
t224 = t142 * t90;
t169 = t220 + t224;
t62 = t90 * qJD(5) - t142 * t157;
t221 = t145 * t62;
t225 = t142 * t61;
t274 = (qJD(5) * (t142 * t88 - t145 * t90) - t221 + t225) * t116 + t169 * t110;
t77 = -pkin(4) * t160 - pkin(9) * t116 + t132;
t85 = -t143 * t123 + t250 * t124;
t242 = t277 * t142 + t145 * t238 + t77 * t198 - t85 * t199;
t14 = qJ(6) * t101 + t19;
t255 = pkin(5) * t99;
t2 = t183 - t255;
t273 = -t101 * t14 + t2;
t18 = -t142 * t43 + t145 * t52;
t3 = t142 * t60 + t145 * t29 + t52 * t198 - t43 * t199;
t272 = -t101 * t18 + t3;
t152 = t116 * t211;
t237 = -qJD(1) * t152 + qJD(3) * t116 + qJD(4) * t85;
t270 = -t90 * t101 + t62;
t239 = -t142 * t62 - t88 * t198;
t267 = (-t220 + t224) * t106 + t239 - t268;
t266 = t114 * t116;
t174 = t138 * (-t118 * t138 + t128) - t140 * t93;
t263 = t174 * t146;
t42 = -qJD(4) * pkin(4) - t269;
t23 = t88 * pkin(5) - t90 * qJ(6) + t42;
t253 = pkin(9) * t99;
t261 = t101 * t23 - t253;
t214 = t110 * t142;
t95 = t142 * t99;
t260 = t101 * (t116 * t198 - t214) + t111 * t88 - t160 * t62 + t116 * t95;
t236 = t142 * t77 + t145 * t85;
t241 = -qJD(5) * t236 - t142 * t238 + t277 * t145;
t258 = t90 ^ 2;
t257 = t101 ^ 2;
t256 = t202 ^ 2;
t254 = pkin(9) * t90;
t252 = pkin(5) * t111 + t241;
t251 = -qJ(6) * t111 + qJD(6) * t160 - t242;
t30 = qJD(4) * t48 + t266;
t5 = t62 * pkin(5) + t61 * qJ(6) - t90 * qJD(6) + t30;
t249 = t142 * t5;
t248 = t145 * t5;
t247 = t23 * t90;
t212 = t139 * t144;
t104 = -t138 * t212 + t140 * t141;
t105 = t138 * t141 + t140 * t212;
t161 = t250 * t104 - t143 * t105;
t246 = t30 * t161;
t245 = t30 * t262;
t244 = t90 * t88;
t177 = pkin(5) * t142 - qJ(6) * t145;
t178 = pkin(5) * t145 + qJ(6) * t142;
t240 = -t177 * t110 + (t178 * qJD(5) - qJD(6) * t145) * t116 + t237;
t75 = pkin(4) * t202 + pkin(9) * t106;
t28 = t142 * t75 + t145 * t269;
t234 = qJ(6) * t99;
t233 = qJD(2) * pkin(2);
t182 = t101 * t88;
t229 = t202 * t88;
t228 = t202 * t90;
t227 = t160 * t99;
t226 = t142 * t30;
t223 = t145 * t30;
t97 = t145 * t99;
t218 = -qJD(6) * t142 + t101 * t177 - t48;
t217 = t101 * t202;
t216 = t101 * t145;
t215 = t202 * t106;
t213 = t110 * t145;
t147 = qJD(2) ^ 2;
t210 = t139 * t147;
t207 = qJD(6) - t18;
t206 = t138 ^ 2 + t140 ^ 2;
t196 = pkin(9) * t199;
t195 = pkin(9) * t198;
t194 = t88 ^ 2 - t258;
t192 = t142 * t211;
t191 = t144 * t210;
t190 = t146 * t210;
t188 = t139 ^ 2 * t203;
t181 = t206 * t114;
t179 = (qJD(5) * t88 - t61) * pkin(9);
t13 = -pkin(5) * t101 + t207;
t175 = t13 * t145 - t14 * t142;
t173 = t142 * t19 + t145 * t18;
t27 = -t142 * t269 + t145 * t75;
t39 = -t142 * t85 + t145 * t77;
t167 = t101 * t198 + t106 * t216 + t95;
t166 = t97 + (-t106 * t142 - t199) * t101;
t69 = t143 * t104 + t250 * t105;
t57 = t142 * t69 + t145 * t211;
t44 = qJD(2) * t151 + qJD(4) * t161;
t16 = -qJD(5) * t57 + t142 * t185 + t145 * t44;
t17 = -qJD(5) * t192 + t142 * t44 - t145 * t185 + t69 * t198;
t58 = t145 * t69 - t192;
t163 = -t16 * t88 + t17 * t90 - t57 * t61 - t58 * t62;
t162 = t101 * t42 - t253;
t45 = qJD(2) * t152 + qJD(4) * t69;
t158 = -t101 * t17 - t161 * t62 + t45 * t88 - t57 * t99;
t156 = t142 * t182 - t221;
t155 = -t166 - t229;
t150 = t101 * t16 - t161 * t61 - t45 * t90 + t58 * t99;
t148 = -t239 * t116 - t88 * t214;
t119 = -pkin(4) - t178;
t117 = t176 - t233;
t103 = t106 ^ 2;
t56 = pkin(5) * t90 + qJ(6) * t88;
t53 = pkin(9) * t221;
t49 = t177 * t116 - t262;
t46 = t101 * t111 - t227;
t37 = pkin(5) * t160 - t39;
t36 = -qJ(6) * t160 + t236;
t34 = -t61 + t182;
t24 = t167 - t228;
t22 = t90 * t216 - t225;
t21 = -pkin(5) * t202 - t27;
t20 = qJ(6) * t202 + t28;
t12 = t268 * t116 - t90 * t213;
t6 = t116 * t97 + t111 * t90 + t160 * t61 + (-t116 * t199 - t213) * t101;
t1 = qJD(6) * t101 + t234 + t3;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, -t190, 0, 0, 0, 0, 0, 0, 0, 0, -t140 * t191, t138 * t191, t206 * t190 (-t104 * t138 + t105 * t140) * t114 + (-t144 * t188 + (t117 * t144 - t263) * t139) * qJD(2), 0, 0, 0, 0, 0, 0, -qJD(4) * t45 + (t106 * t200 - t146 * t99) * t139, -t44 * qJD(4) + (t146 * t157 + t200 * t202) * t139, -t44 * t106 + t157 * t161 + t202 * t45 - t69 * t99, t29 * t69 - t246 + t44 * t48 - t45 * t269 + (t100 * t139 - t188) * t200, 0, 0, 0, 0, 0, 0, t158, -t150, t163, t16 * t19 - t17 * t18 + t183 * t57 + t3 * t58 + t42 * t45 - t246, 0, 0, 0, 0, 0, 0, t158, t163, t150, t1 * t58 + t13 * t17 + t14 * t16 - t161 * t5 + t2 * t57 + t23 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176 * qJD(2) * t206 + t181, -t174 * qJD(3) + qJ(3) * t181 + (t263 + (-t117 - t233) * t144) * t205, -t110 * t202 - t116 * t157, t110 * t106 - t111 * t202 - t116 * t99 - t157 * t160, -t110 * qJD(4), t106 * t111 - t227, -t111 * qJD(4), 0, t100 * t111 + t132 * t99 - t237 * qJD(4) + (-qJD(2) * t160 - t106) * t187, -t100 * t110 - t132 * t126 + (t132 * t129 - t238) * qJD(4), -t238 * t106 + t110 * t269 - t48 * t111 + t30 * t116 + t157 * t262 + t160 * t29 + t202 * t237 - t85 * t99, t29 * t85 - t245 + t238 * t48 - t237 * t269 + (-t100 + t201) * t187, t12, t274, t6, t148, -t260, t46, -t42 * t214 + t111 * t18 + t160 * t183 + t39 * t99 - t62 * t262 + t237 * t88 + (t198 * t42 + t226) * t116 + t241 * t101, -t42 * t213 - t111 * t19 + t160 * t3 - t236 * t99 + t61 * t262 + t237 * t90 + (-t199 * t42 + t223) * t116 - t242 * t101, t39 * t61 - t236 * t62 - t241 * t90 - t242 * t88 + t173 * t110 + (-t142 * t3 + t145 * t183 + (t142 * t18 - t145 * t19) * qJD(5)) * t116, t18 * t241 - t183 * t39 + t19 * t242 + t236 * t3 + t237 * t42 - t245, t12, t6, -t274, t46, t260, t148, -t23 * t214 - t111 * t13 + t160 * t2 - t37 * t99 + t49 * t62 + t240 * t88 + (t198 * t23 + t249) * t116 + t252 * t101, -t36 * t62 - t37 * t61 - t252 * t90 + t251 * t88 - t175 * t110 + (-t1 * t142 + t145 * t2 + (-t13 * t142 - t14 * t145) * qJD(5)) * t116, t23 * t213 - t1 * t160 + t111 * t14 + t36 * t99 + t49 * t61 - t240 * t90 + (t199 * t23 - t248) * t116 - t251 * t101, t1 * t36 - t252 * t13 - t251 * t14 + t2 * t37 + t240 * t23 + t49 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206 * t147, t174 * qJD(2) + t125, 0, 0, 0, 0, 0, 0, 0.2e1 * t275, -t126 + (t129 - t106) * qJD(4), -t103 - t256, t106 * t48 + t202 * t269 + t125, 0, 0, 0, 0, 0, 0, t166 - t229, -t145 * t257 - t228 - t95, t267, t272 * t142 + t145 * t276 - t202 * t42, 0, 0, 0, 0, 0, 0, -t142 * t257 - t229 + t97, t267, t167 + t228, -t202 * t23 - t273 * t145 + (t101 * t13 + t1) * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t215, -t103 + t256, -t126 + (t129 + t106) * qJD(4), -t215, 0, 0, -t100 * t202 - t266, t100 * t106 - t265, 0, 0, t22, -t106 * t169 + t239 + t268, t24, t156, -t155, -t217, -pkin(4) * t62 - t202 * t18 - t223 - t48 * t88 + (-t27 - t195) * t101 + t162 * t142, pkin(4) * t61 + t202 * t19 + t226 - t48 * t90 + (t28 + t196) * t101 + t162 * t145, t27 * t90 + t28 * t88 - t53 + (-t106 * t18 + t3 + (-t18 + t254) * qJD(5)) * t145 + (t179 - t276) * t142, -pkin(4) * t30 - t18 * t27 - t19 * t28 - t42 * t48 + (-qJD(5) * t173 + t142 * t183 + t145 * t3) * pkin(9), t22, t24, t82 + (t106 * t90 + t62) * t142 + (t61 + t182) * t145, -t217, t155, t156, t202 * t13 + t119 * t62 - t248 + t218 * t88 + (t21 - t195) * t101 + t261 * t142, t20 * t88 - t21 * t90 - t53 + (t106 * t13 + t1 + (t13 + t254) * qJD(5)) * t145 + (t179 + t273) * t142, -t202 * t14 + t119 * t61 - t249 - t218 * t90 + (-t20 - t196) * t101 - t261 * t145, t119 * t5 - t13 * t21 - t14 * t20 + t218 * t23 + (qJD(5) * t175 + t1 * t145 + t142 * t2) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, -t194, t34, -t244, -t270, t99, -t42 * t90 + t276, t42 * t88 - t272, 0, 0, t244, t34, t194, t99, t270, -t244, -t56 * t88 - t247 + 0.2e1 * t255 + t276, pkin(5) * t61 - qJ(6) * t62 + (t14 - t19) * t90 + (t13 - t207) * t88, 0.2e1 * t234 - t23 * t88 + t56 * t90 + (0.2e1 * qJD(6) - t18) * t101 + t3, -pkin(5) * t2 + qJ(6) * t1 - t13 * t19 + t14 * t207 - t23 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244 - t275, t34, -t257 - t258, t247 + t273;];
tauc_reg  = t4;
