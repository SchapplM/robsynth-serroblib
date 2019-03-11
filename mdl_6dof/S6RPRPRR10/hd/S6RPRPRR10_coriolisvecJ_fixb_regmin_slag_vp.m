% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:09:20
% EndTime: 2019-03-09 04:09:30
% DurationCPUTime: 3.78s
% Computational Cost: add. (4048->375), mult. (9259->548), div. (0->0), fcn. (6799->8), ass. (0->191)
t179 = sin(qJ(3));
t245 = qJD(1) * t179;
t168 = qJD(5) + t245;
t163 = qJD(6) + t168;
t177 = sin(qJ(6));
t180 = cos(qJ(6));
t175 = sin(pkin(10));
t182 = cos(qJ(3));
t244 = qJD(1) * t182;
t226 = t175 * t244;
t176 = cos(pkin(10));
t233 = t176 * qJD(3);
t139 = t226 - t233;
t225 = t176 * t244;
t243 = qJD(3) * t175;
t141 = t225 + t243;
t178 = sin(qJ(5));
t181 = cos(qJ(5));
t81 = t139 * t178 - t141 * t181;
t82 = t181 * t139 + t141 * t178;
t276 = t177 * t81 - t180 * t82;
t288 = t276 * t163;
t227 = t175 * t245;
t204 = pkin(3) * t182 + qJ(4) * t179;
t149 = t204 * qJD(1);
t183 = -pkin(1) - pkin(7);
t273 = qJD(1) * t183;
t164 = qJD(2) + t273;
t258 = t176 * t182;
t95 = t175 * t149 + t164 * t258;
t72 = pkin(8) * t227 + t95;
t287 = -qJD(4) * t176 + t72;
t274 = -t175 * t178 + t181 * t176;
t286 = qJD(5) * t274;
t285 = t168 * t82;
t199 = t177 * t82 + t180 * t81;
t284 = t199 * t276;
t283 = t81 * t168;
t282 = t163 * t199;
t275 = t274 * t179;
t252 = qJD(1) * t275 + t286;
t147 = t175 * t181 + t176 * t178;
t127 = t147 * qJD(1);
t272 = t147 * qJD(5);
t251 = t179 * t127 + t272;
t281 = t199 ^ 2 - t276 ^ 2;
t154 = pkin(3) * t179 - qJ(4) * t182 + qJ(2);
t132 = t154 * qJD(1);
t153 = t179 * t164;
t133 = qJD(3) * qJ(4) + t153;
t73 = t176 * t132 - t133 * t175;
t45 = pkin(4) * t245 - pkin(8) * t141 + t73;
t74 = t175 * t132 + t176 * t133;
t47 = -pkin(8) * t139 + t74;
t17 = t178 * t45 + t181 * t47;
t13 = -pkin(9) * t82 + t17;
t235 = qJD(6) * t177;
t11 = t13 * t235;
t265 = qJD(3) * pkin(3);
t214 = -qJD(4) + t265;
t261 = t164 * t182;
t124 = -t214 - t261;
t92 = pkin(4) * t139 + t124;
t32 = pkin(5) * t82 + t92;
t280 = -t276 * t32 + t11;
t234 = qJD(6) * t180;
t232 = qJD(1) * qJD(3);
t221 = t179 * t232;
t209 = t181 * t221;
t210 = t178 * t221;
t237 = qJD(5) * t181;
t238 = qJD(5) * t178;
t42 = -t139 * t237 - t141 * t238 + t175 * t210 - t176 * t209;
t43 = -t81 * qJD(5) - t175 * t209 - t176 * t210;
t8 = -t177 * t43 + t180 * t42 - t82 * t234 + t235 * t81;
t279 = t8 - t288;
t16 = -t178 * t47 + t181 * t45;
t12 = pkin(9) * t81 + t16;
t10 = pkin(5) * t168 + t12;
t264 = t13 * t180;
t2 = t10 * t177 + t264;
t171 = t182 * t232;
t230 = pkin(8) * t176 * t179;
t190 = (pkin(4) * t182 + t230) * qJD(1);
t121 = (qJD(4) + t261) * qJD(3);
t119 = t204 * qJD(3) - t182 * qJD(4) + qJD(2);
t98 = t119 * qJD(1);
t52 = -t175 * t121 + t176 * t98;
t33 = qJD(3) * t190 + t52;
t211 = t175 * t221;
t53 = t176 * t121 + t175 * t98;
t46 = pkin(8) * t211 + t53;
t217 = -t178 * t46 + t181 * t33;
t187 = -t17 * qJD(5) + t217;
t4 = pkin(5) * t171 - t42 * pkin(9) + t187;
t195 = t178 * t33 + t181 * t46 + t45 * t237 - t47 * t238;
t5 = -pkin(9) * t43 + t195;
t228 = -t177 * t5 + t180 * t4;
t278 = -t2 * qJD(6) + t32 * t199 + t228;
t186 = t199 * qJD(6) - t177 * t42 - t180 * t43;
t277 = t186 - t282;
t231 = 0.2e1 * qJD(1);
t270 = pkin(8) + qJ(4);
t158 = t270 * t175;
t159 = t270 * t176;
t249 = -t178 * t158 + t181 * t159;
t196 = qJD(4) * t175 + qJD(5) * t159;
t259 = t175 * t182;
t94 = t176 * t149 - t164 * t259;
t58 = t190 + t94;
t271 = t158 * t237 + t287 * t181 + (t196 + t58) * t178;
t87 = t147 * t177 - t180 * t274;
t269 = -t87 * qJD(6) - t251 * t177 + t252 * t180;
t88 = t147 * t180 + t177 * t274;
t268 = t88 * qJD(6) + t252 * t177 + t251 * t180;
t137 = t176 * t154;
t220 = -t175 * t183 + pkin(4);
t80 = -pkin(8) * t258 + t220 * t179 + t137;
t257 = t179 * t183;
t101 = t175 * t154 + t176 * t257;
t93 = -pkin(8) * t259 + t101;
t266 = t178 * t80 + t181 * t93;
t241 = qJD(3) * t182;
t263 = -qJD(1) * t274 - qJD(5) * t275 - t147 * t241;
t118 = t274 * t182;
t262 = -qJD(3) * t118 + t179 * t272 + t127;
t184 = qJD(3) ^ 2;
t255 = t184 * t179;
t254 = t184 * t182;
t185 = qJD(1) ^ 2;
t253 = t185 * qJ(2);
t240 = qJD(3) * t183;
t223 = t182 * t240;
t91 = t175 * t119 + t176 * t223;
t247 = t179 ^ 2 - t182 ^ 2;
t246 = -t184 - t185;
t242 = qJD(3) * t179;
t229 = qJD(2) * t231;
t170 = -pkin(4) * t176 - pkin(3);
t224 = t175 * t242;
t167 = t179 * t240;
t109 = -pkin(4) * t227 + t153;
t222 = t251 * pkin(5) - t109;
t219 = qJD(6) * t10 + t5;
t104 = t176 * t119;
t56 = t104 + (t220 * t182 + t230) * qJD(3);
t70 = pkin(8) * t224 + t91;
t216 = -t178 * t70 + t181 * t56;
t215 = -t178 * t93 + t181 * t80;
t213 = -t181 * t158 - t159 * t178;
t138 = pkin(4) * t259 - t182 * t183;
t212 = -t141 + t243;
t55 = t181 * t58;
t65 = pkin(9) * t274 + t249;
t208 = pkin(5) * t244 + t252 * pkin(9) + t147 * qJD(4) + t249 * qJD(5) + qJD(6) * t65 - t178 * t72 + t55;
t64 = -pkin(9) * t147 + t213;
t207 = t251 * pkin(9) - qJD(6) * t64 + t271;
t115 = t147 * t179;
t206 = qJD(6) * t115 + t262;
t205 = qJD(6) * t275 - t263;
t203 = -t52 * t175 + t53 * t176;
t202 = -t175 * t74 - t176 * t73;
t201 = -t175 * t73 + t176 * t74;
t20 = pkin(5) * t179 - pkin(9) * t118 + t215;
t116 = t147 * t182;
t21 = -pkin(9) * t116 + t266;
t200 = t177 * t20 + t180 * t21;
t60 = t180 * t116 + t118 * t177;
t61 = -t116 * t177 + t118 * t180;
t197 = (-t139 - t233) * t179;
t125 = -pkin(4) * t224 + t167;
t194 = t178 * t56 + t181 * t70 + t80 * t237 - t93 * t238;
t148 = t164 * t242;
t102 = -pkin(4) * t211 + t148;
t191 = -t124 + (t164 + t273) * t182;
t188 = -qJ(4) * t241 + (t124 + t214) * t179;
t161 = t179 * t171;
t108 = -pkin(5) * t274 + t170;
t100 = -t175 * t257 + t137;
t90 = -t175 * t223 + t104;
t89 = pkin(5) * t116 + t138;
t69 = -t178 * t179 * t233 - t181 * t224 + t182 * t286;
t67 = -qJD(3) * t275 - t182 * t272;
t44 = pkin(5) * t69 + t125;
t22 = pkin(5) * t43 + t102;
t15 = t61 * qJD(6) + t177 * t67 + t180 * t69;
t14 = -t60 * qJD(6) - t177 * t69 + t180 * t67;
t7 = -pkin(9) * t69 + t194;
t6 = pkin(5) * t241 - t67 * pkin(9) - qJD(5) * t266 + t216;
t1 = t10 * t180 - t13 * t177;
t3 = [0, 0, 0, 0, t229, qJ(2) * t229, -0.2e1 * t161, 0.2e1 * t247 * t232, -t255, -t254, 0, -t183 * t255 + (qJ(2) * t241 + qJD(2) * t179) * t231, -t183 * t254 + (-qJ(2) * t242 + qJD(2) * t182) * t231 (qJD(1) * t90 + t52) * t179 + ((qJD(1) * t100 + t73) * t182 + (t139 * t183 + t175 * t191) * t179) * qJD(3) (-qJD(1) * t91 - t53) * t179 + ((-qJD(1) * t101 - t74) * t182 + (t141 * t183 + t176 * t191) * t179) * qJD(3), -t91 * t139 - t90 * t141 + (-t175 * t53 - t176 * t52) * t182 + ((t100 * t176 + t101 * t175) * qJD(1) - t202) * t242, t52 * t100 + t53 * t101 + t73 * t90 + t74 * t91 + (t124 - t261) * t167, t118 * t42 - t67 * t81, -t116 * t42 - t118 * t43 - t67 * t82 + t69 * t81, t67 * t168 + t42 * t179 + (qJD(1) * t118 - t81) * t241, -t69 * t168 - t43 * t179 + (-qJD(1) * t116 - t82) * t241, t168 * t241 + t161, t216 * t168 + t217 * t179 + t125 * t82 + t138 * t43 + t102 * t116 + t92 * t69 + (-t168 * t266 - t17 * t179) * qJD(5) + (qJD(1) * t215 + t16) * t241, -t194 * t168 - t195 * t179 - t125 * t81 + t138 * t42 + t102 * t118 + t92 * t67 + (-t266 * qJD(1) - t17) * t241, -t14 * t199 + t61 * t8, t14 * t276 + t15 * t199 + t186 * t61 - t60 * t8, t14 * t163 + t8 * t179 + (qJD(1) * t61 - t199) * t241, -t15 * t163 + t186 * t179 + (-qJD(1) * t60 + t276) * t241, t163 * t241 + t161 (-t177 * t7 + t180 * t6) * t163 + t228 * t179 - t44 * t276 - t89 * t186 + t22 * t60 + t32 * t15 + (-t163 * t200 - t179 * t2) * qJD(6) + ((-t177 * t21 + t180 * t20) * qJD(1) + t1) * t241, t11 * t179 + t32 * t14 + t22 * t61 - t44 * t199 + t89 * t8 + (-(-qJD(6) * t21 + t6) * t163 - t4 * t179) * t177 + (-(qJD(6) * t20 + t7) * t163 - t219 * t179) * t180 + (-qJD(1) * t200 - t2) * t241; 0, 0, 0, 0, -t185, -t253, 0, 0, 0, 0, 0, t246 * t179, t246 * t182 (-t176 * t185 + (t139 - t226) * qJD(3)) * t179 (t175 * t185 + (t141 - t225) * qJD(3)) * t179 (-t139 * t176 + t141 * t175) * t241 + (t139 * t175 + t141 * t176) * qJD(1), t203 * t179 + t202 * qJD(1) + (t124 * t179 + (t201 - t153) * t182) * qJD(3), 0, 0, 0, 0, 0, -t182 * t43 + t263 * t168 + (-t115 * t244 + t179 * t82) * qJD(3), -t182 * t42 + t262 * t168 + (-t179 * t81 - t244 * t275) * qJD(3), 0, 0, 0, 0, 0, t182 * t186 + (t177 * t206 - t180 * t205) * t163 + ((-t115 * t180 - t177 * t275) * t244 - t179 * t276) * qJD(3), -t182 * t8 + (t177 * t205 + t180 * t206) * t163 + (-(-t115 * t177 + t180 * t275) * t244 - t179 * t199) * qJD(3); 0, 0, 0, 0, 0, 0, t182 * t185 * t179, -t247 * t185, 0, 0, 0, -t182 * t253, t179 * t253, t164 * t197 + (t175 * t188 - t179 * t94 - t182 * t73) * qJD(1), t212 * t153 + (t176 * t188 + t179 * t95 + t182 * t74) * qJD(1), t95 * t139 + t94 * t141 + (-qJD(4) * t139 - t73 * t245 + t53) * t176 + (qJD(4) * t141 - t74 * t245 - t52) * t175, -t73 * t94 - t74 * t95 + (-t124 - t265) * t153 + t201 * qJD(4) + t203 * qJ(4), t42 * t147 - t252 * t81, -t147 * t43 + t251 * t81 - t252 * t82 + t274 * t42, t252 * t168 + (qJD(3) * t147 + t81) * t244, -t251 * t168 + (qJD(3) * t274 + t82) * t244, -t168 * t244, -t102 * t274 - t109 * t82 + t170 * t43 + t251 * t92 + (-t55 - t196 * t181 + (qJD(5) * t158 + t287) * t178) * t168 + (qJD(3) * t213 - t16) * t244, t102 * t147 + t109 * t81 + t170 * t42 + t252 * t92 + t271 * t168 + (-qJD(3) * t249 + t17) * t244, -t199 * t269 + t8 * t88, t186 * t88 + t199 * t268 + t269 * t276 - t8 * t87, t269 * t163 + (qJD(3) * t88 + t199) * t244, -t268 * t163 + (-qJD(3) * t87 - t276) * t244, -t163 * t244, -t108 * t186 + t22 * t87 + t268 * t32 - t222 * t276 + (t177 * t207 - t180 * t208) * t163 + ((-t177 * t65 + t180 * t64) * qJD(3) - t1) * t244, t108 * t8 + t22 * t88 + t269 * t32 - t222 * t199 + (t177 * t208 + t180 * t207) * t163 + (-(t177 * t64 + t180 * t65) * qJD(3) + t2) * t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212 * t245, qJD(1) * t197, -t139 ^ 2 - t141 ^ 2, t139 * t74 + t141 * t73 + t148, 0, 0, 0, 0, 0, t43 - t283, t42 - t285, 0, 0, 0, 0, 0, -t186 - t282, t8 + t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81 * t82, t81 ^ 2 - t82 ^ 2, t42 + t285, -t43 - t283, t171, t17 * t168 + t81 * t92 + t187, t16 * t168 + t82 * t92 - t195, t284, t281, t279, t277, t171 -(-t12 * t177 - t264) * t163 + (-t163 * t235 + t171 * t180 - t276 * t81) * pkin(5) + t278 (-t13 * t163 - t4) * t177 + (t12 * t163 - t219) * t180 + (-t163 * t234 - t171 * t177 - t199 * t81) * pkin(5) + t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, t281, t279, t277, t171, t2 * t163 + t278, t1 * t163 - t177 * t4 - t180 * t219 + t280;];
tauc_reg  = t3;
