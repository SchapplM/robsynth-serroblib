% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:07
% EndTime: 2019-03-09 00:11:18
% DurationCPUTime: 3.94s
% Computational Cost: add. (4428->382), mult. (11012->543), div. (0->0), fcn. (8161->10), ass. (0->201)
t177 = sin(qJ(3));
t181 = cos(qJ(3));
t198 = pkin(3) * t177 - pkin(9) * t181;
t136 = t198 * qJD(3);
t176 = sin(qJ(4));
t178 = sin(qJ(2));
t180 = cos(qJ(4));
t239 = qJD(3) * t177;
t173 = sin(pkin(6));
t245 = qJD(1) * t173;
t182 = cos(qJ(2));
t253 = t181 * t182;
t280 = pkin(8) * t176;
t308 = t180 * t136 + t239 * t280 - (-t176 * t253 + t178 * t180) * t245;
t141 = -pkin(3) * t181 - pkin(9) * t177 - pkin(2);
t234 = qJD(4) * t180;
t307 = -(t176 * t178 + t180 * t253) * t245 + t176 * t136 + t141 * t234;
t254 = t180 * t181;
t162 = pkin(8) * t254;
t193 = pkin(4) * t177 - pkin(10) * t254;
t306 = t193 * qJD(3) + (-t162 + (pkin(10) * t177 - t141) * t176) * qJD(4) + t308;
t236 = qJD(4) * t176;
t238 = qJD(3) * t180;
t237 = qJD(3) * t181;
t213 = t176 * t237;
t297 = t177 * t234 + t213;
t305 = -t297 * pkin(10) + (-t177 * t238 - t181 * t236) * pkin(8) + t307;
t241 = qJD(2) * t181;
t221 = t176 * t241;
t284 = pkin(9) + pkin(10);
t225 = qJD(4) * t284;
t133 = t198 * qJD(2);
t215 = t178 * t245;
t139 = qJD(2) * pkin(8) + t215;
t174 = cos(pkin(6));
t244 = qJD(1) * t181;
t99 = -t177 * t139 + t174 * t244;
t267 = t176 * t133 + t180 * t99;
t300 = -pkin(10) * t221 + t176 * t225 + t267;
t203 = t180 * t133 - t176 * t99;
t304 = t193 * qJD(2) + t180 * t225 + t203;
t242 = qJD(2) * t177;
t126 = -t176 * t242 + t238;
t240 = qJD(3) * t176;
t127 = t180 * t242 + t240;
t175 = sin(qJ(5));
t179 = cos(qJ(5));
t195 = t126 * t175 + t179 * t127;
t285 = t195 ^ 2;
t70 = -t179 * t126 + t127 * t175;
t68 = t70 ^ 2;
t303 = -t68 + t285;
t302 = qJ(6) * t70;
t301 = t195 * t70;
t129 = t175 * t180 + t176 * t179;
t229 = qJD(4) + qJD(5);
t79 = t229 * t129;
t269 = t129 * t241 - t79;
t259 = t175 * t176;
t128 = -t179 * t180 + t259;
t232 = qJD(5) * t179;
t299 = -t128 * t241 - t179 * t234 - t180 * t232 + t229 * t259;
t231 = qJD(2) * qJD(3);
t298 = qJD(3) * qJD(4) + t181 * t231;
t160 = -qJD(4) + t241;
t151 = -qJD(5) + t160;
t235 = qJD(4) * t177;
t211 = qJD(2) * t235;
t226 = t298 * t176 + t180 * t211;
t233 = qJD(5) * t175;
t95 = -t176 * t211 + t298 * t180;
t24 = -t126 * t232 + t127 * t233 + t175 * t226 - t179 * t95;
t296 = -t151 * t70 - t24;
t165 = t177 * t231;
t243 = qJD(2) * t173;
t222 = t182 * t243;
t63 = -t139 * t239 + (qJD(3) * t174 + t222) * t244;
t98 = (t136 + t215) * qJD(2);
t205 = t176 * t63 - t180 * t98;
t214 = t182 * t245;
t102 = t141 * qJD(2) - t214;
t258 = t176 * t102;
t260 = t174 * t177;
t156 = qJD(1) * t260;
t100 = t181 * t139 + t156;
t93 = qJD(3) * pkin(9) + t100;
t54 = t180 * t93 + t258;
t186 = -t54 * qJD(4) - t205;
t14 = pkin(4) * t165 - t95 * pkin(10) + t186;
t191 = t102 * t234 + t176 * t98 + t180 * t63 - t93 * t236;
t19 = -t226 * pkin(10) + t191;
t53 = t180 * t102 - t176 * t93;
t40 = -pkin(10) * t127 + t53;
t30 = -pkin(4) * t160 + t40;
t41 = pkin(10) * t126 + t54;
t204 = -t175 * t14 - t179 * t19 - t30 * t232 + t41 * t233;
t92 = -qJD(3) * pkin(3) - t99;
t60 = -pkin(4) * t126 + t92;
t295 = t60 * t70 + t204;
t31 = pkin(5) * t70 + qJD(6) + t60;
t294 = t195 * t31;
t125 = t180 * t141;
t255 = t177 * t180;
t75 = -pkin(10) * t255 + t125 + (-pkin(4) - t280) * t181;
t247 = t176 * t141 + t162;
t257 = t176 * t177;
t84 = -pkin(10) * t257 + t247;
t293 = t306 * t175 + t305 * t179 + t75 * t232 - t84 * t233;
t292 = t305 * t175 - t306 * t179;
t291 = t304 * t179;
t199 = -t100 + (-t221 + t236) * pkin(4);
t290 = qJ(6) * t195;
t145 = t284 * t176;
t146 = t284 * t180;
t248 = -t175 * t145 + t179 * t146;
t289 = -t145 * t232 - t146 * t233 - t304 * t175 - t300 * t179;
t218 = t176 * t235;
t220 = t180 * t237;
t288 = -t218 + t220;
t37 = t179 * t41;
t16 = t175 * t30 + t37;
t210 = t179 * t14 - t175 * t19;
t187 = -t16 * qJD(5) + t210;
t287 = -t60 * t195 + t187;
t25 = t195 * qJD(5) + t175 * t95 + t179 * t226;
t286 = -t151 * t195 - t25;
t35 = t175 * t41;
t15 = t179 * t30 - t35;
t8 = t15 - t290;
t7 = -pkin(5) * t151 + t8;
t283 = t7 - t8;
t106 = t128 * t177;
t275 = t175 * t75 + t179 * t84;
t48 = t175 * t213 + t79 * t177 - t179 * t220;
t282 = pkin(5) * t239 + t48 * qJ(6) - qJD(5) * t275 + t106 * qJD(6) - t292;
t105 = t129 * t177;
t49 = -t233 * t257 + (t229 * t255 + t213) * t179 + t288 * t175;
t281 = -qJ(6) * t49 - qJD(6) * t105 + t293;
t279 = t269 * qJ(6) - qJD(6) * t128 + t289;
t278 = -pkin(5) * t242 + t299 * qJ(6) - t248 * qJD(5) - t129 * qJD(6) + t300 * t175 - t291;
t277 = t179 * t40 - t35;
t274 = qJD(2) * pkin(2);
t273 = t176 * t92;
t201 = t177 * t222;
t64 = qJD(1) * t201 + qJD(3) * t156 + t139 * t237;
t272 = t64 * t176;
t271 = t64 * t180;
t270 = t95 * t176;
t266 = t126 * t160;
t265 = t127 * t160;
t264 = t160 * t180;
t263 = t173 * t178;
t262 = t173 * t182;
t184 = qJD(2) ^ 2;
t261 = t173 * t184;
t256 = t176 * t181;
t183 = qJD(3) ^ 2;
t252 = t183 * t177;
t251 = t183 * t181;
t137 = pkin(4) * t257 + t177 * pkin(8);
t171 = t177 ^ 2;
t246 = -t181 ^ 2 + t171;
t227 = t178 * t261;
t101 = t297 * pkin(4) + pkin(8) * t237;
t168 = -pkin(4) * t180 - pkin(3);
t223 = t178 * t243;
t219 = t160 * t236;
t216 = t160 * t234;
t209 = -t175 * t40 - t37;
t207 = -t175 * t84 + t179 * t75;
t202 = -t179 * t145 - t146 * t175;
t200 = t181 * t222;
t140 = -t214 - t274;
t197 = -t140 - t214;
t112 = t181 * t263 + t260;
t192 = -t112 * t180 + t176 * t262;
t82 = -t112 * t176 - t180 * t262;
t44 = t175 * t192 + t179 * t82;
t45 = t175 * t82 - t179 * t192;
t194 = qJD(2) * t171 - t160 * t181;
t111 = -t174 * t181 + t177 * t263;
t47 = t226 * pkin(4) + t64;
t188 = qJD(3) * (-t197 - t274);
t17 = t25 * pkin(5) + t47;
t167 = pkin(4) * t179 + pkin(5);
t81 = t112 * qJD(3) + t201;
t80 = -t111 * qJD(3) + t200;
t62 = -qJ(6) * t128 + t248;
t61 = -qJ(6) * t129 + t202;
t34 = t82 * qJD(4) + t176 * t223 + t80 * t180;
t33 = t192 * qJD(4) - t80 * t176 + t180 * t223;
t27 = -qJ(6) * t105 + t275;
t26 = -pkin(5) * t181 + qJ(6) * t106 + t207;
t11 = t277 - t290;
t10 = t209 + t302;
t9 = t16 - t302;
t6 = -t45 * qJD(5) - t175 * t34 + t179 * t33;
t5 = t44 * qJD(5) + t175 * t33 + t179 * t34;
t2 = -qJ(6) * t25 - qJD(6) * t70 - t204;
t1 = pkin(5) * t165 + t24 * qJ(6) - qJD(6) * t195 + t187;
t3 = [0, 0, -t227, -t182 * t261, 0, 0, 0, 0, 0, -t181 * t227 + (-t81 - t201) * qJD(3), t177 * t227 + (-t80 - t200) * qJD(3), 0, 0, 0, 0, 0, t111 * t226 - t81 * t126 - t33 * t160 + t82 * t165, t111 * t95 + t127 * t81 + t160 * t34 + t165 * t192, 0, 0, 0, 0, 0, t111 * t25 - t151 * t6 + t165 * t44 + t70 * t81, -t111 * t24 + t151 * t5 - t165 * t45 + t195 * t81, -t195 * t6 + t24 * t44 - t25 * t45 - t5 * t70, t1 * t44 + t111 * t17 + t2 * t45 + t31 * t81 + t5 * t9 + t6 * t7; 0, 0, 0, 0, 0.2e1 * t181 * t165, -0.2e1 * t246 * t231, t251, -t252, 0, -pkin(8) * t251 + t177 * t188, pkin(8) * t252 + t181 * t188, t288 * t127 + t95 * t255 (t126 * t180 - t127 * t176) * t237 + (-t180 * t226 - t270 + (-t126 * t176 - t127 * t180) * qJD(4)) * t177, t160 * t218 - t95 * t181 + (t127 * t177 + t180 * t194) * qJD(3), t177 * t216 + t226 * t181 + (t126 * t177 - t176 * t194) * qJD(3) (-t160 - t241) * t239 (t141 * t236 - t308) * t160 + ((-pkin(8) * t126 + t273) * qJD(3) + (t258 + (pkin(8) * t160 + t93) * t180) * qJD(4) + t205) * t181 + (pkin(8) * t226 + t272 + t92 * t234 + t126 * t214 + ((-pkin(8) * t256 + t125) * qJD(2) + t53) * qJD(3)) * t177, t307 * t160 + (t92 * t238 + (qJD(3) * t127 - t219) * pkin(8) + t191) * t181 + (-t127 * t214 - t92 * t236 + pkin(8) * t95 + t271 + (-pkin(8) * t264 - t247 * qJD(2) - t54) * qJD(3)) * t177, t106 * t24 - t195 * t48, t105 * t24 + t106 * t25 - t195 * t49 + t48 * t70, t48 * t151 + t24 * t181 + (-qJD(2) * t106 + t195) * t239, t49 * t151 + t25 * t181 + (-qJD(2) * t105 - t70) * t239 (-t151 - t241) * t239, -t210 * t181 + t101 * t70 + t137 * t25 + t47 * t105 + t60 * t49 + t292 * t151 + (t151 * t275 + t16 * t181) * qJD(5) + (-t70 * t214 + (qJD(2) * t207 + t15) * qJD(3)) * t177, -t204 * t181 + t101 * t195 - t137 * t24 - t47 * t106 - t60 * t48 + t293 * t151 + (-t195 * t214 + (-t275 * qJD(2) - t16) * qJD(3)) * t177, t1 * t106 - t2 * t105 - t195 * t282 + t26 * t24 - t27 * t25 - t281 * t70 + t7 * t48 - t9 * t49, t2 * t27 + t1 * t26 + t17 * (pkin(5) * t105 + t137) + t281 * t9 + t282 * t7 + (pkin(5) * t49 - t177 * t214 + t101) * t31; 0, 0, 0, 0, -t177 * t184 * t181, t246 * t184, 0, 0, 0, qJD(3) * t100 - t140 * t242 - t64, t197 * t241, -t127 * t264 + t270 (t95 - t266) * t180 + (-t226 + t265) * t176, -t216 + (t160 * t254 + (-t127 + t240) * t177) * qJD(2), t219 + (-t160 * t256 + (-t126 + t238) * t177) * qJD(2), t160 * t242, -pkin(3) * t226 - t271 + t203 * t160 + t100 * t126 + (pkin(9) * t264 + t273) * qJD(4) + (-t53 * t177 + (-pkin(9) * t239 - t181 * t92) * t176) * qJD(2), -pkin(3) * t95 + t272 - t267 * t160 - t100 * t127 + (-pkin(9) * t160 * t176 + t180 * t92) * qJD(4) + (-t92 * t254 + (-pkin(9) * t238 + t54) * t177) * qJD(2), -t24 * t129 - t195 * t299, t24 * t128 - t129 * t25 + t195 * t269 + t299 * t70, t299 * t151 + (qJD(3) * t129 - t195) * t242, -t269 * t151 + (-qJD(3) * t128 + t70) * t242, t151 * t242, t47 * t128 + t168 * t25 + t199 * t70 - t269 * t60 + (t146 * t232 + (-qJD(5) * t145 - t300) * t175 + t291) * t151 + (qJD(3) * t202 - t15) * t242, t47 * t129 - t168 * t24 + t199 * t195 - t299 * t60 + t289 * t151 + (-qJD(3) * t248 + t16) * t242, -t1 * t129 - t2 * t128 - t195 * t278 + t61 * t24 - t62 * t25 + t269 * t9 - t279 * t70 + t299 * t7, t2 * t62 + t1 * t61 + t17 * (pkin(5) * t128 + t168) + t279 * t9 + t278 * t7 + (-t269 * pkin(5) + t199) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127 * t126, -t126 ^ 2 + t127 ^ 2, t95 + t266, -t226 - t265, t165, -t92 * t127 - t54 * t160 + t186, -t126 * t92 - t160 * t53 - t191, t301, t303, t296, t286, t165, t209 * t151 + (-t127 * t70 + t151 * t233 + t165 * t179) * pkin(4) + t287, -t277 * t151 + (-t127 * t195 + t151 * t232 - t165 * t175) * pkin(4) + t295, t10 * t195 + t11 * t70 + t167 * t24 + t9 * t195 - t7 * t70 + (-t175 * t25 + (t175 * t195 - t179 * t70) * qJD(5)) * pkin(4), -pkin(5) * t294 + t1 * t167 - t7 * t10 - t9 * t11 + (-t31 * t127 + t2 * t175 + (-t175 * t7 + t179 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, t303, t296, t286, t165, -t16 * t151 + t287, -t15 * t151 + t295, pkin(5) * t24 - t283 * t70, t283 * t9 + (t1 - t294) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 - t285, t195 * t7 + t9 * t70 + t17;];
tauc_reg  = t3;
