% Calculate vector of inverse dynamics joint torques for
% S5RRPPR3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:33
% DurationCPUTime: 3.46s
% Computational Cost: add. (7580->393), mult. (5531->476), div. (0->0), fcn. (4020->8), ass. (0->221)
t189 = sin(qJ(1));
t191 = cos(qJ(1));
t192 = qJD(1) ^ 2;
t207 = (-qJDD(1) * t189 - t191 * t192) * pkin(1);
t323 = t207 - g(1);
t186 = qJD(1) + qJD(2);
t188 = sin(qJ(5));
t190 = cos(qJ(5));
t273 = Icges(6,4) * t190;
t150 = -Icges(6,2) * t188 + t273;
t221 = Icges(6,1) * t188 + t273;
t257 = t150 + t221;
t274 = Icges(6,4) * t188;
t152 = Icges(6,1) * t190 - t274;
t220 = Icges(6,2) * t190 + t274;
t258 = -t220 + t152;
t322 = (t188 * t257 - t190 * t258) * t186;
t187 = qJ(1) + qJ(2);
t182 = cos(t187);
t175 = pkin(2) * t182;
t180 = pkin(8) + t187;
t173 = sin(t180);
t174 = cos(t180);
t308 = t174 * pkin(3) + t173 * qJ(4);
t244 = t175 + t308;
t278 = pkin(1) * qJD(1);
t249 = t189 * t278;
t181 = sin(t187);
t126 = rSges(3,1) * t181 + rSges(3,2) * t182;
t269 = t126 * t186;
t92 = -t249 - t269;
t230 = rSges(6,1) * t188 + rSges(6,2) * t190;
t252 = qJD(5) * t186;
t105 = qJDD(5) * t173 + t174 * t252;
t131 = t230 * qJD(5);
t281 = rSges(6,2) * t188;
t282 = rSges(6,1) * t190;
t156 = -t281 + t282;
t184 = t186 ^ 2;
t185 = qJDD(1) + qJDD(2);
t255 = qJD(4) * t186;
t198 = qJDD(4) * t173 + t174 * t255 + t207;
t162 = t174 * qJ(4);
t286 = pkin(3) * t173;
t117 = -t162 + t286;
t287 = pkin(2) * t181;
t232 = -pkin(7) * t173 - t287;
t89 = -t173 * rSges(6,3) + t174 * t230;
t213 = -t117 + t232 + t89;
t170 = t174 * pkin(7);
t231 = -t170 - t175;
t254 = qJD(5) * t173;
t168 = t174 * rSges(6,3);
t247 = qJD(5) * t281;
t248 = qJD(5) * t282;
t214 = (t247 - t248) * t174;
t54 = (t173 * t230 + t168) * t186 + t214;
t160 = qJD(4) * t174;
t80 = t186 * t308 - t160;
t11 = -t131 * t254 + t105 * t156 + (-t54 - t80) * t186 + t231 * t184 + t213 * t185 + t198;
t321 = -g(1) + t11;
t263 = t174 * t186;
t142 = rSges(5,2) * t263;
t279 = rSges(5,3) * t174;
t118 = rSges(5,2) * t173 + t279;
t236 = -t117 + t118 - t287;
t261 = t182 * t184;
t266 = t173 * t186;
t23 = -pkin(2) * t261 + (-rSges(5,3) * t266 + t142 - t80) * t186 + t236 * t185 + t198;
t320 = -g(1) + t23;
t119 = rSges(4,1) * t173 + rSges(4,2) * t174;
t140 = rSges(4,2) * t266;
t319 = -t185 * t119 - t186 * (rSges(4,1) * t263 - t140) + (-t181 * t185 - t261) * pkin(2) + t323;
t172 = t182 * rSges(3,1);
t262 = t181 * t186;
t111 = -rSges(3,2) * t262 + t172 * t186;
t318 = -t111 * t186 - t126 * t185 + t323;
t106 = qJDD(5) * t174 - t173 * t252;
t183 = t191 * pkin(1);
t288 = pkin(1) * t189;
t237 = qJDD(1) * t183 - t192 * t288;
t216 = t185 * t175 + t237;
t138 = qJ(4) * t263;
t159 = qJD(4) * t173;
t260 = t138 + t159;
t306 = t185 * t308 + t173 * t255 + t186 * (-pkin(3) * t266 + t260);
t245 = t173 * t248 + t230 * t263;
t55 = (-rSges(6,3) * t186 - t247) * t173 + t245;
t264 = t173 * t190;
t265 = t173 * t188;
t88 = rSges(6,1) * t265 + rSges(6,2) * t264 + t168;
t12 = -t106 * t156 + t185 * t88 + t186 * t55 + t232 * t184 + (pkin(7) * t185 + qJD(5) * t131 - qJDD(4)) * t174 + t216 + t306;
t317 = -g(2) + t12;
t167 = t173 * rSges(5,3);
t121 = -rSges(5,2) * t174 + t167;
t199 = -t184 * t287 + t216;
t259 = rSges(5,2) * t266 + rSges(5,3) * t263;
t24 = -qJDD(4) * t174 + t185 * t121 + t186 * t259 + t199 + t306;
t316 = -g(2) + t24;
t169 = t174 * rSges(4,1);
t122 = -rSges(4,2) * t173 + t169;
t315 = -t119 * t184 + t185 * t122 - g(2) + t199;
t127 = -rSges(3,2) * t181 + t172;
t314 = t127 * t185 - t186 * t269 - g(2) + t237;
t310 = t121 + t244;
t313 = t186 * t310;
t219 = Icges(6,5) * t188 + Icges(6,6) * t190;
t312 = t219 * t186;
t85 = -Icges(6,6) * t173 + t174 * t220;
t87 = -Icges(6,5) * t173 + t174 * t221;
t222 = t188 * t87 + t190 * t85;
t311 = t222 * t174;
t309 = t122 + t175;
t107 = t186 * t117;
t256 = t159 - t107;
t307 = -t186 * t118 - t256 + t259;
t44 = t188 * t85 - t190 * t87;
t209 = t220 * t186;
t302 = -Icges(6,6) * t186 + qJD(5) * t150;
t50 = t173 * t209 - t174 * t302;
t210 = t221 * t186;
t301 = -Icges(6,5) * t186 + qJD(5) * t152;
t52 = t173 * t210 - t174 * t301;
t83 = -Icges(6,3) * t173 + t174 * t219;
t305 = qJD(5) * t44 + t186 * t83 + t188 * t52 + t190 * t50;
t253 = qJD(5) * t174;
t304 = -t156 * t253 + t186 * (t308 - t231 + t88);
t148 = Icges(6,5) * t190 - Icges(6,6) * t188;
t303 = -Icges(6,3) * t186 + qJD(5) * t148;
t129 = t220 * qJD(5);
t130 = t221 * qJD(5);
t217 = t150 * t188 - t152 * t190;
t300 = qJD(5) * t217 + t129 * t190 + t130 * t188 + t148 * t186;
t84 = Icges(6,6) * t174 + t173 * t220;
t143 = Icges(6,4) * t264;
t86 = Icges(6,1) * t265 + Icges(6,5) * t174 + t143;
t223 = t188 * t84 - t190 * t86;
t51 = t173 * t302 + t174 * t209;
t53 = t173 * t301 + t174 * t210;
t82 = Icges(6,3) * t174 + t173 * t219;
t299 = qJD(5) * t223 + t186 * t82 - t188 * t53 - t190 * t51;
t79 = t186 * t89;
t298 = -t173 * t247 - t186 * t232 + t245 + t260 - t79;
t275 = -t152 * t174 + t85;
t284 = t150 * t174 + t87;
t296 = t188 * t275 - t190 * t284;
t276 = t152 * t173 - t84;
t285 = -Icges(6,2) * t265 + t143 + t86;
t295 = -t188 * t276 - t190 * t285;
t294 = -m(5) - m(6);
t293 = -pkin(3) - pkin(7);
t292 = t105 / 0.2e1;
t291 = t106 / 0.2e1;
t290 = t173 / 0.2e1;
t289 = -t174 / 0.2e1;
t218 = t150 * t190 + t152 * t188;
t96 = t173 * t148;
t62 = t174 * t218 - t96;
t277 = t62 * t186;
t267 = t148 * t174;
t251 = -rSges(6,3) + t293;
t28 = t174 * t82 + t84 * t264 + t86 * t265;
t29 = -t174 * t83 - t85 * t264 - t87 * t265;
t250 = t191 * t278;
t243 = -t254 / 0.2e1;
t242 = t254 / 0.2e1;
t241 = -t253 / 0.2e1;
t240 = t253 / 0.2e1;
t94 = -t119 - t287;
t238 = t162 - t287;
t234 = t159 - t249;
t233 = -t160 + t250;
t157 = rSges(2,1) * t191 - rSges(2,2) * t189;
t155 = rSges(2,1) * t189 + rSges(2,2) * t191;
t228 = t173 * t29 + t174 * t28;
t224 = t188 * t86 + t190 * t84;
t73 = t173 * t82;
t30 = -t174 * t224 + t73;
t31 = -t173 * t83 + t311;
t227 = t173 * t31 + t174 * t30;
t112 = t156 * t254;
t215 = t112 + t234;
t32 = t186 * t213 + t215;
t33 = t233 + t304;
t226 = t173 * t32 - t174 * t33;
t225 = -t173 * t88 - t174 * t89;
t211 = t160 - t214;
t206 = -pkin(2) * t262 - t249;
t60 = t170 + t244 + t88;
t202 = t173 * t312 - t174 * t303 - t186 * t222;
t201 = t173 * t303 + t174 * t312 + t186 * t224;
t200 = -t219 * qJD(5) + t186 * t218;
t71 = t279 + (rSges(5,2) - pkin(3)) * t173 + t238;
t59 = t173 * t293 + t238 + t89;
t10 = qJD(5) * t227 - t277;
t16 = -qJD(5) * t224 - t188 * t51 + t190 * t53;
t17 = qJD(5) * t222 - t188 * t50 + t190 * t52;
t20 = t200 * t173 + t174 * t300;
t21 = -t173 * t300 + t200 * t174;
t61 = t173 * t218 + t267;
t58 = t61 * t186;
t9 = qJD(5) * t228 + t58;
t196 = (t58 + ((-t30 + t73 + t29) * t173 + (t31 - t311 + (-t224 + t83) * t173 + t28) * t174) * qJD(5)) * t243 + t44 * t292 - t105 * t62 / 0.2e1 + (-qJD(5) * t218 + t129 * t188 - t130 * t190) * t186 + (-t223 + t61) * t291 + (t277 + (t173 ^ 2 * t83 + (-t73 + t29 + (t224 + t83) * t174) * t174) * qJD(5) + t10) * t241 + (t16 + t21) * t240 + (t17 + t20 + t9) * t242 + (Icges(4,3) + Icges(3,3) + Icges(5,1) - t217) * t185;
t76 = t186 * t94 - t249;
t77 = t186 * t309 + t250;
t195 = (t76 * (-t169 - t175) + t77 * t94) * t186;
t46 = t186 * t236 + t234;
t47 = t233 + t313;
t194 = (t46 * (-t167 - t244) + t47 * (-t286 - t287)) * t186;
t193 = ((-t181 * t33 - t182 * t32) * pkin(2) + t32 * t251 * t174 + (t32 * (-qJ(4) - t230) + t33 * t251) * t173) * t186;
t109 = t186 * t119;
t103 = t156 * t174;
t102 = t156 * t173;
t93 = t127 * t186 + t250;
t40 = qJD(5) * t225 + qJD(3);
t13 = -t105 * t88 - t106 * t89 + qJDD(3) + (-t173 * t55 + t174 * t54) * qJD(5);
t6 = t173 * t305 + t202 * t174;
t5 = -t173 * t299 + t201 * t174;
t4 = t202 * t173 - t174 * t305;
t3 = t201 * t173 + t174 * t299;
t1 = [Icges(2,3) * qJDD(1) + t196 + (t314 * (t127 + t183) + t318 * (-t126 - t288) + (-t111 - t250 + t93) * t92) * m(3) + (g(1) * t155 - g(2) * t157 + (t155 ^ 2 + t157 ^ 2) * qJDD(1)) * m(2) + (t32 * (t211 - t250) + t193 + t317 * (t183 + t60) + t321 * (t59 - t288) + (t107 - t215 - t249 + t298 + t32) * t33) * m(6) + (t46 * (t142 - t233) + t194 + t316 * (t183 + t310) + t320 * (t71 - t288) + (-t206 + t46 + t138 + t234 + t307) * t47) * m(5) + (t76 * (t140 - t250) + t195 + t315 * (t183 + t309) + t319 * (t94 - t288) + (t109 - t206 + t76 - t249) * t77) * m(4); t196 + (t193 + t317 * t60 + t321 * t59 + (-t112 - t256 + t298) * t33 + (-t160 + t211 + t304) * t32) * m(6) + (t194 + t316 * t310 + t320 * t71 + (t142 + t313) * t46 + (t186 * t287 + t260 + t307) * t47) * m(5) + (t76 * t140 + t195 + t77 * t109 - (-t287 * t77 - t309 * t76) * t186 + t315 * t309 + t319 * t94) * m(4) + (-t269 * t93 - t111 * t92 + (t186 * t92 + t314) * t127 + (t186 * t93 - t318) * t126) * m(3); m(6) * t13 + (m(4) + m(5)) * qJDD(3) + (-m(4) + t294) * g(3); t294 * (g(1) * t173 - g(2) * t174) + 0.2e1 * (t11 * t290 + t12 * t289) * m(6) + 0.2e1 * (t23 * t290 + t24 * t289) * m(5); -t9 * t266 / 0.2e1 + t174 * (t105 * t29 + t106 * t28 + t185 * t61 + t186 * t21 + (t173 * t6 + t174 * t5) * qJD(5)) / 0.2e1 + t228 * t291 + ((t186 * t29 + t5) * t174 + (-t186 * t28 + t6) * t173) * t240 + t10 * t263 / 0.2e1 + (t105 * t31 + t106 * t30 - t185 * t62 + t186 * t20 + (t173 * t4 + t174 * t3) * qJD(5)) * t290 + t227 * t292 + ((t186 * t31 + t3) * t174 + (-t186 * t30 + t4) * t173) * t242 + t185 * (t173 * t44 - t174 * t223) / 0.2e1 + t186 * ((t186 * t44 + t16) * t174 + (t186 * t223 + t17) * t173) / 0.2e1 + ((t96 * t253 - t312) * t174 + (-t322 + (t296 * t173 + (-t295 - t267) * t174) * qJD(5)) * t173) * t241 + ((-t254 * t267 - t312) * t173 + (t322 + (t295 * t174 + (-t296 + t96) * t173) * qJD(5)) * t174) * t243 - t186 * ((-t258 * t188 - t257 * t190) * t186 + ((t173 * t275 + t174 * t276) * t190 + (t173 * t284 - t174 * t285) * t188) * qJD(5)) / 0.2e1 + (t13 * t225 + t40 * ((-t186 * t88 + t54) * t174 + (-t55 + t79) * t173) - t226 * t131 + ((t186 * t32 - t12) * t174 + (t186 * t33 + t11) * t173) * t156 - (t102 * t33 + t103 * t32) * t186 - (t40 * (-t102 * t173 - t103 * t174) - t226 * t230) * qJD(5) - g(1) * t102 + g(2) * t103 + g(3) * t230) * m(6);];
tau = t1;
