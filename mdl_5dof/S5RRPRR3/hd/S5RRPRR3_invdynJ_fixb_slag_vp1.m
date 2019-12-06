% Calculate vector of inverse dynamics joint torques for
% S5RRPRR3
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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:18
% EndTime: 2019-12-05 18:30:24
% DurationCPUTime: 4.04s
% Computational Cost: add. (11057->375), mult. (6222->463), div. (0->0), fcn. (4592->10), ass. (0->228)
t182 = sin(qJ(5));
t184 = cos(qJ(5));
t154 = rSges(6,1) * t182 + rSges(6,2) * t184;
t180 = qJD(1) + qJD(2);
t181 = qJ(1) + qJ(2);
t174 = sin(t181);
t175 = cos(t181);
t228 = rSges(3,1) * t174 + rSges(3,2) * t175;
t109 = t228 * t180;
t183 = sin(qJ(1));
t285 = pkin(1) * qJD(1);
t250 = t183 * t285;
t103 = t109 + t250;
t173 = pkin(9) + t181;
t167 = qJ(4) + t173;
t161 = sin(t167);
t162 = cos(t167);
t113 = pkin(4) * t162 + pkin(8) * t161;
t172 = qJD(4) + t180;
t102 = t172 * t113;
t255 = qJD(5) * t161;
t112 = t154 * t255;
t265 = t162 * t182;
t141 = rSges(6,2) * t265;
t327 = t154 * qJD(5);
t244 = t172 * t141 + t327 * t161;
t264 = t162 * t184;
t213 = -rSges(6,1) * t264 - rSges(6,3) * t161;
t81 = -t141 - t213;
t70 = t172 * t81;
t334 = -t244 - t102 - t70 + t112;
t269 = t161 * t172;
t123 = pkin(4) * t269;
t287 = rSges(6,2) * t182;
t288 = rSges(6,1) * t184;
t156 = -t287 + t288;
t133 = t156 * qJD(5);
t179 = qJDD(1) + qJDD(2);
t171 = qJDD(4) + t179;
t178 = t180 ^ 2;
t186 = qJD(1) ^ 2;
t185 = cos(qJ(1));
t303 = pkin(1) * t185;
t304 = pkin(1) * t183;
t233 = -qJDD(1) * t303 + t186 * t304;
t302 = pkin(2) * t174;
t215 = t178 * t302 + t233;
t166 = cos(t173);
t301 = pkin(2) * t175;
t231 = pkin(3) * t166 + t301;
t165 = sin(t173);
t300 = pkin(3) * t165;
t195 = t178 * t300 - t179 * t231 + t215;
t266 = t162 * t172;
t279 = -t113 - t81;
t267 = t161 * t184;
t251 = rSges(6,1) * t267;
t245 = -t327 * t162 - t172 * t251;
t268 = t161 * t182;
t259 = rSges(6,2) * t268 + t162 * rSges(6,3);
t47 = t172 * t259 + t245;
t253 = qJD(5) * t172;
t89 = qJDD(5) * t161 + t162 * t253;
t12 = t133 * t255 + t89 * t154 + (-pkin(8) * t266 + t123 - t47) * t172 + t279 * t171 + t195;
t326 = -g(2) + t12;
t111 = rSges(5,1) * t162 - t161 * rSges(5,2);
t85 = rSges(5,1) * t269 + rSges(5,2) * t266;
t325 = -t111 * t171 + t172 * t85 - g(2) + t195;
t209 = (-qJDD(1) * t183 - t185 * t186) * pkin(1);
t193 = t209 + (-t174 * t179 - t175 * t178) * pkin(2);
t189 = (-t165 * t179 - t166 * t178) * pkin(3) + t193;
t254 = qJD(5) * t162;
t48 = t172 * t213 + t244;
t243 = -pkin(4) - t288;
t299 = t162 * pkin(8);
t61 = t161 * t243 + t259 + t299;
t90 = qJDD(5) * t162 - t161 * t253;
t11 = -t133 * t254 - t90 * t154 + t61 * t171 + (t48 - t102) * t172 + t189;
t324 = -g(3) + t11;
t226 = -rSges(5,1) * t161 - rSges(5,2) * t162;
t86 = -rSges(5,1) * t266 + rSges(5,2) * t269;
t323 = t171 * t226 + t172 * t86 - g(3) + t189;
t263 = t165 * t180;
t142 = rSges(4,2) * t263;
t227 = -rSges(4,1) * t165 - rSges(4,2) * t166;
t262 = t166 * t180;
t333 = t179 * t227 + t180 * (-rSges(4,1) * t262 + t142) + t193 - g(3);
t160 = t165 * rSges(4,2);
t290 = rSges(4,1) * t166;
t119 = -t160 + t290;
t238 = -t119 - t301;
t258 = -rSges(4,1) * t263 - rSges(4,2) * t262;
t332 = t179 * t238 - t180 * t258 - g(2) + t215;
t260 = t175 * t180;
t261 = t174 * t180;
t110 = -rSges(3,1) * t260 + rSges(3,2) * t261;
t331 = t110 * t180 - t179 * t228 - g(3) + t209;
t125 = rSges(3,1) * t175 - t174 * rSges(3,2);
t330 = t109 * t180 - t125 * t179 - g(2) + t233;
t176 = Icges(6,4) * t184;
t220 = -Icges(6,2) * t182 + t176;
t321 = Icges(6,1) * t182 + t176;
t256 = t321 + t220;
t277 = Icges(6,4) * t182;
t148 = Icges(6,2) * t184 + t277;
t151 = Icges(6,1) * t184 - t277;
t257 = t148 - t151;
t329 = (t182 * t256 + t184 * t257) * t172;
t230 = -t290 - t301;
t328 = -pkin(2) * t260 - t142 + (-t119 - t230) * t180;
t249 = t185 * t285;
t197 = t180 * t231 + t249;
t232 = t300 + t302;
t199 = t180 * t232 + t250;
t78 = Icges(6,4) * t264 - Icges(6,2) * t265 + Icges(6,6) * t161;
t139 = Icges(6,4) * t265;
t80 = Icges(6,1) * t264 + Icges(6,5) * t161 - t139;
t222 = t182 * t78 - t184 * t80;
t322 = t222 * t162;
t216 = t251 - t259;
t69 = t172 * t216;
t214 = t154 * t254 - t172 * (-pkin(4) * t161 + t299) + t69;
t30 = t199 + t214;
t31 = t172 * t279 + t112 - t197;
t320 = t161 * t31 + t162 * t30;
t97 = t172 * t226;
t59 = -t97 + t199;
t98 = t172 * t111;
t60 = -t197 - t98;
t319 = -t59 * t86 + t60 * t85;
t236 = t123 - t245;
t318 = (-t214 + t236) * t31 + t334 * t30;
t211 = t220 * t172;
t314 = -Icges(6,6) * t172 + qJD(5) * t148;
t42 = t161 * t314 - t162 * t211;
t212 = t151 * t172;
t312 = -Icges(6,5) * t172 + qJD(5) * t321;
t44 = t161 * t312 - t162 * t212;
t77 = Icges(6,6) * t162 - t161 * t220;
t138 = Icges(6,4) * t268;
t79 = -Icges(6,1) * t267 + Icges(6,5) * t162 + t138;
t45 = t182 * t79 + t184 * t77;
t147 = Icges(6,5) * t184 - Icges(6,6) * t182;
t75 = Icges(6,3) * t162 - t147 * t161;
t317 = qJD(5) * t45 - t172 * t75 + t182 * t42 - t184 * t44;
t41 = -t161 * t211 - t162 * t314;
t43 = -t161 * t212 - t162 * t312;
t46 = t182 * t80 + t184 * t78;
t76 = Icges(6,5) * t264 - Icges(6,6) * t265 + Icges(6,3) * t161;
t316 = qJD(5) * t46 - t172 * t76 + t182 * t41 - t184 * t43;
t146 = Icges(6,5) * t182 + Icges(6,6) * t184;
t315 = -Icges(6,3) * t172 + qJD(5) * t146;
t131 = t220 * qJD(5);
t132 = t151 * qJD(5);
t219 = t148 * t184 + t182 * t321;
t313 = qJD(5) * t219 + t131 * t182 - t132 * t184 - t146 * t172;
t291 = -Icges(6,2) * t264 - t139 + t80;
t293 = t162 * t321 + t78;
t310 = t182 * t291 + t184 * t293;
t292 = Icges(6,2) * t267 + t138 + t79;
t294 = -t161 * t321 + t77;
t309 = -t182 * t292 - t184 * t294;
t308 = t89 / 0.2e1;
t307 = t90 / 0.2e1;
t306 = m(4) + m(5);
t305 = -rSges(6,3) - pkin(8);
t296 = t162 * t75 + t77 * t268;
t295 = t161 * t75 + t79 * t264;
t282 = t182 * t77;
t281 = t184 * t79;
t218 = t182 * t148 - t184 * t321;
t91 = t146 * t161;
t56 = -t162 * t218 + t91;
t280 = t56 * t172;
t272 = t146 * t162;
t271 = t147 * t172;
t99 = t154 * t161;
t270 = t154 * t162;
t242 = -t255 / 0.2e1;
t241 = t255 / 0.2e1;
t240 = -t254 / 0.2e1;
t239 = t254 / 0.2e1;
t237 = -t76 - t281;
t157 = rSges(2,1) * t185 - t183 * rSges(2,2);
t229 = rSges(2,1) * t183 + rSges(2,2) * t185;
t24 = -t267 * t79 + t296;
t25 = t162 * t76 - t267 * t80 + t78 * t268;
t225 = t25 * t161 + t24 * t162;
t26 = -t265 * t77 + t295;
t27 = t161 * t76 - t322;
t224 = t27 * t161 + t26 * t162;
t223 = -t281 + t282;
t106 = t160 + t230;
t104 = -t125 * t180 - t249;
t208 = pkin(2) * t261 + t250;
t105 = t227 - t302;
t203 = -t271 * t161 - t162 * t315 + t172 * t222;
t202 = t161 * t315 - t271 * t162 + t172 * t223;
t84 = -t111 - t231;
t201 = t147 * qJD(5) + t172 * t218;
t83 = t226 - t232;
t200 = t161 * t216 + t162 * t81;
t62 = t161 * t305 + t162 * t243 + t141;
t10 = qJD(5) * t224 + t280;
t16 = -qJD(5) * t223 + t182 * t44 + t184 * t42;
t17 = -qJD(5) * t222 + t182 * t43 + t184 * t41;
t20 = t201 * t161 - t162 * t313;
t21 = t161 * t313 + t201 * t162;
t55 = t161 * t218 + t272;
t52 = t55 * t172;
t9 = qJD(5) * t225 + t52;
t194 = (t52 + ((t27 + t296 + t322) * t162 + (-t26 + (t237 - t282) * t162 + t25 + t295) * t161) * qJD(5)) * t242 + (-qJD(5) * t218 + t131 * t184 + t132 * t182) * t172 + (t46 + t56) * t308 + (t45 + t55) * t307 + (-t280 + ((t25 + (-t76 + t282) * t162 - t295) * t162 + (t161 * t237 - t24 + t296) * t161) * qJD(5) + t10) * t240 + (t16 + t21) * t239 + (Icges(5,3) + t219) * t171 + (t17 + t20 + t9) * t241;
t53 = -t232 + t61;
t54 = -t231 + t62;
t192 = (t231 * t30 + t232 * t31) * t180;
t190 = t194 + (Icges(3,3) + Icges(4,3)) * t179;
t188 = ((-t287 * t31 - t30 * t305) * t161 + (-t243 * t30 + t305 * t31) * t162) * t172;
t187 = t192 + t188;
t107 = t180 * t227;
t74 = t180 * t238 - t249;
t73 = -t107 + t208;
t34 = qJD(5) * t200 + qJD(3);
t13 = qJDD(3) + t90 * t81 + t89 * t216 + (-t161 * t48 + t162 * t47) * qJD(5);
t6 = t161 * t316 + t203 * t162;
t5 = t161 * t317 + t202 * t162;
t4 = t203 * t161 - t162 * t316;
t3 = t202 * t161 - t162 * t317;
t1 = [Icges(2,3) * qJDD(1) + t190 + (t330 * (-t125 - t303) + t331 * (-t228 - t304) + (t104 - t110 + t249) * t103) * m(3) + ((qJDD(1) * t229 + g(3)) * t229 + (qJDD(1) * t157 + g(2)) * t157) * m(2) + (t31 * (t236 + t250) + t187 + t326 * (t54 - t303) + t324 * (t53 - t304) + (-t31 - t197 + t249 + t334) * t30) * m(6) + (t325 * (t84 - t303) + t323 * (t83 - t304) + (t199 + t85) * t60 + (-t60 - t86 - t98) * t59) * m(5) + (t74 * (t208 - t258) + t332 * (t106 - t303) + t333 * (t105 - t304) + (-t74 + t328) * t73) * m(4); t190 + (t324 * t53 + t326 * t54 + t187 - t192 + t318) * m(6) + (t328 * t73 + t332 * t106 + t333 * t105 + (t107 - t258) * t74) * m(4) + (-t103 * t110 + t104 * t109 + (-t103 * t180 - t330) * t125 - (t104 * t180 + t331) * t228) * m(3) + (t323 * t83 + t325 * t84 - t59 * t98 + t60 * t97 + t319) * m(5); m(6) * t13 + t306 * qJDD(3) + (-m(6) - t306) * g(1); t194 + (t324 * t61 + t326 * t62 + t188 + t318) * m(6) + ((t172 * t60 + t323) * t226 + (-t172 * t59 - t325) * t111 + t319) * m(5); t171 * (t46 * t161 + t45 * t162) / 0.2e1 + t172 * ((t172 * t46 + t16) * t162 + (-t172 * t45 + t17) * t161) / 0.2e1 - t9 * t269 / 0.2e1 + t162 * (t55 * t171 + t21 * t172 + t24 * t90 + t25 * t89 + (t161 * t6 + t162 * t5) * qJD(5)) / 0.2e1 + t225 * t307 + ((t172 * t25 + t5) * t162 + (-t172 * t24 + t6) * t161) * t239 + t10 * t266 / 0.2e1 + t161 * (t56 * t171 + t20 * t172 + t26 * t90 + t27 * t89 + (t161 * t4 + t3 * t162) * qJD(5)) / 0.2e1 + t224 * t308 + ((t172 * t27 + t3) * t162 + (-t172 * t26 + t4) * t161) * t241 - t172 * ((-t182 * t257 + t184 * t256) * t172 + ((t161 * t291 + t162 * t292) * t184 + (-t161 * t293 - t162 * t294) * t182) * qJD(5)) / 0.2e1 + ((t91 * t254 + t271) * t162 + (t329 + (t310 * t161 + (-t309 - t272) * t162) * qJD(5)) * t161) * t240 + ((-t255 * t272 + t271) * t161 + (-t329 + (t309 * t162 + (-t310 + t91) * t161) * qJD(5)) * t162) * t242 + (t13 * t200 + t34 * ((-t48 - t70) * t161 + (t47 + t69) * t162) + t12 * t99 - t11 * t270 + (t31 * t266 - t30 * t269) * t154 + t320 * t133 - (t270 * t31 - t30 * t99) * t172 - (t34 * (-t161 * t99 - t162 * t270) + t320 * t156) * qJD(5) - g(1) * t156 - g(2) * t99 + g(3) * t270) * m(6);];
tau = t1;
