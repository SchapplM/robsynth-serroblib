% Calculate vector of inverse dynamics joint torques for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:44
% EndTime: 2019-12-05 18:40:51
% DurationCPUTime: 4.25s
% Computational Cost: add. (11392->396), mult. (6350->466), div. (0->0), fcn. (4656->10), ass. (0->239)
t185 = sin(qJ(5));
t187 = cos(qJ(5));
t156 = rSges(6,1) * t185 + rSges(6,2) * t187;
t184 = qJ(1) + qJ(2);
t176 = sin(t184);
t177 = cos(t184);
t183 = qJD(1) + qJD(2);
t181 = t183 ^ 2;
t182 = qJDD(1) + qJDD(2);
t186 = sin(qJ(1));
t188 = cos(qJ(1));
t189 = qJD(1) ^ 2;
t213 = (-qJDD(1) * t186 - t188 * t189) * pkin(1);
t193 = t213 + (-t176 * t182 - t177 * t181) * pkin(2);
t337 = t193 - g(3);
t179 = qJ(3) + t184;
t168 = pkin(9) + t179;
t163 = sin(t168);
t160 = t163 * rSges(5,2);
t164 = cos(t168);
t292 = rSges(5,1) * t164;
t113 = -t160 + t292;
t175 = qJD(3) + t183;
t272 = t163 * t175;
t122 = rSges(5,2) * t272;
t170 = cos(t179);
t301 = pkin(3) * t170;
t233 = -t292 - t301;
t336 = -t122 + (-t113 - t233) * t175;
t231 = rSges(3,1) * t176 + rSges(3,2) * t177;
t111 = t231 * t183;
t287 = pkin(1) * qJD(1);
t249 = t186 * t287;
t105 = t111 + t249;
t300 = pkin(8) * t163;
t115 = pkin(4) * t164 + t300;
t104 = t175 * t115;
t290 = rSges(6,1) * t187;
t158 = -rSges(6,2) * t185 + t290;
t138 = t158 * qJD(5);
t174 = qJDD(3) + t182;
t256 = qJD(5) * t164;
t173 = t175 ^ 2;
t265 = t170 * t173;
t267 = t164 * t187;
t288 = rSges(6,3) * t163;
t219 = -rSges(6,1) * t267 - t288;
t268 = t164 * t185;
t147 = rSges(6,2) * t268;
t331 = t156 * qJD(5);
t244 = t175 * t147 + t331 * t163;
t48 = t175 * t219 + t244;
t243 = -pkin(4) - t290;
t271 = t163 * t185;
t146 = rSges(6,2) * t271;
t260 = t164 * rSges(6,3) + t146;
t299 = t164 * pkin(8);
t169 = sin(t179);
t302 = pkin(3) * t169;
t57 = t163 * t243 + t260 + t299 - t302;
t255 = qJD(5) * t175;
t92 = qJDD(5) * t164 - t163 * t255;
t11 = -pkin(3) * t265 - t138 * t256 - t92 * t156 + (t48 - t104) * t175 + t57 * t174 + t193;
t330 = t11 - g(3);
t125 = pkin(4) * t272;
t305 = pkin(1) * t188;
t306 = pkin(1) * t186;
t234 = -qJDD(1) * t305 + t189 * t306;
t303 = pkin(2) * t177;
t304 = pkin(2) * t176;
t201 = t181 * t304 - t182 * t303 + t234;
t199 = t173 * t302 + t201;
t81 = -t147 - t219;
t235 = -t115 - t81 - t301;
t257 = qJD(5) * t163;
t269 = t164 * t175;
t270 = t163 * t187;
t250 = rSges(6,1) * t270;
t245 = -t331 * t164 - t175 * t250;
t47 = t175 * t260 + t245;
t91 = qJDD(5) * t163 + t164 * t255;
t12 = t138 * t257 + t91 * t156 + (-pkin(8) * t269 + t125 - t47) * t175 + t235 * t174 + t199;
t329 = t12 - g(2);
t229 = -rSges(5,1) * t163 - rSges(5,2) * t164;
t328 = t174 * t229 + t175 * (-rSges(5,1) * t269 + t122) + (-t169 * t174 - t265) * pkin(3) + t337;
t238 = -t113 - t301;
t261 = -rSges(5,1) * t272 - rSges(5,2) * t269;
t327 = t174 * t238 - t175 * t261 - g(2) + t199;
t230 = -rSges(4,1) * t169 - rSges(4,2) * t170;
t264 = t170 * t175;
t266 = t169 * t175;
t90 = -rSges(4,1) * t264 + rSges(4,2) * t266;
t326 = t174 * t230 + t175 * t90 + t337;
t121 = rSges(4,1) * t170 - t169 * rSges(4,2);
t89 = rSges(4,1) * t266 + rSges(4,2) * t264;
t325 = -t121 * t174 + t175 * t89 - g(2) + t201;
t262 = t177 * t183;
t263 = t176 * t183;
t112 = -rSges(3,1) * t262 + rSges(3,2) * t263;
t335 = t112 * t183 - t182 * t231 - g(3) + t213;
t127 = rSges(3,1) * t177 - t176 * rSges(3,2);
t334 = t111 * t183 - t127 * t182 - g(2) + t234;
t178 = Icges(6,4) * t187;
t223 = -Icges(6,2) * t185 + t178;
t323 = Icges(6,1) * t185 + t178;
t258 = t323 + t223;
t280 = Icges(6,4) * t185;
t150 = Icges(6,2) * t187 + t280;
t153 = Icges(6,1) * t187 - t280;
t259 = t150 - t153;
t333 = (t185 * t258 + t187 * t259) * t175;
t253 = pkin(2) * t262;
t332 = t253 + t336;
t114 = t156 * t257;
t70 = t175 * t81;
t313 = -t244 - t104 - t70 + t114;
t78 = Icges(6,4) * t267 - Icges(6,2) * t268 + Icges(6,6) * t163;
t145 = Icges(6,4) * t268;
t80 = Icges(6,1) * t267 + Icges(6,5) * t163 - t145;
t225 = t185 * t78 - t187 * t80;
t324 = t164 * t225;
t110 = t175 * t121;
t248 = t188 * t287;
t211 = t248 + t253;
t72 = -t211 - t110;
t252 = pkin(3) * t266;
t220 = t250 - t260;
t69 = t175 * t220;
t200 = t156 * t256 - t175 * (-pkin(4) * t163 + t299) + t252 + t69;
t254 = pkin(2) * t263;
t194 = t200 + t254;
t30 = t194 + t249;
t31 = t175 * t235 + t114 - t211;
t322 = t163 * t31 + t164 * t30;
t62 = t175 * t238 - t211;
t99 = t175 * t229;
t320 = (-t261 + t99) * t62;
t215 = t223 * t175;
t316 = -Icges(6,6) * t175 + qJD(5) * t150;
t42 = t163 * t316 - t164 * t215;
t216 = t153 * t175;
t314 = -Icges(6,5) * t175 + qJD(5) * t323;
t44 = t163 * t314 - t164 * t216;
t77 = Icges(6,6) * t164 - t163 * t223;
t144 = Icges(6,4) * t271;
t79 = -Icges(6,1) * t270 + Icges(6,5) * t164 + t144;
t45 = t185 * t79 + t187 * t77;
t149 = Icges(6,5) * t187 - Icges(6,6) * t185;
t75 = Icges(6,3) * t164 - t149 * t163;
t319 = qJD(5) * t45 - t175 * t75 + t185 * t42 - t187 * t44;
t41 = -t163 * t215 - t164 * t316;
t43 = -t163 * t216 - t164 * t314;
t46 = t185 * t80 + t187 * t78;
t76 = Icges(6,5) * t267 - Icges(6,6) * t268 + Icges(6,3) * t163;
t318 = qJD(5) * t46 - t175 * t76 + t185 * t41 - t187 * t43;
t148 = Icges(6,5) * t185 + Icges(6,6) * t187;
t317 = -Icges(6,3) * t175 + qJD(5) * t148;
t133 = t223 * qJD(5);
t134 = t153 * qJD(5);
t222 = t150 * t187 + t185 * t323;
t315 = qJD(5) * t222 + t133 * t185 - t134 * t187 - t148 * t175;
t293 = -Icges(6,2) * t267 - t145 + t80;
t295 = t164 * t323 + t78;
t311 = t185 * t293 + t187 * t295;
t294 = Icges(6,2) * t270 + t144 + t79;
t296 = -t163 * t323 + t77;
t310 = -t185 * t294 - t187 * t296;
t309 = t91 / 0.2e1;
t308 = t92 / 0.2e1;
t307 = -rSges(6,3) - pkin(8);
t298 = t164 * t75 + t77 * t271;
t297 = t163 * t75 + t79 * t267;
t284 = t185 * t77;
t283 = t187 * t79;
t221 = t150 * t185 - t187 * t323;
t93 = t148 * t163;
t56 = -t164 * t221 + t93;
t282 = t56 * t175;
t275 = t148 * t164;
t274 = t149 * t175;
t101 = t156 * t163;
t273 = t156 * t164;
t251 = pkin(3) * t264;
t242 = -t257 / 0.2e1;
t241 = t257 / 0.2e1;
t240 = -t256 / 0.2e1;
t239 = t256 / 0.2e1;
t237 = -t76 - t283;
t236 = t125 - t245;
t159 = rSges(2,1) * t188 - t186 * rSges(2,2);
t232 = rSges(2,1) * t186 + rSges(2,2) * t188;
t24 = -t270 * t79 + t298;
t25 = t164 * t76 - t270 * t80 + t78 * t271;
t228 = t163 * t25 + t24 * t164;
t26 = -t268 * t77 + t297;
t27 = t163 * t76 - t324;
t227 = t27 * t163 + t26 * t164;
t226 = -t283 + t284;
t108 = -t121 - t303;
t86 = t160 + t233;
t106 = -t127 * t183 - t248;
t217 = -t251 - t253;
t212 = t249 + t254;
t107 = t230 - t304;
t85 = t229 - t302;
t208 = t236 + t254;
t206 = t90 - t253;
t204 = -t274 * t163 - t164 * t317 + t175 * t225;
t203 = t163 * t317 - t274 * t164 + t175 * t226;
t84 = t86 - t303;
t202 = t149 * qJD(5) + t175 * t221;
t83 = t85 - t304;
t198 = t163 * t220 + t164 * t81;
t197 = t212 + t252;
t196 = t211 + t251;
t58 = t163 * t307 + t164 * t243 + t147 - t301;
t53 = t57 - t304;
t10 = qJD(5) * t227 + t282;
t16 = -qJD(5) * t226 + t185 * t44 + t187 * t42;
t17 = -qJD(5) * t225 + t185 * t43 + t187 * t41;
t20 = t202 * t163 - t164 * t315;
t21 = t163 * t315 + t202 * t164;
t55 = t163 * t221 + t275;
t52 = t55 * t175;
t9 = qJD(5) * t228 + t52;
t192 = (t52 + ((t27 + t298 + t324) * t164 + (-t26 + (t237 - t284) * t164 + t25 + t297) * t163) * qJD(5)) * t242 + (-qJD(5) * t221 + t133 * t187 + t134 * t185) * t175 + (t46 + t56) * t309 + (t45 + t55) * t308 + (-t282 + ((t25 + (-t76 + t284) * t164 - t297) * t164 + (t163 * t237 - t24 + t298) * t163) * qJD(5) + t10) * t240 + (t16 + t21) * t239 + (t17 + t20 + t9) * t241 + (Icges(5,3) + Icges(4,3) + t222) * t174;
t54 = t58 - t303;
t191 = Icges(3,3) * t182 + t192;
t190 = (t31 * (-t146 + t302) - t30 * (-t288 - t300 - t301) + (-t243 * t30 + t307 * t31) * t164) * t175;
t109 = t175 * t230;
t71 = -t109 + t212;
t61 = t197 - t99;
t34 = qJD(5) * t198 + qJD(4);
t13 = qJDD(4) + t92 * t81 + t91 * t220 + (-t163 * t48 + t164 * t47) * qJD(5);
t6 = t163 * t318 + t204 * t164;
t5 = t163 * t319 + t203 * t164;
t4 = t204 * t163 - t164 * t318;
t3 = t203 * t163 - t164 * t319;
t1 = [Icges(2,3) * qJDD(1) + t191 + (t72 * (t212 + t89) - t71 * (t206 - t248) + t325 * (t108 - t305) + t326 * (t107 - t306)) * m(4) + (t334 * (-t127 - t305) + t335 * (-t231 - t306) + (t106 - t112 + t248) * t105) * m(3) + ((qJDD(1) * t232 + g(3)) * t232 + (qJDD(1) * t159 + g(2)) * t159) * m(2) + (t31 * (t208 + t249) + t190 + t329 * (t54 - t305) + t330 * (t53 - t306) + (t211 - t196 - t31 + t313) * t30) * m(6) + (t62 * (t197 - t261) + t327 * (t84 - t305) + t328 * (t83 - t306) + (t248 - t196 - t62 + t332) * t61) * m(5); t191 + (t190 + t329 * t54 + t330 * t53 + (t208 - t194) * t31 + (t253 + t217 + t313) * t30) * m(6) + (t327 * t84 + t328 * t83 + (t217 + t332) * t61 + t320) * m(5) + ((-t206 - t110 - t253) * t71 + t325 * t108 + t326 * t107 + (t89 + t109) * t72) * m(4) + (-t105 * t112 + t106 * t111 + (-t105 * t183 - t334) * t127 - (t106 * t183 + t335) * t231) * m(3); t192 + (t190 + t329 * t58 + t330 * t57 + (t236 - t200) * t31 + (-t251 + t313) * t30) * m(6) + (t327 * t86 + t328 * t85 + (-t251 + t336) * t61 + t320) * m(5) + (-t71 * t90 + t72 * t89 + (t175 * t72 + t326) * t230 + (-t175 * t71 - t325) * t121) * m(4); (t13 - g(1)) * m(6) + (qJDD(4) - g(1)) * m(5); t174 * (t163 * t46 + t164 * t45) / 0.2e1 + t175 * ((t175 * t46 + t16) * t164 + (-t175 * t45 + t17) * t163) / 0.2e1 - t9 * t272 / 0.2e1 + t164 * (t55 * t174 + t21 * t175 + t24 * t92 + t25 * t91 + (t163 * t6 + t164 * t5) * qJD(5)) / 0.2e1 + t228 * t308 + ((t175 * t25 + t5) * t164 + (-t175 * t24 + t6) * t163) * t239 + t10 * t269 / 0.2e1 + t163 * (t56 * t174 + t20 * t175 + t26 * t92 + t27 * t91 + (t163 * t4 + t164 * t3) * qJD(5)) / 0.2e1 + t227 * t309 + ((t175 * t27 + t3) * t164 + (-t175 * t26 + t4) * t163) * t241 - t175 * ((-t185 * t259 + t187 * t258) * t175 + ((t163 * t293 + t164 * t294) * t187 + (-t163 * t295 - t296 * t164) * t185) * qJD(5)) / 0.2e1 + ((t256 * t93 + t274) * t164 + (t333 + (t311 * t163 + (-t310 - t275) * t164) * qJD(5)) * t163) * t240 + ((-t257 * t275 + t274) * t163 + (-t333 + (t310 * t164 + (-t311 + t93) * t163) * qJD(5)) * t164) * t242 + (t13 * t198 + t34 * ((-t48 - t70) * t163 + (t47 + t69) * t164) + t12 * t101 - t11 * t273 + (t269 * t31 - t272 * t30) * t156 + t322 * t138 - (-t101 * t30 + t273 * t31) * t175 - (t34 * (-t101 * t163 - t164 * t273) + t322 * t158) * qJD(5) - g(1) * t158 - g(2) * t101 + g(3) * t273) * m(6);];
tau = t1;
