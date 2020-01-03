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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:29
% DurationCPUTime: 4.20s
% Computational Cost: add. (11392->381), mult. (6350->471), div. (0->0), fcn. (4656->10), ass. (0->233)
t201 = sin(qJ(1));
t304 = pkin(1) * qJD(1);
t260 = t201 * t304;
t199 = qJ(1) + qJ(2);
t189 = sin(t199);
t198 = qJD(1) + qJD(2);
t280 = t189 * t198;
t263 = pkin(2) * t280;
t222 = t260 + t263;
t188 = qJD(3) + t198;
t192 = qJ(3) + t199;
t179 = sin(t192);
t180 = cos(t192);
t329 = t179 * rSges(4,1) + t180 * rSges(4,2);
t342 = t329 * t188;
t67 = t222 + t342;
t178 = pkin(9) + t192;
t167 = sin(t178);
t168 = cos(t178);
t108 = t167 * rSges(5,1) + t168 * rSges(5,2);
t169 = pkin(3) * t179;
t328 = t169 + t108;
t190 = cos(t199);
t327 = t189 * rSges(3,1) + t190 * rSges(3,2);
t106 = t327 * t198;
t100 = t260 + t106;
t202 = cos(qJ(5));
t285 = t167 * t202;
t140 = rSges(6,1) * t285;
t200 = sin(qJ(5));
t283 = t168 * t200;
t141 = rSges(6,2) * t283;
t153 = rSges(6,1) * t200 + rSges(6,2) * t202;
t267 = qJD(5) * t168;
t257 = t153 * t267;
t221 = -t257 - t263;
t162 = t167 * pkin(4);
t286 = t167 * t200;
t261 = rSges(6,2) * t286;
t248 = t140 - t261;
t77 = -rSges(6,3) * t168 + t248;
t298 = -pkin(8) * t168 + t162 + t77;
t250 = t298 + t169;
t30 = t188 * t250 - t221 + t260;
t281 = t180 * t188;
t136 = pkin(3) * t281;
t268 = qJD(5) * t167;
t247 = -t153 * t268 + t136;
t279 = t190 * t198;
t152 = pkin(2) * t279;
t203 = cos(qJ(1));
t185 = t203 * t304;
t270 = t152 + t185;
t111 = t168 * pkin(4) + t167 * pkin(8);
t282 = t168 * t202;
t262 = rSges(6,1) * t282;
t78 = rSges(6,3) * t167 - t141 + t262;
t297 = t111 + t78;
t31 = t188 * t297 + t247 + t270;
t330 = -t162 - t169;
t97 = t153 * t167;
t98 = t153 * t168;
t341 = -t30 * t97 - t31 * t98;
t205 = (t31 * (-t140 + t330) - t30 * t141) * t188 + t341 * qJD(5);
t303 = t188 * t31;
t284 = t168 * t188;
t287 = t167 * t188;
t274 = pkin(4) * t284 + pkin(8) * t287;
t275 = rSges(6,3) * t287 + t188 * t262;
t66 = t188 * t78;
t324 = -t188 * t111 + t136 - t247 + t274 + t275 - t66;
t345 = t250 * t303 + t324 * t30 + t205;
t191 = Icges(6,4) * t202;
t233 = -Icges(6,2) * t200 + t191;
t326 = Icges(6,1) * t200 + t191;
t271 = t326 + t233;
t295 = Icges(6,4) * t200;
t146 = Icges(6,2) * t202 + t295;
t149 = Icges(6,1) * t202 - t295;
t272 = t146 - t149;
t343 = (t200 * t271 + t202 * t272) * t188;
t122 = pkin(8) * t284;
t155 = rSges(6,1) * t202 - rSges(6,2) * t200;
t132 = t155 * qJD(5);
t197 = qJDD(1) + qJDD(2);
t187 = qJDD(3) + t197;
t194 = t201 * pkin(1);
t204 = qJD(1) ^ 2;
t290 = pkin(1) * qJDD(1);
t249 = -t194 * t204 + t203 * t290;
t312 = pkin(2) * t197;
t313 = pkin(2) * t198 ^ 2;
t211 = -t189 * t313 + t190 * t312 + t249;
t310 = pkin(3) * t187;
t186 = t188 ^ 2;
t311 = pkin(3) * t186;
t209 = -t179 * t311 + t180 * t310 + t211;
t264 = qJD(5) * t202;
t265 = qJD(5) * t200;
t276 = rSges(6,3) * t284 + t188 * t261;
t47 = rSges(6,2) * t168 * t264 + (t168 * t265 + t188 * t285) * rSges(6,1) - t276;
t266 = qJD(5) * t188;
t88 = -qJDD(5) * t167 - t168 * t266;
t12 = -t132 * t268 + t88 * t153 + (-pkin(4) * t287 + t122 - t47) * t188 + t297 * t187 + t209;
t340 = -g(2) + t12;
t195 = t203 * pkin(1);
t269 = t204 * t195 + t201 * t290;
t251 = t189 * t312 + t190 * t313 + t269;
t228 = t179 * t310 + t180 * t311 + t251;
t48 = -rSges(6,1) * t167 * t265 + (-t167 * t264 - t188 * t283) * rSges(6,2) + t275;
t89 = -qJDD(5) * t168 + t167 * t266;
t11 = t132 * t267 - t89 * t153 + (t48 + t274) * t188 + t298 * t187 + t228;
t339 = -g(3) + t11;
t120 = rSges(5,1) * t284;
t338 = t187 * t108 + t188 * (-rSges(5,2) * t287 + t120) + t228 - g(3);
t157 = t167 * rSges(5,2);
t109 = t168 * rSges(5,1) - t157;
t337 = -t108 * t186 + t187 * t109 - g(2) + t209;
t164 = t179 * rSges(4,2);
t87 = rSges(4,1) * t281 - t164 * t188;
t336 = t187 * t329 + t188 * t87 - g(3) + t251;
t117 = rSges(4,1) * t180 - t164;
t335 = t117 * t187 - t188 * t342 - g(2) + t211;
t107 = rSges(3,1) * t279 - rSges(3,2) * t280;
t334 = t107 * t198 + t197 * t327 - g(3) + t269;
t125 = rSges(3,1) * t190 - t189 * rSges(3,2);
t333 = -t106 * t198 + t125 * t197 - g(2) + t249;
t96 = t188 * t109;
t331 = t120 - t96;
t246 = t328 * t188;
t104 = t188 * t117;
t68 = t270 + t104;
t231 = -t200 * t146 + t202 * t326;
t224 = t233 * t188;
t321 = -Icges(6,6) * t188 + qJD(5) * t146;
t42 = -t321 * t167 + t168 * t224;
t225 = t149 * t188;
t318 = -Icges(6,5) * t188 + qJD(5) * t326;
t44 = -t318 * t167 + t168 * t225;
t73 = -Icges(6,6) * t168 + t167 * t233;
t75 = -Icges(6,5) * t168 + t149 * t167;
t45 = t200 * t75 + t202 * t73;
t145 = Icges(6,5) * t202 - Icges(6,6) * t200;
t71 = -Icges(6,3) * t168 + t145 * t167;
t323 = qJD(5) * t45 - t188 * t71 + t200 * t42 - t202 * t44;
t144 = Icges(6,5) * t200 + Icges(6,6) * t202;
t322 = -Icges(6,3) * t188 + qJD(5) * t144;
t129 = t233 * qJD(5);
t130 = t149 * qJD(5);
t232 = t146 * t202 + t200 * t326;
t320 = qJD(5) * t232 + t129 * t200 - t130 * t202 - t144 * t188;
t74 = Icges(6,4) * t282 - Icges(6,2) * t283 + Icges(6,6) * t167;
t139 = Icges(6,4) * t283;
t76 = Icges(6,1) * t282 + Icges(6,5) * t167 - t139;
t237 = t200 * t76 + t202 * t74;
t41 = t167 * t224 + t321 * t168;
t43 = t167 * t225 + t318 * t168;
t72 = Icges(6,5) * t282 - Icges(6,6) * t283 + Icges(6,3) * t167;
t319 = qJD(5) * t237 - t188 * t72 - t200 * t41 + t202 * t43;
t61 = t222 + t246;
t62 = t136 + t270 + t96;
t208 = (-t61 * t157 - t62 * t328) * t188;
t316 = t62 * t246 + t331 * t61 + t208;
t315 = t88 / 0.2e1;
t314 = t89 / 0.2e1;
t308 = t167 * t326 + t73;
t307 = t168 * t326 + t74;
t306 = -t146 * t167 + t75;
t305 = -Icges(6,2) * t282 - t139 + t76;
t302 = t200 * t73;
t301 = t200 * t74;
t300 = t202 * t75;
t288 = t144 * t167;
t56 = t231 * t168 + t288;
t299 = t56 * t188;
t91 = t144 * t168;
t223 = t145 * t188;
t259 = t122 + t276;
t170 = pkin(3) * t180;
t83 = t170 + t109;
t176 = pkin(2) * t189;
t102 = t176 + t329;
t256 = -t268 / 0.2e1;
t255 = t268 / 0.2e1;
t254 = -t267 / 0.2e1;
t253 = t267 / 0.2e1;
t101 = t125 * t198 + t185;
t80 = t176 + t328;
t177 = pkin(2) * t190;
t81 = t177 + t83;
t103 = t117 + t177;
t156 = rSges(2,1) * t203 - t201 * rSges(2,2);
t154 = rSges(2,1) * t201 + rSges(2,2) * t203;
t63 = t75 * t285;
t24 = -t168 * t71 - t286 * t73 + t63;
t64 = t76 * t285;
t25 = t168 * t72 + t286 * t74 - t64;
t242 = -t167 * t25 - t24 * t168;
t65 = t73 * t283;
t26 = -t167 * t71 - t282 * t75 + t65;
t236 = -t202 * t76 + t301;
t27 = t167 * t72 - t236 * t168;
t241 = -t27 * t167 - t26 * t168;
t240 = -t167 * t31 + t168 * t30;
t239 = t167 * t77 + t168 * t78;
t238 = t300 - t302;
t230 = t152 + t87;
t217 = t200 * t306 + t202 * t308;
t216 = t200 * t305 + t202 * t307;
t214 = -t167 * t223 - t322 * t168 + t188 * t236;
t213 = t322 * t167 - t168 * t223 + t188 * t238;
t212 = -t145 * qJD(5) + t188 * t231;
t58 = t170 + t297;
t57 = (-rSges(6,3) - pkin(8)) * t168 + t248 - t330;
t54 = t177 + t58;
t53 = t176 + t57;
t10 = qJD(5) * t241 - t299;
t16 = qJD(5) * t238 + t200 * t44 + t202 * t42;
t17 = qJD(5) * t236 + t200 * t43 + t202 * t41;
t20 = t212 * t167 + t320 * t168;
t21 = -t320 * t167 + t212 * t168;
t55 = t167 * t231 - t91;
t52 = t55 * t188;
t9 = qJD(5) * t242 + t52;
t207 = (t52 + ((t26 + t64 - t65 + (t71 - t301) * t167) * t167 + (-t63 - t27 + (-t236 + t71) * t168 + (t300 + t302) * t167) * t168) * qJD(5)) * t255 - t237 * t315 - t88 * t56 / 0.2e1 + (qJD(5) * t231 + t129 * t202 + t130 * t200) * t188 + (t45 + t55) * t314 + (t299 + ((-t25 + t65 + (t72 - t300) * t168) * t168 + (t24 - t63 + (t72 + t302) * t167) * t167) * qJD(5) + t10) * t253 + (t17 + t20 + t9) * t256 + (t16 + t21) * t254 + (Icges(5,3) + Icges(4,3) + t232) * t187;
t206 = Icges(3,3) * t197 + t207;
t34 = qJD(5) * t239 + qJD(4);
t13 = -t77 * t88 - t78 * t89 + qJDD(4) + (t167 * t48 - t168 * t47) * qJD(5);
t6 = t319 * t167 + t214 * t168;
t5 = -t323 * t167 + t213 * t168;
t4 = t214 * t167 - t319 * t168;
t3 = t213 * t167 + t323 * t168;
t1 = [Icges(2,3) * qJDD(1) + t206 + (t333 * (t125 + t195) + t334 * (t194 + t327) + (-t101 + t107 + t185) * t100) * m(3) + ((t154 ^ 2 + t156 ^ 2) * qJDD(1) - g(2) * t156 - g(3) * t154) * m(2) + (t31 * (-t222 + t259) + t205 + t340 * (t195 + t54) + t339 * (t194 + t53) + (t31 + t324) * t30) * m(6) + (-t62 * t222 + t208 + t337 * (t195 + t81) + t338 * (t194 + t80) + (t62 + t331) * t61) * m(5) + (t335 * (t103 + t195) + t336 * (t194 + t102) + (-t68 + t185 + t230) * t67) * m(4); t206 + (t340 * t54 + t339 * t53 + (t259 - t263 - t221) * t31 + t345) * m(6) + (t337 * t81 + t338 * t80 + t316) * m(5) + ((t230 - t152 - t104) * t67 + t335 * t103 + t336 * t102) * m(4) + (t100 * t107 - t101 * t106 + (-t100 * t198 + t333) * t125 + (t101 * t198 + t334) * t327) * m(3); t207 + (t340 * t58 + t339 * t57 + (t257 + t259) * t31 + t345) * m(6) + (t328 * t338 + t337 * t83 + t316) * m(5) + (t67 * t87 - t68 * t342 + (-t188 * t67 + t335) * t117 + (t188 * t68 + t336) * t329) * m(4); (t13 - g(1)) * m(6) + (qJDD(4) - g(1)) * m(5); t187 * (t167 * t237 - t168 * t45) / 0.2e1 + t188 * ((t188 * t237 - t16) * t168 + (t188 * t45 - t17) * t167) / 0.2e1 + t9 * t287 / 0.2e1 - t168 * (t55 * t187 + t21 * t188 + t24 * t89 + t25 * t88 + (-t167 * t6 - t168 * t5) * qJD(5)) / 0.2e1 + t242 * t314 + ((-t188 * t25 - t5) * t168 + (t188 * t24 - t6) * t167) * t254 - t10 * t284 / 0.2e1 - t167 * (-t56 * t187 + t20 * t188 + t26 * t89 + t27 * t88 + (-t167 * t4 - t168 * t3) * qJD(5)) / 0.2e1 + t241 * t315 + ((-t188 * t27 - t3) * t168 + (t188 * t26 - t4) * t167) * t256 - t188 * ((-t200 * t272 + t202 * t271) * t188 + ((t167 * t305 - t168 * t306) * t202 + (-t167 * t307 + t168 * t308) * t200) * qJD(5)) / 0.2e1 + ((-t267 * t288 - t223) * t168 + (-t343 + (-t216 * t167 + (t91 + t217) * t168) * qJD(5)) * t167) * t253 + ((t268 * t91 - t223) * t167 + (t343 + (-t217 * t168 + (-t288 + t216) * t167) * qJD(5)) * t168) * t255 + (t13 * t239 + t34 * ((t188 * t77 - t47) * t168 + (t48 - t66) * t167) + t240 * t132 + ((t11 - t303) * t168 + (-t188 * t30 - t12) * t167) * t153 - t341 * t188 - (t34 * (-t167 * t97 - t168 * t98) + t240 * t155) * qJD(5) - g(1) * t155 + g(2) * t97 - g(3) * t98) * m(6);];
tau = t1;
