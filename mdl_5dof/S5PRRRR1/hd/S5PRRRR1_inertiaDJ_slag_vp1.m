% Calculate time derivative of joint inertia matrix for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaDJ_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:52
% EndTime: 2019-12-05 17:03:06
% DurationCPUTime: 6.76s
% Computational Cost: add. (15477->529), mult. (22144->819), div. (0->0), fcn. (21028->8), ass. (0->280)
t219 = sin(qJ(2));
t214 = t219 ^ 2;
t222 = cos(qJ(2));
t215 = t222 ^ 2;
t308 = t214 + t215;
t216 = qJ(3) + qJ(4);
t209 = sin(t216);
t305 = qJD(2) * t219;
t291 = t209 * t305;
t210 = cos(t216);
t213 = qJD(3) + qJD(4);
t318 = t213 * t222;
t294 = t210 * t318;
t229 = -t291 + t294;
t218 = sin(qJ(3));
t221 = cos(qJ(3));
t273 = rSges(4,1) * t221 - rSges(4,2) * t218;
t375 = t222 * rSges(4,3) - t219 * t273;
t344 = Icges(4,4) * t221;
t259 = -Icges(4,2) * t218 + t344;
t171 = Icges(4,6) * t219 + t222 * t259;
t345 = Icges(4,4) * t218;
t262 = Icges(4,1) * t221 - t345;
t173 = Icges(4,5) * t219 + t222 * t262;
t244 = t171 * t218 - t173 * t221;
t374 = t219 * t244;
t342 = Icges(5,4) * t210;
t258 = -Icges(5,2) * t209 + t342;
t159 = Icges(5,6) * t219 + t222 * t258;
t343 = Icges(5,4) * t209;
t261 = Icges(5,1) * t210 - t343;
t161 = Icges(5,5) * t219 + t222 * t261;
t246 = t159 * t209 - t161 * t210;
t373 = t219 * t246;
t170 = -Icges(4,6) * t222 + t219 * t259;
t172 = -Icges(4,5) * t222 + t219 * t262;
t245 = t170 * t218 - t172 * t221;
t371 = t222 * t245;
t158 = -Icges(5,6) * t222 + t219 * t258;
t160 = -Icges(5,5) * t222 + t219 * t261;
t247 = t158 * t209 - t160 * t210;
t370 = t222 * t247;
t302 = qJD(3) * t221;
t304 = qJD(2) * t222;
t369 = -t218 * t304 - t219 * t302;
t255 = Icges(5,5) * t210 - Icges(5,6) * t209;
t156 = -Icges(5,3) * t222 + t219 * t255;
t368 = qJD(2) * t156;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t277 = -qJD(5) * t210 + qJD(2);
t241 = t277 * t220;
t276 = qJD(2) * t210 - qJD(5);
t320 = t213 * t219;
t111 = t219 * t241 + (t209 * t320 - t222 * t276) * t217;
t242 = t277 * t217;
t312 = t220 * t222;
t319 = t213 * t220;
t112 = t276 * t312 + (-t209 * t319 + t242) * t219;
t230 = t209 * t304 + t210 * t320;
t69 = t112 * rSges(6,1) + t111 * rSges(6,2) + rSges(6,3) * t230;
t256 = Icges(4,5) * t221 - Icges(4,6) * t218;
t168 = -Icges(4,3) * t222 + t219 * t256;
t295 = t209 * t318;
t367 = t219 * t276 + t295;
t183 = -Icges(5,2) * t210 - t343;
t184 = -Icges(5,1) * t209 - t342;
t243 = t183 * t209 - t184 * t210;
t366 = qJD(2) * t243 - t255 * t213;
t365 = 2 * m(4);
t364 = 2 * m(5);
t363 = 2 * m(6);
t362 = -t219 / 0.2e1;
t361 = t219 / 0.2e1;
t360 = -t222 / 0.2e1;
t359 = t222 / 0.2e1;
t196 = -rSges(4,1) * t218 - rSges(4,2) * t221;
t358 = m(4) * t196;
t268 = -rSges(6,1) * t220 + rSges(6,2) * t217;
t153 = rSges(6,3) * t210 + t209 * t268;
t357 = m(6) * t153;
t356 = pkin(2) * t218;
t355 = pkin(2) * t221;
t353 = rSges(5,1) * t210;
t351 = rSges(5,3) * t222;
t314 = t219 * t217;
t178 = -t210 * t314 - t312;
t313 = t219 * t220;
t317 = t217 * t222;
t179 = t210 * t313 - t317;
t324 = t209 * t219;
t123 = Icges(6,5) * t179 + Icges(6,6) * t178 + Icges(6,3) * t324;
t125 = Icges(6,4) * t179 + Icges(6,2) * t178 + Icges(6,6) * t324;
t127 = Icges(6,1) * t179 + Icges(6,4) * t178 + Icges(6,5) * t324;
t253 = t125 * t217 - t127 * t220;
t63 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t230;
t65 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t230;
t67 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t230;
t19 = (t213 * t253 + t63) * t210 + (-t123 * t213 + t217 * t65 - t220 * t67 + (t125 * t220 + t127 * t217) * qJD(5)) * t209;
t350 = t19 * t222;
t180 = -t210 * t317 + t313;
t181 = t210 * t312 + t314;
t322 = t209 * t222;
t124 = Icges(6,5) * t181 + Icges(6,6) * t180 + Icges(6,3) * t322;
t126 = Icges(6,4) * t181 + Icges(6,2) * t180 + Icges(6,6) * t322;
t128 = Icges(6,1) * t181 + Icges(6,4) * t180 + Icges(6,5) * t322;
t252 = t126 * t217 - t128 * t220;
t109 = t367 * t217 + t222 * t241;
t110 = -t367 * t220 + t222 * t242;
t62 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t229;
t64 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t229;
t66 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t229;
t20 = (t213 * t252 + t62) * t210 + (-t124 * t213 + t217 * t64 - t220 * t66 + (t126 * t220 + t128 * t217) * qJD(5)) * t209;
t349 = t20 * t219;
t254 = -Icges(6,5) * t220 + Icges(6,6) * t217;
t150 = Icges(6,3) * t210 + t209 * t254;
t340 = Icges(6,4) * t220;
t257 = Icges(6,2) * t217 - t340;
t151 = Icges(6,6) * t210 + t209 * t257;
t341 = Icges(6,4) * t217;
t260 = -Icges(6,1) * t220 + t341;
t152 = Icges(6,5) * t210 + t209 * t260;
t323 = t209 * t220;
t325 = t209 * t217;
t348 = t213 * (-t210 * t150 - t151 * t325 + t152 * t323);
t347 = t219 * rSges(4,3);
t211 = t219 * rSges(5,3);
t321 = t210 * t213;
t104 = t268 * t321 + (-rSges(6,3) * t213 + (rSges(6,1) * t217 + rSges(6,2) * t220) * qJD(5)) * t209;
t333 = t104 * t222;
t332 = t170 * t221;
t331 = t171 * t221;
t330 = t172 * t218;
t329 = t173 * t218;
t328 = t375 * t222;
t236 = t183 * t213;
t237 = t184 * t213;
t190 = t273 * qJD(3);
t327 = t190 * t219;
t326 = t209 * t213;
t316 = t219 * t104;
t272 = -rSges(5,2) * t209 + t353;
t167 = t272 * t213;
t315 = t219 * t167;
t311 = t221 * t222;
t310 = rSges(5,2) * t291 + rSges(5,3) * t304;
t303 = qJD(3) * t218;
t299 = pkin(2) * t303;
t309 = t308 * t299;
t157 = Icges(5,3) * t219 + t222 * t255;
t307 = qJD(2) * t157;
t169 = Icges(4,3) * t219 + t222 * t256;
t306 = qJD(2) * t169;
t301 = qJD(3) * t222;
t300 = qJD(5) * t209;
t207 = pkin(2) * t311;
t298 = pkin(2) * t302;
t58 = t123 * t210 + t209 * t253;
t72 = t150 * t324 + t151 * t178 + t152 * t179;
t297 = -t58 / 0.2e1 - t72 / 0.2e1;
t59 = t124 * t210 + t209 * t252;
t73 = t150 * t322 + t180 * t151 + t181 * t152;
t296 = t59 / 0.2e1 + t73 / 0.2e1;
t130 = t181 * rSges(6,1) + t180 * rSges(6,2) + rSges(6,3) * t322;
t292 = t153 * t305;
t290 = t219 * t303;
t287 = -t305 / 0.2e1;
t286 = t305 / 0.2e1;
t284 = t304 / 0.2e1;
t283 = t153 - t356;
t271 = rSges(5,1) * t209 + rSges(5,2) * t210;
t282 = -t271 - t356;
t116 = -t158 * qJD(2) + t222 * t236;
t281 = -t161 * t213 - t116;
t117 = qJD(2) * t159 + t219 * t236;
t280 = t160 * t213 + t117;
t118 = -t160 * qJD(2) + t222 * t237;
t279 = t159 * t213 - t118;
t119 = qJD(2) * t161 + t219 * t237;
t278 = t158 * t213 - t119;
t275 = t308 * t355;
t269 = -rSges(6,1) * t179 - rSges(6,2) * t178;
t48 = t123 * t324 + t125 * t178 + t127 * t179;
t49 = t124 * t324 + t126 * t178 + t128 * t179;
t37 = t49 * t219 - t222 * t48;
t267 = t219 * t48 + t222 * t49;
t50 = t123 * t322 + t180 * t125 + t181 * t127;
t51 = t124 * t322 + t180 * t126 + t181 * t128;
t38 = t51 * t219 - t222 * t50;
t266 = t219 * t50 + t222 * t51;
t265 = t59 * t219 - t58 * t222;
t264 = t58 * t219 + t59 * t222;
t182 = -Icges(5,5) * t209 - Icges(5,6) * t210;
t235 = t213 * t182;
t114 = t222 * t235 - t368;
t115 = t219 * t235 + t307;
t13 = t109 * t125 + t110 * t127 + t123 * t229 + t180 * t65 + t181 * t67 + t322 * t63;
t14 = t109 * t126 + t110 * t128 + t124 * t229 + t180 * t64 + t181 * t66 + t322 * t62;
t8 = qJD(2) * t266 - t13 * t222 + t14 * t219;
t80 = -t156 * t222 - t219 * t247;
t81 = -t157 * t222 - t373;
t82 = t219 * t156 - t370;
t83 = t219 * t157 - t222 * t246;
t263 = t38 * t304 + t37 * t305 + (-t82 * t304 - t80 * t305) * t222 + (t8 + (t83 * qJD(2) + (t117 * t209 - t119 * t210 + t158 * t321 + t160 * t326 - t368) * t222) * t222 + t81 * t305 + t83 * t304 + ((t82 + t373) * qJD(2) + (-t115 - t279 * t210 + t281 * t209 + (t157 - t247) * qJD(2)) * t222 + t219 * t114) * t219) * t219;
t129 = rSges(6,3) * t324 - t269;
t251 = -t129 * t222 + t130 * t219;
t250 = t219 * t129 + t130 * t222;
t162 = t219 * t272 - t351;
t163 = -rSges(5,2) * t322 + t222 * t353 + t211;
t113 = -t219 * t162 - t163 * t222;
t240 = t153 * t304 + t316;
t239 = -t271 * t304 - t315;
t238 = t271 * t213;
t234 = qJD(3) * (-Icges(4,1) * t218 - t344);
t233 = qJD(3) * (-Icges(4,2) * t221 - t345);
t232 = qJD(3) * (-Icges(4,5) * t218 - Icges(4,6) * t221);
t231 = -t272 - t355;
t68 = t110 * rSges(6,1) + t109 * rSges(6,2) + t229 * rSges(6,3);
t228 = t369 * pkin(2);
t12 = (t222 * t115 + (t81 + t370) * qJD(2)) * t222 + (t80 * qJD(2) + (-t116 * t209 + t118 * t210 - t159 * t321 - t161 * t326 + t307) * t219 + (-t114 + t278 * t210 + t280 * t209 + (-t156 - t246) * qJD(2)) * t222) * t219;
t15 = t111 * t125 + t112 * t127 + t123 * t230 + t178 * t65 + t179 * t67 + t324 * t63;
t16 = t111 * t126 + t112 * t128 + t124 * t230 + t178 * t64 + t179 * t66 + t324 * t62;
t9 = qJD(2) * t267 - t15 * t222 + t16 * t219;
t227 = (-t12 - t9) * t222 + t263;
t39 = -t219 * t69 + t130 * t305 + (-qJD(2) * t129 - t68) * t222;
t101 = t254 * t321 + (-Icges(6,3) * t213 + (Icges(6,5) * t217 + Icges(6,6) * t220) * qJD(5)) * t209;
t102 = t257 * t321 + (-Icges(6,6) * t213 + (Icges(6,2) * t220 + t341) * qJD(5)) * t209;
t226 = -t210 * t101 - t102 * t325 + (t210 * t319 - t217 * t300) * t152 + (-t217 * t321 - t220 * t300) * t151;
t57 = -t219 * (-t219 * t238 + (t222 * t272 + t211) * qJD(2)) + t163 * t305 + (-qJD(2) * t162 + rSges(5,2) * t294 - (-t210 * t305 - t295) * rSges(5,1) - t310) * t222;
t23 = t209 * t267 + t72 * t210;
t24 = t209 * t266 + t73 * t210;
t103 = t260 * t321 + (-Icges(6,5) * t213 + (Icges(6,1) * t217 + t340) * qJD(5)) * t209;
t28 = t180 * t102 + t181 * t103 + t109 * t151 + t110 * t152 - t150 * t291 + (t101 * t209 + t150 * t321) * t222;
t3 = (t213 * t266 + t28) * t210 + (-qJD(2) * t38 + t13 * t219 + t14 * t222 - t213 * t73) * t209;
t29 = t101 * t324 + t178 * t102 + t179 * t103 + t111 * t151 + t112 * t152 + t150 * t230;
t4 = (t213 * t267 + t29) * t210 + (-qJD(2) * t37 + t15 * t219 + t16 * t222 - t213 * t72) * t209;
t225 = t3 * t361 + t9 * t324 / 0.2e1 + t4 * t360 + t210 * (qJD(2) * t264 + t349 - t350) / 0.2e1 + t23 * t286 + t24 * t284 - t265 * t326 / 0.2e1 + t8 * t322 / 0.2e1 + (t219 * t37 + t222 * t38) * t321 / 0.2e1 + (t37 * t284 + t38 * t287) * t209;
t165 = t258 * t213;
t166 = t261 * t213;
t223 = -qJD(2) * t182 + (t166 + t236) * t210 + (-t165 + t237) * t209;
t224 = t350 / 0.2e1 - t349 / 0.2e1 + (-t366 * t219 + t223 * t222) * t361 + (t223 * t219 + t366 * t222) * t360 + (t182 * t222 + t219 * t243) * t286 + (-t219 * t182 + t222 * t243) * t284 + (t209 * t279 + t210 * t281 + t28) * t362 + (t209 * t278 - t210 * t280 + t29) * t359 + (-t158 * t210 - t160 * t209 + t58 + t72) * t287 - (-t159 * t210 - t161 * t209 + t59 + t73) * t304 / 0.2e1;
t205 = t305 * t356;
t203 = pkin(2) * t290;
t175 = t222 * t273 + t347;
t155 = t282 * t222;
t154 = t282 * t219;
t148 = t163 + t207;
t147 = t219 * t231 + t351;
t141 = t283 * t222;
t140 = t283 * t219;
t138 = -rSges(4,1) * t290 + (rSges(4,1) * t311 + t347) * qJD(2) + t369 * rSges(4,2);
t137 = t375 * qJD(2) + t196 * t301;
t132 = t219 * t232 + t306;
t131 = -qJD(2) * t168 + t222 * t232;
t108 = t207 + t130;
t107 = (-rSges(6,3) * t209 - t355) * t219 + t269;
t100 = t228 + t239;
t99 = t271 * t305 + t205 + (-t167 - t298) * t222;
t92 = -t275 + t113;
t91 = t203 + t271 * t320 + (t222 * t231 - t211) * qJD(2);
t90 = (-t353 - t355) * t305 + (-t238 - t299) * t222 + t310;
t89 = t219 * t169 - t222 * t244;
t88 = t219 * t168 - t371;
t87 = -t169 * t222 - t374;
t86 = -t168 * t222 - t219 * t245;
t85 = t129 * t210 - t153 * t324;
t84 = -t210 * t130 + t153 * t322;
t77 = t228 + t240;
t76 = -t292 + t205 + (t104 - t298) * t222;
t75 = t251 * t209;
t74 = -t275 - t250;
t61 = -qJD(2) * t207 + t203 - t69;
t60 = (-t218 * t301 - t221 * t305) * pkin(2) + t68;
t54 = t57 + t309;
t43 = (-t153 * t320 + t69) * t210 + (-t129 * t213 - t240) * t209;
t42 = (t153 * t318 - t68) * t210 + (t130 * t213 - t292 + t333) * t209;
t40 = (t103 * t220 + t150 * t213) * t209 + t226;
t32 = t39 + t309;
t25 = t251 * t321 + (qJD(2) * t250 + t219 * t68 - t222 * t69) * t209;
t1 = [0; m(3) * (-t219 * rSges(3,1) - rSges(3,2) * t222) * qJD(2) + m(4) * t137 + m(5) * t90 + m(6) * t60; (t107 * t61 + t108 * t60) * t363 + (t147 * t91 + t148 * t90) * t364 + (t137 * t175 - t138 * t375) * t365 - t226 + t218 * t233 - t221 * t234 + t259 * t302 + t262 * t303 + t209 * t166 + t210 * t165 - t184 * t321 - t103 * t323 + (-t150 + t183) * t326; (t196 * t304 - t327) * m(4) + m(5) * t100 + m(6) * t77; m(6) * (t107 * t76 + t108 * t77 + t140 * t60 + t141 * t61) + m(5) * (t100 * t148 + t147 * t99 + t154 * t90 + t155 * t91) + ((t175 * t358 + t331 / 0.2e1 + t329 / 0.2e1) * t222 + (t332 / 0.2e1 + t330 / 0.2e1 - t375 * t358) * t219) * qJD(2) + t224 + m(4) * (-t175 * t327 - t190 * t328 + (t137 * t219 - t138 * t222) * t196) + (qJD(3) * t245 - (qJD(2) * t171 + t219 * t233) * t221 - (qJD(2) * t173 + t219 * t234) * t218) * t359 + (qJD(3) * t244 - (-t170 * qJD(2) + t222 * t233) * t221 - (-t172 * qJD(2) + t222 * t234) * t218) * t362 - (-t214 / 0.2e1 - t215 / 0.2e1) * t256 * qJD(3); -t222 * t9 + (t140 * t77 + t141 * t76 + t32 * t74) * t363 - t222 * t12 + (t100 * t154 + t155 * t99 + t54 * t92) * t364 + (t87 * t219 - t222 * t86) * t305 - t222 * ((t222 * t132 + (t87 + t371) * qJD(2)) * t222 + (t86 * qJD(2) + (-t171 * t302 - t173 * t303 + t306) * t219 + (-t131 + (t330 + t332) * qJD(3) - t244 * qJD(2)) * t222) * t219) + (t89 * t219 - t222 * t88) * t304 + t219 * ((t219 * t131 + (t88 + t374) * qJD(2)) * t219 + (t89 * qJD(2) + (t170 * t302 + t172 * t303) * t222 + (-t132 + (-t329 - t331) * qJD(3) + (t169 - t245) * qJD(2)) * t219) * t222) + ((-t175 * t222 + t219 * t375) * (-t222 * t137 - t219 * t138 + (t175 * t219 + t328) * qJD(2)) - t308 * t196 * t190) * t365 + t263; m(5) * t239 + m(6) * t240; m(5) * (-(t147 * t222 + t148 * t219) * t167 - (t219 * t90 + t222 * t91 + (-t147 * t219 + t148 * t222) * qJD(2)) * t271) + t224 + (t219 * t60 + t222 * t61 + (-t107 * t219 + t108 * t222) * qJD(2)) * t357 + m(6) * (t107 * t222 + t108 * t219) * t104; m(6) * (t140 * t316 + t141 * t333 - t250 * t32 + t39 * t74) + (t219 * t77 + t222 * t76 + (t140 * t222 - t141 * t219) * qJD(2)) * t357 + t227 + (-t155 * t167 * t222 + t113 * t54 - t154 * t315 + t57 * t92 - (t100 * t219 + t222 * t99 + (t154 * t222 - t155 * t219) * qJD(2)) * t271) * m(5); (t167 * t271 * t308 + t113 * t57) * t364 + (t104 * t153 * t308 - t250 * t39) * t363 + t227; m(6) * t42; m(6) * (t107 * t43 + t108 * t42 + t60 * t84 + t61 * t85) + (t40 + (t219 * t297 - t222 * t296) * t213) * t210 + (-t348 + (-t20 / 0.2e1 - t28 / 0.2e1) * t222 + (-t19 / 0.2e1 - t29 / 0.2e1) * t219 + (t219 * t296 + t222 * t297) * qJD(2)) * t209; t225 + m(6) * (t140 * t42 + t141 * t43 + t25 * t74 + t32 * t75 + t76 * t85 + t77 * t84); t225 + m(6) * (-t25 * t250 + t75 * t39 + (t219 * t84 + t222 * t85) * t104 + (t219 * t42 + t222 * t43 + (-t219 * t85 + t222 * t84) * qJD(2)) * t153); (t25 * t75 + t42 * t84 + t43 * t85) * t363 + (-t210 * t40 + (t210 * t264 + t219 * t23 + t222 * t24) * t213) * t210 + (t222 * t3 + t219 * t4 - t264 * t326 + (t19 * t219 + t20 * t222 + 0.2e1 * t348) * t210 + (-t210 * t265 - t219 * t24 + t222 * t23) * qJD(2)) * t209;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
