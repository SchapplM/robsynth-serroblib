% Calculate joint inertia matrix for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:10:11
% EndTime: 2018-12-10 18:11:06
% DurationCPUTime: 19.11s
% Computational Cost: add. (400374->695), mult. (401676->945), div. (0->0), fcn. (404827->30), ass. (0->321)
t308 = pkin(7) + pkin(14);
t299 = sin(t308) / 0.2e1;
t309 = pkin(7) - pkin(14);
t303 = sin(t309);
t282 = t299 + t303 / 0.2e1;
t300 = cos(t309) / 0.2e1;
t304 = cos(t308);
t284 = t300 + t304 / 0.2e1;
t310 = pkin(6) + qJ(2);
t301 = sin(t310) / 0.2e1;
t311 = pkin(6) - qJ(2);
t305 = sin(t311);
t287 = t301 + t305 / 0.2e1;
t302 = cos(t311) / 0.2e1;
t306 = cos(t310);
t291 = t302 - t306 / 0.2e1;
t312 = sin(pkin(14));
t319 = cos(pkin(6));
t248 = t282 * t319 + t284 * t287 - t291 * t312;
t283 = t299 - t303 / 0.2e1;
t285 = t300 - t304 / 0.2e1;
t316 = cos(pkin(14));
t249 = t283 * t287 + t285 * t319 + t291 * t316;
t314 = sin(pkin(7));
t318 = cos(pkin(7));
t277 = -t287 * t314 + t318 * t319;
t214 = Icges(4,5) * t249 + Icges(4,6) * t248 + Icges(4,3) * t277;
t215 = Icges(4,4) * t249 + Icges(4,2) * t248 + Icges(4,6) * t277;
t216 = Icges(4,1) * t249 + Icges(4,4) * t248 + Icges(4,5) * t277;
t269 = Icges(3,5) * t291 + Icges(3,6) * t287 + Icges(3,3) * t319;
t270 = Icges(3,4) * t291 + Icges(3,2) * t287 + Icges(3,6) * t319;
t271 = Icges(3,1) * t291 + Icges(3,4) * t287 + Icges(3,5) * t319;
t404 = t277 * t214 + t248 * t215 + t249 * t216 + t319 * t269 + t287 * t270 + t291 * t271;
t290 = t302 + t306 / 0.2e1;
t323 = sin(qJ(2));
t324 = sin(qJ(1));
t328 = cos(qJ(1));
t278 = t290 * t328 - t324 * t323;
t315 = sin(pkin(6));
t385 = t315 * t328;
t266 = -t278 * t314 - t318 * t385;
t280 = -t324 * t290 - t323 * t328;
t386 = t315 * t324;
t267 = -t280 * t314 + t318 * t386;
t288 = t301 - t305 / 0.2e1;
t327 = cos(qJ(2));
t281 = -t324 * t288 + t327 * t328;
t238 = t280 * t284 - t281 * t312 + t282 * t386;
t313 = sin(pkin(8));
t317 = cos(pkin(8));
t224 = -t238 * t313 + t267 * t317;
t279 = t288 * t328 + t324 * t327;
t236 = t278 * t284 - t279 * t312 - t282 * t385;
t223 = -t236 * t313 + t266 * t317;
t403 = t404 * t319;
t237 = t278 * t283 + t279 * t316 - t285 * t385;
t376 = pkin(8) + qJ(4);
t345 = sin(t376) / 0.2e1;
t377 = pkin(8) - qJ(4);
t350 = sin(t377);
t286 = t345 - t350 / 0.2e1;
t346 = cos(t377) / 0.2e1;
t351 = cos(t376);
t289 = t346 - t351 / 0.2e1;
t326 = cos(qJ(4));
t188 = t236 * t286 + t237 * t326 + t266 * t289;
t321 = sin(qJ(5));
t396 = cos(qJ(5));
t173 = t188 * t321 - t223 * t396;
t239 = t280 * t283 + t281 * t316 + t285 * t386;
t190 = t238 * t286 + t239 * t326 + t267 * t289;
t175 = t190 * t321 - t224 * t396;
t210 = t248 * t286 + t249 * t326 + t277 * t289;
t231 = -t248 * t313 + t277 * t317;
t180 = t210 * t321 - t231 * t396;
t174 = t188 * t396 + t223 * t321;
t322 = sin(qJ(4));
t333 = t345 + t350 / 0.2e1;
t334 = t346 + t351 / 0.2e1;
t187 = -t236 * t334 + t237 * t322 - t266 * t333;
t320 = sin(qJ(6));
t325 = cos(qJ(6));
t137 = -t174 * t320 + t187 * t325;
t138 = t174 * t325 + t187 * t320;
t92 = Icges(7,5) * t138 + Icges(7,6) * t137 + Icges(7,3) * t173;
t94 = Icges(7,4) * t138 + Icges(7,2) * t137 + Icges(7,6) * t173;
t96 = Icges(7,1) * t138 + Icges(7,4) * t137 + Icges(7,5) * t173;
t28 = t137 * t94 + t138 * t96 + t173 * t92;
t176 = t190 * t396 + t224 * t321;
t189 = -t238 * t334 + t239 * t322 - t267 * t333;
t139 = -t176 * t320 + t189 * t325;
t140 = t176 * t325 + t189 * t320;
t93 = Icges(7,5) * t140 + Icges(7,6) * t139 + Icges(7,3) * t175;
t95 = Icges(7,4) * t140 + Icges(7,2) * t139 + Icges(7,6) * t175;
t97 = Icges(7,1) * t140 + Icges(7,4) * t139 + Icges(7,5) * t175;
t29 = t137 * t95 + t138 * t97 + t173 * t93;
t181 = t210 * t396 + t231 * t321;
t209 = -t248 * t334 + t249 * t322 - t277 * t333;
t163 = -t181 * t320 + t209 * t325;
t164 = t181 * t325 + t209 * t320;
t106 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t180;
t107 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t180;
t108 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t180;
t41 = t106 * t173 + t107 * t137 + t108 * t138;
t1 = t173 * t28 + t175 * t29 + t180 * t41;
t402 = t1 / 0.2e1;
t30 = t139 * t94 + t140 * t96 + t175 * t92;
t31 = t139 * t95 + t140 * t97 + t175 * t93;
t42 = t106 * t175 + t107 * t139 + t108 * t140;
t2 = t173 * t30 + t175 * t31 + t180 * t42;
t401 = t2 / 0.2e1;
t34 = t163 * t94 + t164 * t96 + t180 * t92;
t35 = t163 * t95 + t164 * t97 + t180 * t93;
t49 = t180 * t106 + t163 * t107 + t164 * t108;
t43 = t49 * t180;
t9 = t34 * t173 + t35 * t175 + t43;
t400 = t9 / 0.2e1;
t399 = t173 / 0.2e1;
t398 = t175 / 0.2e1;
t397 = t180 / 0.2e1;
t395 = pkin(5) * t174;
t344 = -rSges(7,1) * t138 - rSges(7,2) * t137;
t98 = rSges(7,3) * t173 - t344;
t394 = pkin(13) * t173 + t395 + t98;
t99 = t140 * rSges(7,1) + t139 * rSges(7,2) + t175 * rSges(7,3);
t393 = t176 * pkin(5) + pkin(13) * t175 + t99;
t109 = rSges(7,1) * t164 + rSges(7,2) * t163 + rSges(7,3) * t180;
t384 = pkin(5) * t181 + pkin(13) * t180 + t109;
t195 = t239 * pkin(3) + pkin(11) * t224;
t247 = t281 * pkin(2) + qJ(3) * t267;
t245 = t319 * t247;
t383 = t319 * t195 + t245;
t194 = pkin(3) * t237 + pkin(11) * t223;
t246 = pkin(2) * t279 + qJ(3) * t266;
t382 = -t194 - t246;
t263 = pkin(2) * t291 + qJ(3) * t277;
t381 = -pkin(3) * t249 - pkin(11) * t231 - t263;
t380 = t246 * t386 + t247 * t385;
t379 = t328 * pkin(1) + pkin(10) * t386;
t112 = Icges(6,5) * t174 - Icges(6,6) * t173 + Icges(6,3) * t187;
t114 = Icges(6,4) * t174 - Icges(6,2) * t173 + Icges(6,6) * t187;
t116 = Icges(6,1) * t174 - Icges(6,4) * t173 + Icges(6,5) * t187;
t51 = t112 * t187 - t114 * t173 + t116 * t174;
t113 = Icges(6,5) * t176 - Icges(6,6) * t175 + Icges(6,3) * t189;
t115 = Icges(6,4) * t176 - Icges(6,2) * t175 + Icges(6,6) * t189;
t117 = Icges(6,1) * t176 - Icges(6,4) * t175 + Icges(6,5) * t189;
t52 = t113 * t187 - t115 * t173 + t117 * t174;
t127 = Icges(6,5) * t181 - Icges(6,6) * t180 + Icges(6,3) * t209;
t128 = Icges(6,4) * t181 - Icges(6,2) * t180 + Icges(6,6) * t209;
t129 = Icges(6,1) * t181 - Icges(6,4) * t180 + Icges(6,5) * t209;
t61 = t127 * t187 - t128 * t173 + t129 * t174;
t13 = t187 * t51 + t189 * t52 + t209 * t61;
t3 = t187 * t28 + t189 * t29 + t209 * t41;
t375 = t3 / 0.2e1 + t13 / 0.2e1;
t53 = t112 * t189 - t114 * t175 + t116 * t176;
t54 = t113 * t189 - t115 * t175 + t117 * t176;
t62 = t127 * t189 - t128 * t175 + t129 * t176;
t14 = t187 * t53 + t189 * t54 + t209 * t62;
t4 = t187 * t30 + t189 * t31 + t209 * t42;
t374 = t4 / 0.2e1 + t14 / 0.2e1;
t15 = t223 * t51 + t224 * t52 + t231 * t61;
t5 = t223 * t28 + t224 * t29 + t231 * t41;
t373 = t5 / 0.2e1 + t15 / 0.2e1;
t16 = t223 * t53 + t224 * t54 + t231 * t62;
t6 = t223 * t30 + t224 * t31 + t231 * t42;
t372 = t6 / 0.2e1 + t16 / 0.2e1;
t17 = t61 * t319 + (t324 * t52 - t328 * t51) * t315;
t7 = t41 * t319 + (-t28 * t328 + t29 * t324) * t315;
t371 = t7 / 0.2e1 + t17 / 0.2e1;
t18 = t62 * t319 + (t324 * t54 - t328 * t53) * t315;
t8 = t42 * t319 + (-t30 * t328 + t31 * t324) * t315;
t370 = t8 / 0.2e1 + t18 / 0.2e1;
t44 = t49 * t209;
t10 = t34 * t187 + t35 * t189 + t44;
t55 = t112 * t209 - t114 * t180 + t116 * t181;
t56 = t113 * t209 - t115 * t180 + t117 * t181;
t72 = t209 * t127 - t180 * t128 + t181 * t129;
t65 = t72 * t209;
t19 = t55 * t187 + t56 * t189 + t65;
t369 = t10 / 0.2e1 + t19 / 0.2e1;
t47 = t49 * t231;
t11 = t34 * t223 + t35 * t224 + t47;
t68 = t72 * t231;
t20 = t55 * t223 + t56 * t224 + t68;
t368 = t11 / 0.2e1 + t20 / 0.2e1;
t48 = t49 * t319;
t12 = t48 + (t35 * t324 - t34 * t328) * t315;
t71 = t72 * t319;
t21 = t71 + (t56 * t324 - t55 * t328) * t315;
t367 = t12 / 0.2e1 + t21 / 0.2e1;
t366 = t35 / 0.2e1 + t42 / 0.2e1;
t365 = t41 / 0.2e1 + t34 / 0.2e1;
t165 = Icges(5,5) * t210 - Icges(5,6) * t209 + Icges(5,3) * t231;
t166 = Icges(5,4) * t210 - Icges(5,2) * t209 + Icges(5,6) * t231;
t167 = Icges(5,1) * t210 - Icges(5,4) * t209 + Icges(5,5) * t231;
t87 = t231 * t165 - t209 * t166 + t210 * t167;
t162 = t190 * pkin(4) + pkin(12) * t189;
t362 = t319 * t162 + t383;
t161 = pkin(4) * t188 + t187 * pkin(12);
t361 = -t161 + t382;
t119 = t176 * rSges(6,1) - t175 * rSges(6,2) + t189 * rSges(6,3);
t177 = pkin(4) * t210 + pkin(12) * t209;
t360 = -t177 + t381;
t148 = t190 * rSges(5,1) - t189 * rSges(5,2) + t224 * rSges(5,3);
t205 = t239 * rSges(4,1) + t238 * rSges(4,2) + t267 * rSges(4,3);
t257 = t281 * rSges(3,1) + t280 * rSges(3,2) + rSges(3,3) * t386;
t357 = -t324 * pkin(1) + pkin(10) * t385;
t349 = t315 * (-rSges(4,1) * t249 - rSges(4,2) * t248 - rSges(4,3) * t277 - t263);
t348 = t194 * t386 + t195 * t385 + t380;
t168 = rSges(5,1) * t210 - rSges(5,2) * t209 + rSges(5,3) * t231;
t347 = t315 * (-t168 + t381);
t130 = rSges(6,1) * t181 - rSges(6,2) * t180 + rSges(6,3) * t209;
t343 = t315 * (-t130 + t360);
t342 = t56 / 0.2e1 + t62 / 0.2e1 + t366;
t341 = t55 / 0.2e1 + t61 / 0.2e1 + t365;
t340 = t161 * t386 + t162 * t385 + t348;
t339 = t315 * (t360 - t384);
t204 = rSges(4,1) * t237 + rSges(4,2) * t236 + rSges(4,3) * t266;
t147 = rSges(5,1) * t188 - rSges(5,2) * t187 + rSges(5,3) * t223;
t118 = rSges(6,1) * t174 - rSges(6,2) * t173 + rSges(6,3) * t187;
t338 = -t246 + t357;
t337 = t247 + t379;
t256 = t279 * rSges(3,1) + t278 * rSges(3,2) - rSges(3,3) * t385;
t142 = Icges(5,5) * t190 - Icges(5,6) * t189 + Icges(5,3) * t224;
t144 = Icges(5,4) * t190 - Icges(5,2) * t189 + Icges(5,6) * t224;
t146 = Icges(5,1) * t190 - Icges(5,4) * t189 + Icges(5,5) * t224;
t79 = t142 * t231 - t144 * t209 + t146 * t210;
t83 = t165 * t224 - t166 * t189 + t167 * t190;
t336 = t83 / 0.2e1 + t79 / 0.2e1 + t342;
t141 = Icges(5,5) * t188 - Icges(5,6) * t187 + Icges(5,3) * t223;
t143 = Icges(5,4) * t188 - Icges(5,2) * t187 + Icges(5,6) * t223;
t145 = Icges(5,1) * t188 - Icges(5,4) * t187 + Icges(5,5) * t223;
t78 = t141 * t231 - t143 * t209 + t145 * t210;
t82 = t165 * t223 - t166 * t187 + t167 * t188;
t335 = t78 / 0.2e1 + t82 / 0.2e1 + t341;
t332 = t195 + t337;
t331 = -t194 + t338;
t330 = t162 + t332;
t329 = -t161 + t331;
t295 = rSges(2,1) * t328 - t324 * rSges(2,2);
t294 = -t324 * rSges(2,1) - rSges(2,2) * t328;
t272 = rSges(3,1) * t291 + rSges(3,2) * t287 + rSges(3,3) * t319;
t255 = Icges(3,1) * t281 + Icges(3,4) * t280 + Icges(3,5) * t386;
t254 = Icges(3,1) * t279 + Icges(3,4) * t278 - Icges(3,5) * t385;
t253 = Icges(3,4) * t281 + Icges(3,2) * t280 + Icges(3,6) * t386;
t252 = Icges(3,4) * t279 + Icges(3,2) * t278 - Icges(3,6) * t385;
t251 = Icges(3,5) * t281 + Icges(3,6) * t280 + Icges(3,3) * t386;
t250 = Icges(3,5) * t279 + Icges(3,6) * t278 - Icges(3,3) * t385;
t244 = t257 + t379;
t243 = -t256 + t357;
t226 = -t319 * t256 - t272 * t385;
t225 = t257 * t319 - t272 * t386;
t222 = (t256 * t324 + t257 * t328) * t315;
t212 = t269 * t386 + t270 * t280 + t271 * t281;
t211 = -t269 * t385 + t278 * t270 + t279 * t271;
t207 = t251 * t319 + t253 * t287 + t255 * t291;
t206 = t250 * t319 + t252 * t287 + t254 * t291;
t203 = Icges(4,1) * t239 + Icges(4,4) * t238 + Icges(4,5) * t267;
t202 = Icges(4,1) * t237 + Icges(4,4) * t236 + Icges(4,5) * t266;
t201 = Icges(4,4) * t239 + Icges(4,2) * t238 + Icges(4,6) * t267;
t200 = Icges(4,4) * t237 + Icges(4,2) * t236 + Icges(4,6) * t266;
t199 = Icges(4,5) * t239 + Icges(4,6) * t238 + Icges(4,3) * t267;
t198 = Icges(4,5) * t237 + Icges(4,6) * t236 + Icges(4,3) * t266;
t179 = t337 + t205;
t178 = -t204 + t338;
t157 = t223 * t177;
t156 = (-t204 - t246) * t319 + t328 * t349;
t155 = t205 * t319 + t324 * t349 + t245;
t152 = t231 * t162;
t151 = (t204 * t324 + t205 * t328) * t315 + t380;
t150 = t224 * t161;
t132 = t214 * t267 + t215 * t238 + t216 * t239;
t131 = t214 * t266 + t215 * t236 + t216 * t237;
t124 = t199 * t277 + t201 * t248 + t203 * t249;
t123 = t198 * t277 + t200 * t248 + t202 * t249;
t122 = t332 + t148;
t121 = -t147 + t331;
t105 = t148 * t231 - t168 * t224;
t104 = -t147 * t231 + t168 * t223;
t103 = t147 * t224 - t148 * t223;
t101 = (-t147 + t382) * t319 + t328 * t347;
t100 = t148 * t319 + t324 * t347 + t383;
t89 = t330 + t119;
t88 = -t118 + t329;
t86 = t87 * t319;
t85 = (t147 * t324 + t148 * t328) * t315 + t348;
t84 = t87 * t231;
t81 = t119 * t209 - t130 * t189;
t80 = -t118 * t209 + t130 * t187;
t77 = t142 * t224 - t144 * t189 + t146 * t190;
t76 = t141 * t224 - t143 * t189 + t145 * t190;
t75 = t142 * t223 - t144 * t187 + t146 * t188;
t74 = t141 * t223 - t143 * t187 + t145 * t188;
t73 = t118 * t189 - t119 * t187;
t70 = t119 * t231 + t152 + (-t130 - t177) * t224;
t69 = t130 * t223 + t157 + (-t118 - t161) * t231;
t67 = (-t118 + t361) * t319 + t328 * t343;
t66 = t119 * t319 + t324 * t343 + t362;
t64 = t330 + t393;
t63 = (-rSges(7,3) - pkin(13)) * t173 - t395 + t329 + t344;
t60 = t118 * t224 + t150 + (-t119 - t162) * t223;
t59 = -t109 * t175 + t180 * t99;
t58 = t109 * t173 - t180 * t98;
t57 = (t118 * t324 + t119 * t328) * t315 + t340;
t50 = -t173 * t99 + t175 * t98;
t46 = -t189 * t384 + t209 * t393;
t45 = t187 * t384 - t209 * t394;
t40 = (t361 - t394) * t319 + t328 * t339;
t39 = t319 * t393 + t324 * t339 + t362;
t38 = t152 + t393 * t231 + (-t177 - t384) * t224;
t37 = t157 + t384 * t223 + (-t161 - t394) * t231;
t36 = -t187 * t393 + t189 * t394;
t33 = (t324 * t394 + t328 * t393) * t315 + t340;
t32 = t86 + (t79 * t324 - t78 * t328) * t315;
t27 = t150 + t394 * t224 + (-t162 - t393) * t223;
t26 = t78 * t223 + t79 * t224 + t84;
t25 = t83 * t319 + (t324 * t77 - t328 * t76) * t315;
t24 = t82 * t319 + (t324 * t75 - t328 * t74) * t315;
t23 = t223 * t76 + t224 * t77 + t231 * t83;
t22 = t223 * t74 + t224 * t75 + t231 * t82;
t90 = [(t63 ^ 2 + t64 ^ 2) * m(7) + Icges(2,3) + (t88 ^ 2 + t89 ^ 2) * m(6) + (t121 ^ 2 + t122 ^ 2) * m(5) + (t178 ^ 2 + t179 ^ 2) * m(4) + (t243 ^ 2 + t244 ^ 2) * m(3) + m(2) * (t294 ^ 2 + t295 ^ 2) + t87 + t72 + t49 + t404; t71 + (t100 * t122 + t101 * t121) * m(5) + t86 + (t66 * t89 + t67 * t88) * m(6) + (t39 * t64 + t40 * t63) * m(7) + t48 + (t155 * t179 + t156 * t178) * m(4) + (t225 * t244 + t226 * t243) * m(3) + ((-t123 / 0.2e1 - t206 / 0.2e1 - t131 / 0.2e1 - t211 / 0.2e1 - t335) * t328 + (t124 / 0.2e1 + t207 / 0.2e1 + t132 / 0.2e1 + t212 / 0.2e1 + t336) * t324) * t315 + t403; (t33 ^ 2 + t39 ^ 2 + t40 ^ 2) * m(7) + (t57 ^ 2 + t66 ^ 2 + t67 ^ 2) * m(6) + (t100 ^ 2 + t101 ^ 2 + t85 ^ 2) * m(5) + (t151 ^ 2 + t155 ^ 2 + t156 ^ 2) * m(4) + (t222 ^ 2 + t225 ^ 2 + t226 ^ 2) * m(3) + (t8 + t18 + t25 + ((t199 * t267 + t201 * t238 + t203 * t239) * t324 - (t198 * t267 + t200 * t238 + t202 * t239) * t328) * t315 + (t251 * t386 + t253 * t280 + t255 * t281) * t386) * t386 + (-t7 - t17 - t24 - ((t199 * t266 + t201 * t236 + t203 * t237) * t324 - (t198 * t266 + t200 * t236 + t202 * t237) * t328) * t315 + (-t250 * t385 + t278 * t252 + t279 * t254) * t385 + (-t250 * t386 + t251 * t385 - t252 * t280 - t278 * t253 - t254 * t281 - t279 * t255) * t386) * t385 + (t12 + t21 + t32 + (t212 + t132) * t386 + (-t131 - t211) * t385 + ((-t123 - t206) * t328 + (t124 + t207) * t324) * t315 + t403) * t319; (m(4) * t178 + m(5) * t121 + m(6) * t88 + m(7) * t63) * t267 + (m(4) * t179 + m(5) * t122 + m(6) * t89 + m(7) * t64) * t266; (m(4) * t151 + m(5) * t85 + m(6) * t57 + m(7) * t33) * t277 + (m(4) * t156 + m(5) * t101 + m(6) * t67 + m(7) * t40) * t267 + (m(4) * t155 + m(5) * t100 + m(6) * t66 + m(7) * t39) * t266; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1 + m(5) / 0.2e1 + m(4) / 0.2e1) * (t266 ^ 2 + t267 ^ 2 + t277 ^ 2); t47 + (t37 * t63 + t38 * t64) * m(7) + (t69 * t88 + t70 * t89) * m(6) + t68 + t84 + (t104 * t121 + t105 * t122) * m(5) + t336 * t224 + t335 * t223; (t27 * t33 + t37 * t40 + t38 * t39) * m(7) + (t60 * t57 + t66 * t70 + t67 * t69) * m(6) + (t100 * t105 + t101 * t104 + t103 * t85) * m(5) + (t26 / 0.2e1 + t368) * t319 + (t32 / 0.2e1 + t367) * t231 + (t25 / 0.2e1 + t370) * t224 + (t24 / 0.2e1 + t371) * t223 + ((-t22 / 0.2e1 - t373) * t328 + (t23 / 0.2e1 + t372) * t324) * t315; (m(5) * t103 + m(6) * t60 + m(7) * t27) * t277 + (m(5) * t104 + m(6) * t69 + m(7) * t37) * t267 + (m(5) * t105 + m(6) * t70 + m(7) * t38) * t266; (t27 ^ 2 + t37 ^ 2 + t38 ^ 2) * m(7) + (t60 ^ 2 + t69 ^ 2 + t70 ^ 2) * m(6) + (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) * m(5) + (t11 + t20 + t26) * t231 + (t6 + t16 + t23) * t224 + (t5 + t15 + t22) * t223; t44 + (t45 * t63 + t46 * t64) * m(7) + (t80 * t88 + t81 * t89) * m(6) + t65 + t342 * t189 + t341 * t187; (t33 * t36 + t39 * t46 + t40 * t45) * m(7) + (t57 * t73 + t66 * t81 + t67 * t80) * m(6) + t369 * t319 + t367 * t209 + t370 * t189 + t371 * t187 + (t324 * t374 - t328 * t375) * t315; (m(6) * t73 + m(7) * t36) * t277 + (m(6) * t80 + m(7) * t45) * t267 + (m(6) * t81 + m(7) * t46) * t266; (t27 * t36 + t37 * t45 + t38 * t46) * m(7) + (t60 * t73 + t69 * t80 + t70 * t81) * m(6) + t369 * t231 + t374 * t224 + t375 * t223 + t368 * t209 + t372 * t189 + t373 * t187; (t36 ^ 2 + t45 ^ 2 + t46 ^ 2) * m(7) + (t73 ^ 2 + t80 ^ 2 + t81 ^ 2) * m(6) + (t10 + t19) * t209 + (t4 + t14) * t189 + (t3 + t13) * t187; t43 + (t58 * t63 + t59 * t64) * m(7) + t366 * t175 + t365 * t173; t8 * t398 + t319 * t400 + (t33 * t50 + t39 * t59 + t40 * t58) * m(7) + t12 * t397 + t7 * t399 + (-t328 * t1 / 0.2e1 + t324 * t401) * t315; (t266 * t59 + t267 * t58 + t277 * t50) * m(7); (t27 * t50 + t37 * t58 + t38 * t59) * m(7) + t6 * t398 + t11 * t397 + t231 * t400 + t224 * t401 + t223 * t402 + t5 * t399; (t36 * t50 + t45 * t58 + t46 * t59) * m(7) + t187 * t402 + t189 * t401 + t4 * t398 + t3 * t399 + t10 * t397 + t209 * t400; (t50 ^ 2 + t58 ^ 2 + t59 ^ 2) * m(7) + t173 * t1 + t180 * t9 + t175 * t2;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t90(1) t90(2) t90(4) t90(7) t90(11) t90(16); t90(2) t90(3) t90(5) t90(8) t90(12) t90(17); t90(4) t90(5) t90(6) t90(9) t90(13) t90(18); t90(7) t90(8) t90(9) t90(10) t90(14) t90(19); t90(11) t90(12) t90(13) t90(14) t90(15) t90(20); t90(16) t90(17) t90(18) t90(19) t90(20) t90(21);];
Mq  = res;
