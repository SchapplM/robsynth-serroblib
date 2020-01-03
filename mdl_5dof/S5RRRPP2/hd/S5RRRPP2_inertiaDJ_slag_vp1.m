% Calculate time derivative of joint inertia matrix for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:30
% EndTime: 2019-12-31 20:51:42
% DurationCPUTime: 7.78s
% Computational Cost: add. (6360->433), mult. (7114->582), div. (0->0), fcn. (5302->6), ass. (0->234)
t213 = cos(qJ(3));
t211 = sin(qJ(3));
t322 = Icges(5,5) * t211;
t324 = Icges(6,4) * t211;
t326 = Icges(4,4) * t211;
t370 = t322 + t324 - t326 + (-Icges(4,2) - Icges(6,2) - Icges(5,3)) * t213;
t321 = Icges(5,5) * t213;
t323 = Icges(6,4) * t213;
t325 = Icges(4,4) * t213;
t369 = -t321 - t323 + t325 + (Icges(4,1) + Icges(5,1) + Icges(6,1)) * t211;
t259 = Icges(5,3) * t211 + t321;
t261 = Icges(6,2) * t211 + t323;
t263 = -Icges(4,2) * t211 + t325;
t264 = Icges(6,1) * t213 + t324;
t265 = Icges(5,1) * t213 + t322;
t266 = Icges(4,1) * t213 - t326;
t368 = (t264 + t265 + t266) * t211 - (t259 + t261 - t263) * t213;
t210 = qJ(1) + qJ(2);
t207 = cos(t210);
t363 = rSges(6,3) + qJ(5);
t281 = t363 * t207;
t359 = t370 * t211 + t369 * t213;
t206 = sin(t210);
t316 = t207 * t213;
t317 = t207 * t211;
t364 = rSges(6,1) + pkin(4);
t313 = rSges(6,2) * t317 - t363 * t206 + t364 * t316;
t209 = qJD(1) + qJD(2);
t318 = t207 * t209;
t299 = qJD(3) * t213;
t288 = t206 * t299;
t315 = t209 * t211;
t362 = t207 * t315 + t288;
t168 = Icges(6,5) * t211 - Icges(6,6) * t213;
t170 = Icges(4,5) * t211 + Icges(4,6) * t213;
t172 = Icges(5,4) * t211 - Icges(5,6) * t213;
t361 = t168 - t170 - t172;
t258 = Icges(6,5) * t213 + Icges(6,6) * t211;
t260 = Icges(4,5) * t213 - Icges(4,6) * t211;
t262 = Icges(5,4) * t213 + Icges(5,6) * t211;
t358 = t359 * t209 + (t258 - t260 - t262) * qJD(3);
t240 = t263 * t207;
t108 = Icges(4,6) * t206 + t240;
t243 = t266 * t207;
t114 = Icges(4,5) * t206 + t243;
t248 = t108 * t211 - t114 * t213;
t356 = t206 * t248;
t238 = t261 * t207;
t104 = -Icges(6,6) * t206 + t238;
t241 = t264 * t207;
t110 = -Icges(6,5) * t206 + t241;
t252 = t104 * t211 + t110 * t213;
t355 = t206 * t252;
t236 = t259 * t207;
t100 = Icges(5,6) * t206 + t236;
t242 = t265 * t207;
t112 = Icges(5,4) * t206 + t242;
t256 = t100 * t211 + t112 * t213;
t354 = t206 * t256;
t107 = -Icges(4,6) * t207 + t206 * t263;
t113 = -Icges(4,5) * t207 + t206 * t266;
t250 = t107 * t211 - t113 * t213;
t353 = t207 * t250;
t103 = Icges(6,6) * t207 + t206 * t261;
t109 = Icges(6,5) * t207 + t206 * t264;
t254 = t103 * t211 + t109 * t213;
t352 = t207 * t254;
t111 = -Icges(5,4) * t207 + t206 * t265;
t99 = -Icges(5,6) * t207 + t206 * t259;
t268 = t111 * t213 + t211 * t99;
t351 = t207 * t268;
t300 = qJD(3) * t211;
t216 = t368 * qJD(3) + t369 * t299 + t370 * t300;
t346 = 2 * m(3);
t345 = 2 * m(4);
t344 = 2 * m(5);
t343 = 2 * m(6);
t342 = m(5) / 0.2e1;
t341 = m(6) / 0.2e1;
t338 = -rSges(5,1) - pkin(3);
t332 = rSges(4,2) * t211;
t333 = rSges(4,1) * t213;
t157 = (-t332 + t333) * qJD(3);
t336 = m(4) * t157;
t185 = rSges(4,1) * t211 + rSges(4,2) * t213;
t335 = m(4) * t185;
t212 = sin(qJ(1));
t334 = pkin(1) * t212;
t331 = pkin(1) * qJD(1);
t197 = t206 * rSges(5,2);
t196 = t206 * rSges(4,3);
t182 = pkin(3) * t211 - qJ(4) * t213;
t282 = -rSges(6,2) * t213 + t364 * t211;
t278 = -t182 - t282;
t88 = t278 * t207;
t330 = t209 * t88;
t329 = -rSges(6,2) - qJ(4);
t328 = -rSges(5,3) - qJ(4);
t320 = t206 * t209;
t319 = t206 * t213;
t276 = rSges(6,1) * t213 + rSges(6,2) * t211;
t314 = pkin(4) * t319 + t206 * t276 + t281;
t275 = pkin(3) * t213 + qJ(4) * t211;
t125 = t275 * t206;
t195 = pkin(3) * t316;
t126 = qJ(4) * t317 + t195;
t312 = t206 * t125 + t207 * t126;
t130 = qJD(3) * t275 - qJD(4) * t213;
t277 = rSges(5,1) * t213 + rSges(5,3) * t211;
t311 = -t277 * qJD(3) - t130;
t293 = t206 * t315;
t310 = rSges(4,2) * t293 + rSges(4,3) * t318;
t286 = t207 * t299;
t298 = qJD(4) * t211;
t309 = qJ(4) * t286 + t207 * t298;
t308 = rSges(5,2) * t318 + rSges(5,3) * t286;
t184 = rSges(5,1) * t211 - rSges(5,3) * t213;
t307 = -t182 - t184;
t306 = t207 * rSges(4,3) + t206 * t332;
t304 = t207 * pkin(2) + t206 * pkin(7);
t303 = t206 ^ 2 + t207 ^ 2;
t302 = qJD(3) * t206;
t301 = qJD(3) * t207;
t297 = -pkin(3) - t364;
t289 = t206 * t300;
t166 = pkin(3) * t289;
t287 = t207 * t300;
t234 = -t209 * t319 - t287;
t296 = t206 * (t362 * qJ(4) + t209 * t195 + t206 * t298 - t166) + t207 * (pkin(3) * t234 - qJ(4) * t293 + t309) + t125 * t318;
t295 = t212 * t331;
t214 = cos(qJ(1));
t294 = t214 * t331;
t291 = -rSges(4,1) * t289 - t362 * rSges(4,2);
t181 = pkin(7) * t318;
t290 = t181 + t309;
t119 = rSges(5,1) * t316 + rSges(5,3) * t317 + t197;
t283 = -pkin(2) - t333;
t143 = t207 * rSges(3,1) - rSges(3,2) * t206;
t93 = t307 * t207;
t280 = -t364 * t289 - t363 * t320;
t279 = t304 + t126;
t124 = -rSges(3,1) * t318 + rSges(3,2) * t320;
t142 = -rSges(3,1) * t206 - rSges(3,2) * t207;
t267 = t111 * t211 - t213 * t99;
t257 = t100 * t213 - t112 * t211;
t255 = -t103 * t213 + t109 * t211;
t253 = t104 * t213 - t110 * t211;
t251 = t107 * t213 + t113 * t211;
t249 = t108 * t213 + t114 * t211;
t120 = rSges(4,1) * t316 - rSges(4,2) * t317 + t196;
t244 = -pkin(4) * t299 - t276 * qJD(3) - t130;
t123 = t142 * t209;
t239 = t262 * t207;
t237 = t260 * t207;
t235 = t258 * t207;
t84 = t279 + t119;
t90 = t120 + t304;
t202 = t207 * pkin(7);
t89 = t206 * t283 + t202 + t306;
t233 = t211 * t328 + t213 * t338 - pkin(2);
t232 = t211 * t329 + t213 * t297 - pkin(2);
t231 = t233 * t206;
t227 = Icges(5,2) * t209 - qJD(3) * t172;
t223 = Icges(4,3) * t209 - qJD(3) * t170;
t222 = -Icges(6,3) * t209 - qJD(3) * t168;
t48 = t279 + t313;
t199 = t207 * rSges(5,2);
t83 = t199 + t202 + t231;
t44 = (t283 * t207 + (-rSges(4,3) - pkin(7)) * t206) * t209 - t291;
t217 = t206 * t232 - t281;
t47 = t202 + t217;
t215 = (-t361 * t206 + t359 * t207 + t249) * t318 / 0.2e1 + (-t358 * t206 + (-t248 + t252 + t256) * qJD(3) - t368 * t320) * t206 / 0.2e1 - (t358 * t207 + (-t250 + t254 + t268) * qJD(3)) * t207 / 0.2e1 + (t359 * t206 + t361 * t207 + t251 + t255 + t267) * t320 / 0.2e1 - (t253 + t257 + (-t236 - t238 + t240) * t213 + (t241 + t242 + t243) * t211) * t318 / 0.2e1;
t43 = rSges(4,1) * t234 - rSges(4,2) * t286 - pkin(2) * t320 + t181 + t310;
t22 = t209 * t231 + t287 * t338 + t290 + t308;
t160 = rSges(5,1) * t289;
t23 = t160 + t166 + (t328 * t299 - t298) * t206 + ((-rSges(5,2) - pkin(7)) * t206 + t233 * t207) * t209;
t164 = rSges(6,2) * t286;
t11 = -qJD(5) * t206 + t209 * t217 + t287 * t297 + t164 + t290;
t12 = t166 + (-pkin(7) * t209 + t299 * t329 - t298) * t206 + (t209 * t232 - qJD(5)) * t207 - t280;
t208 = t214 * pkin(1);
t129 = t143 + t208;
t128 = t142 - t334;
t127 = t182 * t320;
t117 = rSges(4,1) * t319 - t306;
t116 = t206 * t277 - t199;
t106 = Icges(5,2) * t206 + t239;
t105 = -Icges(5,2) * t207 + t206 * t262;
t102 = Icges(4,3) * t206 + t237;
t101 = -Icges(4,3) * t207 + t206 * t260;
t98 = -Icges(6,3) * t206 + t235;
t97 = Icges(6,3) * t207 + t206 * t258;
t96 = t124 - t294;
t95 = t123 - t295;
t92 = t307 * t206;
t87 = t278 * t206;
t86 = t208 + t90;
t85 = t89 - t334;
t74 = t206 * t227 + t209 * t239;
t73 = t207 * t227 - t262 * t320;
t70 = t206 * t223 + t209 * t237;
t69 = t207 * t223 - t260 * t320;
t66 = t206 * t222 + t209 * t235;
t65 = t207 * t222 - t258 * t320;
t64 = t208 + t84;
t63 = t83 - t334;
t46 = t208 + t48;
t45 = t47 - t334;
t42 = t206 * t311 + t209 * t93;
t41 = t184 * t320 + t207 * t311 + t127;
t40 = t44 - t294;
t39 = t43 - t295;
t38 = t102 * t206 - t248 * t207;
t37 = t101 * t206 - t353;
t36 = t106 * t206 + t256 * t207;
t35 = t105 * t206 + t351;
t34 = -t206 * t98 + t252 * t207;
t33 = -t206 * t97 + t352;
t32 = -t102 * t207 - t356;
t31 = -t101 * t207 - t250 * t206;
t30 = -t106 * t207 + t354;
t29 = -t105 * t207 + t268 * t206;
t28 = t207 * t98 + t355;
t27 = t254 * t206 + t207 * t97;
t26 = t206 * t244 + t330;
t25 = t207 * t244 + t282 * t320 + t127;
t24 = t116 * t206 + t119 * t207 + t312;
t21 = t23 - t294;
t20 = t22 - t295;
t19 = t206 * t314 + t207 * t313 + t312;
t10 = t12 - t294;
t9 = t11 - t295;
t2 = t206 * (rSges(5,3) * t288 - t160) + t207 * (-rSges(5,1) * t287 + t308) + (t207 * t116 + (-t119 - t126 + t197) * t206) * t209 + t296;
t1 = (rSges(6,2) * t288 + (-t126 - t313) * t209 + t280) * t206 + (t164 - t364 * t287 + (-t281 + t314) * t209) * t207 + t296;
t3 = [(t10 * t45 + t46 * t9) * t343 + (t20 * t64 + t21 * t63) * t344 + (t39 * t86 + t40 * t85) * t345 + (t128 * t96 + t129 * t95) * t346 + t216; m(6) * (t10 * t47 + t11 * t46 + t12 * t45 + t48 * t9) + m(5) * (t20 * t84 + t21 * t83 + t22 * t64 + t23 * t63) + m(4) * (t39 * t90 + t40 * t89 + t43 * t86 + t44 * t85) + m(3) * (t123 * t129 + t124 * t128 + t142 * t96 + t143 * t95) + t216; (t22 * t84 + t23 * t83) * t344 + (t11 * t48 + t12 * t47) * t343 + (t43 * t90 + t44 * t89) * t345 + (t123 * t143 + t124 * t142) * t346 + t216; t215 + (-t206 * t86 - t207 * t85) * t336 + ((-t209 * t86 - t40) * t207 + (t209 * t85 - t39) * t206) * t335 + m(6) * (t10 * t88 + t25 * t45 + t26 * t46 + t87 * t9) + m(5) * (t20 * t92 + t21 * t93 + t41 * t63 + t42 * t64); t215 + m(6) * (t11 * t87 + t12 * t88 + t25 * t47 + t26 * t48) + m(5) * (t22 * t92 + t23 * t93 + t41 * t83 + t42 * t84) + (-t206 * t90 - t207 * t89) * t336 + ((-t209 * t90 - t44) * t207 + (t209 * t89 - t43) * t206) * t335; (t2 * t24 + t41 * t93 + t42 * t92) * t344 + (t1 * t19 + t25 * t88 + t26 * t87) * t343 + t206 * ((t206 * t69 + (t37 + t356) * t209) * t206 + (t38 * t209 + (t107 * t299 + t113 * t300) * t207 + (-t249 * qJD(3) - t209 * t250 - t70) * t206) * t207) - t207 * ((-t207 * t66 + (t28 - t352) * t209) * t207 + (t27 * t209 + (t104 * t299 - t110 * t300) * t206 + (t255 * qJD(3) + t209 * t252 + t65) * t207) * t206) - t207 * ((t207 * t70 + (t32 + t353) * t209) * t207 + (t31 * t209 + (-t108 * t299 - t114 * t300) * t206 + (t251 * qJD(3) - t209 * t248 - t69) * t207) * t206) - t207 * ((t207 * t74 + (t30 - t351) * t209) * t207 + (t29 * t209 + (t100 * t299 - t112 * t300) * t206 + (t267 * qJD(3) + t209 * t256 - t73) * t207) * t206) + t206 * ((-t206 * t65 + (t33 - t355) * t209) * t206 + (t34 * t209 + (-t103 * t299 + t109 * t300) * t207 + (t253 * qJD(3) + t209 * t254 + t66) * t206) * t207) + t206 * ((t206 * t73 + (t35 - t354) * t209) * t206 + (t36 * t209 + (t111 * t300 - t299 * t99) * t207 + (t257 * qJD(3) + t209 * t268 - t74) * t206) * t207) + ((t117 * t206 + t120 * t207) * (((-t120 + t196) * t209 + t291) * t206 + (t209 * t117 - t185 * t301 + t310) * t207) + t303 * t185 * t157) * t345 + ((-t27 - t29 - t31) * t207 + (t28 + t30 + t32) * t206) * t320 + ((-t33 - t35 - t37) * t207 + (t34 + t36 + t38) * t206) * t318; 0.2e1 * ((t206 * t46 + t207 * t45) * t341 + (t206 * t64 + t207 * t63) * t342) * t299 + 0.2e1 * ((t10 * t207 + t206 * t9 + t318 * t46 - t320 * t45) * t341 + (t20 * t206 + t207 * t21 + t318 * t64 - t320 * t63) * t342) * t211; 0.2e1 * ((t206 * t84 + t207 * t83) * t342 + (t206 * t48 + t207 * t47) * t341) * t299 + 0.2e1 * ((t206 * t22 + t207 * t23 + t318 * t84 - t320 * t83) * t342 + (t11 * t206 + t12 * t207 + t318 * t48 - t320 * t47) * t341) * t211; 0.2e1 * ((t301 * t93 + t302 * t92 - t2) * t342 + (t301 * t88 + t302 * t87 - t1) * t341) * t213 + 0.2e1 * ((qJD(3) * t24 + t206 * t42 + t207 * t41 + t318 * t92 - t320 * t93) * t342 + (qJD(3) * t19 + t206 * t26 + t207 * t25 + t318 * t87 - t320 * t88) * t341) * t211; 0.4e1 * (t342 + t341) * (-0.1e1 + t303) * t211 * t299; m(6) * ((-t209 * t45 + t9) * t207 + (-t209 * t46 - t10) * t206); m(6) * ((-t209 * t47 + t11) * t207 + (-t209 * t48 - t12) * t206); m(6) * ((t26 - t330) * t207 + (-t209 * t87 - t25) * t206); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
