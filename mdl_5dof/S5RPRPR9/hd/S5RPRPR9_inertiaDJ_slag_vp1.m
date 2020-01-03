% Calculate time derivative of joint inertia matrix for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR9_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR9_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:43
% EndTime: 2019-12-31 18:23:58
% DurationCPUTime: 8.63s
% Computational Cost: add. (10982->512), mult. (14908->735), div. (0->0), fcn. (13855->8), ass. (0->269)
t386 = Icges(5,1) + Icges(4,3);
t198 = cos(qJ(3));
t195 = sin(qJ(3));
t314 = Icges(5,6) * t195;
t322 = Icges(4,4) * t195;
t385 = -t314 - t322 + (-Icges(4,2) - Icges(5,3)) * t198;
t313 = Icges(5,6) * t198;
t321 = Icges(4,4) * t198;
t384 = -t313 - t321 + (-Icges(4,1) - Icges(5,2)) * t195;
t231 = Icges(4,5) * t198 - Icges(4,6) * t195;
t233 = Icges(5,4) * t198 - Icges(5,5) * t195;
t383 = t231 - t233;
t193 = qJ(1) + pkin(8);
t191 = cos(t193);
t190 = sin(t193);
t238 = Icges(4,1) * t198 - t322;
t117 = Icges(4,5) * t190 + t191 * t238;
t229 = Icges(5,2) * t198 - t314;
t120 = Icges(5,4) * t190 - t191 * t229;
t355 = t117 - t120;
t235 = -Icges(4,2) * t195 + t321;
t115 = Icges(4,6) * t190 + t191 * t235;
t227 = -Icges(5,3) * t195 + t313;
t118 = Icges(5,5) * t190 - t191 * t227;
t357 = t115 - t118;
t375 = -t195 * t357 + t198 * t355;
t382 = t375 * t191;
t114 = -Icges(4,6) * t191 + t190 * t235;
t350 = Icges(5,5) * t191 + t190 * t227;
t358 = t114 + t350;
t116 = -Icges(4,5) * t191 + t190 * t238;
t349 = Icges(5,4) * t191 + t190 * t229;
t356 = t116 + t349;
t377 = t190 * t386 + t191 * t383;
t381 = t385 * qJD(3);
t380 = t384 * qJD(3);
t352 = -t195 * t358 + t198 * t356;
t379 = t352 * t191;
t378 = t190 * t383 - t191 * t386;
t376 = ((Icges(5,5) - Icges(4,6)) * t198 + (Icges(5,4) - Icges(4,5)) * t195) * qJD(3);
t288 = qJD(3) * t195;
t270 = t191 * t288;
t291 = qJD(1) * t198;
t273 = t190 * t291;
t374 = t270 + t273;
t373 = t377 * qJD(1);
t197 = cos(qJ(5));
t194 = sin(qJ(5));
t319 = Icges(6,4) * t194;
t232 = Icges(6,2) * t197 + t319;
t139 = Icges(6,6) * t195 - t198 * t232;
t318 = Icges(6,4) * t197;
t236 = Icges(6,1) * t194 + t318;
t140 = Icges(6,5) * t195 - t198 * t236;
t372 = t139 * t197 + t140 * t194;
t371 = -t377 * t191 + t379 + (t375 + t378) * t190;
t370 = t190 * t355 + t191 * t356;
t369 = t190 * t357 + t191 * t358;
t368 = t190 ^ 2;
t367 = t191 ^ 2;
t366 = qJD(3) / 0.2e1;
t365 = -t190 * t352 + t191 * t378;
t362 = t190 * t377 + t382;
t361 = -qJD(1) * t378 + t191 * t376;
t360 = -t190 * t376 - t373;
t330 = rSges(4,1) * t198;
t256 = -rSges(4,2) * t195 + t330;
t327 = rSges(4,3) * t191;
t125 = t190 * t256 - t327;
t305 = t191 * t198;
t306 = t191 * t195;
t354 = -rSges(4,2) * t306 + t190 * rSges(4,3);
t126 = rSges(4,1) * t305 + t354;
t171 = rSges(4,1) * t195 + rSges(4,2) * t198;
t215 = qJD(3) * t171;
t292 = qJD(1) * t195;
t274 = t190 * t292;
t293 = qJD(1) * t191;
t201 = rSges(4,2) * t274 + rSges(4,3) * t293 - t191 * t215;
t30 = (qJD(1) * t125 + t201) * t191 + (-t190 * t215 + (-t126 + t354) * qJD(1)) * t190;
t347 = 2 * m(4);
t359 = t30 * t347;
t353 = t190 * rSges(5,1) - rSges(5,2) * t305;
t294 = qJD(1) * t190;
t346 = 2 * m(5);
t345 = 2 * m(6);
t344 = m(5) / 0.2e1;
t343 = m(6) / 0.2e1;
t342 = t190 / 0.2e1;
t341 = -t191 / 0.2e1;
t340 = t191 / 0.2e1;
t339 = t195 / 0.2e1;
t338 = rSges(6,3) + pkin(7);
t337 = m(4) * t171;
t336 = sin(qJ(1)) * pkin(1);
t335 = pkin(3) * t198;
t334 = pkin(4) * t191;
t185 = t190 * pkin(4);
t192 = cos(qJ(1)) * pkin(1);
t333 = qJD(1) / 0.2e1;
t285 = qJD(5) * t198;
t103 = (Icges(6,2) * t194 - t318) * t285 + (Icges(6,6) * t198 + t195 * t232) * qJD(3);
t230 = Icges(6,5) * t194 + Icges(6,6) * t197;
t102 = (-Icges(6,5) * t197 + Icges(6,6) * t194) * t285 + (Icges(6,3) * t198 + t195 * t230) * qJD(3);
t138 = Icges(6,3) * t195 - t198 * t230;
t287 = qJD(3) * t198;
t258 = t194 * t139 * t285 + t195 * t102 + t138 * t287 + t288 * t372;
t104 = (-Icges(6,1) * t197 + t319) * t285 + (Icges(6,5) * t198 + t195 * t236) * qJD(3);
t304 = t194 * t104;
t62 = t138 * t195 - t198 * t372;
t332 = ((-t304 + (-qJD(5) * t140 - t103) * t197) * t198 + t258) * t195 + t62 * t287;
t264 = qJD(5) * t195 + qJD(1);
t202 = -t194 * t264 + t197 * t287;
t263 = qJD(5) + t292;
t221 = t190 * t263;
t75 = t191 * t202 - t197 * t221;
t203 = t194 * t287 + t197 * t264;
t76 = t191 * t203 - t194 * t221;
t331 = rSges(6,1) * t76 + rSges(6,2) * t75;
t329 = rSges(5,1) * t191;
t328 = rSges(5,2) * t195;
t326 = -rSges(5,3) - qJ(4);
t302 = t195 * t197;
t134 = -t190 * t194 + t191 * t302;
t303 = t194 * t195;
t135 = t190 * t197 + t191 * t303;
t83 = rSges(6,1) * t135 + rSges(6,2) * t134 + rSges(6,3) * t305;
t325 = pkin(7) * t305 + t185 + t83;
t307 = t190 * t198;
t136 = t190 * t302 + t191 * t194;
t137 = t190 * t303 - t191 * t197;
t255 = -rSges(6,1) * t137 - rSges(6,2) * t136;
t84 = rSges(6,3) * t307 - t255;
t324 = pkin(7) * t307 - t334 + t84;
t311 = qJ(4) * t195;
t310 = qJ(4) * t198;
t251 = t311 + t335;
t141 = t251 * t190;
t142 = pkin(3) * t305 + qJ(4) * t306;
t301 = t141 * t190 + t142 * t191;
t144 = qJD(3) * t251 - qJD(4) * t198;
t253 = -rSges(5,2) * t198 + rSges(5,3) * t195;
t300 = -qJD(3) * t253 - t144;
t269 = t191 * t287;
t286 = qJD(4) * t195;
t299 = qJ(4) * t269 + t191 * t286;
t169 = pkin(3) * t195 - t310;
t252 = rSges(5,3) * t198 + t328;
t298 = -t169 + t252;
t297 = t367 + t368;
t290 = qJD(3) * t190;
t289 = qJD(3) * t191;
t284 = -pkin(3) - t338;
t18 = t102 * t305 + t103 * t134 + t104 * t135 - t138 * t374 + t139 * t75 + t140 * t76;
t79 = Icges(6,4) * t135 + Icges(6,2) * t134 + Icges(6,6) * t305;
t81 = Icges(6,1) * t135 + Icges(6,4) * t134 + Icges(6,5) * t305;
t245 = t194 * t81 + t197 * t79;
t34 = Icges(6,5) * t76 + Icges(6,6) * t75 - Icges(6,3) * t374;
t36 = Icges(6,4) * t76 + Icges(6,2) * t75 - Icges(6,6) * t374;
t38 = Icges(6,1) * t76 + Icges(6,4) * t75 - Icges(6,5) * t374;
t77 = Icges(6,5) * t135 + Icges(6,6) * t134 + Icges(6,3) * t305;
t9 = (qJD(3) * t245 + t34) * t195 + (qJD(3) * t77 - t194 * t38 - t197 * t36 + (t194 * t79 - t197 * t81) * qJD(5)) * t198;
t283 = t18 / 0.2e1 + t9 / 0.2e1;
t271 = t190 * t288;
t161 = pkin(3) * t271;
t272 = t191 * t291;
t282 = t141 * t293 + t190 * (pkin(3) * t272 + t190 * t286 - t161 + (t190 * t287 + t191 * t292) * qJ(4)) + t191 * (-pkin(3) * t374 - qJ(4) * t274 + t299);
t80 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t307;
t82 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t307;
t244 = t194 * t82 + t197 * t80;
t207 = -t271 + t272;
t220 = t191 * t263;
t73 = t190 * t202 + t197 * t220;
t74 = t190 * t203 + t194 * t220;
t33 = Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t207;
t35 = Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t207;
t37 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t207;
t78 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t307;
t10 = (qJD(3) * t244 + t33) * t195 + (qJD(3) * t78 - t194 * t37 - t197 * t35 + (t194 * t80 - t197 * t82) * qJD(5)) * t198;
t17 = t102 * t307 + t103 * t136 + t104 * t137 + t138 * t207 + t139 * t73 + t140 * t74;
t279 = t10 / 0.2e1 + t17 / 0.2e1;
t32 = t195 * t78 - t198 * t244;
t48 = t136 * t139 + t137 * t140 + t138 * t307;
t278 = t32 / 0.2e1 + t48 / 0.2e1;
t31 = t195 * t77 - t198 * t245;
t47 = t134 * t139 + t135 * t140 + t138 * t305;
t277 = -t47 / 0.2e1 - t31 / 0.2e1;
t180 = pkin(6) * t293;
t276 = t180 + t299;
t275 = pkin(2) * t191 + t190 * pkin(6) + t192;
t268 = -t233 * qJD(3) / 0.2e1 + t231 * t366;
t267 = -t288 / 0.2e1;
t186 = t191 * pkin(6);
t266 = t186 - t336;
t254 = rSges(6,1) * t194 + rSges(6,2) * t197;
t143 = rSges(6,3) * t195 - t198 * t254;
t265 = pkin(7) * t195 + t143;
t108 = t298 * t191;
t262 = rSges(5,1) * t293 + rSges(5,2) * t374 + rSges(5,3) * t269;
t261 = -t169 - t265;
t260 = t74 * rSges(6,1) + t73 * rSges(6,2);
t259 = t195 * t326 - pkin(2);
t26 = t134 * t79 + t135 * t81 + t305 * t77;
t27 = t134 * t80 + t135 * t82 + t305 * t78;
t250 = t190 * t27 + t191 * t26;
t15 = t190 * t26 - t191 * t27;
t28 = t136 * t79 + t137 * t81 + t307 * t77;
t29 = t136 * t80 + t137 * t82 + t307 * t78;
t249 = t190 * t29 + t191 * t28;
t16 = t190 * t28 - t191 * t29;
t248 = t190 * t32 + t191 * t31;
t247 = t190 * t31 - t191 * t32;
t246 = t190 * t83 - t191 * t84;
t239 = t275 + t142;
t127 = rSges(5,3) * t306 + t353;
t105 = (-rSges(6,1) * t197 + rSges(6,2) * t194) * t285 + (rSges(6,3) * t198 + t195 * t254) * qJD(3);
t219 = -pkin(7) * t287 - t105 - t144;
t218 = -pkin(2) - t256;
t98 = t261 * t191;
t217 = t198 * t33 - t288 * t78;
t216 = t198 * t34 - t288 * t77;
t205 = t198 * t284 - pkin(2) - t311;
t204 = (rSges(5,2) - pkin(3)) * t198 + t259;
t200 = t190 * t205 - t336;
t181 = pkin(4) * t293;
t155 = t256 * qJD(3);
t145 = t169 * t294;
t128 = t190 * t253 - t329;
t107 = t298 * t190;
t101 = t126 + t275;
t100 = t190 * t218 + t266 + t327;
t97 = t261 * t190;
t68 = t127 + t239;
t67 = t190 * t204 + t266 + t329;
t66 = t171 * t290 + (-t192 + (-rSges(4,3) - pkin(6)) * t190 + t218 * t191) * qJD(1);
t65 = t180 + (-t336 + (-pkin(2) - t330) * t190) * qJD(1) + t201;
t64 = qJD(1) * t108 + t190 * t300;
t63 = t191 * t300 - t252 * t294 + t145;
t60 = -t143 * t305 + t195 * t83;
t59 = t143 * t307 - t195 * t84;
t50 = t239 + t325;
t49 = t186 + t200 + t255 + t334;
t46 = t127 * t191 + t128 * t190 + t301;
t45 = t246 * t198;
t44 = t161 + (-t286 + (t198 * t326 - t328) * qJD(3)) * t190 + (-t192 + (-rSges(5,1) - pkin(6)) * t190 + t204 * t191) * qJD(1);
t43 = -pkin(3) * t270 + (-t336 + (t259 - t335) * t190) * qJD(1) + t262 + t276;
t42 = qJD(1) * t98 + t190 * t219;
t41 = t191 * t219 + t265 * t294 + t145;
t40 = -rSges(6,3) * t374 + t331;
t39 = rSges(6,3) * t207 + t260;
t25 = t190 * t324 + t191 * t325 + t301;
t24 = t161 + (-t286 + (t195 * t338 - t310) * qJD(3)) * t190 + (-t192 + (-pkin(4) - pkin(6)) * t190 + t205 * t191) * qJD(1) - t260;
t23 = qJD(1) * t200 + t270 * t284 + t181 + t276 + t331;
t21 = (-t143 * t290 - t39) * t195 + (-qJD(3) * t84 + t105 * t190 + t143 * t293) * t198;
t20 = (t143 * t289 + t40) * t195 + (qJD(3) * t83 - t105 * t191 + t143 * t294) * t198;
t19 = (qJD(1) * t128 + t262) * t191 + (t252 * t290 + (-t127 - t142 + t353) * qJD(1)) * t190 + t282;
t14 = t246 * t288 + (-t190 * t40 + t191 * t39 + (-t190 * t84 - t191 * t83) * qJD(1)) * t198;
t13 = t195 * t48 + t198 * t249;
t12 = t195 * t47 + t198 * t250;
t11 = (-pkin(7) * t270 + qJD(1) * t324 + t181 + t40) * t191 + (-pkin(7) * t271 + t39 + (-t142 - t325 + t185) * qJD(1)) * t190 + t282;
t8 = t134 * t35 + t135 * t37 + t191 * t217 - t273 * t78 + t75 * t80 + t76 * t82;
t7 = t134 * t36 + t135 * t38 + t191 * t216 - t273 * t77 + t75 * t79 + t76 * t81;
t6 = t136 * t35 + t137 * t37 + t190 * t217 + t272 * t78 + t73 * t80 + t74 * t82;
t5 = t136 * t36 + t137 * t38 + t190 * t216 + t272 * t77 + t73 * t79 + t74 * t81;
t4 = qJD(1) * t250 + t190 * t7 - t191 * t8;
t3 = qJD(1) * t249 + t190 * t5 - t191 * t6;
t2 = (-qJD(3) * t250 + t18) * t195 + (-qJD(1) * t15 + qJD(3) * t47 + t190 * t8 + t191 * t7) * t198;
t1 = (-qJD(3) * t249 + t17) * t195 + (-qJD(1) * t16 + qJD(3) * t48 + t190 * t6 + t191 * t5) * t198;
t22 = [(t23 * t50 + t24 * t49) * t345 + (t100 * t66 + t101 * t65) * t347 + (t43 * t68 + t44 * t67) * t346 + t258 - t198 * t304 + (-t103 * t198 - t140 * t285) * t197 + (t238 + t229 + t385) * t288 + (t227 + t235 - t384) * t287; 0; 0; (t268 * t191 - t279) * t191 + (t268 * t190 + t283) * t190 + m(4) * ((-t190 * t65 - t191 * t66) * t171 + (-t100 * t191 - t101 * t190) * t155) + m(5) * (t107 * t43 + t108 * t44 + t63 * t67 + t64 * t68) + m(6) * (t23 * t97 + t24 * t98 + t41 * t49 + t42 * t50) + ((-t357 * qJD(3) + t191 * t380) * t342 + (-t358 * qJD(3) + t190 * t380) * t341 + (t341 * t355 - t342 * t356) * qJD(1)) * t195 + ((t355 * qJD(3) + t191 * t381) * t342 + (t356 * qJD(3) + t190 * t381) * t341 + (t341 * t357 - t342 * t358) * qJD(1)) * t198 + ((-t101 * t337 + (t115 / 0.2e1 - t118 / 0.2e1) * t198 + (t117 / 0.2e1 - t120 / 0.2e1) * t195 - t277) * t191 + (t100 * t337 + (t350 / 0.2e1 + t114 / 0.2e1) * t198 + (t349 / 0.2e1 + t116 / 0.2e1) * t195 + t278) * t190) * qJD(1); m(4) * t30 + m(5) * t19 + m(6) * t11; t16 * t294 + t15 * t293 + (t11 * t25 + t41 * t98 + t42 * t97) * t345 + (t107 * t64 + t108 * t63 + t19 * t46) * t346 + t297 * t171 * t155 * t347 + (t126 * t359 + t365 * t294 + t360 * t367 - t3 + (-t371 + t379) * t293) * t191 + (t4 + t125 * t359 + t361 * t368 + t362 * t293 + ((t360 - t373) * t190 + t361 * t191 + t370 * t288 + t369 * t287 + (-t195 * t370 - t198 * t369) * qJD(3) + ((t352 + t377) * t190 - t382 + t362 + t365) * qJD(1)) * t191 + (-t190 * t375 + t371) * t294) * t190; 0.2e1 * ((t190 * t50 + t191 * t49) * t343 + (t190 * t68 + t191 * t67) * t344) * t287 + 0.2e1 * ((t190 * t23 + t191 * t24 + t293 * t50 - t294 * t49) * t343 + (t190 * t43 + t191 * t44 + t293 * t68 - t294 * t67) * t344) * t195; (m(5) + m(6)) * t288; 0.2e1 * ((t289 * t98 + t290 * t97 - t11) * t343 + (t107 * t290 + t108 * t289 - t19) * t344) * t198 + 0.2e1 * ((qJD(3) * t25 + t190 * t42 + t191 * t41 + t293 * t97 - t294 * t98) * t343 + (qJD(3) * t46 + t107 * t293 - t108 * t294 + t190 * t64 + t191 * t63) * t344) * t195; 0.4e1 * (t344 + t343) * (-0.1e1 + t297) * t195 * t287; m(6) * (t20 * t50 + t21 * t49 + t23 * t60 + t24 * t59) + (-t190 * t278 + t191 * t277) * t288 + (t283 * t191 + t279 * t190 + (t190 * t277 + t191 * t278) * qJD(1)) * t198 + t332; m(6) * t14; m(6) * (-t11 * t45 + t14 * t25 + t20 * t97 + t21 * t98 + t41 * t59 + t42 * t60) + (-t1 / 0.2e1 + t15 * t267 + t12 * t333 + (qJD(1) * t31 - t10) * t339) * t191 + (t13 * t333 + t16 * t267 + t2 / 0.2e1 + (qJD(1) * t32 + t9) * t339) * t190 + (t3 * t342 + t4 * t340 + t247 * t366 + (t16 * t340 - t190 * t15 / 0.2e1) * qJD(1)) * t198; m(6) * ((-t14 + (t190 * t60 + t191 * t59) * qJD(3)) * t198 + (-qJD(3) * t45 + t190 * t20 + t191 * t21 + (-t190 * t59 + t191 * t60) * qJD(1)) * t195); (-t14 * t45 + t20 * t60 + t21 * t59) * t345 + ((-t12 * t191 - t13 * t190 - t195 * t248) * qJD(3) + t332) * t195 + (t191 * t2 + t190 * t1 + t195 * (t10 * t190 + t191 * t9) + (t195 * t62 + t198 * t248) * qJD(3) + (-t12 * t190 + t13 * t191 - t195 * t247) * qJD(1)) * t198;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t22(1), t22(2), t22(4), t22(7), t22(11); t22(2), t22(3), t22(5), t22(8), t22(12); t22(4), t22(5), t22(6), t22(9), t22(13); t22(7), t22(8), t22(9), t22(10), t22(14); t22(11), t22(12), t22(13), t22(14), t22(15);];
Mq = res;
