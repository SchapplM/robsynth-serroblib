% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:44
% DurationCPUTime: 3.69s
% Computational Cost: add. (13167->315), mult. (12356->444), div. (0->0), fcn. (11074->8), ass. (0->203)
t248 = cos(qJ(1));
t242 = pkin(8) + qJ(4);
t231 = qJ(5) + t242;
t226 = sin(t231);
t227 = cos(t231);
t201 = rSges(6,1) * t227 - rSges(6,2) * t226;
t230 = cos(t242);
t322 = pkin(4) * t230;
t272 = (t201 + t322) * t248;
t134 = t272 * t248;
t247 = sin(qJ(1));
t295 = t230 * t247;
t222 = pkin(4) * t295;
t145 = t201 * t247 + t222;
t229 = sin(t242);
t209 = rSges(5,1) * t230 - rSges(5,2) * t229;
t243 = t247 ^ 2;
t244 = t248 ^ 2;
t223 = t243 + t244;
t354 = m(6) / 0.2e1;
t380 = -m(5) / 0.2e1;
t310 = t223 * t209 * t380 + (-t145 * t247 - t134) * t354;
t191 = t209 * t247;
t192 = t209 * t248;
t114 = t247 * t191 + t192 * t248;
t297 = t227 * t247;
t300 = t226 * t247;
t182 = rSges(6,1) * t297 - rSges(6,2) * t300;
t139 = t182 + t222;
t356 = m(5) / 0.2e1;
t311 = t114 * t356 + (t247 * t139 + t134) * t354;
t26 = t311 - t310;
t383 = t26 * qJD(1);
t307 = Icges(6,4) * t226;
t265 = Icges(6,2) * t227 + t307;
t154 = Icges(6,6) * t248 + t265 * t247;
t217 = Icges(6,4) * t297;
t156 = Icges(6,1) * t300 + Icges(6,5) * t248 + t217;
t259 = -t154 * t227 - t156 * t226;
t382 = t259 * t248;
t306 = Icges(6,4) * t227;
t197 = -Icges(6,2) * t226 + t306;
t267 = Icges(6,1) * t226 + t306;
t381 = t197 + t267;
t378 = -t247 / 0.2e1;
t330 = t247 / 0.2e1;
t328 = t248 / 0.2e1;
t269 = rSges(6,1) * t226 + rSges(6,2) * t227;
t119 = t223 * t269;
t261 = Icges(6,5) * t226 + Icges(6,6) * t227;
t377 = t261 * t247;
t376 = t261 * t248;
t303 = t269 * t247;
t308 = Icges(5,4) * t230;
t205 = -Icges(5,2) * t229 + t308;
t268 = Icges(5,1) * t229 + t308;
t375 = (t268 / 0.2e1 + t205 / 0.2e1) * t230;
t157 = -Icges(6,5) * t247 + t267 * t248;
t288 = t197 * t248 + t157;
t289 = -Icges(6,2) * t300 + t156 + t217;
t155 = -Icges(6,6) * t247 + t265 * t248;
t199 = Icges(6,1) * t227 - t307;
t290 = -t199 * t248 + t155;
t291 = -t199 * t247 + t154;
t374 = -(t290 * t247 - t291 * t248) * t226 + (t288 * t247 - t289 * t248) * t227;
t174 = -Icges(5,5) * t247 + t268 * t248;
t284 = t205 * t248 + t174;
t221 = Icges(5,4) * t295;
t296 = t229 * t247;
t173 = Icges(5,1) * t296 + Icges(5,5) * t248 + t221;
t285 = -Icges(5,2) * t296 + t173 + t221;
t309 = Icges(5,4) * t229;
t266 = Icges(5,2) * t230 + t309;
t172 = -Icges(5,6) * t247 + t266 * t248;
t207 = Icges(5,1) * t230 - t309;
t286 = -t207 * t248 + t172;
t171 = Icges(5,6) * t248 + t266 * t247;
t287 = -t207 * t247 + t171;
t373 = -(t286 * t247 - t287 * t248) * t229 + (t284 * t247 - t285 * t248) * t230;
t262 = Icges(6,5) * t227 - Icges(6,6) * t226;
t176 = t247 * t262;
t177 = t262 * t248;
t319 = (t244 * t176 + (-t248 * t177 - t374) * t247) * t328 + (-t243 * t177 + (t247 * t176 + t374) * t248) * t330;
t183 = t201 * t248;
t366 = t247 * t182 + t248 * t183;
t249 = -t247 * rSges(6,3) + t269 * t248;
t148 = t248 * t249;
t158 = rSges(6,3) * t248 + t303;
t245 = sin(pkin(8));
t321 = t245 * pkin(3);
t323 = pkin(4) * t229;
t216 = t321 + t323;
t320 = -pkin(6) - qJ(3);
t280 = -pkin(7) + t320;
t225 = t248 * t280;
t228 = t248 * t320;
t281 = -t248 * t216 - t247 * t280;
t50 = t248 * (t247 * t320 + t248 * t321 + t281) - t148 + (t225 - t228 - t158 + (-t216 + t321) * t247) * t247;
t6 = t319 + m(6) * (-t134 * t269 - t145 * t303 - t366 * t50);
t371 = t6 * qJD(5);
t370 = t247 * t272;
t369 = (t230 * t172 + t229 * t174) * t248;
t368 = (t155 * t227 + t157 * t226) * t248;
t232 = t248 * qJ(2);
t240 = t247 * rSges(5,3);
t270 = t229 * rSges(5,1) + t230 * rSges(5,2);
t251 = t270 + t321;
t115 = t232 - t240 + (-pkin(1) + t320) * t247 + t251 * t248;
t116 = -t228 + (rSges(5,3) + pkin(1)) * t248 + (qJ(2) + t251) * t247;
t367 = -t115 * t247 + t248 * t116;
t113 = -t191 * t248 + t247 * t192;
t88 = -t139 * t248 + t370;
t312 = t113 * t356 + t88 * t354;
t364 = t381 * t226 - (-t265 + t199) * t227;
t275 = -t381 * t227 / 0.2e1 + (-t199 / 0.2e1 + t265 / 0.2e1) * t226;
t153 = -Icges(6,3) * t247 + t376;
t59 = t248 * (Icges(6,3) * t248 + t377) + t154 * t297 + t156 * t300;
t60 = -t248 * t153 - t155 * t297 - t157 * t300;
t62 = -t247 * t153 + t368;
t5 = (t60 * t247 + t248 * t59) * t378 + ((t60 - t382) * t247 + (t62 - t368 + (t153 + t259) * t247 + t59) * t248) * t330 + (t243 * t153 + t62 * t247 + (t60 + (t153 - t259) * t248 + t382) * t248) * t328;
t359 = 0.2e1 * t223;
t358 = 0.4e1 * qJD(1);
t357 = 2 * qJD(4);
t355 = -m(6) / 0.2e1;
t271 = rSges(4,1) * t245 + rSges(4,2) * cos(pkin(8));
t279 = rSges(4,3) + pkin(1) + qJ(3);
t131 = -t279 * t247 + t271 * t248 + t232;
t132 = t279 * t248 + (qJ(2) + t271) * t247;
t353 = m(4) * (-t131 * t247 + t248 * t132);
t352 = m(4) * (t131 * t248 + t247 * t132);
t351 = m(5) * (t115 * t192 + t116 * t191);
t350 = m(5) * t367;
t349 = m(5) * (t115 * t248 + t247 * t116);
t100 = -t247 * pkin(1) + t232 + t249 - t281;
t101 = -t225 + (rSges(6,3) + pkin(1)) * t248 + (qJ(2) + t216 + t269) * t247;
t96 = t248 * t101;
t273 = -t100 * t303 + t269 * t96;
t345 = m(6) * (t145 * t183 - t182 * t272 + t273);
t344 = m(6) * (t88 * t201 + t273);
t343 = m(6) * (t100 * t272 + t101 * t139);
t342 = m(6) * (t100 * t183 + t101 * t182);
t341 = m(6) * (-t100 * t247 + t96);
t340 = m(6) * (t100 * t248 + t247 * t101);
t336 = m(6) * (t145 * t248 - t370);
t105 = -t182 * t248 + t247 * t183;
t333 = t105 / 0.2e1;
t329 = -t248 / 0.2e1;
t327 = m(3) * ((rSges(3,3) * t248 + t232) * t248 + (rSges(3,3) + qJ(2)) * t243);
t318 = m(6) * qJD(5);
t253 = m(6) * t366;
t254 = t223 * t201 * t355;
t63 = -t253 / 0.2e1 + t254;
t293 = t63 * qJD(1);
t86 = (t355 - m(4) / 0.2e1 + t380) * t359;
t292 = t86 * qJD(1);
t278 = t105 * t355;
t263 = Icges(5,5) * t229 + Icges(5,6) * t230;
t169 = Icges(5,3) * t248 + t263 * t247;
t67 = t248 * t169 + t171 * t295 + t173 * t296;
t170 = -Icges(5,3) * t247 + t263 * t248;
t68 = -t248 * t170 - t172 * t295 - t174 * t296;
t277 = t269 + t323;
t276 = t223 * t270;
t264 = Icges(5,5) * t230 - Icges(5,6) * t229;
t144 = t277 * t247;
t146 = t277 * t248;
t260 = -t144 * t247 - t146 * t248;
t256 = -t230 * t171 - t229 * t173;
t252 = -t5 + (t288 * t226 + t290 * t227 + t364 * t248 - t377) * t330 + (-t289 * t226 - t291 * t227 - t364 * t247 - t376) * t328;
t250 = -t275 + (t328 + t329) * (t155 * t226 - t157 * t227);
t186 = t264 * t248;
t185 = t247 * t264;
t149 = t247 * t169;
t112 = t119 * t318;
t99 = t318 * t333;
t98 = -t247 * t158 - t148;
t87 = -t223 * t322 - t366;
t85 = (-m(6) / 0.4e1 - m(4) / 0.4e1 - m(5) / 0.4e1) * t359 + (m(4) + m(5) + m(6)) * t223 / 0.2e1;
t82 = t336 / 0.2e1;
t70 = -t247 * t170 + t369;
t69 = t256 * t248 + t149;
t64 = t253 / 0.2e1 + t254;
t44 = t70 * t247 + t248 * t69;
t43 = t68 * t247 + t248 * t67;
t34 = t344 / 0.2e1;
t33 = t275 + t342;
t32 = t345 / 0.2e1;
t31 = t82 - t312;
t30 = -t336 / 0.2e1 + t312;
t29 = t82 + t312;
t25 = t310 + t311;
t24 = t341 + t350 + t353;
t21 = t327 + t340 + t349 + t352;
t18 = -t375 + (-t207 / 0.2e1 + t266 / 0.2e1) * t229 + t351 + t343 + t275;
t16 = t243 * t170 + (-t149 + t68 + (t170 - t256) * t248) * t248;
t15 = (-t69 + t149 + t68) * t247 + (t70 - t369 + (t170 + t256) * t247 + t67) * t248;
t8 = m(6) * (-t201 * t119 - t366 * t98) + t319;
t7 = t8 * qJD(5);
t4 = t32 - t344 / 0.2e1 + t5;
t3 = t34 - t345 / 0.2e1 + t5;
t2 = t32 + t34 + t252;
t1 = (t16 / 0.2e1 + t44 / 0.2e1) * t248 + (-t43 / 0.2e1 + t15 / 0.2e1) * t247 + t5;
t9 = [qJD(2) * t21 + qJD(3) * t24 + qJD(4) * t18 + qJD(5) * t33, qJD(1) * t21 + qJD(3) * t85 + qJD(4) * t29 + t99, qJD(1) * t24 + qJD(2) * t85 + qJD(4) * t25 + qJD(5) * t64, t18 * qJD(1) + t29 * qJD(2) + t25 * qJD(3) + t2 * qJD(5) + ((-t100 * t144 + t101 * t146 + (-t139 + t145) * t272) * t354 + (t113 * t209 + t270 * t367) * t356) * t357 + (t15 * t378 + (-t285 * t229 - t287 * t230) * t328 + t252 - (t244 / 0.2e1 + t243 / 0.2e1) * t263 + (t284 * t229 + t286 * t230 + t43) * t330 + (t16 + t44) * t329) * qJD(4), t33 * qJD(1) + t64 * qJD(3) + t2 * qJD(4) + t252 * qJD(5) + (qJD(2) * t333 + (t105 * t201 + t273) * qJD(5)) * m(6); t86 * qJD(3) + t30 * qJD(4) + t99 + (-t340 / 0.4e1 - t327 / 0.4e1 - t352 / 0.4e1 - t349 / 0.4e1) * t358, 0, t292, t30 * qJD(1) + (t260 * t354 - t276 * t356) * t357 - t112, -t112 + 0.2e1 * (t105 * qJD(1) / 0.4e1 - t119 * qJD(4) / 0.2e1) * m(6); -t86 * qJD(2) + t26 * qJD(4) - t63 * qJD(5) + (-t341 / 0.4e1 - t353 / 0.4e1 - t350 / 0.4e1) * t358, -t292, 0, t383 + m(6) * (-t144 * t248 + t146 * t247) * qJD(4), -t293; (t250 + (-t266 + t207) * t229 / 0.2e1 + t375) * qJD(1) + t31 * qJD(2) - t26 * qJD(3) + t1 * qJD(4) + t4 * qJD(5) + (-t343 / 0.4e1 - t351 / 0.4e1) * t358, t31 * qJD(1), -t383, t1 * qJD(1) + (m(5) * (-(-t248 * (t270 * t248 - t240) + (-t248 * rSges(5,3) - t270 * t247) * t247) * t114 - t209 * t276) + (t244 * t185 + (-t248 * t186 - t373) * t247) * t328 + (-t243 * t186 + (t247 * t185 + t373) * t248) * t330 + m(6) * (-t144 * t145 - t146 * t272 + t50 * t87) + t319) * qJD(4) + t371, t4 * qJD(1) + t6 * qJD(4) + t371; (t250 - t342) * qJD(1) + qJD(2) * t278 + t63 * qJD(3) + t3 * qJD(4) + t5 * qJD(5), qJD(1) * t278, t293, t3 * qJD(1) + ((t260 * t201 + t98 * t87) * m(6) + t319) * qJD(4) + t7, qJD(1) * t5 + qJD(4) * t8 + t7;];
Cq = t9;
