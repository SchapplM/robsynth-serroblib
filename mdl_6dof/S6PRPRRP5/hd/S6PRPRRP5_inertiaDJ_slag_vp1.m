% Calculate time derivative of joint inertia matrix for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:05
% EndTime: 2019-03-08 20:14:29
% DurationCPUTime: 16.68s
% Computational Cost: add. (53289->1119), mult. (149990->1560), div. (0->0), fcn. (170506->10), ass. (0->405)
t469 = rSges(7,3) + qJ(6);
t356 = sin(pkin(6));
t467 = t356 / 0.2e1;
t358 = cos(pkin(6));
t465 = t358 / 0.2e1;
t355 = sin(pkin(10));
t442 = t355 * t356;
t357 = cos(pkin(10));
t441 = t356 * t357;
t361 = sin(qJ(2));
t363 = cos(qJ(2));
t436 = t363 * t358;
t383 = t355 * t361 - t357 * t436;
t456 = cos(qJ(4));
t393 = t356 * t456;
t455 = sin(qJ(4));
t310 = -t357 * t393 + t383 * t455;
t392 = t356 * t455;
t344 = t358 * t456 - t363 * t392;
t370 = -t358 * t455 - t363 * t393;
t414 = qJD(2) * t361;
t388 = t356 * t414;
t311 = qJD(4) * t370 + t388 * t455;
t312 = qJD(4) * t344 - t388 * t456;
t413 = qJD(2) * t363;
t387 = t356 * t413;
t239 = Icges(5,5) * t311 - Icges(5,6) * t312 + Icges(5,3) * t387;
t240 = Icges(5,4) * t311 - Icges(5,2) * t312 + Icges(5,6) * t387;
t241 = Icges(5,1) * t311 - Icges(5,4) * t312 + Icges(5,5) * t387;
t341 = t355 * t436 + t357 * t361;
t307 = -t341 * t456 + t355 * t392;
t386 = t357 * t413;
t389 = t355 * t414;
t331 = -t358 * t389 + t386;
t257 = -qJD(4) * t307 + t331 * t455;
t308 = t341 * t455 + t355 * t393;
t258 = qJD(4) * t308 - t331 * t456;
t440 = t356 * t361;
t278 = Icges(5,5) * t344 + Icges(5,6) * t370 + Icges(5,3) * t440;
t279 = Icges(5,4) * t344 + Icges(5,2) * t370 + Icges(5,6) * t440;
t280 = Icges(5,1) * t344 + Icges(5,4) * t370 + Icges(5,5) * t440;
t330 = t341 * qJD(2);
t437 = t358 * t361;
t342 = -t355 * t437 + t357 * t363;
t103 = t239 * t342 - t240 * t307 + t241 * t308 + t257 * t280 - t258 * t279 - t278 * t330;
t219 = Icges(5,5) * t308 - Icges(5,6) * t307 + Icges(5,3) * t342;
t221 = Icges(5,4) * t308 - Icges(5,2) * t307 + Icges(5,6) * t342;
t223 = Icges(5,1) * t308 - Icges(5,4) * t307 + Icges(5,5) * t342;
t145 = t219 * t342 - t221 * t307 + t223 * t308;
t309 = t357 * t392 + t383 * t456;
t340 = t355 * t363 + t357 * t437;
t220 = Icges(5,5) * t310 + Icges(5,6) * t309 + Icges(5,3) * t340;
t222 = Icges(5,4) * t310 + Icges(5,2) * t309 + Icges(5,6) * t340;
t224 = Icges(5,1) * t310 + Icges(5,4) * t309 + Icges(5,5) * t340;
t146 = t220 * t342 - t222 * t307 + t224 * t308;
t158 = t278 * t342 - t279 * t307 + t280 * t308;
t328 = -t358 * t386 + t389;
t360 = sin(qJ(5));
t362 = cos(qJ(5));
t313 = -t344 * t360 + t362 * t440;
t402 = t360 * t440;
t314 = t344 * t362 + t402;
t227 = Icges(7,5) * t314 + Icges(7,6) * t313 - Icges(7,3) * t370;
t229 = Icges(7,4) * t314 + Icges(7,2) * t313 - Icges(7,6) * t370;
t231 = Icges(7,1) * t314 + Icges(7,4) * t313 - Icges(7,5) * t370;
t261 = -t308 * t360 + t342 * t362;
t443 = t342 * t360;
t262 = t308 * t362 + t443;
t117 = t227 * t307 + t229 * t261 + t231 * t262;
t196 = -qJD(5) * t262 - t257 * t360 - t330 * t362;
t366 = qJD(5) * t261 - t330 * t360;
t197 = t257 * t362 + t366;
t128 = Icges(7,5) * t197 + Icges(7,6) * t196 + Icges(7,3) * t258;
t132 = Icges(7,4) * t197 + Icges(7,2) * t196 + Icges(7,6) * t258;
t136 = Icges(7,1) * t197 + Icges(7,4) * t196 + Icges(7,5) * t258;
t178 = Icges(7,5) * t262 + Icges(7,6) * t261 + Icges(7,3) * t307;
t182 = Icges(7,4) * t262 + Icges(7,2) * t261 + Icges(7,6) * t307;
t186 = Icges(7,1) * t262 + Icges(7,4) * t261 + Icges(7,5) * t307;
t29 = t128 * t307 + t132 * t261 + t136 * t262 + t178 * t258 + t182 * t196 + t186 * t197;
t329 = t340 * qJD(2);
t259 = qJD(4) * t309 + t329 * t455;
t444 = t340 * t360;
t264 = t310 * t362 + t444;
t198 = -qJD(5) * t264 - t259 * t360 - t328 * t362;
t263 = -t310 * t360 + t340 * t362;
t367 = qJD(5) * t263 - t328 * t360;
t199 = t259 * t362 + t367;
t260 = -qJD(4) * t310 + t329 * t456;
t129 = Icges(7,5) * t199 + Icges(7,6) * t198 - Icges(7,3) * t260;
t133 = Icges(7,4) * t199 + Icges(7,2) * t198 - Icges(7,6) * t260;
t137 = Icges(7,1) * t199 + Icges(7,4) * t198 - Icges(7,5) * t260;
t179 = Icges(7,5) * t264 + Icges(7,6) * t263 - Icges(7,3) * t309;
t183 = Icges(7,4) * t264 + Icges(7,2) * t263 - Icges(7,6) * t309;
t187 = Icges(7,1) * t264 + Icges(7,4) * t263 - Icges(7,5) * t309;
t30 = t129 * t307 + t133 * t261 + t137 * t262 + t179 * t258 + t183 * t196 + t187 * t197;
t237 = -qJD(5) * t314 - t311 * t360 + t362 * t387;
t364 = qJD(5) * t313 + t360 * t387;
t238 = t311 * t362 + t364;
t163 = Icges(7,5) * t238 + Icges(7,6) * t237 + Icges(7,3) * t312;
t165 = Icges(7,4) * t238 + Icges(7,2) * t237 + Icges(7,6) * t312;
t167 = Icges(7,1) * t238 + Icges(7,4) * t237 + Icges(7,5) * t312;
t65 = t163 * t307 + t165 * t261 + t167 * t262 + t196 * t229 + t197 * t231 + t227 * t258;
t93 = t178 * t307 + t182 * t261 + t186 * t262;
t94 = t179 * t307 + t183 * t261 + t187 * t262;
t5 = t29 * t342 + t30 * t340 - t94 * t328 - t93 * t330 + (t117 * t413 + t361 * t65) * t356;
t228 = Icges(6,5) * t314 + Icges(6,6) * t313 - Icges(6,3) * t370;
t230 = Icges(6,4) * t314 + Icges(6,2) * t313 - Icges(6,6) * t370;
t232 = Icges(6,1) * t314 + Icges(6,4) * t313 - Icges(6,5) * t370;
t118 = t228 * t307 + t230 * t261 + t232 * t262;
t130 = Icges(6,5) * t197 + Icges(6,6) * t196 + Icges(6,3) * t258;
t134 = Icges(6,4) * t197 + Icges(6,2) * t196 + Icges(6,6) * t258;
t138 = Icges(6,1) * t197 + Icges(6,4) * t196 + Icges(6,5) * t258;
t180 = Icges(6,5) * t262 + Icges(6,6) * t261 + Icges(6,3) * t307;
t184 = Icges(6,4) * t262 + Icges(6,2) * t261 + Icges(6,6) * t307;
t188 = Icges(6,1) * t262 + Icges(6,4) * t261 + Icges(6,5) * t307;
t31 = t130 * t307 + t134 * t261 + t138 * t262 + t180 * t258 + t184 * t196 + t188 * t197;
t131 = Icges(6,5) * t199 + Icges(6,6) * t198 - Icges(6,3) * t260;
t135 = Icges(6,4) * t199 + Icges(6,2) * t198 - Icges(6,6) * t260;
t139 = Icges(6,1) * t199 + Icges(6,4) * t198 - Icges(6,5) * t260;
t181 = Icges(6,5) * t264 + Icges(6,6) * t263 - Icges(6,3) * t309;
t185 = Icges(6,4) * t264 + Icges(6,2) * t263 - Icges(6,6) * t309;
t189 = Icges(6,1) * t264 + Icges(6,4) * t263 - Icges(6,5) * t309;
t32 = t131 * t307 + t135 * t261 + t139 * t262 + t181 * t258 + t185 * t196 + t189 * t197;
t164 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t312;
t166 = Icges(6,4) * t238 + Icges(6,2) * t237 + Icges(6,6) * t312;
t168 = Icges(6,1) * t238 + Icges(6,4) * t237 + Icges(6,5) * t312;
t66 = t164 * t307 + t166 * t261 + t168 * t262 + t196 * t230 + t197 * t232 + t228 * t258;
t95 = t180 * t307 + t184 * t261 + t188 * t262;
t96 = t181 * t307 + t185 * t261 + t189 * t262;
t6 = t31 * t342 + t32 * t340 - t96 * t328 - t95 * t330 + (t118 * t413 + t361 * t66) * t356;
t200 = Icges(5,5) * t257 - Icges(5,6) * t258 - Icges(5,3) * t330;
t202 = Icges(5,4) * t257 - Icges(5,2) * t258 - Icges(5,6) * t330;
t204 = Icges(5,1) * t257 - Icges(5,4) * t258 - Icges(5,5) * t330;
t78 = t200 * t342 - t202 * t307 + t204 * t308 - t219 * t330 - t221 * t258 + t223 * t257;
t201 = Icges(5,5) * t259 + Icges(5,6) * t260 - Icges(5,3) * t328;
t203 = Icges(5,4) * t259 + Icges(5,2) * t260 - Icges(5,6) * t328;
t205 = Icges(5,1) * t259 + Icges(5,4) * t260 - Icges(5,5) * t328;
t79 = t201 * t342 - t203 * t307 + t205 * t308 - t220 * t330 - t222 * t258 + t224 * t257;
t464 = -t145 * t330 - t146 * t328 + t79 * t340 + t78 * t342 + (t103 * t361 + t158 * t413) * t356 + t5 + t6;
t104 = t239 * t340 + t240 * t309 + t241 * t310 + t259 * t280 + t260 * t279 - t278 * t328;
t147 = t219 * t340 + t221 * t309 + t223 * t310;
t148 = t220 * t340 + t222 * t309 + t224 * t310;
t159 = t278 * t340 + t279 * t309 + t280 * t310;
t119 = -t227 * t309 + t229 * t263 + t231 * t264;
t33 = -t128 * t309 + t132 * t263 + t136 * t264 - t178 * t260 + t182 * t198 + t186 * t199;
t34 = -t129 * t309 + t133 * t263 + t137 * t264 - t179 * t260 + t183 * t198 + t187 * t199;
t67 = -t163 * t309 + t165 * t263 + t167 * t264 + t198 * t229 + t199 * t231 - t227 * t260;
t97 = -t178 * t309 + t182 * t263 + t186 * t264;
t98 = -t179 * t309 + t183 * t263 + t187 * t264;
t7 = -t98 * t328 + t33 * t342 - t97 * t330 + t34 * t340 + (t119 * t413 + t361 * t67) * t356;
t100 = -t181 * t309 + t185 * t263 + t189 * t264;
t120 = -t228 * t309 + t230 * t263 + t232 * t264;
t35 = -t130 * t309 + t134 * t263 + t138 * t264 - t180 * t260 + t184 * t198 + t188 * t199;
t36 = -t131 * t309 + t135 * t263 + t139 * t264 - t181 * t260 + t185 * t198 + t189 * t199;
t68 = -t164 * t309 + t166 * t263 + t168 * t264 + t198 * t230 + t199 * t232 - t228 * t260;
t99 = -t180 * t309 + t184 * t263 + t188 * t264;
t8 = -t100 * t328 - t99 * t330 + t36 * t340 + t35 * t342 + (t120 * t413 + t361 * t68) * t356;
t80 = t200 * t340 + t202 * t309 + t204 * t310 - t219 * t328 + t221 * t260 + t223 * t259;
t81 = t201 * t340 + t203 * t309 + t205 * t310 - t220 * t328 + t222 * t260 + t224 * t259;
t463 = -t147 * t330 - t148 * t328 + t81 * t340 + t80 * t342 + (t104 * t361 + t159 * t413) * t356 + t7 + t8;
t462 = 2 * m(5);
t461 = 2 * m(6);
t460 = 2 * m(7);
t459 = t358 ^ 2;
t454 = pkin(8) * t328;
t453 = pkin(8) * t330;
t452 = pkin(5) * t362;
t450 = Icges(3,4) * t361;
t449 = Icges(3,4) * t363;
t448 = Icges(4,6) * t361;
t447 = Icges(4,6) * t363;
t446 = t328 * t159;
t445 = t330 * t158;
t288 = -Icges(3,5) * t328 - Icges(3,6) * t329;
t439 = t357 * t288;
t290 = Icges(4,4) * t328 + Icges(4,5) * t329;
t438 = t357 * t290;
t435 = rSges(7,1) * t197 + rSges(7,2) * t196 + pkin(5) * t366 + qJD(6) * t307 + t257 * t452 + t258 * t469;
t434 = rSges(7,1) * t199 + rSges(7,2) * t198 + pkin(5) * t367 - qJD(6) * t309 + t259 * t452 - t260 * t469;
t143 = rSges(6,1) * t199 + rSges(6,2) * t198 - rSges(6,3) * t260;
t214 = pkin(4) * t259 - pkin(9) * t260;
t433 = -t143 - t214;
t432 = rSges(7,1) * t238 + rSges(7,2) * t237 + pkin(5) * t364 - t370 * qJD(6) + t311 * t452 + t312 * t469;
t170 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t312;
t253 = t311 * pkin(4) + t312 * pkin(9);
t431 = -t170 - t253;
t430 = rSges(7,1) * t262 + rSges(7,2) * t261 + pkin(5) * t443 + t307 * t469 + t308 * t452;
t429 = rSges(7,1) * t264 + rSges(7,2) * t263 + pkin(5) * t444 - t309 * t469 + t310 * t452;
t191 = rSges(6,1) * t262 + rSges(6,2) * t261 + rSges(6,3) * t307;
t251 = pkin(4) * t308 + pkin(9) * t307;
t428 = t191 + t251;
t193 = rSges(6,1) * t264 + rSges(6,2) * t263 - rSges(6,3) * t309;
t252 = pkin(4) * t310 - pkin(9) * t309;
t427 = -t193 - t252;
t426 = t342 * t214 - t330 * t252;
t213 = pkin(4) * t257 + pkin(9) * t258;
t425 = t213 * t440 + t251 * t387;
t255 = -pkin(2) * t330 + qJ(3) * t331 + qJD(3) * t341;
t250 = t358 * t255;
t424 = t358 * t213 + t250;
t423 = rSges(7,1) * t314 + rSges(7,2) * t313 + pkin(5) * t402 + t344 * t452 - t370 * t469;
t235 = rSges(6,1) * t314 + rSges(6,2) * t313 - rSges(6,3) * t370;
t303 = pkin(4) * t344 - pkin(9) * t370;
t422 = t235 + t303;
t421 = t340 * t253 - t328 * t303;
t254 = -t328 * pkin(2) + t329 * qJ(3) + qJD(3) * t383;
t420 = t254 * t442 + t255 * t441;
t301 = t340 * pkin(2) + qJ(3) * t383;
t302 = pkin(2) * t342 + qJ(3) * t341;
t419 = t301 * t442 + t302 * t441;
t300 = t358 * t302;
t315 = pkin(3) * t442 + pkin(8) * t342;
t418 = t358 * t315 + t300;
t316 = -pkin(3) * t441 + pkin(8) * t340;
t417 = -t301 - t316;
t345 = (pkin(2) * t361 - qJ(3) * t363) * t356;
t416 = -pkin(3) * t358 - pkin(8) * t440 - t345;
t415 = qJD(2) * t356;
t1 = t117 * t312 + t258 * t93 - t260 * t94 + t29 * t307 - t30 * t309 - t370 * t65;
t2 = t118 * t312 + t258 * t95 - t260 * t96 + t307 * t31 - t309 * t32 - t370 * t66;
t412 = t1 / 0.2e1 + t2 / 0.2e1;
t3 = t119 * t312 + t258 * t97 - t260 * t98 + t307 * t33 - t309 * t34 - t370 * t67;
t4 = -t100 * t260 + t120 * t312 + t258 * t99 + t307 * t35 - t309 * t36 - t370 * t68;
t411 = t4 / 0.2e1 + t3 / 0.2e1;
t108 = -t180 * t370 + t184 * t313 + t188 * t314;
t109 = -t181 * t370 + t185 * t313 + t189 * t314;
t150 = -t228 * t370 + t230 * t313 + t232 * t314;
t42 = -t130 * t370 + t134 * t313 + t138 * t314 + t180 * t312 + t184 * t237 + t188 * t238;
t43 = -t131 * t370 + t135 * t313 + t139 * t314 + t181 * t312 + t185 * t237 + t189 * t238;
t74 = -t164 * t370 + t166 * t313 + t168 * t314 + t228 * t312 + t230 * t237 + t232 * t238;
t10 = t108 * t258 - t109 * t260 + t150 * t312 + t307 * t42 - t309 * t43 - t370 * t74;
t106 = -t178 * t370 + t182 * t313 + t186 * t314;
t107 = -t179 * t370 + t183 * t313 + t187 * t314;
t149 = -t227 * t370 + t229 * t313 + t231 * t314;
t40 = -t128 * t370 + t132 * t313 + t136 * t314 + t178 * t312 + t182 * t237 + t186 * t238;
t41 = -t129 * t370 + t133 * t313 + t137 * t314 + t179 * t312 + t183 * t237 + t187 * t238;
t73 = -t163 * t370 + t165 * t313 + t167 * t314 + t227 * t312 + t229 * t237 + t231 * t238;
t9 = t106 * t258 - t107 * t260 + t149 * t312 + t307 * t40 - t309 * t41 - t370 * t73;
t410 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = -t106 * t330 - t107 * t328 + t41 * t340 + t40 * t342 + (t149 * t413 + t361 * t73) * t356;
t12 = -t108 * t330 - t109 * t328 + t43 * t340 + t42 * t342 + (t150 * t413 + t361 * t74) * t356;
t409 = t12 / 0.2e1 + t11 / 0.2e1;
t13 = t358 * t65 + (t29 * t355 - t30 * t357) * t356;
t14 = t358 * t66 + (t31 * t355 - t32 * t357) * t356;
t408 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t358 * t67 + (t33 * t355 - t34 * t357) * t356;
t16 = t358 * t68 + (t35 * t355 - t357 * t36) * t356;
t407 = -t15 / 0.2e1 - t16 / 0.2e1;
t17 = t358 * t73 + (t355 * t40 - t357 * t41) * t356;
t18 = t358 * t74 + (t355 * t42 - t357 * t43) * t356;
t406 = t17 / 0.2e1 + t18 / 0.2e1;
t405 = ((-t94 - t96) * t357 + (t93 + t95) * t355) * t467 + (t118 + t117) * t465;
t404 = -((-t100 - t98) * t357 + (t97 + t99) * t355) * t356 / 0.2e1 - (t119 + t120) * t358 / 0.2e1;
t403 = ((-t107 - t109) * t357 + (t106 + t108) * t355) * t467 + (t149 + t150) * t465;
t401 = -t214 - t434;
t400 = -t253 - t432;
t399 = t251 + t430;
t398 = -t252 - t429;
t397 = t303 + t423;
t396 = t358 * t251 + t418;
t395 = -t252 + t417;
t394 = -t303 + t416;
t390 = t356 ^ 2 * t413;
t317 = (-qJD(3) * t363 + (pkin(2) * t363 + qJ(3) * t361) * qJD(2)) * t356;
t385 = (-t317 - (-rSges(4,2) * t363 + rSges(4,3) * t361) * t415) * t356;
t384 = (-t358 * rSges(4,1) - (-rSges(4,2) * t361 - rSges(4,3) * t363) * t356 - t345) * t356;
t382 = pkin(8) * t390;
t381 = t213 * t441 + t214 * t442 + t420;
t380 = t315 * t441 + t316 * t442 + t419;
t377 = t435 - t453;
t281 = rSges(5,1) * t344 + rSges(5,2) * t370 + rSges(5,3) * t440;
t376 = (-t281 + t416) * t356;
t375 = (-t235 + t394) * t356;
t374 = t251 * t441 + t252 * t442 + t380;
t372 = (-t328 * t355 - t330 * t357) * pkin(8);
t371 = (t394 - t423) * t356;
t242 = t311 * rSges(5,1) - t312 * rSges(5,2) + rSges(5,3) * t387;
t369 = -t382 + (-t242 - t317) * t356;
t368 = -t382 + (-t317 + t431) * t356;
t365 = -t382 + (-t317 + t400) * t356;
t339 = (rSges(3,1) * t363 - rSges(3,2) * t361) * t415;
t337 = (Icges(3,1) * t363 - t450) * t415;
t336 = (-Icges(3,2) * t361 + t449) * t415;
t335 = (-Icges(4,4) * t363 + Icges(4,5) * t361) * t415;
t334 = (Icges(3,5) * t363 - Icges(3,6) * t361) * t415;
t333 = (-Icges(4,2) * t363 + t448) * t415;
t332 = (Icges(4,3) * t361 - t447) * t415;
t325 = t358 * rSges(3,3) + (rSges(3,1) * t361 + rSges(3,2) * t363) * t356;
t324 = Icges(4,4) * t358 + (-Icges(4,2) * t361 - t447) * t356;
t323 = Icges(4,5) * t358 + (-Icges(4,3) * t363 - t448) * t356;
t322 = Icges(3,5) * t358 + (Icges(3,1) * t361 + t449) * t356;
t321 = Icges(3,6) * t358 + (Icges(3,2) * t363 + t450) * t356;
t318 = t358 * t454;
t299 = -rSges(3,1) * t330 - rSges(3,2) * t331;
t298 = rSges(4,2) * t330 + rSges(4,3) * t331;
t297 = -rSges(3,1) * t328 - rSges(3,2) * t329;
t296 = rSges(4,2) * t328 + rSges(4,3) * t329;
t295 = -Icges(3,1) * t330 - Icges(3,4) * t331;
t294 = -Icges(3,1) * t328 - Icges(3,4) * t329;
t293 = -Icges(3,4) * t330 - Icges(3,2) * t331;
t292 = -Icges(3,4) * t328 - Icges(3,2) * t329;
t291 = Icges(4,4) * t330 + Icges(4,5) * t331;
t289 = -Icges(3,5) * t330 - Icges(3,6) * t331;
t287 = Icges(4,2) * t330 + Icges(4,6) * t331;
t286 = Icges(4,2) * t328 + Icges(4,6) * t329;
t285 = Icges(4,6) * t330 + Icges(4,3) * t331;
t284 = Icges(4,6) * t328 + Icges(4,3) * t329;
t277 = rSges(3,1) * t342 - rSges(3,2) * t341 + rSges(3,3) * t442;
t276 = t340 * rSges(3,1) - rSges(3,2) * t383 - rSges(3,3) * t441;
t275 = -rSges(4,1) * t441 - t340 * rSges(4,2) + rSges(4,3) * t383;
t274 = rSges(4,1) * t442 - rSges(4,2) * t342 + rSges(4,3) * t341;
t273 = Icges(3,1) * t342 - Icges(3,4) * t341 + Icges(3,5) * t442;
t272 = Icges(3,1) * t340 - Icges(3,4) * t383 - Icges(3,5) * t441;
t271 = Icges(3,4) * t342 - Icges(3,2) * t341 + Icges(3,6) * t442;
t270 = Icges(3,4) * t340 - Icges(3,2) * t383 - Icges(3,6) * t441;
t269 = -Icges(4,4) * t441 - Icges(4,2) * t340 + Icges(4,6) * t383;
t268 = Icges(4,4) * t442 - Icges(4,2) * t342 + Icges(4,6) * t341;
t267 = -Icges(4,5) * t441 - Icges(4,6) * t340 + Icges(4,3) * t383;
t266 = Icges(4,5) * t442 - Icges(4,6) * t342 + Icges(4,3) * t341;
t265 = t340 * t303;
t246 = t251 * t440;
t233 = t342 * t252;
t226 = rSges(5,1) * t310 + rSges(5,2) * t309 + rSges(5,3) * t340;
t225 = rSges(5,1) * t308 - rSges(5,2) * t307 + rSges(5,3) * t342;
t215 = (t297 * t355 + t299 * t357) * t356;
t207 = rSges(5,1) * t259 + rSges(5,2) * t260 - rSges(5,3) * t328;
t206 = rSges(5,1) * t257 - rSges(5,2) * t258 - rSges(5,3) * t330;
t195 = (-t275 - t301) * t358 + t357 * t384;
t194 = t274 * t358 + t355 * t384 + t300;
t177 = t225 * t440 - t281 * t342;
t176 = -t226 * t440 + t281 * t340;
t173 = (-t254 - t296) * t358 + t357 * t385;
t172 = t298 * t358 + t355 * t385 + t250;
t171 = t278 * t440 + t279 * t370 + t280 * t344;
t162 = (t274 * t357 + t275 * t355) * t356 + t419;
t160 = -t225 * t340 + t226 * t342;
t157 = (t296 * t355 + t298 * t357) * t356 + t420;
t156 = (-t226 + t417) * t358 + t357 * t376;
t155 = t225 * t358 + t355 * t376 + t418;
t154 = t193 * t370 - t235 * t309;
t153 = -t191 * t370 - t235 * t307;
t152 = t220 * t440 + t222 * t370 + t224 * t344;
t151 = t219 * t440 + t221 * t370 + t223 * t344;
t144 = (t225 * t357 + t226 * t355) * t356 + t380;
t141 = rSges(6,1) * t197 + rSges(6,2) * t196 + rSges(6,3) * t258;
t125 = t191 * t309 + t193 * t307;
t124 = t318 + (-t207 - t254) * t358 + t369 * t357;
t123 = t250 + (t206 - t453) * t358 + t369 * t355;
t122 = t191 * t440 - t342 * t422 + t246;
t121 = t235 * t340 + t427 * t440 + t265;
t116 = -t342 * t242 + t330 * t281 + (t206 * t361 + t225 * t413) * t356;
t115 = t340 * t242 - t328 * t281 + (-t207 * t361 - t226 * t413) * t356;
t114 = (-t193 + t395) * t358 + t357 * t375;
t113 = t191 * t358 + t355 * t375 + t396;
t112 = t370 * t240 + t344 * t241 - t312 * t279 + t311 * t280 + (t239 * t361 + t278 * t413) * t356;
t111 = (t206 * t357 + t207 * t355 + t372) * t356 + t420;
t110 = t193 * t342 - t340 * t428 + t233;
t105 = -t206 * t340 + t207 * t342 + t225 * t328 - t226 * t330;
t102 = -t309 * t423 + t370 * t429;
t101 = -t307 * t423 - t370 * t430;
t92 = (t191 * t357 + t193 * t355) * t356 + t374;
t91 = -t342 * t397 + t430 * t440 + t246;
t90 = t340 * t423 + t398 * t440 + t265;
t89 = (t395 - t429) * t358 + t357 * t371;
t88 = t355 * t371 + t358 * t430 + t396;
t87 = t307 * t429 + t309 * t430;
t86 = t370 * t203 + t344 * t205 - t312 * t222 + t311 * t224 + (t201 * t361 + t220 * t413) * t356;
t85 = t370 * t202 + t344 * t204 - t312 * t221 + t311 * t223 + (t200 * t361 + t219 * t413) * t356;
t84 = -t340 * t399 + t342 * t429 + t233;
t83 = t318 + (-t254 + t433) * t358 + t368 * t357;
t82 = (t141 - t453) * t358 + t368 * t355 + t424;
t77 = t143 * t370 - t170 * t309 - t193 * t312 - t235 * t260;
t76 = -t141 * t370 - t170 * t307 + t191 * t312 - t235 * t258;
t75 = (t355 * t429 + t357 * t430) * t356 + t374;
t72 = (t141 * t357 + t143 * t355 + t372) * t356 + t381;
t71 = (t141 * t361 + t191 * t413) * t356 + t431 * t342 + t422 * t330 + t425;
t70 = t340 * t170 - t328 * t235 + (t361 * t433 + t413 * t427) * t356 + t421;
t69 = t141 * t309 + t143 * t307 + t191 * t260 + t193 * t258;
t62 = t108 * t342 + t109 * t340 + t150 * t440;
t61 = t106 * t342 + t107 * t340 + t149 * t440;
t60 = t318 + (-t254 + t401) * t358 + t365 * t357;
t59 = t355 * t365 + t358 * t377 + t424;
t58 = t108 * t307 - t109 * t309 - t150 * t370;
t57 = t106 * t307 - t107 * t309 - t149 * t370;
t56 = t143 * t342 - t193 * t330 + (-t141 - t213) * t340 + t428 * t328 + t426;
t51 = t100 * t340 + t120 * t440 + t342 * t99;
t50 = t119 * t440 + t340 * t98 + t342 * t97;
t49 = t118 * t440 + t340 * t96 + t342 * t95;
t48 = t117 * t440 + t340 * t94 + t342 * t93;
t47 = -t100 * t309 - t120 * t370 + t307 * t99;
t46 = -t119 * t370 + t307 * t97 - t309 * t98;
t45 = -t118 * t370 + t307 * t95 - t309 * t96;
t44 = -t117 * t370 + t307 * t93 - t309 * t94;
t39 = -t260 * t423 - t309 * t432 - t312 * t429 + t370 * t434;
t38 = -t258 * t423 - t307 * t432 + t312 * t430 - t370 * t435;
t37 = (t377 * t357 + (t434 - t454) * t355) * t356 + t381;
t28 = t400 * t342 + t397 * t330 + (t361 * t435 + t413 * t430) * t356 + t425;
t27 = t432 * t340 - t423 * t328 + (t361 * t401 + t398 * t413) * t356 + t421;
t26 = t112 * t358 + (t355 * t85 - t357 * t86) * t356;
t25 = t258 * t429 + t260 * t430 + t307 * t434 + t309 * t435;
t24 = t104 * t358 + (t355 * t80 - t357 * t81) * t356;
t23 = t103 * t358 + (t355 * t78 - t357 * t79) * t356;
t22 = t434 * t342 - t429 * t330 + (-t213 - t435) * t340 + t399 * t328 + t426;
t21 = -t151 * t330 - t152 * t328 + t86 * t340 + t85 * t342 + (t112 * t361 + t171 * t413) * t356;
t19 = [0; m(3) * t215 + m(4) * t157 + m(5) * t111 + m(6) * t72 + m(7) * t37; -t15 * t441 - t16 * t441 + t14 * t442 + t13 * t442 - t24 * t441 + t23 * t442 - ((-t329 * t271 - t328 * t273 - t289 * t441 - t293 * t383 + t340 * t295) * t442 - (-t329 * t270 - t328 * t272 - t292 * t383 + t340 * t294 - t356 * t439) * t441 + (-t329 * t321 - t328 * t322 - t334 * t441 - t336 * t383 + t340 * t337) * t358) * t441 - ((t329 * t266 + t328 * t268 + t285 * t383 - t340 * t287 - t291 * t441) * t442 - (t329 * t267 + t328 * t269 + t284 * t383 - t340 * t286 - t356 * t438) * t441 + (t329 * t323 + t328 * t324 + t332 * t383 - t340 * t333 - t335 * t441) * t358) * t441 + ((-t271 * t331 - t273 * t330 + t289 * t442 - t293 * t341 + t295 * t342) * t442 - (-t270 * t331 - t272 * t330 + t288 * t442 - t292 * t341 + t294 * t342) * t441 + (-t321 * t331 - t322 * t330 + t334 * t442 - t336 * t341 + t337 * t342) * t358) * t442 + ((t266 * t331 + t268 * t330 + t285 * t341 - t287 * t342 + t291 * t442) * t442 - (t267 * t331 + t269 * t330 + t284 * t341 - t286 * t342 + t290 * t442) * t441 + (t323 * t331 + t324 * t330 + t332 * t341 - t333 * t342 + t335 * t442) * t358) * t442 + (t37 * t75 + t59 * t88 + t60 * t89) * t460 + (t113 * t82 + t114 * t83 + t72 * t92) * t461 + t358 * t17 + t358 * t18 + t358 * t26 + (t111 * t144 + t123 * t155 + t124 * t156) * t462 + 0.2e1 * m(4) * (t157 * t162 + t172 * t194 + t173 * t195) + 0.2e1 * m(3) * ((-t276 * t358 - t325 * t441) * (-t297 * t358 - t339 * t441) + (t277 * t358 - t325 * t442) * (t299 * t358 - t339 * t442) + (t276 * t355 + t277 * t357) * t356 * t215) + t358 * (t459 * t335 + (((-t285 * t363 - t287 * t361) * t355 - (-t284 * t363 - t286 * t361) * t357 + ((t266 * t361 - t268 * t363) * t355 - (t267 * t361 - t269 * t363) * t357) * qJD(2)) * t356 + (-t438 + t291 * t355 - t332 * t363 - t333 * t361 + (t323 * t361 - t324 * t363) * qJD(2)) * t358) * t356) + t358 * (t459 * t334 + (((t293 * t363 + t295 * t361) * t355 - (t292 * t363 + t294 * t361) * t357 + ((-t271 * t361 + t273 * t363) * t355 - (-t270 * t361 + t272 * t363) * t357) * qJD(2)) * t356 + (-t439 + t289 * t355 + t336 * t363 + t337 * t361 + (-t321 * t361 + t322 * t363) * qJD(2)) * t358) * t356); (m(4) + m(5) + m(6) + m(7)) * t388; m(7) * (t331 * t89 + t341 * t60 + t329 * t88 + t383 * t59 + (-t363 * t37 + t414 * t75) * t356) + m(6) * (t331 * t114 + t341 * t83 + t329 * t113 + t383 * t82 + (-t363 * t72 + t414 * t92) * t356) + m(5) * (t331 * t156 + t341 * t124 + t329 * t155 + t383 * t123 + (-t111 * t363 + t144 * t414) * t356) + m(4) * (t331 * t195 + t341 * t173 + t329 * t194 + t383 * t172 + (-t157 * t363 + t162 * t414) * t356); 0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t329 * t383 + t341 * t331 - t361 * t390); m(5) * t105 + m(6) * t56 + m(7) * t22; -t405 * t330 + t404 * t328 + (t23 / 0.2e1 + t408) * t342 + (t24 / 0.2e1 - t407) * t340 + m(5) * (t105 * t144 + t111 * t160 + t115 * t156 + t116 * t155 + t123 * t177 + t124 * t176) + m(6) * (t110 * t72 + t113 * t71 + t114 * t70 + t121 * t83 + t122 * t82 + t56 * t92) + m(7) * (t22 * t75 + t27 * t89 + t28 * t88 + t37 * t84 + t59 * t91 + t60 * t90) + (-t445 / 0.2e1 - t446 / 0.2e1 + t21 / 0.2e1 + t409) * t358 + (-t330 * (t145 * t355 - t146 * t357) / 0.2e1 - t328 * (t147 * t355 - t148 * t357) / 0.2e1 + (t26 / 0.2e1 + t406) * t361 + (t171 * t465 + (t151 * t355 - t152 * t357) * t467 + t403) * t413) * t356 + t464 * t442 / 0.2e1 - t463 * t441 / 0.2e1; m(5) * (t115 * t341 + t176 * t331 + t116 * t383 + t177 * t329 + (-t105 * t363 + t160 * t414) * t356) + m(6) * (t70 * t341 + t121 * t331 + t71 * t383 + t122 * t329 + (t110 * t414 - t363 * t56) * t356) + m(7) * (t27 * t341 + t90 * t331 + t28 * t383 + t91 * t329 + (-t22 * t363 + t414 * t84) * t356); t464 * t342 + t463 * t340 + (-t145 * t342 - t146 * t340 - t48 - t49) * t330 + (-t147 * t342 - t148 * t340 - t50 - t51) * t328 + (t22 * t84 + t27 * t90 + t28 * t91) * t460 + (t110 * t56 + t121 * t70 + t122 * t71) * t461 + (t105 * t160 + t115 * t176 + t116 * t177) * t462 + ((t151 * t342 + t152 * t340 + t61 + t62) * t413 + (t171 * t387 + t11 + t12 + t21 - t445 - t446) * t361) * t356; m(6) * t69 + m(7) * t25; t410 * t358 - t406 * t370 + t403 * t312 + t407 * t309 + t408 * t307 + t404 * t260 + t405 * t258 + m(7) * (t101 * t59 + t102 * t60 + t25 * t75 + t37 * t87 + t38 * t88 + t39 * t89) + m(6) * (t113 * t76 + t114 * t77 + t125 * t72 + t153 * t82 + t154 * t83 + t69 * t92) + (t355 * t412 - t357 * t411) * t356; m(6) * (t77 * t341 + t154 * t331 + t76 * t383 + t153 * t329 + (t125 * t414 - t363 * t69) * t356) + m(7) * (t39 * t341 + t102 * t331 + t38 * t383 + t101 * t329 + (-t25 * t363 + t414 * t87) * t356); -t409 * t370 + t412 * t342 + t411 * t340 + (-t45 / 0.2e1 - t44 / 0.2e1) * t330 + (-t47 / 0.2e1 - t46 / 0.2e1) * t328 + (t62 / 0.2e1 + t61 / 0.2e1) * t312 + (-t8 / 0.2e1 - t7 / 0.2e1) * t309 + (t6 / 0.2e1 + t5 / 0.2e1) * t307 + (-t50 / 0.2e1 - t51 / 0.2e1) * t260 + (t48 / 0.2e1 + t49 / 0.2e1) * t258 + m(6) * (t110 * t69 + t121 * t77 + t122 * t76 + t125 * t56 + t153 * t71 + t154 * t70) + m(7) * (t101 * t28 + t102 * t27 + t22 * t87 + t25 * t84 + t38 * t91 + t39 * t90) + (t410 * t361 + (t57 / 0.2e1 + t58 / 0.2e1) * t413) * t356; -(t9 + t10) * t370 + (t57 + t58) * t312 + (-t4 - t3) * t309 + (t1 + t2) * t307 + (-t47 - t46) * t260 + (t44 + t45) * t258 + (t101 * t38 + t102 * t39 + t25 * t87) * t460 + (t125 * t69 + t153 * t76 + t154 * t77) * t461; m(7) * t312; m(7) * (t258 * t89 - t260 * t88 + t307 * t60 - t309 * t59 + t312 * t75 - t37 * t370); m(7) * (t258 * t341 + t307 * t331 - t260 * t383 - t309 * t329 + (-t312 * t363 - t370 * t414) * t356); m(7) * (-t22 * t370 + t258 * t90 - t260 * t91 + t27 * t307 - t28 * t309 + t312 * t84); m(7) * (-t101 * t260 + t102 * t258 - t25 * t370 + t307 * t39 - t309 * t38 + t312 * t87); (t258 * t307 + t260 * t309 - t312 * t370) * t460;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
