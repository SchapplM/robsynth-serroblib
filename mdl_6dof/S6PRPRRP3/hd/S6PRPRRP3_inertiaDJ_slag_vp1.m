% Calculate time derivative of joint inertia matrix for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:45
% EndTime: 2019-03-08 20:05:14
% DurationCPUTime: 18.19s
% Computational Cost: add. (83749->1133), mult. (152883->1602), div. (0->0), fcn. (173423->12), ass. (0->410)
t482 = rSges(7,3) + qJ(6);
t369 = sin(pkin(6));
t481 = t369 / 0.2e1;
t372 = cos(pkin(6));
t480 = t372 / 0.2e1;
t368 = sin(pkin(10));
t462 = t368 * t369;
t371 = cos(pkin(10));
t461 = t369 * t371;
t378 = cos(qJ(2));
t376 = sin(qJ(2));
t458 = t372 * t376;
t355 = t368 * t378 + t371 * t458;
t433 = pkin(11) + qJ(4);
t366 = sin(t433);
t406 = cos(t433);
t325 = t355 * t406 - t366 * t461;
t397 = t369 * t406;
t383 = -t355 * t366 - t371 * t397;
t457 = t372 * t378;
t387 = -t368 * t376 + t371 * t457;
t233 = Icges(5,5) * t325 + Icges(5,6) * t383 - Icges(5,3) * t387;
t235 = Icges(5,4) * t325 + Icges(5,2) * t383 - Icges(5,6) * t387;
t237 = Icges(5,1) * t325 + Icges(5,4) * t383 - Icges(5,5) * t387;
t149 = -t233 * t387 + t235 * t383 + t237 * t325;
t420 = t368 * t458;
t357 = t371 * t378 - t420;
t327 = t357 * t406 + t366 * t462;
t356 = t368 * t457 + t371 * t376;
t384 = -t357 * t366 + t368 * t397;
t234 = Icges(5,5) * t327 + Icges(5,6) * t384 + Icges(5,3) * t356;
t236 = Icges(5,4) * t327 + Icges(5,2) * t384 + Icges(5,6) * t356;
t238 = Icges(5,1) * t327 + Icges(5,4) * t384 + Icges(5,5) * t356;
t150 = -t234 * t387 + t236 * t383 + t238 * t325;
t460 = t369 * t376;
t341 = t366 * t460 - t372 * t406;
t342 = t372 * t366 + t376 * t397;
t459 = t369 * t378;
t292 = Icges(5,5) * t342 - Icges(5,6) * t341 - Icges(5,3) * t459;
t293 = Icges(5,4) * t342 - Icges(5,2) * t341 - Icges(5,6) * t459;
t294 = Icges(5,1) * t342 - Icges(5,4) * t341 - Icges(5,5) * t459;
t159 = -t292 * t387 + t293 * t383 + t294 * t325;
t345 = t355 * qJD(2);
t434 = qJD(2) * t378;
t347 = -qJD(2) * t420 + t371 * t434;
t435 = qJD(2) * t376;
t375 = sin(qJ(5));
t377 = cos(qJ(5));
t331 = -t342 * t375 - t377 * t459;
t419 = t375 * t459;
t388 = -t342 * t377 + t419;
t222 = -Icges(7,5) * t388 + Icges(7,6) * t331 + Icges(7,3) * t341;
t224 = -Icges(7,4) * t388 + Icges(7,2) * t331 + Icges(7,6) * t341;
t226 = -Icges(7,1) * t388 + Icges(7,4) * t331 + Icges(7,5) * t341;
t288 = -t325 * t375 - t377 * t387;
t465 = t387 * t375;
t289 = t325 * t377 - t465;
t118 = -t222 * t383 + t224 * t288 + t226 * t289;
t344 = t387 * qJD(2);
t283 = qJD(4) * t383 + t344 * t406;
t204 = -qJD(5) * t289 - t283 * t375 + t345 * t377;
t382 = qJD(5) * t288 + t345 * t375;
t205 = t283 * t377 + t382;
t282 = qJD(4) * t325 + t344 * t366;
t131 = Icges(7,5) * t205 + Icges(7,6) * t204 + Icges(7,3) * t282;
t135 = Icges(7,4) * t205 + Icges(7,2) * t204 + Icges(7,6) * t282;
t139 = Icges(7,1) * t205 + Icges(7,4) * t204 + Icges(7,5) * t282;
t180 = Icges(7,5) * t289 + Icges(7,6) * t288 - Icges(7,3) * t383;
t184 = Icges(7,4) * t289 + Icges(7,2) * t288 - Icges(7,6) * t383;
t188 = Icges(7,1) * t289 + Icges(7,4) * t288 - Icges(7,5) * t383;
t30 = -t131 * t383 + t135 * t288 + t139 * t289 + t180 * t282 + t184 * t204 + t188 * t205;
t346 = t356 * qJD(2);
t285 = qJD(4) * t384 - t346 * t406;
t464 = t356 * t375;
t291 = t327 * t377 + t464;
t206 = -qJD(5) * t291 - t285 * t375 + t347 * t377;
t290 = -t327 * t375 + t356 * t377;
t381 = qJD(5) * t290 + t347 * t375;
t207 = t285 * t377 + t381;
t284 = qJD(4) * t327 - t346 * t366;
t132 = Icges(7,5) * t207 + Icges(7,6) * t206 + Icges(7,3) * t284;
t136 = Icges(7,4) * t207 + Icges(7,2) * t206 + Icges(7,6) * t284;
t140 = Icges(7,1) * t207 + Icges(7,4) * t206 + Icges(7,5) * t284;
t181 = Icges(7,5) * t291 + Icges(7,6) * t290 - Icges(7,3) * t384;
t185 = Icges(7,4) * t291 + Icges(7,2) * t290 - Icges(7,6) * t384;
t189 = Icges(7,1) * t291 + Icges(7,4) * t290 - Icges(7,5) * t384;
t31 = -t132 * t383 + t136 * t288 + t140 * t289 + t181 * t282 + t185 * t204 + t189 * t205;
t323 = -qJD(4) * t341 + t397 * t434;
t407 = t369 * t435;
t246 = qJD(5) * t388 - t323 * t375 + t377 * t407;
t380 = qJD(5) * t331 + t375 * t407;
t247 = t323 * t377 + t380;
t322 = t366 * t369 * t434 + qJD(4) * t342;
t163 = Icges(7,5) * t247 + Icges(7,6) * t246 + Icges(7,3) * t322;
t165 = Icges(7,4) * t247 + Icges(7,2) * t246 + Icges(7,6) * t322;
t167 = Icges(7,1) * t247 + Icges(7,4) * t246 + Icges(7,5) * t322;
t65 = -t163 * t383 + t165 * t288 + t167 * t289 + t204 * t224 + t205 * t226 + t222 * t282;
t98 = -t180 * t383 + t184 * t288 + t188 * t289;
t99 = -t181 * t383 + t185 * t288 + t189 * t289;
t5 = -t30 * t387 + t31 * t356 + t98 * t345 + t99 * t347 + (t118 * t435 - t378 * t65) * t369;
t182 = Icges(6,5) * t289 + Icges(6,6) * t288 - Icges(6,3) * t383;
t186 = Icges(6,4) * t289 + Icges(6,2) * t288 - Icges(6,6) * t383;
t190 = Icges(6,1) * t289 + Icges(6,4) * t288 - Icges(6,5) * t383;
t100 = -t182 * t383 + t186 * t288 + t190 * t289;
t183 = Icges(6,5) * t291 + Icges(6,6) * t290 - Icges(6,3) * t384;
t187 = Icges(6,4) * t291 + Icges(6,2) * t290 - Icges(6,6) * t384;
t191 = Icges(6,1) * t291 + Icges(6,4) * t290 - Icges(6,5) * t384;
t101 = -t183 * t383 + t187 * t288 + t191 * t289;
t223 = -Icges(6,5) * t388 + Icges(6,6) * t331 + Icges(6,3) * t341;
t225 = -Icges(6,4) * t388 + Icges(6,2) * t331 + Icges(6,6) * t341;
t227 = -Icges(6,1) * t388 + Icges(6,4) * t331 + Icges(6,5) * t341;
t119 = -t223 * t383 + t225 * t288 + t227 * t289;
t133 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t282;
t137 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t282;
t141 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t282;
t32 = -t133 * t383 + t137 * t288 + t141 * t289 + t182 * t282 + t186 * t204 + t190 * t205;
t134 = Icges(6,5) * t207 + Icges(6,6) * t206 + Icges(6,3) * t284;
t138 = Icges(6,4) * t207 + Icges(6,2) * t206 + Icges(6,6) * t284;
t142 = Icges(6,1) * t207 + Icges(6,4) * t206 + Icges(6,5) * t284;
t33 = -t134 * t383 + t138 * t288 + t142 * t289 + t183 * t282 + t187 * t204 + t191 * t205;
t164 = Icges(6,5) * t247 + Icges(6,6) * t246 + Icges(6,3) * t322;
t166 = Icges(6,4) * t247 + Icges(6,2) * t246 + Icges(6,6) * t322;
t168 = Icges(6,1) * t247 + Icges(6,4) * t246 + Icges(6,5) * t322;
t66 = -t164 * t383 + t166 * t288 + t168 * t289 + t204 * t225 + t205 * t227 + t223 * t282;
t6 = t100 * t345 + t101 * t347 - t32 * t387 + t33 * t356 + (t119 * t435 - t378 * t66) * t369;
t196 = Icges(5,5) * t283 - Icges(5,6) * t282 + Icges(5,3) * t345;
t198 = Icges(5,4) * t283 - Icges(5,2) * t282 + Icges(5,6) * t345;
t200 = Icges(5,1) * t283 - Icges(5,4) * t282 + Icges(5,5) * t345;
t80 = -t196 * t387 + t198 * t383 + t200 * t325 + t233 * t345 - t235 * t282 + t237 * t283;
t197 = Icges(5,5) * t285 - Icges(5,6) * t284 + Icges(5,3) * t347;
t199 = Icges(5,4) * t285 - Icges(5,2) * t284 + Icges(5,6) * t347;
t201 = Icges(5,1) * t285 - Icges(5,4) * t284 + Icges(5,5) * t347;
t81 = -t197 * t387 + t199 * t383 + t201 * t325 + t234 * t345 - t236 * t282 + t238 * t283;
t248 = Icges(5,5) * t323 - Icges(5,6) * t322 + Icges(5,3) * t407;
t249 = Icges(5,4) * t323 - Icges(5,2) * t322 + Icges(5,6) * t407;
t250 = Icges(5,1) * t323 - Icges(5,4) * t322 + Icges(5,5) * t407;
t93 = -t248 * t387 + t249 * t383 + t250 * t325 - t282 * t293 + t283 * t294 + t292 * t345;
t479 = t149 * t345 + t150 * t347 - t80 * t387 + t81 * t356 + (t159 * t435 - t378 * t93) * t369 + t5 + t6;
t151 = t233 * t356 + t235 * t384 + t237 * t327;
t152 = t234 * t356 + t236 * t384 + t238 * t327;
t160 = t292 * t356 + t293 * t384 + t294 * t327;
t102 = -t180 * t384 + t184 * t290 + t188 * t291;
t103 = -t181 * t384 + t185 * t290 + t189 * t291;
t120 = -t222 * t384 + t224 * t290 + t226 * t291;
t34 = -t131 * t384 + t135 * t290 + t139 * t291 + t180 * t284 + t184 * t206 + t188 * t207;
t35 = -t132 * t384 + t136 * t290 + t140 * t291 + t181 * t284 + t185 * t206 + t189 * t207;
t67 = -t163 * t384 + t165 * t290 + t167 * t291 + t206 * t224 + t207 * t226 + t222 * t284;
t7 = t102 * t345 + t103 * t347 - t34 * t387 + t35 * t356 + (t120 * t435 - t378 * t67) * t369;
t104 = -t182 * t384 + t186 * t290 + t190 * t291;
t105 = -t183 * t384 + t187 * t290 + t191 * t291;
t121 = -t223 * t384 + t225 * t290 + t227 * t291;
t36 = -t133 * t384 + t137 * t290 + t141 * t291 + t182 * t284 + t186 * t206 + t190 * t207;
t37 = -t134 * t384 + t138 * t290 + t142 * t291 + t183 * t284 + t187 * t206 + t191 * t207;
t68 = -t164 * t384 + t166 * t290 + t168 * t291 + t206 * t225 + t207 * t227 + t223 * t284;
t8 = t104 * t345 + t105 * t347 - t36 * t387 + t37 * t356 + (t121 * t435 - t378 * t68) * t369;
t82 = t196 * t356 + t198 * t384 + t200 * t327 + t233 * t347 - t235 * t284 + t237 * t285;
t83 = t197 * t356 + t199 * t384 + t201 * t327 + t234 * t347 - t236 * t284 + t238 * t285;
t94 = t248 * t356 + t249 * t384 + t250 * t327 - t284 * t293 + t285 * t294 + t292 * t347;
t478 = t151 * t345 + t152 * t347 - t82 * t387 + t83 * t356 + (t160 * t435 - t378 * t94) * t369 + t7 + t8;
t477 = 2 * m(5);
t476 = 2 * m(6);
t475 = 2 * m(7);
t370 = cos(pkin(11));
t472 = pkin(3) * t370;
t471 = pkin(5) * t377;
t469 = Icges(3,4) * t376;
t468 = Icges(3,4) * t378;
t467 = t345 * t159;
t466 = t347 * t160;
t367 = sin(pkin(11));
t463 = t367 * t372;
t455 = rSges(7,1) * t205 + rSges(7,2) * t204 + pkin(5) * t382 - qJD(6) * t383 + t282 * t482 + t283 * t471;
t454 = rSges(7,1) * t207 + rSges(7,2) * t206 + pkin(5) * t381 - qJD(6) * t384 + t284 * t482 + t285 * t471;
t146 = rSges(6,1) * t207 + rSges(6,2) * t206 + rSges(6,3) * t284;
t214 = pkin(4) * t285 + pkin(9) * t284;
t453 = -t146 - t214;
t452 = rSges(7,1) * t247 + rSges(7,2) * t246 + pkin(5) * t380 + t341 * qJD(6) + t322 * t482 + t323 * t471;
t451 = rSges(7,1) * t289 + rSges(7,2) * t288 - pkin(5) * t465 + t325 * t471 - t383 * t482;
t450 = rSges(7,1) * t291 + rSges(7,2) * t290 + pkin(5) * t464 + t327 * t471 - t384 * t482;
t195 = rSges(6,1) * t291 + rSges(6,2) * t290 - rSges(6,3) * t384;
t268 = pkin(4) * t327 - pkin(9) * t384;
t449 = -t195 - t268;
t213 = pkin(4) * t283 + pkin(9) * t282;
t267 = pkin(4) * t325 - pkin(9) * t383;
t448 = t356 * t213 + t347 * t267;
t447 = -rSges(7,1) * t388 + rSges(7,2) * t331 - pkin(5) * t419 + t341 * t482 + t342 * t471;
t240 = pkin(8) * t347 - t346 * t472;
t287 = -pkin(2) * t346 + qJ(3) * t347 + qJD(3) * t356;
t280 = t372 * t287;
t446 = t372 * t240 + t280;
t422 = t367 * t462;
t244 = pkin(3) * t422 + pkin(8) * t356 + t357 * t472;
t321 = pkin(2) * t357 + qJ(3) * t356;
t319 = t372 * t321;
t445 = t372 * t244 + t319;
t230 = -rSges(6,1) * t388 + rSges(6,2) * t331 + rSges(6,3) * t341;
t304 = t342 * pkin(4) + t341 * pkin(9);
t444 = -t230 - t304;
t239 = pkin(8) * t345 + t344 * t472;
t286 = pkin(2) * t344 + qJ(3) * t345 - qJD(3) * t387;
t443 = -t239 - t286;
t421 = t367 * t461;
t243 = -pkin(3) * t421 - pkin(8) * t387 + t355 * t472;
t320 = pkin(2) * t355 - qJ(3) * t387;
t442 = -t243 - t320;
t441 = t267 * t459 - t304 * t387;
t440 = t286 * t462 + t287 * t461;
t358 = (pkin(2) * t376 - qJ(3) * t378) * t369;
t439 = -pkin(3) * t463 - (-pkin(8) * t378 + t376 * t472) * t369 - t358;
t438 = t320 * t462 + t321 * t461;
t338 = (-qJD(3) * t378 + (pkin(2) * t378 + qJ(3) * t376) * qJD(2)) * t369;
t436 = qJD(2) * t369;
t437 = -(pkin(8) * t376 + t378 * t472) * t436 - t338;
t1 = t118 * t322 + t282 * t98 + t284 * t99 - t30 * t383 - t31 * t384 + t341 * t65;
t2 = t100 * t282 + t101 * t284 + t119 * t322 - t32 * t383 - t33 * t384 + t341 * t66;
t432 = t2 / 0.2e1 + t1 / 0.2e1;
t3 = t102 * t282 + t103 * t284 + t120 * t322 - t34 * t383 + t341 * t67 - t35 * t384;
t4 = t104 * t282 + t105 * t284 + t121 * t322 + t341 * t68 - t36 * t383 - t37 * t384;
t431 = t4 / 0.2e1 + t3 / 0.2e1;
t110 = t182 * t341 + t186 * t331 - t190 * t388;
t111 = t183 * t341 + t187 * t331 - t191 * t388;
t130 = t223 * t341 + t225 * t331 - t227 * t388;
t42 = t133 * t341 + t137 * t331 - t141 * t388 + t182 * t322 + t186 * t246 + t190 * t247;
t43 = t134 * t341 + t138 * t331 - t142 * t388 + t183 * t322 + t187 * t246 + t191 * t247;
t74 = t164 * t341 + t166 * t331 - t168 * t388 + t223 * t322 + t225 * t246 + t227 * t247;
t10 = t110 * t282 + t111 * t284 + t130 * t322 + t341 * t74 - t383 * t42 - t384 * t43;
t108 = t180 * t341 + t184 * t331 - t188 * t388;
t109 = t181 * t341 + t185 * t331 - t189 * t388;
t129 = t222 * t341 + t224 * t331 - t226 * t388;
t40 = t131 * t341 + t135 * t331 - t139 * t388 + t180 * t322 + t184 * t246 + t188 * t247;
t41 = t132 * t341 + t136 * t331 - t140 * t388 + t181 * t322 + t185 * t246 + t189 * t247;
t73 = t163 * t341 + t165 * t331 - t167 * t388 + t222 * t322 + t224 * t246 + t226 * t247;
t9 = t108 * t282 + t109 * t284 + t129 * t322 + t341 * t73 - t383 * t40 - t384 * t41;
t430 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t108 * t345 + t109 * t347 - t40 * t387 + t41 * t356 + (t129 * t435 - t378 * t73) * t369;
t12 = t110 * t345 + t111 * t347 - t42 * t387 + t43 * t356 + (t130 * t435 - t378 * t74) * t369;
t429 = t12 / 0.2e1 + t11 / 0.2e1;
t13 = t372 * t65 + (-t30 * t371 + t31 * t368) * t369;
t14 = t372 * t66 + (-t32 * t371 + t33 * t368) * t369;
t428 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t372 * t67 + (-t34 * t371 + t35 * t368) * t369;
t16 = t372 * t68 + (-t36 * t371 + t368 * t37) * t369;
t427 = t15 / 0.2e1 + t16 / 0.2e1;
t17 = t372 * t73 + (t368 * t41 - t371 * t40) * t369;
t18 = t372 * t74 + (t368 * t43 - t371 * t42) * t369;
t426 = t18 / 0.2e1 + t17 / 0.2e1;
t425 = ((-t100 - t98) * t371 + (t101 + t99) * t368) * t481 + (t119 + t118) * t480;
t424 = ((-t102 - t104) * t371 + (t103 + t105) * t368) * t481 + (t120 + t121) * t480;
t423 = ((-t108 - t110) * t371 + (t109 + t111) * t368) * t481 + (t130 + t129) * t480;
t417 = -t214 - t454;
t416 = -t268 - t450;
t266 = t323 * pkin(4) + t322 * pkin(9);
t415 = t213 * t459 - t266 * t387 + t345 * t304;
t414 = t372 * t214 + t446;
t413 = -t213 + t443;
t412 = -t304 - t447;
t411 = t372 * t268 + t445;
t410 = -t267 + t442;
t409 = -t266 + t437;
t408 = -t304 + t439;
t352 = -t367 * t460 + t370 * t372;
t353 = t370 * t460 + t463;
t405 = (-t353 * rSges(4,1) - t352 * rSges(4,2) + rSges(4,3) * t459 - t358) * t369;
t399 = rSges(4,1) * t370 - rSges(4,2) * t367;
t404 = (-(rSges(4,3) * t376 + t378 * t399) * t436 - t338) * t369;
t403 = t239 * t462 + t240 * t461 + t440;
t402 = t243 * t462 + t244 * t461 + t438;
t251 = rSges(5,1) * t323 - rSges(5,2) * t322 + rSges(5,3) * t407;
t401 = (-t251 + t437) * t369;
t295 = t342 * rSges(5,1) - t341 * rSges(5,2) - rSges(5,3) * t459;
t400 = (-t295 + t439) * t369;
t396 = Icges(4,1) * t370 - Icges(4,4) * t367;
t395 = Icges(4,4) * t370 - Icges(4,2) * t367;
t394 = Icges(4,5) * t370 - Icges(4,6) * t367;
t393 = -(Icges(4,4) * t353 + Icges(4,2) * t352 - Icges(4,6) * t459) * t367 + (Icges(4,1) * t353 + Icges(4,4) * t352 - Icges(4,5) * t459) * t370;
t170 = rSges(6,1) * t247 + rSges(6,2) * t246 + rSges(6,3) * t322;
t392 = (-t170 + t409) * t369;
t391 = (-t230 + t408) * t369;
t390 = t213 * t462 + t214 * t461 + t403;
t389 = t267 * t462 + t268 * t461 + t402;
t386 = (t409 - t452) * t369;
t385 = (t408 - t447) * t369;
t334 = -t355 * t367 - t370 * t461;
t335 = t355 * t370 - t421;
t336 = -t357 * t367 + t370 * t462;
t337 = t357 * t370 + t422;
t379 = (-(Icges(4,4) * t337 + Icges(4,2) * t336 + Icges(4,6) * t356) * t367 + (Icges(4,1) * t337 + Icges(4,4) * t336 + Icges(4,5) * t356) * t370) * t368 - (-(Icges(4,4) * t335 + Icges(4,2) * t334 - Icges(4,6) * t387) * t367 + (Icges(4,1) * t335 + Icges(4,4) * t334 - Icges(4,5) * t387) * t370) * t371;
t351 = (rSges(3,1) * t378 - rSges(3,2) * t376) * t436;
t350 = (Icges(3,1) * t378 - t469) * t436;
t349 = (-Icges(3,2) * t376 + t468) * t436;
t348 = (Icges(3,5) * t378 - Icges(3,6) * t376) * t436;
t343 = t372 * rSges(3,3) + (rSges(3,1) * t376 + rSges(3,2) * t378) * t369;
t340 = Icges(3,5) * t372 + (Icges(3,1) * t376 + t468) * t369;
t339 = Icges(3,6) * t372 + (Icges(3,2) * t378 + t469) * t369;
t330 = (Icges(4,5) * t376 + t378 * t396) * t436;
t329 = (Icges(4,6) * t376 + t378 * t395) * t436;
t328 = (Icges(4,3) * t376 + t378 * t394) * t436;
t318 = -rSges(3,1) * t346 - rSges(3,2) * t347;
t317 = rSges(3,1) * t344 - rSges(3,2) * t345;
t315 = -Icges(3,1) * t346 - Icges(3,4) * t347;
t314 = Icges(3,1) * t344 - Icges(3,4) * t345;
t313 = -Icges(3,4) * t346 - Icges(3,2) * t347;
t312 = Icges(3,4) * t344 - Icges(3,2) * t345;
t311 = -Icges(3,5) * t346 - Icges(3,6) * t347;
t310 = Icges(3,5) * t344 - Icges(3,6) * t345;
t307 = rSges(3,1) * t357 - rSges(3,2) * t356 + rSges(3,3) * t462;
t306 = rSges(3,1) * t355 + rSges(3,2) * t387 - rSges(3,3) * t461;
t302 = Icges(3,1) * t357 - Icges(3,4) * t356 + Icges(3,5) * t462;
t301 = Icges(3,1) * t355 + Icges(3,4) * t387 - Icges(3,5) * t461;
t300 = Icges(3,4) * t357 - Icges(3,2) * t356 + Icges(3,6) * t462;
t299 = Icges(3,4) * t355 + Icges(3,2) * t387 - Icges(3,6) * t461;
t296 = Icges(4,5) * t353 + Icges(4,6) * t352 - Icges(4,3) * t459;
t276 = rSges(4,3) * t347 - t346 * t399;
t275 = rSges(4,3) * t345 + t344 * t399;
t274 = Icges(4,5) * t347 - t346 * t396;
t273 = Icges(4,5) * t345 + t344 * t396;
t272 = Icges(4,6) * t347 - t346 * t395;
t271 = Icges(4,6) * t345 + t344 * t395;
t270 = Icges(4,3) * t347 - t346 * t394;
t269 = Icges(4,3) * t345 + t344 * t394;
t261 = t268 * t407;
t260 = rSges(4,1) * t337 + rSges(4,2) * t336 + rSges(4,3) * t356;
t259 = rSges(4,1) * t335 + rSges(4,2) * t334 - rSges(4,3) * t387;
t253 = Icges(4,5) * t337 + Icges(4,6) * t336 + Icges(4,3) * t356;
t252 = Icges(4,5) * t335 + Icges(4,6) * t334 - Icges(4,3) * t387;
t245 = t356 * t267;
t242 = rSges(5,1) * t327 + rSges(5,2) * t384 + rSges(5,3) * t356;
t241 = rSges(5,1) * t325 + rSges(5,2) * t383 - rSges(5,3) * t387;
t216 = (t317 * t368 + t318 * t371) * t369;
t203 = rSges(5,1) * t285 - rSges(5,2) * t284 + rSges(5,3) * t347;
t202 = rSges(5,1) * t283 - rSges(5,2) * t282 + rSges(5,3) * t345;
t193 = rSges(6,1) * t289 + rSges(6,2) * t288 - rSges(6,3) * t383;
t179 = -t242 * t459 - t356 * t295;
t178 = t241 * t459 - t295 * t387;
t175 = (-t259 - t320) * t372 + t371 * t405;
t174 = t260 * t372 + t368 * t405 + t319;
t173 = -t292 * t459 - t341 * t293 + t342 * t294;
t172 = (-t275 - t286) * t372 + t371 * t404;
t171 = t276 * t372 + t368 * t404 + t280;
t162 = t241 * t356 + t242 * t387;
t158 = (t259 * t368 + t260 * t371) * t369 + t438;
t157 = (t275 * t368 + t276 * t371) * t369 + t440;
t156 = t195 * t341 + t230 * t384;
t155 = -t193 * t341 - t230 * t383;
t154 = -t234 * t459 - t341 * t236 + t342 * t238;
t153 = -t233 * t459 - t341 * t235 + t342 * t237;
t148 = (-t241 + t442) * t372 + t371 * t400;
t147 = t242 * t372 + t368 * t400 + t445;
t144 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t282;
t126 = -t193 * t384 + t195 * t383;
t125 = t356 * t444 + t449 * t459;
t124 = t193 * t459 - t230 * t387 + t441;
t123 = -t356 * t251 - t347 * t295 + (-t203 * t378 + t242 * t435) * t369;
t122 = -t387 * t251 + t345 * t295 + (t202 * t378 - t241 * t435) * t369;
t117 = (t241 * t368 + t242 * t371) * t369 + t402;
t116 = (-t202 + t443) * t372 + t371 * t401;
t115 = t203 * t372 + t368 * t401 + t446;
t114 = t193 * t356 - t387 * t449 + t245;
t113 = -t341 * t249 + t342 * t250 - t322 * t293 + t323 * t294 + (-t248 * t378 + t292 * t435) * t369;
t112 = t202 * t356 + t203 * t387 + t241 * t347 - t242 * t345;
t107 = (-t193 + t410) * t372 + t371 * t391;
t106 = t195 * t372 + t368 * t391 + t411;
t97 = t341 * t450 + t384 * t447;
t96 = -t341 * t451 - t383 * t447;
t95 = (t202 * t368 + t203 * t371) * t369 + t403;
t92 = t356 * t412 + t416 * t459;
t91 = -t387 * t447 + t451 * t459 + t441;
t90 = (t193 * t368 + t195 * t371) * t369 + t389;
t89 = t383 * t450 - t384 * t451;
t88 = -t341 * t199 + t342 * t201 - t322 * t236 + t323 * t238 + (-t197 * t378 + t234 * t435) * t369;
t87 = -t341 * t198 + t342 * t200 - t322 * t235 + t323 * t237 + (-t196 * t378 + t233 * t435) * t369;
t86 = (t410 - t451) * t372 + t371 * t385;
t85 = t368 * t385 + t372 * t450 + t411;
t84 = t356 * t451 - t387 * t416 + t245;
t79 = t146 * t341 + t170 * t384 + t195 * t322 - t230 * t284;
t78 = -t144 * t341 - t170 * t383 - t193 * t322 + t230 * t282;
t77 = (-t144 + t413) * t372 + t371 * t392;
t76 = t146 * t372 + t368 * t392 + t414;
t75 = (t368 * t451 + t371 * t450) * t369 + t389;
t72 = -t144 * t384 + t146 * t383 + t193 * t284 - t195 * t282;
t71 = t261 + (-t170 - t266) * t356 + t444 * t347 + (t195 * t435 + t378 * t453) * t369;
t70 = -t387 * t170 + t345 * t230 + (t144 * t378 + (-t193 - t267) * t435) * t369 + t415;
t69 = (t144 * t368 + t146 * t371) * t369 + t390;
t62 = -t110 * t387 + t111 * t356 - t130 * t459;
t61 = -t108 * t387 + t109 * t356 - t129 * t459;
t60 = -t110 * t383 - t111 * t384 + t130 * t341;
t59 = -t108 * t383 - t109 * t384 + t129 * t341;
t58 = t144 * t356 + t193 * t347 + t345 * t449 - t387 * t453 + t448;
t57 = (t413 - t455) * t372 + t371 * t386;
t56 = t368 * t386 + t372 * t454 + t414;
t51 = -t104 * t387 + t105 * t356 - t121 * t459;
t50 = -t102 * t387 + t103 * t356 - t120 * t459;
t49 = -t100 * t387 + t101 * t356 - t119 * t459;
t48 = -t118 * t459 + t99 * t356 - t387 * t98;
t47 = -t104 * t383 - t105 * t384 + t121 * t341;
t46 = -t102 * t383 - t103 * t384 + t120 * t341;
t45 = -t100 * t383 - t101 * t384 + t119 * t341;
t44 = t118 * t341 - t383 * t98 - t384 * t99;
t39 = -t284 * t447 + t322 * t450 + t341 * t454 + t384 * t452;
t38 = t282 * t447 - t322 * t451 - t341 * t455 - t383 * t452;
t29 = (t368 * t455 + t371 * t454) * t369 + t390;
t28 = t261 + (-t266 - t452) * t356 + t412 * t347 + (t378 * t417 + t435 * t450) * t369;
t27 = -t452 * t387 + t447 * t345 + (t455 * t378 + (-t267 - t451) * t435) * t369 + t415;
t26 = t113 * t372 + (t368 * t88 - t371 * t87) * t369;
t25 = -t282 * t450 + t284 * t451 + t383 * t454 - t384 * t455;
t24 = t372 * t94 + (t368 * t83 - t371 * t82) * t369;
t23 = t372 * t93 + (t368 * t81 - t371 * t80) * t369;
t22 = t345 * t416 + t347 * t451 + t356 * t455 - t387 * t417 + t448;
t21 = t153 * t345 + t154 * t347 - t87 * t387 + t88 * t356 + (-t113 * t378 + t173 * t435) * t369;
t19 = [0; m(3) * t216 + m(4) * t157 + m(5) * t95 + m(6) * t69 + m(7) * t29; -t14 * t461 - t13 * t461 + t15 * t462 + t16 * t462 - t23 * t461 + t24 * t462 - ((-t300 * t345 + t302 * t344 - t311 * t461 + t313 * t387 + t315 * t355) * t462 - (-t299 * t345 + t301 * t344 - t310 * t461 + t312 * t387 + t314 * t355) * t461 + (-t339 * t345 + t340 * t344 - t348 * t461 + t349 * t387 + t350 * t355) * t372) * t461 - ((t345 * t296 - t328 * t387 + t334 * t329 + t335 * t330 + t344 * t393) * t372 + ((t253 * t345 - t270 * t387 + t272 * t334 + t274 * t335) * t368 - (t252 * t345 - t269 * t387 + t271 * t334 + t335 * t273) * t371 + t379 * t344) * t369) * t461 + ((-t300 * t347 - t302 * t346 + t311 * t462 - t313 * t356 + t315 * t357) * t462 - (-t299 * t347 - t301 * t346 + t310 * t462 - t312 * t356 + t314 * t357) * t461 + (-t339 * t347 - t340 * t346 + t348 * t462 - t349 * t356 + t350 * t357) * t372) * t462 + ((t347 * t296 + t356 * t328 + t336 * t329 + t337 * t330 - t346 * t393) * t372 + ((t253 * t347 + t270 * t356 + t272 * t336 + t274 * t337) * t368 - (t252 * t347 + t269 * t356 + t271 * t336 + t273 * t337) * t371 - t379 * t346) * t369) * t462 + t372 * t17 + (t29 * t75 + t56 * t85 + t57 * t86) * t475 + (t106 * t76 + t107 * t77 + t69 * t90) * t476 + t372 * t18 + t372 * t26 + (t115 * t147 + t116 * t148 + t117 * t95) * t477 + 0.2e1 * m(4) * (t157 * t158 + t171 * t174 + t172 * t175) + t372 * (t372 ^ 2 * t348 + (((t313 * t378 + t315 * t376) * t368 - (t312 * t378 + t314 * t376) * t371 + ((-t300 * t376 + t302 * t378) * t368 - (-t299 * t376 + t301 * t378) * t371) * qJD(2)) * t369 + (-t310 * t371 + t311 * t368 + t349 * t378 + t350 * t376 + (-t339 * t376 + t340 * t378) * qJD(2)) * t372) * t369) + t372 * ((t352 * t329 + t353 * t330) * t372 + ((t352 * t272 + t353 * t274) * t368 - (t352 * t271 + t353 * t273) * t371 + (-t328 * t372 + (t269 * t371 - t270 * t368) * t369) * t378 + ((t296 * t372 + (-t252 * t371 + t253 * t368) * t369) * t376 + (t369 * t379 + t372 * t393) * t378) * qJD(2)) * t369) + 0.2e1 * m(3) * ((-t306 * t372 - t343 * t461) * (-t317 * t372 - t351 * t461) + (t307 * t372 - t343 * t462) * (t318 * t372 - t351 * t462) + (t306 * t368 + t307 * t371) * t369 * t216); (m(4) + m(5) + m(6) + m(7)) * t407; m(7) * (t345 * t85 + t347 * t86 - t387 * t56 + t356 * t57 + (-t29 * t378 + t435 * t75) * t369) + m(6) * (t345 * t106 + t347 * t107 - t387 * t76 + t356 * t77 + (-t378 * t69 + t435 * t90) * t369) + m(5) * (-t387 * t115 + t356 * t116 + t345 * t147 + t347 * t148 + (t117 * t435 - t378 * t95) * t369) + m(4) * (-t387 * t171 + t356 * t172 + t345 * t174 + t347 * t175 + (-t157 * t378 + t158 * t435) * t369); 0.4e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (-t369 ^ 2 * t376 * t434 - t345 * t387 + t356 * t347); m(5) * t112 + m(6) * t58 + m(7) * t22; t424 * t347 + t425 * t345 + (t24 / 0.2e1 + t427) * t356 - (t23 / 0.2e1 + t428) * t387 + m(5) * (t112 * t117 + t115 * t179 + t116 * t178 + t122 * t148 + t123 * t147 + t162 * t95) + m(6) * (t106 * t71 + t107 * t70 + t114 * t69 + t124 * t77 + t125 * t76 + t58 * t90) + m(7) * (t22 * t75 + t27 * t86 + t28 * t85 + t29 * t84 + t56 * t92 + t57 * t91) + (t21 / 0.2e1 + t466 / 0.2e1 + t467 / 0.2e1 + t429) * t372 + (t347 * (-t151 * t371 + t152 * t368) / 0.2e1 + t345 * (-t149 * t371 + t150 * t368) / 0.2e1 + (-t26 / 0.2e1 - t426) * t378 + (t173 * t480 + (-t153 * t371 + t154 * t368) * t481 + t423) * t435) * t369 + t478 * t462 / 0.2e1 - t479 * t461 / 0.2e1; m(5) * (t122 * t356 - t123 * t387 + t178 * t347 + t179 * t345 + (-t112 * t378 + t162 * t435) * t369) + m(6) * (t124 * t347 + t125 * t345 - t71 * t387 + t70 * t356 + (t114 * t435 - t378 * t58) * t369) + m(7) * (t27 * t356 - t28 * t387 + t92 * t345 + t91 * t347 + (-t22 * t378 + t435 * t84) * t369); t478 * t356 - t479 * t387 + (-t151 * t387 + t152 * t356 + t50 + t51) * t347 + (-t149 * t387 + t150 * t356 + t48 + t49) * t345 + (t22 * t84 + t27 * t91 + t28 * t92) * t475 + (t114 * t58 + t124 * t70 + t125 * t71) * t476 + (t112 * t162 + t122 * t178 + t123 * t179) * t477 + ((-t153 * t387 + t154 * t356 + t61 + t62) * t435 + (-t173 * t407 - t11 - t12 - t21 - t466 - t467) * t378) * t369; m(6) * t72 + m(7) * t25; t430 * t372 + t426 * t341 - t427 * t384 - t428 * t383 + t423 * t322 + t424 * t284 + t425 * t282 + m(7) * (t25 * t75 + t29 * t89 + t38 * t86 + t39 * t85 + t56 * t97 + t57 * t96) + m(6) * (t106 * t79 + t107 * t78 + t126 * t69 + t155 * t77 + t156 * t76 + t72 * t90) + (t368 * t431 - t371 * t432) * t369; m(6) * (t155 * t347 + t156 * t345 - t79 * t387 + t78 * t356 + (t126 * t435 - t378 * t72) * t369) + m(7) * (t97 * t345 + t96 * t347 - t39 * t387 + t38 * t356 + (-t25 * t378 + t435 * t89) * t369); t431 * t356 - t432 * t387 + (t47 / 0.2e1 + t46 / 0.2e1) * t347 + (t45 / 0.2e1 + t44 / 0.2e1) * t345 + t429 * t341 - (t8 / 0.2e1 + t7 / 0.2e1) * t384 - (t6 / 0.2e1 + t5 / 0.2e1) * t383 + (t62 / 0.2e1 + t61 / 0.2e1) * t322 + (t51 / 0.2e1 + t50 / 0.2e1) * t284 + (t49 / 0.2e1 + t48 / 0.2e1) * t282 + m(6) * (t114 * t72 + t124 * t78 + t125 * t79 + t126 * t58 + t155 * t70 + t156 * t71) + m(7) * (t22 * t89 + t25 * t84 + t27 * t96 + t28 * t97 + t38 * t91 + t39 * t92) + (-t430 * t378 + (t59 / 0.2e1 + t60 / 0.2e1) * t435) * t369; (t9 + t10) * t341 - (t3 + t4) * t384 - (t2 + t1) * t383 + (t59 + t60) * t322 + (t46 + t47) * t284 + (t45 + t44) * t282 + (t25 * t89 + t38 * t96 + t39 * t97) * t475 + (t126 * t72 + t155 * t78 + t156 * t79) * t476; m(7) * t322; m(7) * (t282 * t85 + t284 * t86 + t29 * t341 + t322 * t75 - t383 * t56 - t384 * t57); m(7) * (-t282 * t387 + t284 * t356 - t383 * t345 - t384 * t347 + (-t322 * t378 + t341 * t435) * t369); m(7) * (t22 * t341 - t27 * t384 - t28 * t383 + t282 * t92 + t284 * t91 + t322 * t84); m(7) * (t25 * t341 + t282 * t97 + t284 * t96 + t322 * t89 - t38 * t384 - t383 * t39); (-t282 * t383 - t284 * t384 + t322 * t341) * t475;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
