% Calculate time derivative of joint inertia matrix for
% S6PRPRRP4
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:18
% EndTime: 2019-03-08 20:09:49
% DurationCPUTime: 18.40s
% Computational Cost: add. (81738->1122), mult. (151173->1589), div. (0->0), fcn. (172407->12), ass. (0->404)
t483 = rSges(7,1) + pkin(5);
t482 = rSges(7,3) + qJ(6);
t374 = sin(pkin(6));
t481 = t374 / 0.2e1;
t377 = cos(pkin(6));
t480 = t377 / 0.2e1;
t373 = sin(pkin(10));
t465 = t373 * t374;
t376 = cos(pkin(10));
t464 = t374 * t376;
t381 = cos(qJ(2));
t380 = sin(qJ(2));
t461 = t377 * t380;
t359 = t373 * t381 + t376 * t461;
t436 = pkin(11) + qJ(4);
t371 = sin(t436);
t408 = cos(t436);
t325 = t359 * t408 - t371 * t464;
t399 = t374 * t408;
t344 = t377 * t371 + t380 * t399;
t379 = sin(qJ(5));
t462 = t374 * t381;
t472 = cos(qJ(5));
t332 = t344 * t472 - t379 * t462;
t383 = -t359 * t371 - t376 * t399;
t460 = t377 * t381;
t390 = -t373 * t380 + t376 * t460;
t232 = Icges(5,5) * t325 + Icges(5,6) * t383 - Icges(5,3) * t390;
t234 = Icges(5,4) * t325 + Icges(5,2) * t383 - Icges(5,6) * t390;
t236 = Icges(5,1) * t325 + Icges(5,4) * t383 - Icges(5,5) * t390;
t147 = -t232 * t390 + t234 * t383 + t236 * t325;
t423 = t373 * t461;
t361 = t376 * t381 - t423;
t327 = t361 * t408 + t371 * t465;
t360 = t373 * t460 + t376 * t380;
t384 = -t361 * t371 + t373 * t399;
t233 = Icges(5,5) * t327 + Icges(5,6) * t384 + Icges(5,3) * t360;
t235 = Icges(5,4) * t327 + Icges(5,2) * t384 + Icges(5,6) * t360;
t237 = Icges(5,1) * t327 + Icges(5,4) * t384 + Icges(5,5) * t360;
t148 = -t233 * t390 + t235 * t383 + t237 * t325;
t463 = t374 * t380;
t343 = t371 * t463 - t377 * t408;
t292 = Icges(5,5) * t344 - Icges(5,6) * t343 - Icges(5,3) * t462;
t293 = Icges(5,4) * t344 - Icges(5,2) * t343 - Icges(5,6) * t462;
t294 = Icges(5,1) * t344 - Icges(5,4) * t343 - Icges(5,5) * t462;
t159 = -t292 * t390 + t293 * t383 + t294 * t325;
t349 = t359 * qJD(2);
t437 = qJD(2) * t381;
t351 = -qJD(2) * t423 + t376 * t437;
t438 = qJD(2) * t380;
t385 = -t344 * t379 - t462 * t472;
t221 = Icges(7,5) * t332 + Icges(7,6) * t343 - Icges(7,3) * t385;
t223 = Icges(7,4) * t332 + Icges(7,2) * t343 - Icges(7,6) * t385;
t225 = Icges(7,1) * t332 + Icges(7,4) * t343 - Icges(7,5) * t385;
t289 = t325 * t472 - t379 * t390;
t389 = -t325 * t379 - t390 * t472;
t118 = -t221 * t389 - t223 * t383 + t225 * t289;
t348 = t390 * qJD(2);
t283 = qJD(4) * t383 + t348 * t408;
t202 = qJD(5) * t289 + t283 * t379 - t349 * t472;
t203 = qJD(5) * t389 + t283 * t472 + t349 * t379;
t282 = qJD(4) * t325 + t348 * t371;
t129 = Icges(7,5) * t203 + Icges(7,6) * t282 + Icges(7,3) * t202;
t133 = Icges(7,4) * t203 + Icges(7,2) * t282 + Icges(7,6) * t202;
t137 = Icges(7,1) * t203 + Icges(7,4) * t282 + Icges(7,5) * t202;
t178 = Icges(7,5) * t289 - Icges(7,6) * t383 - Icges(7,3) * t389;
t182 = Icges(7,4) * t289 - Icges(7,2) * t383 - Icges(7,6) * t389;
t186 = Icges(7,1) * t289 - Icges(7,4) * t383 - Icges(7,5) * t389;
t27 = -t129 * t389 - t133 * t383 + t137 * t289 + t178 * t202 + t182 * t282 + t186 * t203;
t350 = t360 * qJD(2);
t285 = qJD(4) * t384 - t350 * t408;
t291 = t327 * t472 + t360 * t379;
t204 = qJD(5) * t291 + t285 * t379 - t351 * t472;
t388 = -t327 * t379 + t360 * t472;
t205 = qJD(5) * t388 + t285 * t472 + t351 * t379;
t284 = qJD(4) * t327 - t350 * t371;
t130 = Icges(7,5) * t205 + Icges(7,6) * t284 + Icges(7,3) * t204;
t134 = Icges(7,4) * t205 + Icges(7,2) * t284 + Icges(7,6) * t204;
t138 = Icges(7,1) * t205 + Icges(7,4) * t284 + Icges(7,5) * t204;
t179 = Icges(7,5) * t291 - Icges(7,6) * t384 - Icges(7,3) * t388;
t183 = Icges(7,4) * t291 - Icges(7,2) * t384 - Icges(7,6) * t388;
t187 = Icges(7,1) * t291 - Icges(7,4) * t384 - Icges(7,5) * t388;
t28 = -t130 * t389 - t134 * t383 + t138 * t289 + t179 * t202 + t183 * t282 + t187 * t203;
t323 = -qJD(4) * t343 + t399 * t437;
t409 = t374 * t438;
t245 = qJD(5) * t332 + t323 * t379 - t409 * t472;
t246 = qJD(5) * t385 + t323 * t472 + t379 * t409;
t322 = t371 * t374 * t437 + qJD(4) * t344;
t162 = Icges(7,5) * t246 + Icges(7,6) * t322 + Icges(7,3) * t245;
t164 = Icges(7,4) * t246 + Icges(7,2) * t322 + Icges(7,6) * t245;
t166 = Icges(7,1) * t246 + Icges(7,4) * t322 + Icges(7,5) * t245;
t65 = -t162 * t389 - t164 * t383 + t166 * t289 + t202 * t221 + t203 * t225 + t223 * t282;
t96 = -t178 * t389 - t182 * t383 + t186 * t289;
t97 = -t179 * t389 - t183 * t383 + t187 * t289;
t5 = -t27 * t390 + t28 * t360 + t96 * t349 + t97 * t351 + (t118 * t438 - t381 * t65) * t374;
t222 = Icges(6,5) * t332 + Icges(6,6) * t385 + Icges(6,3) * t343;
t224 = Icges(6,4) * t332 + Icges(6,2) * t385 + Icges(6,6) * t343;
t226 = Icges(6,1) * t332 + Icges(6,4) * t385 + Icges(6,5) * t343;
t119 = -t222 * t383 + t224 * t389 + t226 * t289;
t131 = Icges(6,5) * t203 - Icges(6,6) * t202 + Icges(6,3) * t282;
t135 = Icges(6,4) * t203 - Icges(6,2) * t202 + Icges(6,6) * t282;
t139 = Icges(6,1) * t203 - Icges(6,4) * t202 + Icges(6,5) * t282;
t180 = Icges(6,5) * t289 + Icges(6,6) * t389 - Icges(6,3) * t383;
t184 = Icges(6,4) * t289 + Icges(6,2) * t389 - Icges(6,6) * t383;
t188 = Icges(6,1) * t289 + Icges(6,4) * t389 - Icges(6,5) * t383;
t29 = -t131 * t383 + t135 * t389 + t139 * t289 + t180 * t282 - t184 * t202 + t188 * t203;
t132 = Icges(6,5) * t205 - Icges(6,6) * t204 + Icges(6,3) * t284;
t136 = Icges(6,4) * t205 - Icges(6,2) * t204 + Icges(6,6) * t284;
t140 = Icges(6,1) * t205 - Icges(6,4) * t204 + Icges(6,5) * t284;
t181 = Icges(6,5) * t291 + Icges(6,6) * t388 - Icges(6,3) * t384;
t185 = Icges(6,4) * t291 + Icges(6,2) * t388 - Icges(6,6) * t384;
t189 = Icges(6,1) * t291 + Icges(6,4) * t388 - Icges(6,5) * t384;
t30 = -t132 * t383 + t136 * t389 + t140 * t289 + t181 * t282 - t185 * t202 + t189 * t203;
t163 = Icges(6,5) * t246 - Icges(6,6) * t245 + Icges(6,3) * t322;
t165 = Icges(6,4) * t246 - Icges(6,2) * t245 + Icges(6,6) * t322;
t167 = Icges(6,1) * t246 - Icges(6,4) * t245 + Icges(6,5) * t322;
t66 = -t163 * t383 + t165 * t389 + t167 * t289 - t202 * t224 + t203 * t226 + t222 * t282;
t98 = -t180 * t383 + t184 * t389 + t188 * t289;
t99 = -t181 * t383 + t185 * t389 + t189 * t289;
t6 = -t29 * t390 + t30 * t360 + t98 * t349 + t99 * t351 + (t119 * t438 - t381 * t66) * t374;
t194 = Icges(5,5) * t283 - Icges(5,6) * t282 + Icges(5,3) * t349;
t196 = Icges(5,4) * t283 - Icges(5,2) * t282 + Icges(5,6) * t349;
t198 = Icges(5,1) * t283 - Icges(5,4) * t282 + Icges(5,5) * t349;
t80 = -t194 * t390 + t196 * t383 + t198 * t325 + t232 * t349 - t234 * t282 + t236 * t283;
t195 = Icges(5,5) * t285 - Icges(5,6) * t284 + Icges(5,3) * t351;
t197 = Icges(5,4) * t285 - Icges(5,2) * t284 + Icges(5,6) * t351;
t199 = Icges(5,1) * t285 - Icges(5,4) * t284 + Icges(5,5) * t351;
t81 = -t195 * t390 + t197 * t383 + t199 * t325 + t233 * t349 - t235 * t282 + t237 * t283;
t247 = Icges(5,5) * t323 - Icges(5,6) * t322 + Icges(5,3) * t409;
t248 = Icges(5,4) * t323 - Icges(5,2) * t322 + Icges(5,6) * t409;
t249 = Icges(5,1) * t323 - Icges(5,4) * t322 + Icges(5,5) * t409;
t91 = -t247 * t390 + t248 * t383 + t249 * t325 - t282 * t293 + t283 * t294 + t292 * t349;
t479 = t147 * t349 + t148 * t351 - t80 * t390 + t81 * t360 + (t159 * t438 - t381 * t91) * t374 + t5 + t6;
t149 = t232 * t360 + t234 * t384 + t236 * t327;
t150 = t233 * t360 + t235 * t384 + t237 * t327;
t160 = t292 * t360 + t293 * t384 + t294 * t327;
t100 = -t178 * t388 - t182 * t384 + t186 * t291;
t101 = -t179 * t388 - t183 * t384 + t187 * t291;
t120 = -t221 * t388 - t223 * t384 + t225 * t291;
t31 = -t129 * t388 - t133 * t384 + t137 * t291 + t178 * t204 + t182 * t284 + t186 * t205;
t32 = -t130 * t388 - t134 * t384 + t138 * t291 + t179 * t204 + t183 * t284 + t187 * t205;
t67 = -t162 * t388 - t164 * t384 + t166 * t291 + t204 * t221 + t205 * t225 + t223 * t284;
t7 = t100 * t349 + t101 * t351 - t31 * t390 + t32 * t360 + (t120 * t438 - t381 * t67) * t374;
t102 = -t180 * t384 + t184 * t388 + t188 * t291;
t103 = -t181 * t384 + t185 * t388 + t189 * t291;
t121 = -t222 * t384 + t224 * t388 + t226 * t291;
t33 = -t131 * t384 + t135 * t388 + t139 * t291 + t180 * t284 - t184 * t204 + t188 * t205;
t34 = -t132 * t384 + t136 * t388 + t140 * t291 + t181 * t284 - t185 * t204 + t189 * t205;
t68 = -t163 * t384 + t165 * t388 + t167 * t291 - t204 * t224 + t205 * t226 + t222 * t284;
t8 = t102 * t349 + t103 * t351 - t33 * t390 + t34 * t360 + (t121 * t438 - t381 * t68) * t374;
t82 = t194 * t360 + t196 * t384 + t198 * t327 + t232 * t351 - t234 * t284 + t236 * t285;
t83 = t195 * t360 + t197 * t384 + t199 * t327 + t233 * t351 - t235 * t284 + t237 * t285;
t92 = t247 * t360 + t248 * t384 + t249 * t327 - t284 * t293 + t285 * t294 + t292 * t351;
t478 = t149 * t349 + t150 * t351 - t82 * t390 + t83 * t360 + (t160 * t438 - t381 * t92) * t374 + t7 + t8;
t477 = 2 * m(5);
t476 = 2 * m(6);
t475 = 2 * m(7);
t375 = cos(pkin(11));
t471 = pkin(3) * t375;
t470 = Icges(3,4) * t380;
t469 = Icges(3,4) * t381;
t468 = t349 * t159;
t467 = t351 * t160;
t372 = sin(pkin(11));
t466 = t372 * t377;
t458 = rSges(7,2) * t282 - qJD(6) * t389 + t482 * t202 + t203 * t483;
t457 = rSges(7,2) * t284 - qJD(6) * t388 + t482 * t204 + t205 * t483;
t144 = rSges(6,1) * t205 - rSges(6,2) * t204 + rSges(6,3) * t284;
t212 = pkin(4) * t285 + pkin(9) * t284;
t456 = -t144 - t212;
t455 = rSges(7,2) * t322 - qJD(6) * t385 + t482 * t245 + t246 * t483;
t454 = -rSges(7,2) * t383 + t289 * t483 - t482 * t389;
t453 = -rSges(7,2) * t384 + t291 * t483 - t482 * t388;
t193 = rSges(6,1) * t291 + rSges(6,2) * t388 - rSges(6,3) * t384;
t267 = pkin(4) * t327 - pkin(9) * t384;
t452 = -t193 - t267;
t211 = pkin(4) * t283 + pkin(9) * t282;
t266 = pkin(4) * t325 - pkin(9) * t383;
t451 = t360 * t211 + t351 * t266;
t239 = pkin(8) * t351 - t350 * t471;
t287 = -pkin(2) * t350 + qJ(3) * t351 + qJD(3) * t360;
t280 = t377 * t287;
t450 = t377 * t239 + t280;
t425 = t372 * t465;
t243 = pkin(3) * t425 + pkin(8) * t360 + t361 * t471;
t321 = pkin(2) * t361 + qJ(3) * t360;
t319 = t377 * t321;
t449 = t377 * t243 + t319;
t448 = rSges(7,2) * t343 + t332 * t483 - t482 * t385;
t229 = rSges(6,1) * t332 + rSges(6,2) * t385 + rSges(6,3) * t343;
t304 = pkin(4) * t344 + pkin(9) * t343;
t447 = -t229 - t304;
t238 = pkin(8) * t349 + t348 * t471;
t286 = pkin(2) * t348 + qJ(3) * t349 - qJD(3) * t390;
t446 = -t238 - t286;
t424 = t372 * t464;
t242 = -pkin(3) * t424 - pkin(8) * t390 + t359 * t471;
t320 = pkin(2) * t359 - qJ(3) * t390;
t445 = -t242 - t320;
t444 = t266 * t462 - t304 * t390;
t443 = t286 * t465 + t287 * t464;
t362 = (pkin(2) * t380 - qJ(3) * t381) * t374;
t442 = -pkin(3) * t466 - (-pkin(8) * t381 + t380 * t471) * t374 - t362;
t441 = t320 * t465 + t321 * t464;
t338 = (-qJD(3) * t381 + (pkin(2) * t381 + qJ(3) * t380) * qJD(2)) * t374;
t439 = qJD(2) * t374;
t440 = -(pkin(8) * t380 + t381 * t471) * t439 - t338;
t1 = t118 * t322 - t27 * t383 - t28 * t384 + t282 * t96 + t284 * t97 + t343 * t65;
t2 = t119 * t322 + t282 * t98 + t284 * t99 - t29 * t383 - t30 * t384 + t343 * t66;
t435 = t1 / 0.2e1 + t2 / 0.2e1;
t3 = t100 * t282 + t101 * t284 + t120 * t322 - t31 * t383 - t32 * t384 + t343 * t67;
t4 = t102 * t282 + t103 * t284 + t121 * t322 - t33 * t383 - t34 * t384 + t343 * t68;
t434 = t3 / 0.2e1 + t4 / 0.2e1;
t108 = t180 * t343 + t184 * t385 + t188 * t332;
t109 = t181 * t343 + t185 * t385 + t189 * t332;
t128 = t222 * t343 + t224 * t385 + t226 * t332;
t42 = t131 * t343 + t135 * t385 + t139 * t332 + t180 * t322 - t184 * t245 + t188 * t246;
t43 = t132 * t343 + t136 * t385 + t140 * t332 + t181 * t322 - t185 * t245 + t189 * t246;
t74 = t163 * t343 + t165 * t385 + t167 * t332 + t222 * t322 - t224 * t245 + t226 * t246;
t10 = t108 * t282 + t109 * t284 + t128 * t322 + t343 * t74 - t383 * t42 - t384 * t43;
t106 = -t178 * t385 + t182 * t343 + t186 * t332;
t107 = -t179 * t385 + t183 * t343 + t187 * t332;
t127 = -t221 * t385 + t223 * t343 + t225 * t332;
t40 = -t129 * t385 + t133 * t343 + t137 * t332 + t178 * t245 + t182 * t322 + t186 * t246;
t41 = -t130 * t385 + t134 * t343 + t138 * t332 + t179 * t245 + t183 * t322 + t187 * t246;
t73 = -t162 * t385 + t164 * t343 + t166 * t332 + t221 * t245 + t223 * t322 + t225 * t246;
t9 = t106 * t282 + t107 * t284 + t127 * t322 + t343 * t73 - t383 * t40 - t384 * t41;
t433 = -t9 / 0.2e1 - t10 / 0.2e1;
t11 = t106 * t349 + t107 * t351 - t40 * t390 + t41 * t360 + (t127 * t438 - t381 * t73) * t374;
t12 = t108 * t349 + t109 * t351 - t42 * t390 + t43 * t360 + (t128 * t438 - t381 * t74) * t374;
t432 = t11 / 0.2e1 + t12 / 0.2e1;
t13 = t377 * t65 + (-t27 * t376 + t28 * t373) * t374;
t14 = t377 * t66 + (-t29 * t376 + t30 * t373) * t374;
t431 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = t377 * t67 + (-t31 * t376 + t32 * t373) * t374;
t16 = t377 * t68 + (-t33 * t376 + t34 * t373) * t374;
t430 = t16 / 0.2e1 + t15 / 0.2e1;
t17 = t377 * t73 + (t373 * t41 - t376 * t40) * t374;
t18 = t377 * t74 + (t373 * t43 - t376 * t42) * t374;
t429 = t18 / 0.2e1 + t17 / 0.2e1;
t428 = ((-t96 - t98) * t376 + (t97 + t99) * t373) * t481 + (t118 + t119) * t480;
t427 = ((-t100 - t102) * t376 + (t101 + t103) * t373) * t481 + (t121 + t120) * t480;
t426 = ((-t106 - t108) * t376 + (t107 + t109) * t373) * t481 + (t128 + t127) * t480;
t420 = -t212 - t457;
t419 = -t267 - t453;
t265 = pkin(4) * t323 + pkin(9) * t322;
t418 = t211 * t462 - t265 * t390 + t349 * t304;
t417 = t377 * t212 + t450;
t416 = -t211 + t446;
t415 = t377 * t267 + t449;
t414 = -t304 - t448;
t413 = -t266 + t445;
t412 = -t265 + t440;
t411 = -t304 + t442;
t356 = -t372 * t463 + t375 * t377;
t357 = t375 * t463 + t466;
t407 = (-t357 * rSges(4,1) - t356 * rSges(4,2) + rSges(4,3) * t462 - t362) * t374;
t401 = rSges(4,1) * t375 - rSges(4,2) * t372;
t406 = (-(rSges(4,3) * t380 + t381 * t401) * t439 - t338) * t374;
t405 = t238 * t465 + t239 * t464 + t443;
t404 = t242 * t465 + t243 * t464 + t441;
t250 = rSges(5,1) * t323 - rSges(5,2) * t322 + rSges(5,3) * t409;
t403 = (-t250 + t440) * t374;
t295 = t344 * rSges(5,1) - t343 * rSges(5,2) - rSges(5,3) * t462;
t402 = (-t295 + t442) * t374;
t398 = Icges(4,1) * t375 - Icges(4,4) * t372;
t397 = Icges(4,4) * t375 - Icges(4,2) * t372;
t396 = Icges(4,5) * t375 - Icges(4,6) * t372;
t395 = -(Icges(4,4) * t357 + Icges(4,2) * t356 - Icges(4,6) * t462) * t372 + (Icges(4,1) * t357 + Icges(4,4) * t356 - Icges(4,5) * t462) * t375;
t169 = rSges(6,1) * t246 - rSges(6,2) * t245 + rSges(6,3) * t322;
t394 = (-t169 + t412) * t374;
t393 = (-t229 + t411) * t374;
t392 = t211 * t465 + t212 * t464 + t405;
t391 = t266 * t465 + t267 * t464 + t404;
t387 = (t412 - t455) * t374;
t386 = (t411 - t448) * t374;
t334 = -t359 * t372 - t375 * t464;
t335 = t359 * t375 - t424;
t336 = -t361 * t372 + t375 * t465;
t337 = t361 * t375 + t425;
t382 = (-(Icges(4,4) * t337 + Icges(4,2) * t336 + Icges(4,6) * t360) * t372 + (Icges(4,1) * t337 + Icges(4,4) * t336 + Icges(4,5) * t360) * t375) * t373 - (-(Icges(4,4) * t335 + Icges(4,2) * t334 - Icges(4,6) * t390) * t372 + (Icges(4,1) * t335 + Icges(4,4) * t334 - Icges(4,5) * t390) * t375) * t376;
t355 = (rSges(3,1) * t381 - rSges(3,2) * t380) * t439;
t354 = (Icges(3,1) * t381 - t470) * t439;
t353 = (-Icges(3,2) * t380 + t469) * t439;
t352 = (Icges(3,5) * t381 - Icges(3,6) * t380) * t439;
t345 = t377 * rSges(3,3) + (rSges(3,1) * t380 + rSges(3,2) * t381) * t374;
t342 = Icges(3,5) * t377 + (Icges(3,1) * t380 + t469) * t374;
t341 = Icges(3,6) * t377 + (Icges(3,2) * t381 + t470) * t374;
t330 = (Icges(4,5) * t380 + t381 * t398) * t439;
t329 = (Icges(4,6) * t380 + t381 * t397) * t439;
t328 = (Icges(4,3) * t380 + t381 * t396) * t439;
t318 = -rSges(3,1) * t350 - rSges(3,2) * t351;
t317 = rSges(3,1) * t348 - rSges(3,2) * t349;
t315 = -Icges(3,1) * t350 - Icges(3,4) * t351;
t314 = Icges(3,1) * t348 - Icges(3,4) * t349;
t313 = -Icges(3,4) * t350 - Icges(3,2) * t351;
t312 = Icges(3,4) * t348 - Icges(3,2) * t349;
t311 = -Icges(3,5) * t350 - Icges(3,6) * t351;
t310 = Icges(3,5) * t348 - Icges(3,6) * t349;
t307 = rSges(3,1) * t361 - rSges(3,2) * t360 + rSges(3,3) * t465;
t306 = rSges(3,1) * t359 + rSges(3,2) * t390 - rSges(3,3) * t464;
t302 = Icges(3,1) * t361 - Icges(3,4) * t360 + Icges(3,5) * t465;
t301 = Icges(3,1) * t359 + Icges(3,4) * t390 - Icges(3,5) * t464;
t300 = Icges(3,4) * t361 - Icges(3,2) * t360 + Icges(3,6) * t465;
t299 = Icges(3,4) * t359 + Icges(3,2) * t390 - Icges(3,6) * t464;
t296 = Icges(4,5) * t357 + Icges(4,6) * t356 - Icges(4,3) * t462;
t276 = rSges(4,3) * t351 - t350 * t401;
t275 = rSges(4,3) * t349 + t348 * t401;
t274 = Icges(4,5) * t351 - t350 * t398;
t273 = Icges(4,5) * t349 + t348 * t398;
t272 = Icges(4,6) * t351 - t350 * t397;
t271 = Icges(4,6) * t349 + t348 * t397;
t270 = Icges(4,3) * t351 - t350 * t396;
t269 = Icges(4,3) * t349 + t348 * t396;
t260 = t267 * t409;
t259 = rSges(4,1) * t337 + rSges(4,2) * t336 + rSges(4,3) * t360;
t258 = rSges(4,1) * t335 + rSges(4,2) * t334 - rSges(4,3) * t390;
t252 = Icges(4,5) * t337 + Icges(4,6) * t336 + Icges(4,3) * t360;
t251 = Icges(4,5) * t335 + Icges(4,6) * t334 - Icges(4,3) * t390;
t244 = t360 * t266;
t241 = rSges(5,1) * t327 + rSges(5,2) * t384 + rSges(5,3) * t360;
t240 = rSges(5,1) * t325 + rSges(5,2) * t383 - rSges(5,3) * t390;
t215 = (t317 * t373 + t318 * t376) * t374;
t201 = rSges(5,1) * t285 - rSges(5,2) * t284 + rSges(5,3) * t351;
t200 = rSges(5,1) * t283 - rSges(5,2) * t282 + rSges(5,3) * t349;
t191 = rSges(6,1) * t289 + rSges(6,2) * t389 - rSges(6,3) * t383;
t177 = -t241 * t462 - t360 * t295;
t176 = t240 * t462 - t295 * t390;
t175 = (-t258 - t320) * t377 + t376 * t407;
t174 = t259 * t377 + t373 * t407 + t319;
t173 = -t292 * t462 - t343 * t293 + t344 * t294;
t172 = (-t275 - t286) * t377 + t376 * t406;
t171 = t276 * t377 + t373 * t406 + t280;
t161 = t240 * t360 + t241 * t390;
t158 = (t258 * t373 + t259 * t376) * t374 + t441;
t157 = (t275 * t373 + t276 * t376) * t374 + t443;
t156 = t193 * t343 + t229 * t384;
t155 = -t191 * t343 - t229 * t383;
t154 = -t233 * t462 - t343 * t235 + t344 * t237;
t153 = -t232 * t462 - t343 * t234 + t344 * t236;
t146 = (-t240 + t445) * t377 + t376 * t402;
t145 = t241 * t377 + t373 * t402 + t449;
t142 = rSges(6,1) * t203 - rSges(6,2) * t202 + rSges(6,3) * t282;
t126 = -t191 * t384 + t193 * t383;
t125 = t360 * t447 + t452 * t462;
t124 = t191 * t462 - t229 * t390 + t444;
t123 = -t360 * t250 - t351 * t295 + (-t201 * t381 + t241 * t438) * t374;
t122 = -t390 * t250 + t349 * t295 + (t200 * t381 - t240 * t438) * t374;
t117 = (t240 * t373 + t241 * t376) * t374 + t404;
t116 = (-t200 + t446) * t377 + t376 * t403;
t115 = t201 * t377 + t373 * t403 + t450;
t114 = t191 * t360 - t390 * t452 + t244;
t113 = -t343 * t248 + t344 * t249 - t322 * t293 + t323 * t294 + (-t247 * t381 + t292 * t438) * t374;
t112 = t200 * t360 + t201 * t390 + t240 * t351 - t241 * t349;
t111 = t343 * t453 + t384 * t448;
t110 = -t343 * t454 - t383 * t448;
t105 = (-t191 + t413) * t377 + t376 * t393;
t104 = t193 * t377 + t373 * t393 + t415;
t95 = t360 * t414 + t419 * t462;
t94 = -t390 * t448 + t454 * t462 + t444;
t93 = (t200 * t373 + t201 * t376) * t374 + t405;
t90 = t383 * t453 - t384 * t454;
t89 = (t191 * t373 + t193 * t376) * t374 + t391;
t88 = (t413 - t454) * t377 + t376 * t386;
t87 = t373 * t386 + t377 * t453 + t415;
t86 = t360 * t454 - t390 * t419 + t244;
t85 = -t343 * t197 + t344 * t199 - t322 * t235 + t323 * t237 + (-t195 * t381 + t233 * t438) * t374;
t84 = -t343 * t196 + t344 * t198 - t322 * t234 + t323 * t236 + (-t194 * t381 + t232 * t438) * t374;
t79 = t144 * t343 + t169 * t384 + t193 * t322 - t229 * t284;
t78 = -t142 * t343 - t169 * t383 - t191 * t322 + t229 * t282;
t77 = (t373 * t454 + t376 * t453) * t374 + t391;
t76 = (-t142 + t416) * t377 + t376 * t394;
t75 = t144 * t377 + t373 * t394 + t417;
t72 = -t142 * t384 + t144 * t383 + t191 * t284 - t193 * t282;
t71 = t260 + (-t169 - t265) * t360 + t447 * t351 + (t193 * t438 + t381 * t456) * t374;
t70 = -t390 * t169 + t349 * t229 + (t142 * t381 + (-t191 - t266) * t438) * t374 + t418;
t69 = (t142 * t373 + t144 * t376) * t374 + t392;
t62 = -t108 * t390 + t109 * t360 - t128 * t462;
t61 = -t106 * t390 + t107 * t360 - t127 * t462;
t60 = -t108 * t383 - t109 * t384 + t128 * t343;
t59 = -t106 * t383 - t107 * t384 + t127 * t343;
t58 = t142 * t360 + t191 * t351 + t349 * t452 - t390 * t456 + t451;
t57 = (t416 - t458) * t377 + t376 * t387;
t56 = t373 * t387 + t377 * t457 + t417;
t51 = -t102 * t390 + t103 * t360 - t121 * t462;
t50 = -t100 * t390 + t101 * t360 - t120 * t462;
t49 = -t119 * t462 + t99 * t360 - t390 * t98;
t48 = -t118 * t462 + t97 * t360 - t390 * t96;
t47 = -t102 * t383 - t103 * t384 + t121 * t343;
t46 = -t100 * t383 - t101 * t384 + t120 * t343;
t45 = t119 * t343 - t383 * t98 - t384 * t99;
t44 = t118 * t343 - t383 * t96 - t384 * t97;
t39 = -t284 * t448 + t322 * t453 + t343 * t457 + t384 * t455;
t38 = t282 * t448 - t322 * t454 - t343 * t458 - t383 * t455;
t37 = t260 + (-t265 - t455) * t360 + t414 * t351 + (t381 * t420 + t438 * t453) * t374;
t36 = -t455 * t390 + t448 * t349 + (t458 * t381 + (-t266 - t454) * t438) * t374 + t418;
t35 = (t373 * t458 + t376 * t457) * t374 + t392;
t26 = t113 * t377 + (t373 * t85 - t376 * t84) * t374;
t25 = -t282 * t453 + t284 * t454 + t383 * t457 - t384 * t458;
t24 = t377 * t92 + (t373 * t83 - t376 * t82) * t374;
t23 = t377 * t91 + (t373 * t81 - t376 * t80) * t374;
t22 = t349 * t419 + t351 * t454 + t360 * t458 - t390 * t420 + t451;
t21 = t153 * t349 + t154 * t351 - t84 * t390 + t85 * t360 + (-t113 * t381 + t173 * t438) * t374;
t19 = [0; m(3) * t215 + m(4) * t157 + m(5) * t93 + m(6) * t69 + m(7) * t35; t16 * t465 + t15 * t465 - t13 * t464 - t14 * t464 - t23 * t464 + t24 * t465 - ((-t300 * t349 + t302 * t348 - t311 * t464 + t313 * t390 + t315 * t359) * t465 - (-t299 * t349 + t301 * t348 - t310 * t464 + t312 * t390 + t314 * t359) * t464 + (-t341 * t349 + t342 * t348 - t352 * t464 + t353 * t390 + t354 * t359) * t377) * t464 - ((t349 * t296 - t328 * t390 + t334 * t329 + t335 * t330 + t348 * t395) * t377 + ((t252 * t349 - t270 * t390 + t272 * t334 + t274 * t335) * t373 - (t251 * t349 - t269 * t390 + t271 * t334 + t273 * t335) * t376 + t382 * t348) * t374) * t464 + ((-t300 * t351 - t302 * t350 + t311 * t465 - t313 * t360 + t315 * t361) * t465 - (-t299 * t351 - t301 * t350 + t310 * t465 - t312 * t360 + t314 * t361) * t464 + (-t341 * t351 - t342 * t350 + t352 * t465 - t353 * t360 + t354 * t361) * t377) * t465 + ((t351 * t296 + t360 * t328 + t336 * t329 + t337 * t330 - t350 * t395) * t377 + ((t252 * t351 + t270 * t360 + t272 * t336 + t274 * t337) * t373 - (t251 * t351 + t269 * t360 + t271 * t336 + t273 * t337) * t376 - t382 * t350) * t374) * t465 + (t35 * t77 + t56 * t87 + t57 * t88) * t475 + t377 * t18 + (t104 * t75 + t105 * t76 + t69 * t89) * t476 + t377 * t17 + t377 * t26 + (t115 * t145 + t116 * t146 + t117 * t93) * t477 + 0.2e1 * m(4) * (t157 * t158 + t171 * t174 + t172 * t175) + t377 * (t377 ^ 2 * t352 + (((t313 * t381 + t315 * t380) * t373 - (t312 * t381 + t314 * t380) * t376 + ((-t300 * t380 + t302 * t381) * t373 - (-t299 * t380 + t301 * t381) * t376) * qJD(2)) * t374 + (-t310 * t376 + t311 * t373 + t353 * t381 + t354 * t380 + (-t341 * t380 + t342 * t381) * qJD(2)) * t377) * t374) + t377 * ((t356 * t329 + t357 * t330) * t377 + ((t356 * t272 + t357 * t274) * t373 - (t356 * t271 + t357 * t273) * t376 + (-t328 * t377 + (t269 * t376 - t270 * t373) * t374) * t381 + ((t296 * t377 + (-t251 * t376 + t252 * t373) * t374) * t380 + (t374 * t382 + t377 * t395) * t381) * qJD(2)) * t374) + 0.2e1 * m(3) * ((-t306 * t377 - t345 * t464) * (-t317 * t377 - t355 * t464) + (t307 * t377 - t345 * t465) * (t318 * t377 - t355 * t465) + (t306 * t373 + t307 * t376) * t374 * t215); (m(4) + m(5) + m(6) + m(7)) * t409; m(7) * (t349 * t87 + t351 * t88 - t390 * t56 + t360 * t57 + (-t35 * t381 + t438 * t77) * t374) + m(6) * (t349 * t104 + t351 * t105 - t390 * t75 + t360 * t76 + (-t381 * t69 + t438 * t89) * t374) + m(5) * (-t390 * t115 + t360 * t116 + t349 * t145 + t351 * t146 + (t117 * t438 - t381 * t93) * t374) + m(4) * (-t390 * t171 + t360 * t172 + t349 * t174 + t351 * t175 + (-t157 * t381 + t158 * t438) * t374); 0.4e1 * (m(7) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * (-t374 ^ 2 * t380 * t437 - t349 * t390 + t360 * t351); m(5) * t112 + m(6) * t58 + m(7) * t22; t427 * t351 + t428 * t349 + (t24 / 0.2e1 + t430) * t360 - (t23 / 0.2e1 + t431) * t390 + m(5) * (t112 * t117 + t115 * t177 + t116 * t176 + t122 * t146 + t123 * t145 + t161 * t93) + m(6) * (t104 * t71 + t105 * t70 + t114 * t69 + t124 * t76 + t125 * t75 + t58 * t89) + m(7) * (t22 * t77 + t35 * t86 + t36 * t88 + t37 * t87 + t56 * t95 + t57 * t94) + (t467 / 0.2e1 + t468 / 0.2e1 + t21 / 0.2e1 + t432) * t377 + (t351 * (-t149 * t376 + t150 * t373) / 0.2e1 + t349 * (-t147 * t376 + t148 * t373) / 0.2e1 + (-t26 / 0.2e1 - t429) * t381 + (t173 * t480 + (-t153 * t376 + t154 * t373) * t481 + t426) * t438) * t374 + t478 * t465 / 0.2e1 - t479 * t464 / 0.2e1; m(7) * (t95 * t349 + t94 * t351 - t37 * t390 + t36 * t360 + (-t22 * t381 + t438 * t86) * t374) + m(5) * (t122 * t360 - t123 * t390 + t176 * t351 + t177 * t349 + (-t112 * t381 + t161 * t438) * t374) + m(6) * (t124 * t351 + t125 * t349 - t71 * t390 + t70 * t360 + (t114 * t438 - t381 * t58) * t374); t478 * t360 - t479 * t390 + (-t149 * t390 + t150 * t360 + t50 + t51) * t351 + (-t147 * t390 + t148 * t360 + t48 + t49) * t349 + (t22 * t86 + t36 * t94 + t37 * t95) * t475 + (t114 * t58 + t124 * t70 + t125 * t71) * t476 + (t112 * t161 + t122 * t176 + t123 * t177) * t477 + ((-t153 * t390 + t154 * t360 + t61 + t62) * t438 + (-t173 * t409 - t11 - t12 - t21 - t467 - t468) * t381) * t374; m(6) * t72 + m(7) * t25; -t433 * t377 + t429 * t343 - t430 * t384 - t431 * t383 + t426 * t322 + t427 * t284 + t428 * t282 + m(7) * (t110 * t57 + t111 * t56 + t25 * t77 + t35 * t90 + t38 * t88 + t39 * t87) + m(6) * (t104 * t79 + t105 * t78 + t126 * t69 + t155 * t76 + t156 * t75 + t72 * t89) + (t373 * t434 - t376 * t435) * t374; m(7) * (t110 * t351 + t111 * t349 - t39 * t390 + t38 * t360 + (-t25 * t381 + t438 * t90) * t374) + m(6) * (t155 * t351 + t156 * t349 - t79 * t390 + t78 * t360 + (t126 * t438 - t381 * t72) * t374); t434 * t360 - t435 * t390 + (t47 / 0.2e1 + t46 / 0.2e1) * t351 + (t45 / 0.2e1 + t44 / 0.2e1) * t349 + t432 * t343 - (t7 / 0.2e1 + t8 / 0.2e1) * t384 - (t5 / 0.2e1 + t6 / 0.2e1) * t383 + (t62 / 0.2e1 + t61 / 0.2e1) * t322 + (t51 / 0.2e1 + t50 / 0.2e1) * t284 + (t48 / 0.2e1 + t49 / 0.2e1) * t282 + m(6) * (t114 * t72 + t124 * t78 + t125 * t79 + t126 * t58 + t155 * t70 + t156 * t71) + m(7) * (t110 * t36 + t111 * t37 + t22 * t90 + t25 * t86 + t38 * t94 + t39 * t95) + (t433 * t381 + (t59 / 0.2e1 + t60 / 0.2e1) * t438) * t374; (t9 + t10) * t343 - (t4 + t3) * t384 - (t1 + t2) * t383 + (t59 + t60) * t322 + (t47 + t46) * t284 + (t44 + t45) * t282 + (t110 * t38 + t111 * t39 + t25 * t90) * t475 + (t126 * t72 + t155 * t78 + t156 * t79) * t476; m(7) * t245; m(7) * (t202 * t87 + t204 * t88 + t245 * t77 - t35 * t385 - t388 * t57 - t389 * t56); m(7) * (-t202 * t390 + t204 * t360 - t389 * t349 - t388 * t351 + (-t245 * t381 - t385 * t438) * t374); m(7) * (t202 * t95 + t204 * t94 - t22 * t385 + t245 * t86 - t36 * t388 - t37 * t389); m(7) * (t110 * t204 + t111 * t202 + t245 * t90 - t25 * t385 - t38 * t388 - t389 * t39); (-t202 * t389 - t204 * t388 - t245 * t385) * t475;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
