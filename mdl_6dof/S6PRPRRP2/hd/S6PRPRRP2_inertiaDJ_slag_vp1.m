% Calculate time derivative of joint inertia matrix for
% S6PRPRRP2
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:00:01
% EndTime: 2019-03-08 20:00:34
% DurationCPUTime: 22.74s
% Computational Cost: add. (92724->1090), mult. (263149->1528), div. (0->0), fcn. (315061->12), ass. (0->428)
t530 = rSges(7,1) + pkin(5);
t529 = rSges(7,3) + qJ(6);
t511 = m(7) / 0.2e1;
t531 = 0.2e1 * t511;
t496 = sin(pkin(11));
t497 = cos(pkin(11));
t501 = sin(qJ(2));
t504 = cos(qJ(2));
t398 = t501 * t496 - t504 * t497;
t376 = t398 * qJD(2);
t388 = sin(pkin(10));
t390 = cos(pkin(10));
t498 = cos(pkin(6));
t418 = t498 * t496;
t419 = t498 * t497;
t362 = t418 * t504 + t419 * t501;
t393 = qJD(2) * t362;
t314 = t388 * t376 - t390 * t393;
t394 = -t418 * t501 + t419 * t504;
t356 = t394 * qJD(2);
t399 = t496 * t504 + t497 * t501;
t377 = t399 * qJD(2);
t315 = t356 * t390 - t377 * t388;
t424 = t498 * t504;
t372 = -t388 * t501 + t390 * t424;
t363 = t372 * qJD(2);
t423 = t498 * t501;
t400 = -t388 * t504 - t390 * t423;
t364 = t400 * qJD(2);
t528 = -Icges(3,5) * t363 - Icges(4,5) * t315 - Icges(3,6) * t364 - Icges(4,6) * t314;
t316 = t376 * t390 + t388 * t393;
t317 = -t356 * t388 - t377 * t390;
t374 = -t388 * t424 - t390 * t501;
t365 = t374 * qJD(2);
t401 = t388 * t423 - t390 * t504;
t366 = t401 * qJD(2);
t527 = Icges(3,5) * t365 + Icges(4,5) * t317 + Icges(3,6) * t366 + Icges(4,6) * t316;
t389 = sin(pkin(6));
t361 = t399 * t389;
t354 = qJD(2) * t361;
t471 = qJD(2) * t389;
t355 = t398 * t471;
t526 = -Icges(4,5) * t355 - Icges(4,6) * t354 + (Icges(3,5) * t504 - Icges(3,6) * t501) * t471;
t338 = t362 * t390 - t388 * t398;
t392 = sin(qJ(4));
t503 = cos(qJ(4));
t452 = t389 * t503;
t406 = -t338 * t392 - t390 * t452;
t254 = qJD(4) * t406 + t315 * t503;
t491 = t389 * t392;
t298 = t338 * t503 - t390 * t491;
t337 = -t388 * t399 + t390 * t394;
t391 = sin(qJ(5));
t502 = cos(qJ(5));
t263 = t298 * t502 - t337 * t391;
t193 = qJD(5) * t263 + t254 * t391 + t314 * t502;
t412 = -t298 * t391 - t337 * t502;
t194 = qJD(5) * t412 + t254 * t502 - t314 * t391;
t253 = qJD(4) * t298 + t315 * t392;
t123 = Icges(7,5) * t194 + Icges(7,6) * t253 + Icges(7,3) * t193;
t127 = Icges(7,4) * t194 + Icges(7,2) * t253 + Icges(7,6) * t193;
t131 = Icges(7,1) * t194 + Icges(7,4) * t253 + Icges(7,5) * t193;
t181 = Icges(7,5) * t263 - Icges(7,6) * t406 - Icges(7,3) * t412;
t185 = Icges(7,4) * t263 - Icges(7,2) * t406 - Icges(7,6) * t412;
t189 = Icges(7,1) * t263 - Icges(7,4) * t406 - Icges(7,5) * t412;
t340 = -t362 * t388 - t390 * t398;
t407 = -t340 * t392 + t388 * t452;
t256 = qJD(4) * t407 + t317 * t503;
t300 = t340 * t503 + t388 * t491;
t339 = -t388 * t394 - t390 * t399;
t265 = t300 * t502 - t339 * t391;
t195 = qJD(5) * t265 + t256 * t391 + t316 * t502;
t411 = -t300 * t391 - t339 * t502;
t196 = qJD(5) * t411 + t256 * t502 - t316 * t391;
t255 = qJD(4) * t300 + t317 * t392;
t33 = -t123 * t411 - t127 * t407 + t131 * t265 + t181 * t195 + t185 * t255 + t189 * t196;
t124 = Icges(7,5) * t196 + Icges(7,6) * t255 + Icges(7,3) * t195;
t128 = Icges(7,4) * t196 + Icges(7,2) * t255 + Icges(7,6) * t195;
t132 = Icges(7,1) * t196 + Icges(7,4) * t255 + Icges(7,5) * t195;
t182 = Icges(7,5) * t265 - Icges(7,6) * t407 - Icges(7,3) * t411;
t186 = Icges(7,4) * t265 - Icges(7,2) * t407 - Icges(7,6) * t411;
t190 = Icges(7,1) * t265 - Icges(7,4) * t407 - Icges(7,5) * t411;
t34 = -t124 * t411 - t128 * t407 + t132 * t265 + t182 * t195 + t186 * t255 + t190 * t196;
t405 = -t361 * t392 + t498 * t503;
t294 = qJD(4) * t405 - t355 * t503;
t347 = t361 * t503 + t392 * t498;
t360 = t398 * t389;
t296 = t347 * t502 + t360 * t391;
t237 = t296 * qJD(5) + t294 * t391 - t354 * t502;
t410 = -t347 * t391 + t360 * t502;
t238 = qJD(5) * t410 + t294 * t502 + t354 * t391;
t293 = qJD(4) * t347 - t355 * t392;
t160 = Icges(7,5) * t238 + Icges(7,6) * t293 + Icges(7,3) * t237;
t162 = Icges(7,4) * t238 + Icges(7,2) * t293 + Icges(7,6) * t237;
t164 = Icges(7,1) * t238 + Icges(7,4) * t293 + Icges(7,5) * t237;
t229 = Icges(7,5) * t296 - Icges(7,6) * t405 - Icges(7,3) * t410;
t231 = Icges(7,4) * t296 - Icges(7,2) * t405 - Icges(7,6) * t410;
t233 = Icges(7,1) * t296 - Icges(7,4) * t405 - Icges(7,5) * t410;
t67 = -t160 * t411 - t162 * t407 + t164 * t265 + t195 * t229 + t196 * t233 + t231 * t255;
t15 = t67 * t498 + (-t33 * t390 + t34 * t388) * t389;
t125 = Icges(6,5) * t194 - Icges(6,6) * t193 + Icges(6,3) * t253;
t129 = Icges(6,4) * t194 - Icges(6,2) * t193 + Icges(6,6) * t253;
t133 = Icges(6,1) * t194 - Icges(6,4) * t193 + Icges(6,5) * t253;
t183 = Icges(6,5) * t263 + Icges(6,6) * t412 - Icges(6,3) * t406;
t187 = Icges(6,4) * t263 + Icges(6,2) * t412 - Icges(6,6) * t406;
t191 = Icges(6,1) * t263 + Icges(6,4) * t412 - Icges(6,5) * t406;
t35 = -t125 * t407 + t129 * t411 + t133 * t265 + t183 * t255 - t187 * t195 + t191 * t196;
t126 = Icges(6,5) * t196 - Icges(6,6) * t195 + Icges(6,3) * t255;
t130 = Icges(6,4) * t196 - Icges(6,2) * t195 + Icges(6,6) * t255;
t134 = Icges(6,1) * t196 - Icges(6,4) * t195 + Icges(6,5) * t255;
t184 = Icges(6,5) * t265 + Icges(6,6) * t411 - Icges(6,3) * t407;
t188 = Icges(6,4) * t265 + Icges(6,2) * t411 - Icges(6,6) * t407;
t192 = Icges(6,1) * t265 + Icges(6,4) * t411 - Icges(6,5) * t407;
t36 = -t126 * t407 + t130 * t411 + t134 * t265 + t184 * t255 - t188 * t195 + t192 * t196;
t161 = Icges(6,5) * t238 - Icges(6,6) * t237 + Icges(6,3) * t293;
t163 = Icges(6,4) * t238 - Icges(6,2) * t237 + Icges(6,6) * t293;
t165 = Icges(6,1) * t238 - Icges(6,4) * t237 + Icges(6,5) * t293;
t230 = Icges(6,5) * t296 + Icges(6,6) * t410 - Icges(6,3) * t405;
t232 = Icges(6,4) * t296 + Icges(6,2) * t410 - Icges(6,6) * t405;
t234 = Icges(6,1) * t296 + Icges(6,4) * t410 - Icges(6,5) * t405;
t68 = -t161 * t407 + t163 * t411 + t165 * t265 - t195 * t232 + t196 * t234 + t230 * t255;
t16 = t68 * t498 + (-t35 * t390 + t36 * t388) * t389;
t173 = Icges(5,5) * t254 - Icges(5,6) * t253 - Icges(5,3) * t314;
t175 = Icges(5,4) * t254 - Icges(5,2) * t253 - Icges(5,6) * t314;
t177 = Icges(5,1) * t254 - Icges(5,4) * t253 - Icges(5,5) * t314;
t220 = Icges(5,5) * t298 + Icges(5,6) * t406 - Icges(5,3) * t337;
t222 = Icges(5,4) * t298 + Icges(5,2) * t406 - Icges(5,6) * t337;
t224 = Icges(5,1) * t298 + Icges(5,4) * t406 - Icges(5,5) * t337;
t77 = -t173 * t339 + t175 * t407 + t177 * t300 - t220 * t316 - t222 * t255 + t224 * t256;
t174 = Icges(5,5) * t256 - Icges(5,6) * t255 - Icges(5,3) * t316;
t176 = Icges(5,4) * t256 - Icges(5,2) * t255 - Icges(5,6) * t316;
t178 = Icges(5,1) * t256 - Icges(5,4) * t255 - Icges(5,5) * t316;
t221 = Icges(5,5) * t300 + Icges(5,6) * t407 - Icges(5,3) * t339;
t223 = Icges(5,4) * t300 + Icges(5,2) * t407 - Icges(5,6) * t339;
t225 = Icges(5,1) * t300 + Icges(5,4) * t407 - Icges(5,5) * t339;
t78 = -t174 * t339 + t176 * t407 + t178 * t300 - t221 * t316 - t223 * t255 + t225 * t256;
t239 = Icges(5,5) * t294 - Icges(5,6) * t293 + Icges(5,3) * t354;
t240 = Icges(5,4) * t294 - Icges(5,2) * t293 + Icges(5,6) * t354;
t241 = Icges(5,1) * t294 - Icges(5,4) * t293 + Icges(5,5) * t354;
t288 = Icges(5,5) * t347 + Icges(5,6) * t405 + Icges(5,3) * t360;
t289 = Icges(5,4) * t347 + Icges(5,2) * t405 + Icges(5,6) * t360;
t290 = Icges(5,1) * t347 + Icges(5,4) * t405 + Icges(5,5) * t360;
t91 = -t239 * t339 + t240 * t407 + t241 * t300 - t255 * t289 + t256 * t290 - t288 * t316;
t525 = t16 + t15 + t91 * t498 + (t388 * t78 - t390 * t77) * t389;
t109 = t239 * t360 + t240 * t405 + t241 * t347 + t288 * t354 - t289 * t293 + t290 * t294;
t40 = -t123 * t410 - t127 * t405 + t131 * t296 + t181 * t237 + t185 * t293 + t189 * t238;
t41 = -t124 * t410 - t128 * t405 + t132 * t296 + t182 * t237 + t186 * t293 + t190 * t238;
t73 = -t160 * t410 - t162 * t405 + t164 * t296 + t229 * t237 + t231 * t293 + t233 * t238;
t17 = t73 * t498 + (t388 * t41 - t390 * t40) * t389;
t42 = -t125 * t405 + t129 * t410 + t133 * t296 + t183 * t293 - t187 * t237 + t191 * t238;
t43 = -t126 * t405 + t130 * t410 + t134 * t296 + t184 * t293 - t188 * t237 + t192 * t238;
t74 = -t161 * t405 + t163 * t410 + t165 * t296 + t230 * t293 - t232 * t237 + t234 * t238;
t18 = t74 * t498 + (t388 * t43 - t390 * t42) * t389;
t84 = t173 * t360 + t175 * t405 + t177 * t347 + t220 * t354 - t222 * t293 + t224 * t294;
t85 = t174 * t360 + t176 * t405 + t178 * t347 + t221 * t354 - t223 * t293 + t225 * t294;
t524 = t17 + t18 + t109 * t498 + (t388 * t85 - t390 * t84) * t389;
t29 = -t123 * t412 - t127 * t406 + t131 * t263 + t181 * t193 + t185 * t253 + t189 * t194;
t30 = -t124 * t412 - t128 * t406 + t132 * t263 + t182 * t193 + t186 * t253 + t190 * t194;
t65 = -t160 * t412 - t162 * t406 + t164 * t263 + t193 * t229 + t194 * t233 + t231 * t253;
t13 = t65 * t498 + (-t29 * t390 + t30 * t388) * t389;
t31 = -t125 * t406 + t129 * t412 + t133 * t263 + t183 * t253 - t187 * t193 + t191 * t194;
t32 = -t126 * t406 + t130 * t412 + t134 * t263 + t184 * t253 - t188 * t193 + t192 * t194;
t66 = -t161 * t406 + t163 * t412 + t165 * t263 - t193 * t232 + t194 * t234 + t230 * t253;
t14 = t66 * t498 + (-t31 * t390 + t32 * t388) * t389;
t75 = -t173 * t337 + t175 * t406 + t177 * t298 - t220 * t314 - t222 * t253 + t224 * t254;
t76 = -t174 * t337 + t176 * t406 + t178 * t298 - t221 * t314 - t223 * t253 + t225 * t254;
t90 = -t239 * t337 + t240 * t406 + t241 * t298 - t253 * t289 + t254 * t290 - t288 * t314;
t523 = t90 * t498 + (t388 * t76 - t390 * t75) * t389 + t14 + t13;
t522 = t389 ^ 2;
t106 = -t183 * t405 + t187 * t410 + t191 * t296;
t107 = -t184 * t405 + t188 * t410 + t192 * t296;
t150 = -t230 * t405 + t232 * t410 + t234 * t296;
t10 = t106 * t253 + t107 * t255 + t150 * t293 - t405 * t74 - t406 * t42 - t407 * t43;
t104 = -t181 * t410 - t185 * t405 + t189 * t296;
t105 = -t182 * t410 - t186 * t405 + t190 * t296;
t149 = -t229 * t410 - t231 * t405 + t233 * t296;
t9 = t104 * t253 + t105 * t255 + t149 * t293 - t40 * t406 - t405 * t73 - t407 * t41;
t521 = t9 + t10;
t490 = rSges(7,2) * t253 - qJD(6) * t412 + t529 * t193 + t194 * t530;
t520 = t490 * t531;
t489 = rSges(7,2) * t255 - qJD(6) * t411 + t529 * t195 + t196 * t530;
t487 = -rSges(7,2) * t406 + t263 * t530 - t529 * t412;
t485 = -rSges(7,2) * t407 + t265 * t530 - t529 * t411;
t210 = pkin(4) * t254 + pkin(9) * t253;
t250 = pkin(4) * t298 - pkin(9) * t406;
t483 = -t339 * t210 - t316 * t250;
t379 = pkin(2) * qJD(2) * t424 - t389 * qJD(3);
t446 = qJD(2) * t501;
t432 = pkin(2) * t446;
t348 = t390 * t379 - t388 * t432;
t349 = -t388 * t379 - t390 * t432;
t492 = t389 * t390;
t493 = t388 * t389;
t472 = t348 * t493 + t349 * t492;
t512 = m(6) / 0.2e1;
t519 = t512 + t511;
t144 = -t220 * t337 + t222 * t406 + t224 * t298;
t145 = -t221 * t337 + t223 * t406 + t225 * t298;
t157 = -t288 * t337 + t289 * t406 + t290 * t298;
t119 = -t229 * t412 - t231 * t406 + t233 * t263;
t95 = -t181 * t412 - t185 * t406 + t189 * t263;
t96 = -t182 * t412 - t186 * t406 + t190 * t263;
t5 = t119 * t354 - t29 * t337 - t30 * t339 - t314 * t95 - t316 * t96 + t360 * t65;
t120 = -t230 * t406 + t232 * t412 + t234 * t263;
t97 = -t183 * t406 + t187 * t412 + t191 * t263;
t98 = -t184 * t406 + t188 * t412 + t192 * t263;
t6 = t120 * t354 - t31 * t337 - t314 * t97 - t316 * t98 - t32 * t339 + t360 * t66;
t518 = -t314 * t144 - t145 * t316 + t157 * t354 - t337 * t75 - t339 * t76 + t360 * t90 + t5 + t6;
t146 = -t220 * t339 + t222 * t407 + t224 * t300;
t147 = -t221 * t339 + t223 * t407 + t225 * t300;
t158 = -t288 * t339 + t289 * t407 + t290 * t300;
t100 = -t182 * t411 - t186 * t407 + t190 * t265;
t121 = -t229 * t411 - t231 * t407 + t233 * t265;
t99 = -t181 * t411 - t185 * t407 + t189 * t265;
t7 = -t100 * t316 + t121 * t354 - t314 * t99 - t33 * t337 - t339 * t34 + t360 * t67;
t101 = -t183 * t407 + t187 * t411 + t191 * t265;
t102 = -t184 * t407 + t188 * t411 + t192 * t265;
t122 = -t230 * t407 + t232 * t411 + t234 * t265;
t8 = -t101 * t314 - t102 * t316 + t122 * t354 - t337 * t35 - t339 * t36 + t360 * t68;
t517 = -t146 * t314 - t316 * t147 + t158 * t354 - t337 * t77 - t339 * t78 + t360 * t91 + t7 + t8;
t11 = -t104 * t314 - t105 * t316 + t149 * t354 - t337 * t40 - t339 * t41 + t360 * t73;
t12 = -t106 * t314 - t107 * t316 + t150 * t354 - t337 * t42 - t339 * t43 + t360 * t74;
t151 = t220 * t360 + t222 * t405 + t224 * t347;
t152 = t221 * t360 + t223 * t405 + t225 * t347;
t169 = t288 * t360 + t289 * t405 + t290 * t347;
t495 = t169 * t354;
t516 = t109 * t360 - t151 * t314 - t152 * t316 - t337 * t84 - t339 * t85 + t11 + t12 + t495;
t515 = -t485 * t253 + t487 * t255 + t406 * t489;
t514 = m(4) / 0.2e1;
t513 = m(5) / 0.2e1;
t500 = pkin(2) * t389;
t499 = pkin(2) * t504;
t488 = rSges(7,2) * t293 - qJD(6) * t410 + t237 * t529 + t238 * t530;
t198 = rSges(6,1) * t263 + rSges(6,2) * t412 - rSges(6,3) * t406;
t486 = -t198 - t250;
t200 = rSges(6,1) * t265 + rSges(6,2) * t411 - rSges(6,3) * t407;
t251 = pkin(4) * t300 - pkin(9) * t407;
t484 = t200 + t251;
t211 = pkin(4) * t256 + pkin(9) * t255;
t482 = t360 * t211 + t354 * t251;
t245 = pkin(4) * t294 + pkin(9) * t293;
t292 = pkin(4) * t347 - pkin(9) * t405;
t481 = -t337 * t245 - t314 * t292;
t480 = -rSges(7,2) * t405 + t296 * t530 - t410 * t529;
t236 = rSges(6,1) * t296 + rSges(6,2) * t410 - rSges(6,3) * t405;
t479 = t236 + t292;
t276 = pkin(3) * t317 - pkin(8) * t316;
t345 = t498 * t349;
t478 = t498 * t276 + t345;
t287 = pkin(3) * t340 - pkin(8) * t339;
t445 = pkin(2) * t423 - qJ(3) * t389;
t327 = -t388 * t445 + t390 * t499;
t313 = t498 * t327;
t477 = t498 * t287 + t313;
t326 = t388 * t499 + t390 * t445;
t476 = t326 * t493 + t327 * t492;
t447 = qJD(2) * t504;
t378 = qJD(3) * t498 + t447 * t500;
t475 = pkin(3) * t355 - pkin(8) * t354 - t378;
t382 = qJ(3) * t498 + t500 * t501;
t474 = -pkin(3) * t361 - pkin(8) * t360 - t382;
t473 = 0.2e1 * t472;
t468 = 0.2e1 * t498;
t467 = 0.2e1 * t389;
t1 = t119 * t293 + t253 * t95 + t255 * t96 - t29 * t406 - t30 * t407 - t405 * t65;
t2 = t120 * t293 + t253 * t97 + t255 * t98 - t31 * t406 - t32 * t407 - t405 * t66;
t466 = -t2 / 0.2e1 - t1 / 0.2e1;
t3 = t100 * t255 + t121 * t293 + t253 * t99 - t33 * t406 - t34 * t407 - t405 * t67;
t4 = t101 * t253 + t102 * t255 + t122 * t293 - t35 * t406 - t36 * t407 - t405 * t68;
t465 = -t4 / 0.2e1 - t3 / 0.2e1;
t463 = t498 / 0.2e1;
t462 = -t250 - t487;
t461 = t251 + t485;
t460 = t498 * t211 + t478;
t459 = t292 + t480;
t458 = -t245 + t475;
t457 = t498 * t251 + t477;
t456 = -t292 + t474;
t455 = m(6) * t498;
t454 = m(7) * t498;
t451 = t504 * Icges(3,4);
t450 = t501 * Icges(3,4);
t444 = t498 * t326;
t443 = t498 * t348;
t442 = 0.2e1 * m(5);
t440 = 0.2e1 * m(6);
t438 = 0.2e1 * m(7);
t437 = t389 * (rSges(4,1) * t355 + rSges(4,2) * t354 - t378);
t436 = t389 * (-t361 * rSges(4,1) + t360 * rSges(4,2) - rSges(4,3) * t498 - t382);
t275 = pkin(3) * t315 - pkin(8) * t314;
t260 = t275 * t493;
t261 = t276 * t492;
t430 = 0.2e1 * t260 + 0.2e1 * t261 + t473;
t429 = t260 + t261 + t472;
t286 = pkin(3) * t338 - pkin(8) * t337;
t428 = t286 * t493 + t287 * t492 + t476;
t242 = rSges(5,1) * t294 - rSges(5,2) * t293 + rSges(5,3) * t354;
t426 = (-t242 + t475) * t389;
t291 = rSges(5,1) * t347 + rSges(5,2) * t405 + rSges(5,3) * t360;
t425 = (-t291 + t474) * t389;
t136 = rSges(6,1) * t194 - rSges(6,2) * t193 + rSges(6,3) * t253;
t179 = rSges(5,1) * t254 - rSges(5,2) * t253 - rSges(5,3) * t314;
t420 = -m(5) * t179 - m(6) * t136;
t167 = rSges(6,1) * t238 - rSges(6,2) * t237 + rSges(6,3) * t293;
t417 = (-t167 + t458) * t389;
t416 = (-t236 + t456) * t389;
t207 = t210 * t493;
t208 = t211 * t492;
t414 = t207 + t208 + t429;
t413 = t250 * t493 + t251 * t492 + t428;
t409 = t389 * (t458 - t488);
t408 = (t456 - t480) * t389;
t403 = -t275 * t498 - t443;
t402 = -t286 * t498 - t444;
t138 = rSges(6,1) * t196 - rSges(6,2) * t195 + rSges(6,3) * t255;
t71 = -t136 * t407 + t138 * t406 + t198 * t255 - t200 * t253;
t180 = rSges(5,1) * t256 - rSges(5,2) * t255 - rSges(5,3) * t316;
t397 = m(5) * t180 + m(6) * t138 + m(7) * t489;
t396 = -t210 * t498 + t403;
t395 = -t250 * t498 + t402;
t370 = (rSges(3,1) * t504 - rSges(3,2) * t501) * t471;
t369 = (Icges(3,1) * t504 - t450) * t471;
t368 = (-Icges(3,2) * t501 + t451) * t471;
t359 = t498 * rSges(3,3) + (rSges(3,1) * t501 + rSges(3,2) * t504) * t389;
t358 = Icges(3,5) * t498 + (Icges(3,1) * t501 + t451) * t389;
t357 = Icges(3,6) * t498 + (Icges(3,2) * t504 + t450) * t389;
t336 = rSges(3,1) * t365 + rSges(3,2) * t366;
t334 = rSges(3,1) * t363 + rSges(3,2) * t364;
t333 = Icges(3,1) * t365 + Icges(3,4) * t366;
t332 = Icges(3,1) * t363 + Icges(3,4) * t364;
t331 = Icges(3,4) * t365 + Icges(3,2) * t366;
t330 = Icges(3,4) * t363 + Icges(3,2) * t364;
t323 = -rSges(3,1) * t401 + rSges(3,2) * t374 + rSges(3,3) * t493;
t322 = -rSges(3,1) * t400 + rSges(3,2) * t372 - rSges(3,3) * t492;
t321 = -Icges(3,1) * t401 + Icges(3,4) * t374 + Icges(3,5) * t493;
t320 = -Icges(3,1) * t400 + Icges(3,4) * t372 - Icges(3,5) * t492;
t319 = -Icges(3,4) * t401 + Icges(3,2) * t374 + Icges(3,6) * t493;
t318 = -Icges(3,4) * t400 + Icges(3,2) * t372 - Icges(3,6) * t492;
t310 = Icges(4,1) * t361 - Icges(4,4) * t360 + Icges(4,5) * t498;
t309 = Icges(4,4) * t361 - Icges(4,2) * t360 + Icges(4,6) * t498;
t307 = -Icges(4,1) * t355 - Icges(4,4) * t354;
t306 = -Icges(4,4) * t355 - Icges(4,2) * t354;
t282 = rSges(4,1) * t340 + rSges(4,2) * t339 + rSges(4,3) * t493;
t281 = rSges(4,1) * t338 + rSges(4,2) * t337 - rSges(4,3) * t492;
t280 = Icges(4,1) * t340 + Icges(4,4) * t339 + Icges(4,5) * t493;
t279 = Icges(4,1) * t338 + Icges(4,4) * t337 - Icges(4,5) * t492;
t278 = Icges(4,4) * t340 + Icges(4,2) * t339 + Icges(4,6) * t493;
t277 = Icges(4,4) * t338 + Icges(4,2) * t337 - Icges(4,6) * t492;
t274 = rSges(4,1) * t317 + rSges(4,2) * t316;
t273 = rSges(4,1) * t315 + rSges(4,2) * t314;
t272 = Icges(4,1) * t317 + Icges(4,4) * t316;
t271 = Icges(4,1) * t315 + Icges(4,4) * t314;
t270 = Icges(4,4) * t317 + Icges(4,2) * t316;
t269 = Icges(4,4) * t315 + Icges(4,2) * t314;
t257 = t337 * t292;
t244 = t360 * t251;
t228 = t339 * t250;
t227 = rSges(5,1) * t300 + rSges(5,2) * t407 - rSges(5,3) * t339;
t226 = rSges(5,1) * t298 + rSges(5,2) * t406 - rSges(5,3) * t337;
t213 = -t273 * t498 + t390 * t437 - t443;
t212 = t274 * t498 + t388 * t437 + t345;
t172 = (t273 * t388 + t274 * t390) * t389 + t472;
t171 = t227 * t360 + t291 * t339;
t170 = -t226 * t360 - t291 * t337;
t159 = -t226 * t339 + t227 * t337;
t156 = -t226 * t498 + t390 * t425 + t402;
t155 = t227 * t498 + t388 * t425 + t477;
t154 = -t200 * t405 + t236 * t407;
t153 = t198 * t405 - t236 * t406;
t148 = (t226 * t388 + t227 * t390) * t389 + t428;
t143 = -t179 * t498 + t390 * t426 + t403;
t142 = t180 * t498 + t388 * t426 + t478;
t139 = -t198 * t407 + t200 * t406;
t118 = t200 * t360 + t339 * t479 + t244;
t117 = -t236 * t337 + t360 * t486 - t257;
t116 = t180 * t360 + t227 * t354 + t242 * t339 + t291 * t316;
t115 = -t179 * t360 - t226 * t354 - t242 * t337 - t291 * t314;
t114 = -t198 * t498 + t390 * t416 + t395;
t113 = t200 * t498 + t388 * t416 + t457;
t112 = (t179 * t388 + t180 * t390) * t389 + t429;
t111 = -t405 * t485 + t407 * t480;
t110 = t405 * t487 - t406 * t480;
t108 = -t198 * t339 + t337 * t484 - t228;
t103 = -t179 * t339 + t180 * t337 - t226 * t316 + t227 * t314;
t94 = (t198 * t388 + t200 * t390) * t389 + t413;
t93 = t339 * t459 + t360 * t485 + t244;
t92 = -t337 * t480 + t360 * t462 - t257;
t89 = t390 * t408 - t487 * t498 + t395;
t88 = t388 * t408 + t485 * t498 + t457;
t87 = t406 * t485 - t407 * t487;
t86 = t337 * t461 - t339 * t487 - t228;
t83 = (t388 * t487 + t390 * t485) * t389 + t413;
t82 = -t136 * t498 + t390 * t417 + t396;
t81 = t138 * t498 + t388 * t417 + t460;
t80 = -t138 * t405 + t167 * t407 + t200 * t293 - t236 * t255;
t79 = t136 * t405 - t167 * t406 - t198 * t293 + t236 * t253;
t72 = (t136 * t388 + t138 * t390) * t389 + t414;
t70 = t138 * t360 + t200 * t354 + (t167 + t245) * t339 + t479 * t316 + t482;
t69 = -t167 * t337 - t236 * t314 + (-t136 - t210) * t360 + t486 * t354 + t481;
t64 = t150 * t498 + (-t106 * t390 + t107 * t388) * t389;
t63 = t149 * t498 + (-t104 * t390 + t105 * t388) * t389;
t62 = t390 * t409 - t490 * t498 + t396;
t61 = t388 * t409 + t489 * t498 + t460;
t60 = -t106 * t337 - t107 * t339 + t150 * t360;
t59 = -t104 * t337 - t105 * t339 + t149 * t360;
t58 = -t106 * t406 - t107 * t407 - t150 * t405;
t57 = -t104 * t406 - t105 * t407 - t149 * t405;
t56 = -t136 * t339 - t198 * t316 + (t138 + t211) * t337 + t484 * t314 + t483;
t55 = t122 * t498 + (-t101 * t390 + t102 * t388) * t389;
t54 = t121 * t498 + (t100 * t388 - t390 * t99) * t389;
t53 = t120 * t498 + (t388 * t98 - t390 * t97) * t389;
t52 = t119 * t498 + (t388 * t96 - t390 * t95) * t389;
t51 = -t101 * t337 - t102 * t339 + t122 * t360;
t50 = -t100 * t339 + t121 * t360 - t337 * t99;
t49 = t120 * t360 - t337 * t97 - t339 * t98;
t48 = t119 * t360 - t337 * t95 - t339 * t96;
t47 = -t101 * t406 - t102 * t407 - t122 * t405;
t46 = -t100 * t407 - t121 * t405 - t406 * t99;
t45 = -t120 * t405 - t406 * t97 - t407 * t98;
t44 = -t119 * t405 - t406 * t95 - t407 * t96;
t39 = -t255 * t480 + t293 * t485 - t405 * t489 + t407 * t488;
t38 = t253 * t480 - t293 * t487 + t405 * t490 - t406 * t488;
t37 = (t388 * t490 + t390 * t489) * t389 + t414;
t28 = t489 * t360 + t485 * t354 + (t245 + t488) * t339 + t459 * t316 + t482;
t27 = -t488 * t337 - t480 * t314 + (-t210 - t490) * t360 + t462 * t354 + t481;
t25 = -t407 * t490 + t515;
t22 = -t490 * t339 - t487 * t316 + (t211 + t489) * t337 + t461 * t314 + t483;
t19 = [0; t473 * t514 + t430 * t513 + (m(3) * t336 + m(4) * t274 + t397) * t492 + (m(3) * t334 + m(4) * t273 - t420 + t520) * t493 + t519 * (0.2e1 * t207 + 0.2e1 * t208 + t430); (t37 * t83 + t61 * t88 + t62 * t89) * t438 + (t113 * t81 + t114 * t82 + t72 * t94) * t440 + (t112 * t148 + t142 * t155 + t143 * t156) * t442 + (((-t319 * t446 + t321 * t447 + t331 * t504 + t333 * t501) * t388 - (-t318 * t446 + t320 * t447 + t330 * t504 + t332 * t501) * t390) * t522 + (-t360 * t306 + t361 * t307 - t354 * t309 - t355 * t310 + t526 * t498) * t498 + ((-t357 * t446 + t358 * t447 + t368 * t504 + t369 * t501) * t498 + (t360 * t269 - t361 * t271 + t354 * t277 + t355 * t279 + t498 * t528) * t390 + (-t360 * t270 + t361 * t272 - t354 * t278 - t355 * t280 + t498 * t527) * t388) * t389 + t524) * t498 + ((t270 * t339 + t272 * t340 + t278 * t316 + t280 * t317 + t319 * t366 + t321 * t365 + t331 * t374 - t333 * t401 + t493 * t527) * t493 + (t306 * t339 + t307 * t340 + t309 * t316 + t310 * t317 + t357 * t366 + t358 * t365 + t368 * t374 - t369 * t401 + t493 * t526) * t498 + t525) * t493 + ((t269 * t337 + t271 * t338 + t277 * t314 + t279 * t315 + t318 * t364 + t320 * t363 + t330 * t372 - t332 * t400 + t492 * t528) * t492 + (-t306 * t337 - t307 * t338 - t309 * t314 - t310 * t315 - t357 * t364 - t358 * t363 - t368 * t372 + t369 * t400 + t492 * t526) * t498 + (-t319 * t364 - t321 * t363 - t331 * t372 + t333 * t400 - t270 * t337 - t272 * t338 - t278 * t314 - t280 * t315 - t318 * t366 - t320 * t365 - t330 * t374 + t332 * t401 - t269 * t339 - t271 * t340 - t277 * t316 - t279 * t317 + t527 * t492 + t528 * t493) * t493 - t523) * t492 + 0.2e1 * m(4) * ((-t281 * t498 + t390 * t436 - t444) * t213 + (t282 * t498 + t388 * t436 + t313) * t212 + ((t281 * t388 + t282 * t390) * t389 + t476) * t172) + 0.2e1 * m(3) * ((-t322 * t498 - t359 * t492) * (-t334 * t498 - t370 * t492) + (t323 * t498 - t359 * t493) * (t336 * t498 - t370 * t493) + (t322 * t388 + t323 * t390) * t522 * (t334 * t388 + t336 * t390)); 0; (t37 * t468 + (t388 * t62 - t390 * t61) * t467) * t511 + (t72 * t468 + (t388 * t82 - t390 * t81) * t467) * t512 + (t112 * t468 + (-t142 * t390 + t143 * t388) * t467) * t513 + (t172 * t468 + (-t212 * t390 + t213 * t388) * t467) * t514; 0; (-m(7) * t490 + t420) * t339 + t397 * t337 + (-m(5) * t226 - m(6) * t198 - m(7) * t487) * t316 + (m(5) * t227 + m(6) * t200 + m(7) * t485) * t314 + 0.2e1 * t519 * (t337 * t211 + t314 * t251 + t483); m(5) * (t103 * t148 + t112 * t159 + t115 * t156 + t116 * t155 + t142 * t171 + t143 * t170) + (t108 * t72 + t113 * t70 + t114 * t69 + t117 * t82 + t118 * t81 + t56 * t94) * m(6) + (t22 * t83 + t27 * t89 + t28 * t88 + t37 * t86 + t61 * t93 + t62 * t92) * m(7) - (t52 + t157 * t498 + (-t144 * t390 + t145 * t388) * t389 + t53) * t314 / 0.2e1 - (t158 * t498 + (-t146 * t390 + t147 * t388) * t389 + t55 + t54) * t316 / 0.2e1 - t523 * t337 / 0.2e1 - t525 * t339 / 0.2e1 + (t169 * t498 + (-t151 * t390 + t152 * t388) * t389 + t64 + t63) * t354 / 0.2e1 + t524 * t360 / 0.2e1 + t516 * t463 + t517 * t493 / 0.2e1 - t518 * t492 / 0.2e1; m(5) * t103 * t498 + t56 * t455 + t22 * t454 + ((-m(5) * t116 - m(6) * t70 - m(7) * t28) * t390 + (m(5) * t115 + m(6) * t69 + m(7) * t27) * t388) * t389; (t22 * t86 + t27 * t92 + t28 * t93) * t438 + (t108 * t56 + t117 * t69 + t118 * t70) * t440 + (t103 * t159 + t115 * t170 + t116 * t171) * t442 + (t60 + t59) * t354 + (t495 + t516) * t360 + (-t354 * t152 - t517) * t339 + (-t354 * t151 - t518) * t337 + (t146 * t337 + t147 * t339 - t158 * t360 - t50 - t51) * t316 + (t144 * t337 + t145 * t339 - t157 * t360 - t48 - t49) * t314; m(6) * t71 - t407 * t520 + t515 * t531; (t110 * t62 + t111 * t61 + t25 * t83 + t37 * t87 + t38 * t89 + t39 * t88) * m(7) + (t113 * t80 + t114 * t79 + t139 * t72 + t153 * t82 + t154 * t81 + t71 * t94) * m(6) - (t18 / 0.2e1 + t17 / 0.2e1) * t405 - (t16 / 0.2e1 + t15 / 0.2e1) * t407 - (t14 / 0.2e1 + t13 / 0.2e1) * t406 + (t64 / 0.2e1 + t63 / 0.2e1) * t293 + (t55 / 0.2e1 + t54 / 0.2e1) * t255 + (t53 / 0.2e1 + t52 / 0.2e1) * t253 + (-t388 * t465 + t390 * t466) * t389 + t521 * t463; t71 * t455 + t25 * t454 + ((-m(6) * t80 - m(7) * t39) * t390 + (m(6) * t79 + m(7) * t38) * t388) * t389; (t108 * t71 + t117 * t79 + t118 * t80 + t139 * t56 + t153 * t69 + t154 * t70) * m(6) + (t110 * t27 + t111 * t28 + t22 * t87 + t25 * t86 + t38 * t92 + t39 * t93) * m(7) + (t10 / 0.2e1 + t9 / 0.2e1) * t360 + (t57 / 0.2e1 + t58 / 0.2e1) * t354 - (t12 / 0.2e1 + t11 / 0.2e1) * t405 + t465 * t339 + t466 * t337 + (-t46 / 0.2e1 - t47 / 0.2e1) * t316 + (-t44 / 0.2e1 - t45 / 0.2e1) * t314 - (t8 / 0.2e1 + t7 / 0.2e1) * t407 - (t6 / 0.2e1 + t5 / 0.2e1) * t406 + (t59 / 0.2e1 + t60 / 0.2e1) * t293 + (t50 / 0.2e1 + t51 / 0.2e1) * t255 + (t48 / 0.2e1 + t49 / 0.2e1) * t253; (t110 * t38 + t111 * t39 + t25 * t87) * t438 + (t139 * t71 + t153 * t79 + t154 * t80) * t440 - t521 * t405 - (t3 + t4) * t407 - (t1 + t2) * t406 + (t57 + t58) * t293 + (t46 + t47) * t255 + (t44 + t45) * t253; t237 * t531; (t193 * t88 + t195 * t89 + t237 * t83 - t37 * t410 - t411 * t62 - t412 * t61) * m(7); (t237 * t468 + (-t193 * t390 + t195 * t388) * t467) * t511; (t193 * t93 + t195 * t92 - t22 * t410 + t237 * t86 - t27 * t411 - t28 * t412) * m(7); (t110 * t195 + t111 * t193 + t237 * t87 - t25 * t410 - t38 * t411 - t39 * t412) * m(7); (-t193 * t412 - t195 * t411 - t237 * t410) * t438;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
