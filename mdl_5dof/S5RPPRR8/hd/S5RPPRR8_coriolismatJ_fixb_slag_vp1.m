% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR8
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR8_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR8_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:10
% DurationCPUTime: 3.52s
% Computational Cost: add. (14626->256), mult. (15475->352), div. (0->0), fcn. (17970->8), ass. (0->182)
t259 = pkin(8) + qJ(4);
t244 = sin(t259);
t245 = cos(t259);
t319 = sin(qJ(1));
t320 = cos(qJ(1));
t175 = t320 * t244 - t319 * t245;
t174 = -t319 * t244 - t320 * t245;
t207 = sin(qJ(5));
t208 = cos(qJ(5));
t294 = Icges(6,4) * t208;
t230 = -Icges(6,2) * t207 + t294;
t130 = Icges(6,6) * t174 + t230 * t175;
t295 = Icges(6,4) * t207;
t232 = Icges(6,1) * t208 - t295;
t134 = Icges(6,5) * t174 + t232 * t175;
t377 = t208 * t130 + t207 * t134;
t299 = t377 * t175;
t129 = -Icges(6,6) * t175 + t230 * t174;
t133 = -Icges(6,5) * t175 + t232 * t174;
t378 = t208 * t129 + t207 * t133;
t393 = t378 * t174;
t394 = t299 / 0.2e1 - t393 / 0.2e1;
t392 = -t175 / 0.4e1;
t193 = -t207 * rSges(6,1) - t208 * rSges(6,2);
t152 = t175 * t193;
t153 = t193 * t174;
t388 = m(6) * (t319 * t152 + t153 * t320);
t257 = -t388 / 0.2e1;
t95 = t388 / 0.2e1;
t343 = m(6) / 0.2e1;
t346 = m(5) / 0.2e1;
t157 = t174 * rSges(5,1) + t175 * rSges(5,2);
t234 = -t175 * rSges(5,1) + t174 * rSges(5,2);
t355 = t157 * t319 + t320 * t234;
t301 = t207 * rSges(6,2);
t303 = t174 * rSges(6,3);
t246 = -t175 * t301 + t303;
t300 = t208 * rSges(6,1);
t247 = pkin(4) + t300;
t93 = -t174 * pkin(7) - t247 * t175 - t246;
t89 = t320 * t93;
t373 = t175 * rSges(6,3) + t174 * t301;
t91 = -t175 * pkin(7) + t247 * t174 - t373;
t363 = t91 * t319 + t89;
t297 = t363 * t343 + t355 * t346;
t386 = t95 + t257;
t390 = qJD(5) * t386;
t233 = t300 - t301;
t92 = -(rSges(6,3) + pkin(7)) * t174 - (pkin(4) + t233) * t175;
t305 = t92 - t93;
t364 = m(6) * t305;
t330 = t91 * t364;
t389 = t330 * qJD(4);
t273 = t386 * qJD(2);
t304 = m(6) * qJD(5);
t137 = -t175 * t300 - t246;
t56 = (t233 * t175 + t137 + t303) * t174;
t387 = t56 * t304;
t279 = t129 * t207;
t280 = t130 * t207;
t276 = t133 * t208;
t277 = t134 * t208;
t228 = Icges(6,5) * t208 - Icges(6,6) * t207;
t127 = -Icges(6,3) * t174 - t228 * t175;
t385 = t175 * t127;
t128 = Icges(6,3) * t175 - t228 * t174;
t384 = t175 * t128;
t381 = t56 * qJD(3);
t296 = cos(pkin(8));
t204 = t296 * pkin(3) + pkin(2);
t237 = -t319 * pkin(1) + t320 * qJ(2);
t206 = sin(pkin(8));
t251 = t320 * t206;
t213 = pkin(3) * t251 - t319 * t204 + t237;
t77 = t213 - t93;
t218 = t320 * pkin(1) + t319 * qJ(2);
t250 = t319 * t206;
t209 = pkin(3) * t250 + t320 * t204 + t218;
t78 = t209 - t91;
t380 = -t77 * t91 + t78 * t93;
t272 = t174 * t127 - t175 * t277;
t358 = -t128 - t280;
t379 = -t358 * t175 + t272;
t344 = -m(6) / 0.2e1;
t376 = -t207 / 0.2e1;
t365 = -t208 / 0.2e1;
t227 = Icges(6,5) * t207 + Icges(6,6) * t208;
t144 = t227 * t174;
t375 = t227 * t175;
t374 = 0.2e1 * t343 * qJD(1);
t309 = t363 * t344 - m(5) * t355 / 0.2e1;
t310 = (-t92 * t320 + t89) * t343;
t372 = t309 - t310;
t370 = 0.2e1 * t257;
t124 = -t157 + t209;
t210 = t213 - t234;
t369 = t124 * t234 - t210 * t157;
t291 = Icges(6,2) * t208;
t229 = t291 + t295;
t220 = t207 * t229;
t231 = Icges(6,1) * t207 + t294;
t274 = t208 * t231;
t357 = t220 - t274;
t101 = t357 * t174 + t375;
t368 = -t174 / 0.2e1;
t367 = t175 / 0.2e1;
t366 = -t228 / 0.2e1;
t313 = m(6) * t193;
t361 = t175 * t174;
t359 = -t127 + t279;
t102 = t357 * t175 - t144;
t282 = t175 * t102;
t283 = t174 * t101;
t356 = t282 / 0.4e1 + t283 / 0.4e1;
t354 = t393 / 0.4e1 + t152 * t364 / 0.2e1;
t353 = t299 / 0.4e1 - ((-t77 - t93) * t175 + (t78 + t91) * t174) * t313 / 0.2e1;
t240 = t220 / 0.2e1 + t232 * t376 + t230 * t365 - t274 / 0.2e1;
t352 = t174 ^ 2;
t351 = t175 ^ 2;
t349 = 0.4e1 * qJD(1);
t348 = 0.2e1 * qJD(4);
t36 = -t175 * t280 - t272;
t271 = -t174 * t128 + t175 * t276;
t37 = -t175 * t279 + t271;
t20 = -t36 * t174 + t37 * t175;
t342 = -t20 / 0.2e1;
t270 = t174 * t277 + t385;
t38 = -t174 * t280 + t270;
t269 = t174 * t276 + t384;
t39 = -t174 * t279 + t269;
t21 = -t38 * t174 + t39 * t175;
t341 = -t21 / 0.2e1;
t185 = -t320 * t296 - t250;
t186 = -t319 * t296 + t251;
t340 = m(4) * ((-t185 * rSges(4,1) - t186 * rSges(4,2) + t218) * t319 + (t186 * rSges(4,1) - t185 * rSges(4,2) + t237) * t320);
t339 = m(5) * t369;
t336 = m(5) * (t124 * t319 + t210 * t320);
t335 = m(6) * ((t78 - t91) * t153 + (-t77 + t93) * t152);
t332 = m(6) * t380;
t329 = m(6) * (t152 * t77 - t153 * t78);
t328 = m(6) * (-t152 * t93 + t153 * t91);
t75 = t77 * t320;
t327 = m(6) * (t78 * t319 + t75);
t318 = m(3) * ((t320 * rSges(3,3) + t237) * t320 + (t319 * rSges(3,3) + t218) * t319);
t308 = t56 * qJD(1) * t344;
t79 = t213 - t92;
t307 = t77 - t79;
t285 = t152 * t193;
t284 = t153 * t193;
t261 = qJD(5) * t174;
t260 = qJD(5) * t175;
t67 = t393 / 0.2e1;
t212 = (-t232 + t291) * t208 + (t230 + t231 + t294) * t207;
t239 = t342 + (-t231 * t174 - t129) * t207 / 0.2e1 + (t229 * t174 - t133) * t365 + t212 * t368 + t175 * t366;
t238 = t341 + (-t231 * t175 - t130) * t376 + (-t229 * t175 + t134) * t365 + t174 * t366 + t212 * t367;
t225 = t240 + t394;
t224 = -t299 / 0.2e1 + t67 + t240;
t222 = -t282 / 0.4e1 - t353 + t356 - (-t378 + t101) * t174 / 0.4e1;
t221 = -t283 / 0.4e1 - t354 + t356 + (-t377 + t102) * t392;
t80 = -t175 * t152 - t174 * t153;
t58 = t175 * t137 + t174 * (-t174 * t300 + t373);
t41 = 0.2e1 * t95;
t32 = -t240 + t328;
t31 = -t240 + t329;
t28 = t58 * t56;
t24 = (-t129 * t175 - t130 * t174) * t207 + t270 + t271;
t19 = t318 + t327 + t336 + t340;
t14 = t335 / 0.2e1;
t13 = t332 + t339;
t12 = t297 - t372;
t11 = t309 + t310 - t297;
t10 = t297 + t372;
t9 = (t36 + (-t277 + t280) * t175 + t269) * t175 + (-t175 * t359 - t24 + t37) * t174;
t8 = (t24 - t38) * t175 + (-t39 + (t276 - t279) * t174 + t379) * t174;
t7 = (t38 + (t359 - t276) * t175) * t175 + (t359 * t174 - t269 + t384 + t39) * t174;
t6 = (-t36 - t379) * t175 + (-t37 + (t358 + t277) * t174 + t385) * t174;
t5 = -t335 / 0.2e1 + t221 + t222 - t240;
t4 = t377 * t392 + t14 + t222 + t225 + t354;
t3 = t14 + t221 - t393 / 0.4e1 + t224 + t353;
t2 = (t8 / 0.2e1 + t342) * t175 + (t341 - t6 / 0.2e1) * t174;
t1 = m(6) * t28 + (t7 / 0.2e1 + t20 / 0.2e1) * t175 + (t21 / 0.2e1 - t9 / 0.2e1) * t174;
t15 = [-m(6) * t307 * t78 * t349 / 0.4e1 + t19 * qJD(2) + t13 * qJD(4) + t31 * qJD(5), qJD(1) * t19 + qJD(4) * t10 + t390, t387 / 0.2e1, t13 * qJD(1) + t10 * qJD(2) + t4 * qJD(5) - t389 + (-t380 * t343 - t369 * t346) * t348, t31 * qJD(1) + t273 + t4 * qJD(4) + (-t7 / 0.2e1 + t239) * t260 + (t9 / 0.2e1 + t238) * t261 + (t381 / 0.2e1 + ((t233 * t78 + t284) * t175 + (t233 * t77 - t285) * t174 - t28) * qJD(5)) * m(6); t11 * qJD(4) + t41 * qJD(5) + (-t327 / 0.4e1 - t336 / 0.4e1 - t318 / 0.4e1 - t340 / 0.4e1) * t349 + (-t320 * t79 + t75) * t374, 0, 0, t11 * qJD(1) + t370 * qJD(5) + t297 * t348, t41 * qJD(1) + t370 * qJD(4) - (-t319 * t174 + t320 * t175) * t233 * t304; -t387 / 0.2e1, 0, 0, 0, -t80 * t304 + t308; t12 * qJD(2) + t389 + t3 * qJD(5) + (-t332 / 0.4e1 - t339 / 0.4e1) * t349 + (t305 * t78 - t307 * t91) * t374, qJD(1) * t12 - t390, 0, qJD(1) * t330 + t32 * qJD(5), t3 * qJD(1) - t273 + t32 * qJD(4) + (-t8 / 0.2e1 - t239) * t260 + (t6 / 0.2e1 - t238) * t261 + ((t233 * t91 - t284) * t175 + (t233 * t93 + t285) * t174) * t304; t370 * qJD(2) + t5 * qJD(4) + t1 * qJD(5) + t344 * t381 + (t225 + t67 - t329 + (-t377 / 0.2e1 + t307 * t313) * t175) * qJD(1), qJD(1) * t370 + qJD(4) * t386, t308, t5 * qJD(1) + t2 * qJD(5) + t273 + (t224 - t328 + t394) * qJD(4), t1 * qJD(1) + t2 * qJD(4) + (m(6) * (t58 * t80 - (t351 + t352) * t233 * t193) + (t351 * t144 - t361 * t375) * t367 + (-t144 * t361 + t352 * t375) * t368) * qJD(5);];
Cq = t15;
