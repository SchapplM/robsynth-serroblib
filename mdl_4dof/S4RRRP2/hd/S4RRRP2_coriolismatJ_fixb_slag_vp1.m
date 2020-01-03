% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:13:00
% DurationCPUTime: 4.13s
% Computational Cost: add. (11523->257), mult. (11578->335), div. (0->0), fcn. (10284->6), ass. (0->168)
t255 = qJ(1) + qJ(2);
t251 = sin(t255);
t252 = cos(t255);
t256 = sin(qJ(3));
t258 = cos(qJ(3));
t400 = rSges(5,2) * t258 + (rSges(5,1) + pkin(3)) * t256;
t407 = t400 * t252;
t408 = t400 * t251;
t411 = m(5) * (t251 * t408 + t252 * t407);
t107 = -t411 / 0.2e1;
t98 = t411 / 0.2e1;
t410 = Icges(4,2) + Icges(5,2);
t330 = Icges(5,4) * t256;
t331 = Icges(4,4) * t256;
t403 = t258 * t410 + t330 + t331;
t253 = Icges(5,4) * t258;
t254 = Icges(4,4) * t258;
t406 = t253 + t254 + (Icges(4,1) + Icges(5,1)) * t256;
t314 = t251 * t258;
t315 = t251 * t256;
t161 = Icges(5,4) * t314 - Icges(5,2) * t315 - Icges(5,6) * t252;
t163 = Icges(4,4) * t314 - Icges(4,2) * t315 - Icges(4,6) * t252;
t405 = t161 + t163;
t230 = Icges(5,4) * t315;
t165 = Icges(5,1) * t314 - Icges(5,5) * t252 - t230;
t231 = Icges(4,4) * t315;
t167 = Icges(4,1) * t314 - Icges(4,5) * t252 - t231;
t404 = t165 + t167;
t217 = Icges(5,1) * t258 - t330;
t219 = Icges(4,1) * t258 - t331;
t402 = t217 + t219;
t213 = -Icges(5,2) * t256 + t253;
t215 = -Icges(4,2) * t256 + t254;
t401 = -t215 - t213 - t406;
t397 = (-Icges(4,6) - Icges(5,6)) * t258 + (-Icges(4,5) - Icges(5,5)) * t256;
t166 = Icges(5,5) * t251 + t217 * t252;
t168 = Icges(4,5) * t251 + t219 * t252;
t396 = -t252 * t403 + t166 + t168;
t395 = -t314 * t410 - t230 - t231 + t404;
t162 = Icges(5,6) * t251 + t213 * t252;
t164 = Icges(4,6) * t251 + t215 * t252;
t394 = -t252 * t406 - t162 - t164;
t393 = t251 * t406 + t405;
t344 = pkin(3) * t258;
t248 = pkin(2) + t344;
t342 = -qJ(4) - pkin(6);
t133 = -rSges(5,1) * t314 + rSges(5,2) * t315 - t251 * t248 + (rSges(5,3) - t342) * t252;
t335 = rSges(5,1) * t258;
t283 = t248 + t335;
t313 = t252 * t256;
t289 = -rSges(5,2) * t313 - t251 * t342;
t134 = t251 * rSges(5,3) + t283 * t252 + t289;
t58 = t133 * t408 - t134 * t407;
t247 = t252 * pkin(6);
t336 = rSges(4,1) * t258;
t286 = pkin(2) + t336;
t290 = rSges(4,2) * t315 + rSges(4,3) * t252;
t137 = -t251 * t286 + t247 + t290;
t235 = rSges(4,2) * t313;
t138 = -t235 + t286 * t252 + (rSges(4,3) + pkin(6)) * t251;
t221 = rSges(4,1) * t256 + rSges(4,2) * t258;
t187 = t221 * t251;
t188 = t221 * t252;
t61 = t137 * t187 - t138 * t188;
t392 = -m(4) * t61 - m(5) * t58;
t375 = m(4) / 0.2e1;
t374 = m(5) / 0.2e1;
t71 = t133 * t252 + t134 * t251;
t388 = t71 * m(5) * qJD(2);
t387 = (t402 - t403) * t258 + t401 * t256;
t345 = cos(qJ(1)) * pkin(1);
t346 = sin(qJ(1)) * pkin(1);
t351 = m(3) * (t345 * (-rSges(3,1) * t251 - rSges(3,2) * t252) + (rSges(3,1) * t252 - rSges(3,2) * t251) * t346);
t131 = t133 - t346;
t132 = t134 + t345;
t67 = t131 * t252 + t132 * t251;
t385 = t67 * m(5) * qJD(1);
t384 = t397 * t251;
t383 = t397 * t252;
t249 = t251 ^ 2;
t250 = t252 ^ 2;
t288 = t249 + t250;
t382 = -t256 * t396 + t258 * t394;
t381 = t256 * t395 + t258 * t393;
t287 = qJD(1) + qJD(2);
t135 = t137 - t346;
t136 = t138 + t345;
t341 = ((-t132 + t134) * t407 + (t131 - t133) * t408) * t374 + ((-t136 + t138) * t252 + (t135 - t137) * t251) * t221 * t375;
t56 = t131 * t408 - t132 * t407;
t60 = t135 * t187 - t136 * t188;
t380 = (t58 + t56) * t374 + (t61 + t60) * t375;
t269 = -t401 * t258 / 0.2e1 + (-t403 / 0.2e1 + t402 / 0.2e1) * t256;
t378 = 0.4e1 * qJD(1);
t376 = 2 * qJD(3);
t55 = -t135 * t138 + t136 * t137;
t371 = m(4) * t55;
t369 = m(4) * t60;
t364 = m(5) * (t67 - t71);
t363 = m(5) * (t71 + t67);
t53 = -t131 * t134 + t132 * t133;
t362 = m(5) * t53;
t360 = m(5) * t56;
t357 = -t251 / 0.2e1;
t356 = t251 / 0.2e1;
t355 = -t252 / 0.2e1;
t157 = Icges(5,5) * t314 - Icges(5,6) * t315 - Icges(5,3) * t252;
t280 = t162 * t256 - t157;
t139 = t166 * t314;
t210 = Icges(5,5) * t258 - Icges(5,6) * t256;
t319 = t210 * t252;
t158 = Icges(5,3) * t251 + t319;
t282 = t158 * t252 - t139;
t312 = t252 * t258;
t301 = t158 * t251 + t166 * t312;
t302 = -t251 * t157 - t165 * t312;
t76 = -t161 * t313 - t302;
t77 = -t162 * t313 + t301;
t10 = (t252 * t280 - t301 + t77) * t252 + (t251 * t280 + t282 + t76) * t251;
t159 = Icges(4,5) * t314 - Icges(4,6) * t315 - Icges(4,3) * t252;
t279 = t164 * t256 - t159;
t140 = t168 * t314;
t211 = Icges(4,5) * t258 - Icges(4,6) * t256;
t318 = t211 * t252;
t160 = Icges(4,3) * t251 + t318;
t281 = t160 * t252 - t140;
t299 = t160 * t251 + t168 * t312;
t300 = -t251 * t159 - t167 * t312;
t78 = -t163 * t313 - t300;
t79 = -t164 * t313 + t299;
t11 = (t252 * t279 - t299 + t79) * t252 + (t279 * t251 + t281 + t78) * t251;
t322 = t161 * t256;
t73 = -t162 * t315 - t282;
t12 = (t73 - t139 + (t158 + t322) * t252 + t302) * t252 + t301 * t251;
t321 = t163 * t256;
t75 = -t164 * t315 - t281;
t13 = (t75 - t140 + (t160 + t321) * t252 + t300) * t252 + t299 * t251;
t40 = t251 * t73 - t252 * (-(-t165 * t258 + t322) * t251 - t157 * t252);
t41 = t251 * t75 - t252 * (-(-t167 * t258 + t321) * t251 - t159 * t252);
t42 = t251 * t77 - t252 * t76;
t43 = t251 * t79 - t252 * t78;
t2 = (-t12 / 0.2e1 + t43 / 0.2e1 - t13 / 0.2e1 + t42 / 0.2e1) * t252 + (t40 / 0.2e1 + t11 / 0.2e1 + t41 / 0.2e1 + t10 / 0.2e1) * t251;
t52 = 0.2e1 * t107;
t352 = qJD(3) * t2 + qJD(4) * t52;
t337 = m(5) * qJD(4);
t284 = rSges(5,2) * t256 - t335 - t344;
t270 = (-t187 * t252 + t188 * t251) * t221;
t262 = t269 + t380;
t261 = -t269 + (t256 * t404 + t258 * t405) * (t356 + t357);
t51 = t107 + t98;
t260 = t51 * qJD(4) + ((t12 + t13) * t252 / 0.2e1 + (t10 + t11 + t40 + t41) * t357 + (t396 * t258 + t394 * t256 + t387 * t252 + (t210 + t211) * t251) * t356 + (t387 * t251 - t256 * t393 + t258 * t395 - t318 - t319 + t42 + t43) * t355) * qJD(3);
t224 = -rSges(4,2) * t256 + t336;
t174 = t284 * t252;
t172 = t284 * t251;
t50 = 0.2e1 * t98;
t48 = t51 * qJD(3);
t46 = t50 * qJD(3);
t39 = t363 / 0.2e1;
t38 = t364 / 0.2e1;
t21 = t269 - t392;
t20 = t269 + t360 + t369;
t17 = t351 + t362 + t371;
t16 = t39 - t364 / 0.2e1;
t15 = t39 + t38;
t14 = t38 - t363 / 0.2e1;
t5 = t262 + t341;
t4 = t262 - t341;
t3 = t261 + t341 - t380;
t1 = [qJD(2) * t17 + qJD(3) * t20 + t337 * t67, t17 * qJD(1) + t5 * qJD(3) + t15 * qJD(4) + 0.2e1 * (t351 / 0.2e1 + t55 * t375 + t53 * t374) * qJD(2), t20 * qJD(1) + t5 * qJD(2) + (((-t135 * t252 - t136 * t251) * t224 + t270) * t375 + (t131 * t174 + t132 * t172) * t374) * t376 + t260, qJD(2) * t15 + t385 + t48; t4 * qJD(3) + t16 * qJD(4) + (-t351 / 0.4e1 - t371 / 0.4e1 - t362 / 0.4e1) * t378, qJD(3) * t21 + t337 * t71, t4 * qJD(1) + t21 * qJD(2) + (((-t137 * t252 - t138 * t251) * t224 + t270) * t375 + (t133 * t174 + t134 * t172) * t374) * t376 + t260, qJD(1) * t16 + t388 + t48; t261 * qJD(1) + t3 * qJD(2) + (-t369 / 0.4e1 - t360 / 0.4e1) * t378 + t352, t3 * qJD(1) + t352 + (t261 + t392) * qJD(2), (m(4) * ((t251 * (rSges(4,1) * t314 - t290) + t252 * (rSges(4,1) * t312 + t251 * rSges(4,3) - t235)) * (-t187 * t251 - t188 * t252) + t288 * t224 * t221) + m(5) * (-t408 * t172 - t407 * t174 - ((-pkin(2) * t251 - t133 + t247) * t251 + ((-pkin(2) + t283) * t252 + (rSges(5,3) - pkin(6)) * t251 + t289) * t252) * t400 * t288) + (t383 * t249 + (t381 * t252 + (t382 - t384) * t251) * t252) * t356 + (t384 * t250 + (t382 * t251 + (t381 - t383) * t252) * t251) * t355) * qJD(3) + t287 * t2, t287 * t52; t14 * qJD(2) - t385 + t46, t14 * qJD(1) - t388 + t46, m(5) * (-t172 * t252 + t174 * t251) * qJD(3) + t287 * t50, 0;];
Cq = t1;
