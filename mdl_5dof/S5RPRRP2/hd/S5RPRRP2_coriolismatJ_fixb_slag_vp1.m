% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:08
% EndTime: 2020-01-03 11:45:17
% DurationCPUTime: 4.74s
% Computational Cost: add. (18685->269), mult. (12347->341), div. (0->0), fcn. (10902->8), ass. (0->180)
t260 = qJ(1) + pkin(8);
t256 = qJ(3) + t260;
t251 = sin(t256);
t252 = cos(t256);
t262 = sin(qJ(4));
t264 = cos(qJ(4));
t286 = rSges(6,2) * t264 + (rSges(6,1) + pkin(4)) * t262;
t418 = t286 * t252;
t419 = t286 * t251;
t422 = m(6) * (t251 * t419 + t252 * t418);
t120 = -t422 / 0.2e1;
t107 = t422 / 0.2e1;
t421 = Icges(5,2) + Icges(6,2);
t339 = Icges(6,4) * t262;
t340 = Icges(5,4) * t262;
t413 = t421 * t264 + t339 + t340;
t257 = Icges(6,4) * t264;
t258 = Icges(5,4) * t264;
t417 = t257 + t258 + (Icges(5,1) + Icges(6,1)) * t262;
t320 = t252 * t264;
t321 = t252 * t262;
t166 = Icges(6,4) * t320 - Icges(6,2) * t321 + Icges(6,6) * t251;
t168 = Icges(5,4) * t320 - Icges(5,2) * t321 + Icges(5,6) * t251;
t416 = t166 + t168;
t220 = Icges(6,4) * t321;
t170 = Icges(6,1) * t320 + Icges(6,5) * t251 - t220;
t221 = Icges(5,4) * t321;
t172 = Icges(5,1) * t320 + Icges(5,5) * t251 - t221;
t415 = t170 + t172;
t233 = Icges(6,1) * t264 - t339;
t235 = Icges(5,1) * t264 - t340;
t412 = t233 + t235;
t411 = (Icges(5,6) + Icges(6,6)) * t264 + (Icges(5,5) + Icges(6,5)) * t262;
t410 = t421 * t320 + t220 + t221 - t415;
t169 = -Icges(6,5) * t252 + t233 * t251;
t171 = -Icges(5,5) * t252 + t235 * t251;
t409 = -t413 * t251 + t169 + t171;
t408 = t417 * t252 + t416;
t229 = -Icges(6,2) * t262 + t257;
t165 = -Icges(6,6) * t252 + t229 * t251;
t231 = -Icges(5,2) * t262 + t258;
t167 = -Icges(5,6) * t252 + t231 * t251;
t407 = -t417 * t251 - t165 - t167;
t406 = t229 + t231 + t417;
t390 = rSges(6,3) + qJ(5) + pkin(7);
t405 = t390 * t252;
t325 = t251 * t262;
t281 = -rSges(6,2) * t325 - t405;
t351 = pkin(4) * t264;
t253 = pkin(3) + t351;
t344 = rSges(6,1) * t264;
t284 = t253 + t344;
t143 = t251 * t284 + t281;
t283 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t260);
t133 = t143 + t283;
t404 = -t133 + t143;
t272 = rSges(6,1) * t320 - rSges(6,2) * t321 + t252 * t253;
t144 = t390 * t251 + t272;
t65 = -t143 * t419 - t144 * t418;
t246 = t252 * pkin(7);
t298 = -rSges(5,2) * t325 - t252 * rSges(5,3);
t345 = rSges(5,1) * t264;
t145 = -t246 + (pkin(3) + t345) * t251 + t298;
t270 = rSges(5,1) * t320 - rSges(5,2) * t321 + rSges(5,3) * t251;
t293 = t252 * pkin(3) + t251 * pkin(7);
t146 = t270 + t293;
t237 = rSges(5,1) * t262 + rSges(5,2) * t264;
t193 = t237 * t251;
t194 = t237 * t252;
t70 = -t145 * t193 - t146 * t194;
t401 = -m(5) * t70 - m(6) * t65;
t377 = m(5) / 0.2e1;
t376 = m(6) / 0.2e1;
t240 = -rSges(5,2) * t262 + t345;
t397 = t240 * t377;
t396 = t411 * t251;
t395 = t411 * t252;
t393 = -t409 * t262 + t407 * t264;
t392 = t410 * t262 - t408 * t264;
t391 = (-t412 + t413) * t264 + t406 * t262;
t290 = qJD(1) + qJD(3);
t291 = pkin(2) * cos(t260) + cos(qJ(1)) * pkin(1);
t356 = m(4) * (-t291 * (rSges(4,1) * t251 + rSges(4,2) * t252) + t283 * (t252 * rSges(4,1) - rSges(4,2) * t251));
t248 = t251 ^ 2;
t249 = t252 ^ 2;
t292 = t248 + t249;
t134 = t144 + t291;
t141 = t145 + t283;
t142 = t146 + t291;
t350 = ((-t134 + t144) * t418 + t404 * t419) * t376 + ((-t142 + t146) * t252 + (-t141 + t145) * t251) * t237 * t377;
t63 = -t133 * t419 - t134 * t418;
t68 = -t141 * t193 - t142 * t194;
t387 = (t65 + t63) * t376 + (t70 + t68) * t377;
t269 = t406 * t264 / 0.2e1 + (-t413 / 0.2e1 + t412 / 0.2e1) * t262;
t380 = 4 * qJD(1);
t378 = 2 * qJD(4);
t62 = t141 * t146 - t145 * t142;
t373 = m(5) * t62;
t371 = m(5) * t68;
t123 = t134 * t251;
t137 = t144 * t251;
t366 = m(6) * (t404 * t252 + t123 - t137);
t365 = m(6) * (t123 + t137 + (-t133 - t143) * t252);
t57 = t133 * t144 - t143 * t134;
t364 = m(6) * t57;
t362 = m(6) * t63;
t360 = -t251 / 0.2e1;
t358 = -t252 / 0.2e1;
t357 = t252 / 0.2e1;
t346 = m(6) * qJD(5);
t334 = t165 * t262;
t333 = t166 * t262;
t332 = t167 * t262;
t331 = t168 * t262;
t330 = t169 * t264;
t329 = t171 * t264;
t226 = Icges(6,5) * t264 - Icges(6,6) * t262;
t327 = t226 * t251;
t227 = Icges(5,5) * t264 - Icges(5,6) * t262;
t326 = t227 * t251;
t324 = t251 * t264;
t149 = t169 * t324;
t150 = t170 * t324;
t153 = t165 * t321;
t161 = -Icges(6,3) * t252 + t327;
t275 = t170 * t264 - t333;
t84 = -t161 * t251 - t169 * t320 + t153;
t162 = Icges(6,5) * t320 - Icges(6,6) * t321 + Icges(6,3) * t251;
t85 = t162 * t251 + t252 * t275;
t10 = (t84 + t150 - t153 + (t161 - t333) * t251) * t251 + (-t149 - t85 + (t161 + t275) * t252 + (t330 + t334) * t251) * t252;
t151 = t171 * t324;
t152 = t172 * t324;
t154 = t167 * t321;
t163 = -Icges(5,3) * t252 + t326;
t273 = t172 * t264 - t331;
t86 = -t163 * t251 - t171 * t320 + t154;
t164 = Icges(5,5) * t320 - Icges(5,6) * t321 + Icges(5,3) * t251;
t87 = t164 * t251 + t252 * t273;
t11 = (t86 + t152 - t154 + (t163 - t331) * t251) * t251 + (-t151 - t87 + (t163 + t273) * t252 + (t329 + t332) * t251) * t252;
t80 = -t161 * t252 - t165 * t325 + t149;
t81 = t162 * t252 + t166 * t325 - t150;
t12 = (t153 - t81 + (t162 - t330) * t252) * t252 + (-t149 + t80 + (t162 + t334) * t251) * t251;
t82 = -t163 * t252 - t167 * t325 + t151;
t83 = t164 * t252 + t168 * t325 - t152;
t13 = (t154 - t83 + (t164 - t329) * t252) * t252 + (-t151 + t82 + (t164 + t332) * t251) * t251;
t48 = -t251 * t81 - t252 * t80;
t49 = -t251 * t83 - t252 * t82;
t50 = -t251 * t85 - t252 * t84;
t51 = -t251 * t87 - t252 * t86;
t2 = (-t13 / 0.2e1 - t51 / 0.2e1 - t12 / 0.2e1 - t50 / 0.2e1) * t252 + (t49 / 0.2e1 - t11 / 0.2e1 + t48 / 0.2e1 - t10 / 0.2e1) * t251;
t60 = 0.2e1 * t120;
t289 = t2 * qJD(4) + t60 * qJD(5);
t69 = -t133 * t252 + t123;
t288 = m(6) * t69 * qJD(1);
t79 = -t143 * t252 + t137;
t287 = t79 * m(6) * qJD(3);
t285 = -rSges(6,2) * t262 + t344 + t351;
t268 = t269 + t387;
t267 = -t269 + (t415 * t262 + t416 * t264) * (t357 + t358);
t59 = t120 + t107;
t266 = t59 * qJD(5) + ((t10 + t11) * t251 / 0.2e1 + (t409 * t264 + t407 * t262 + (-t226 - t227) * t252 - t391 * t251) * t358 + (t12 + t13 + t50 + t51) * t357 + (t391 * t252 + t408 * t262 + t410 * t264 - t326 - t327 + t48 + t49) * t360) * qJD(4);
t178 = t285 * t252;
t176 = t285 * t251;
t138 = -t193 * t251 - t194 * t252;
t110 = t286 * t292;
t58 = 0.2e1 * t107;
t55 = t59 * qJD(4);
t53 = t58 * qJD(4);
t35 = t365 / 0.2e1;
t34 = t366 / 0.2e1;
t23 = t269 - t401;
t20 = t269 + t362 + t371;
t17 = t356 + t364 + t373;
t16 = t35 - t366 / 0.2e1;
t15 = t35 + t34;
t14 = t34 - t365 / 0.2e1;
t5 = t268 + t350;
t4 = t268 - t350;
t3 = t267 + t350 - t387;
t1 = [t17 * qJD(3) + t20 * qJD(4) + t346 * t69, 0, t17 * qJD(1) + t5 * qJD(4) + t15 * qJD(5) + 0.2e1 * (t356 / 0.2e1 + t62 * t377 + t57 * t376) * qJD(3), t20 * qJD(1) + t5 * qJD(3) + ((t141 * t252 - t142 * t251) * t397 + (t133 * t178 - t134 * t176) * t376) * t378 + t266, t15 * qJD(3) + t288 + t55; 0, 0, 0, (-t110 * t376 + t138 * t377) * t378, 0; t4 * qJD(4) + t16 * qJD(5) + (-t356 / 0.4e1 - t373 / 0.4e1 - t364 / 0.4e1) * t380, 0, t23 * qJD(4) + t346 * t79, t4 * qJD(1) + t23 * qJD(3) + ((t143 * t178 - t144 * t176) * t376 + (t145 * t252 - t146 * t251) * t397) * t378 + t266, t16 * qJD(1) + t287 + t55; t267 * qJD(1) + t3 * qJD(3) + (-t371 / 0.4e1 - t362 / 0.4e1) * t380 + t289, 0, t3 * qJD(1) + t289 + (t267 + t401) * qJD(3), (m(5) * (t237 * t240 * t292 + (t252 * t270 + t251 * (rSges(5,1) * t324 + t298)) * t138) + m(6) * (-t110 * ((t272 - t293) * t252 + (t405 + t246 + t281 + (-pkin(3) + t284) * t251) * t251) + t419 * t176 + t418 * t178) + (t395 * t248 + (t393 * t252 + (-t392 - t396) * t251) * t252) * t360 + (-t396 * t249 + (t392 * t251 + (-t393 + t395) * t252) * t251) * t358) * qJD(4) + t290 * t2, t290 * t60; t14 * qJD(3) - t288 + t53, 0, t14 * qJD(1) - t287 + t53, m(6) * (t176 * t252 - t178 * t251) * qJD(4) + t290 * t58, 0;];
Cq = t1;
