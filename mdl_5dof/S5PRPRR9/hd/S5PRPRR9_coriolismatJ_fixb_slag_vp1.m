% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR9_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR9_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:38
% EndTime: 2019-12-31 17:39:42
% DurationCPUTime: 3.19s
% Computational Cost: add. (14864->246), mult. (15051->338), div. (0->0), fcn. (17542->6), ass. (0->180)
t242 = pkin(8) + qJ(2);
t230 = sin(t242);
t231 = cos(t242);
t301 = sin(qJ(4));
t302 = cos(qJ(4));
t170 = -t230 * t302 + t231 * t301;
t169 = -t230 * t301 - t231 * t302;
t194 = sin(qJ(5));
t195 = cos(qJ(5));
t260 = t195 * Icges(6,4);
t215 = -Icges(6,2) * t194 + t260;
t127 = Icges(6,6) * t169 + t215 * t170;
t279 = Icges(6,4) * t194;
t217 = Icges(6,1) * t195 - t279;
t131 = Icges(6,5) * t169 + t217 * t170;
t356 = t127 * t195 + t194 * t131;
t281 = t356 * t170;
t126 = -Icges(6,6) * t170 + t215 * t169;
t130 = -Icges(6,5) * t170 + t217 * t169;
t357 = t126 * t195 + t194 * t130;
t375 = t357 * t169;
t376 = t281 / 0.2e1 - t375 / 0.2e1;
t374 = -t170 / 0.4e1;
t188 = -t194 * rSges(6,1) - rSges(6,2) * t195;
t148 = t170 * t188;
t149 = t188 * t169;
t370 = m(6) * (t230 * t148 + t149 * t231);
t240 = -t370 / 0.2e1;
t91 = t370 / 0.2e1;
t323 = m(6) / 0.2e1;
t325 = m(5) / 0.2e1;
t150 = t169 * rSges(5,1) + t170 * rSges(5,2);
t219 = -t170 * rSges(5,1) + t169 * rSges(5,2);
t335 = t150 * t230 + t231 * t219;
t166 = t169 * pkin(7);
t282 = t194 * rSges(6,2);
t163 = t170 * t282;
t165 = t169 * rSges(6,3);
t232 = -t163 + t165;
t284 = rSges(6,1) * t195;
t234 = pkin(4) + t284;
t95 = -t234 * t170 - t166 - t232;
t89 = t231 * t95;
t352 = t170 * rSges(6,3) + t169 * t282;
t93 = -t170 * pkin(7) + t234 * t169 - t352;
t342 = t93 * t230 + t89;
t290 = t342 * t323 + t335 * t325;
t368 = t91 + t240;
t372 = qJD(5) * t368;
t264 = t170 * t195;
t162 = rSges(6,1) * t264;
t245 = t162 + t165;
t94 = (-pkin(4) + t282) * t170 - t166 - t245;
t286 = t94 - t95;
t343 = m(6) * t286;
t311 = t93 * t343;
t371 = t311 * qJD(4);
t257 = t368 * qJD(3);
t285 = m(6) * qJD(5);
t134 = -t232 - t162;
t56 = (t134 - t163 + t245) * t169;
t369 = t56 * t285;
t262 = t126 * t194;
t263 = t127 * t194;
t367 = t130 * t195;
t366 = t131 * t195;
t213 = Icges(6,5) * t195 - Icges(6,6) * t194;
t124 = -Icges(6,3) * t169 - t213 * t170;
t365 = t170 * t124;
t125 = Icges(6,3) * t170 - t213 * t169;
t364 = t170 * t125;
t361 = t56 * qJD(1);
t360 = t323 * qJD(2);
t206 = -t230 * pkin(2) + t231 * qJ(3);
t200 = -t230 * pkin(3) + t206;
t78 = t200 - t95;
t201 = t231 * pkin(2) + t230 * qJ(3);
t197 = t231 * pkin(3) + t201;
t79 = t197 - t93;
t359 = -t78 * t93 + t79 * t95;
t256 = t169 * t124 - t131 * t264;
t337 = -t125 - t263;
t358 = -t170 * t337 + t256;
t355 = -t194 / 0.2e1;
t344 = -t195 / 0.2e1;
t212 = Icges(6,5) * t194 + Icges(6,6) * t195;
t141 = t212 * t169;
t354 = t212 * t170;
t353 = 0.2e1 * t360;
t291 = -m(6) * t342 / 0.2e1 - m(5) * t335 / 0.2e1;
t292 = (-t94 * t231 + t89) * t323;
t351 = t291 - t292;
t349 = 0.2e1 * t240;
t138 = -t150 + t197;
t196 = t200 - t219;
t348 = t138 * t219 - t196 * t150;
t276 = Icges(6,2) * t195;
t214 = t276 + t279;
t205 = t194 * t214;
t216 = Icges(6,1) * t194 + t260;
t258 = t195 * t216;
t336 = t205 - t258;
t100 = t336 * t169 + t354;
t347 = -t169 / 0.2e1;
t346 = t170 / 0.2e1;
t345 = -t213 / 0.2e1;
t295 = m(6) * t188;
t340 = t170 * t169;
t338 = -t124 + t262;
t101 = t170 * t336 - t141;
t266 = t170 * t101;
t268 = t169 * t100;
t334 = t266 / 0.4e1 + t268 / 0.4e1;
t333 = t375 / 0.4e1 + t148 * t343 / 0.2e1;
t332 = t281 / 0.4e1 - ((-t78 - t95) * t170 + (t79 + t93) * t169) * t295 / 0.2e1;
t226 = t205 / 0.2e1 + t217 * t355 + t215 * t344 - t258 / 0.2e1;
t331 = t169 ^ 2;
t330 = t170 ^ 2;
t328 = 0.4e1 * qJD(2);
t327 = 0.2e1 * qJD(4);
t36 = -t170 * t263 - t256;
t255 = -t169 * t125 + t130 * t264;
t37 = -t170 * t262 + t255;
t19 = -t169 * t36 + t170 * t37;
t322 = -t19 / 0.2e1;
t267 = t169 * t195;
t254 = t131 * t267 + t365;
t38 = -t169 * t263 + t254;
t253 = t130 * t267 + t364;
t39 = -t169 * t262 + t253;
t20 = -t169 * t38 + t170 * t39;
t321 = -t20 / 0.2e1;
t320 = m(5) * t348;
t317 = m(5) * (t138 * t230 + t196 * t231);
t316 = m(6) * ((t79 - t93) * t149 + (-t78 + t95) * t148);
t313 = m(6) * t359;
t310 = m(6) * (t148 * t78 - t149 * t79);
t309 = m(6) * (-t148 * t95 + t149 * t93);
t75 = t78 * t231;
t308 = m(6) * (t79 * t230 + t75);
t300 = m(4) * ((t231 * rSges(4,3) + t206) * t231 + (t230 * rSges(4,3) + t201) * t230);
t289 = t56 * t360;
t80 = t200 - t94;
t288 = t78 - t80;
t270 = t148 * t188;
t269 = t149 * t188;
t244 = qJD(5) * t169;
t243 = qJD(5) * t170;
t67 = t375 / 0.2e1;
t198 = (-t217 + t276) * t195 + (t215 + t216 + t260) * t194;
t225 = t322 + (-t216 * t169 - t126) * t194 / 0.2e1 + (t214 * t169 - t130) * t344 + t198 * t347 + t170 * t345;
t224 = t321 + (-t216 * t170 - t127) * t355 + (-t214 * t170 + t131) * t344 + t169 * t345 + t198 * t346;
t211 = t226 + t376;
t210 = -t281 / 0.2e1 + t67 + t226;
t208 = -t266 / 0.4e1 - t332 + t334 - (-t357 + t100) * t169 / 0.4e1;
t207 = -t268 / 0.4e1 - t333 + t334 + (-t356 + t101) * t374;
t189 = t282 - t284;
t77 = -t148 * t170 - t149 * t169;
t58 = t134 * t170 + (-rSges(6,1) * t267 + t352) * t169;
t41 = 0.2e1 * t91;
t32 = -t226 + t309;
t31 = -t226 + t310;
t28 = t56 * t58;
t24 = t300 + t308 + t317;
t23 = (-t126 * t170 - t127 * t169) * t194 + t254 + t255;
t14 = t316 / 0.2e1;
t13 = t313 + t320;
t12 = t290 - t351;
t11 = t290 + t351;
t10 = t291 + t292 - t290;
t9 = (t36 + (-t366 + t263) * t170 + t253) * t170 + (-t170 * t338 - t23 + t37) * t169;
t8 = (t23 - t38) * t170 + (-t39 + (t367 - t262) * t169 + t358) * t169;
t7 = (t38 + (t338 - t367) * t170) * t170 + (t338 * t169 - t253 + t364 + t39) * t169;
t6 = (-t358 - t36) * t170 + (-t37 + (t337 + t366) * t169 + t365) * t169;
t5 = -t316 / 0.2e1 + t207 + t208 - t226;
t4 = t356 * t374 + t14 + t208 + t211 + t333;
t3 = -t375 / 0.4e1 + t14 + t207 + t210 + t332;
t2 = (t8 / 0.2e1 + t322) * t170 + (t321 - t6 / 0.2e1) * t169;
t1 = m(6) * t28 + (t7 / 0.2e1 + t19 / 0.2e1) * t170 + (t20 / 0.2e1 - t9 / 0.2e1) * t169;
t15 = [0, t369 / 0.2e1, 0, 0, t77 * t285 + t289; -t369 / 0.2e1, -m(6) * t288 * t79 * t328 / 0.4e1 + t24 * qJD(3) + t13 * qJD(4) + t31 * qJD(5), qJD(2) * t24 + qJD(4) * t11 + t372, t13 * qJD(2) + t11 * qJD(3) + t4 * qJD(5) - t371 + (-t359 * t323 - t348 * t325) * t327, t31 * qJD(2) + t257 + t4 * qJD(4) + (-t7 / 0.2e1 + t225) * t243 + (t9 / 0.2e1 + t224) * t244 + (-t361 / 0.2e1 + ((-t189 * t79 + t269) * t170 + (-t189 * t78 - t270) * t169 - t28) * qJD(5)) * m(6); 0, t10 * qJD(4) + t41 * qJD(5) + (-t308 / 0.4e1 - t317 / 0.4e1 - t300 / 0.4e1) * t328 + (-t231 * t80 + t75) * t353, 0, t10 * qJD(2) + t349 * qJD(5) + t290 * t327, t41 * qJD(2) + t349 * qJD(4) + (-t169 * t230 + t170 * t231) * t189 * t285; 0, t12 * qJD(3) + t371 + t3 * qJD(5) + (-t313 / 0.4e1 - t320 / 0.4e1) * t328 + (t286 * t79 - t288 * t93) * t353, qJD(2) * t12 - t372, qJD(2) * t311 + t32 * qJD(5), t3 * qJD(2) - t257 + t32 * qJD(4) + (-t8 / 0.2e1 - t225) * t243 + (t6 / 0.2e1 - t224) * t244 + ((-t189 * t93 - t269) * t170 + (-t189 * t95 + t270) * t169) * t285; t289, t349 * qJD(3) + t5 * qJD(4) + t1 * qJD(5) + t323 * t361 + (t211 + t67 - t310 + (-t356 / 0.2e1 + t288 * t295) * t170) * qJD(2), qJD(2) * t349 + qJD(4) * t368, t5 * qJD(2) + t2 * qJD(5) + t257 + (t210 - t309 + t376) * qJD(4), t1 * qJD(2) + t2 * qJD(4) + (m(6) * (t58 * t77 + (t330 + t331) * t189 * t188) + (t330 * t141 - t340 * t354) * t346 + (-t141 * t340 + t331 * t354) * t347) * qJD(5);];
Cq = t15;
