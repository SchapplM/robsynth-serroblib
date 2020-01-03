% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:22
% EndTime: 2020-01-03 12:07:27
% DurationCPUTime: 2.70s
% Computational Cost: add. (27928->239), mult. (13707->313), div. (0->0), fcn. (11584->10), ass. (0->163)
t245 = qJ(1) + qJ(2);
t243 = qJ(3) + t245;
t238 = sin(t243);
t239 = cos(t243);
t198 = rSges(4,1) * t238 + rSges(4,2) * t239;
t240 = sin(t245);
t311 = pkin(2) * t240;
t194 = t198 + t311;
t312 = sin(qJ(1)) * pkin(1);
t184 = t194 + t312;
t199 = t239 * rSges(4,1) - rSges(4,2) * t238;
t241 = cos(t245);
t236 = pkin(2) * t241;
t195 = t199 + t236;
t244 = cos(qJ(1)) * pkin(1);
t185 = t195 + t244;
t116 = t184 * t199 - t198 * t185;
t275 = -t194 * t199 + t198 * t195;
t237 = pkin(9) + t243;
t232 = sin(t237);
t233 = cos(t237);
t310 = pkin(3) * t238;
t182 = rSges(5,1) * t232 + rSges(5,2) * t233 + t310;
t174 = t182 + t311;
t234 = pkin(3) * t239;
t183 = t233 * rSges(5,1) - rSges(5,2) * t232 + t234;
t175 = t183 + t236;
t276 = -t174 * t183 + t182 * t175;
t246 = sin(qJ(5));
t285 = t232 * t246;
t270 = -rSges(6,2) * t285 - t233 * rSges(6,3);
t248 = cos(qJ(5));
t307 = rSges(6,1) * t248;
t146 = t310 - t233 * pkin(8) + (pkin(4) + t307) * t232 + t270;
t136 = t146 + t311;
t281 = t233 * t248;
t282 = t233 * t246;
t264 = rSges(6,1) * t281 - rSges(6,2) * t282 + t232 * rSges(6,3);
t147 = t233 * pkin(4) + t232 * pkin(8) + t234 + t264;
t137 = t236 + t147;
t308 = -t136 * t147 + t146 * t137;
t348 = m(6) / 0.2e1;
t349 = m(5) / 0.2e1;
t350 = m(4) / 0.2e1;
t130 = t136 + t312;
t131 = t244 + t137;
t50 = t130 * t147 - t146 * t131;
t165 = t174 + t312;
t166 = t175 + t244;
t91 = t165 * t183 - t182 * t166;
t266 = (t308 + t50) * t348 + (t276 + t91) * t349 + (t275 + t116) * t350;
t267 = (-t308 + t50) * t348 + (-t276 + t91) * t349 + (-t275 + t116) * t350;
t3 = t266 - t267;
t370 = t3 * qJD(1);
t242 = Icges(6,4) * t248;
t218 = -Icges(6,2) * t246 + t242;
t363 = Icges(6,1) * t246 + t242;
t369 = t218 + t363;
t320 = m(3) * (-t244 * (rSges(3,1) * t240 + rSges(3,2) * t241) + t312 * (t241 * rSges(3,1) - rSges(3,2) * t240));
t223 = -rSges(6,2) * t246 + t307;
t367 = m(6) * t223;
t305 = Icges(6,4) * t246;
t217 = Icges(6,2) * t248 + t305;
t220 = Icges(6,1) * t248 - t305;
t263 = t369 * t248 / 0.2e1 + (-t217 / 0.2e1 + t220 / 0.2e1) * t246;
t221 = rSges(6,1) * t246 + rSges(6,2) * t248;
t192 = t221 * t232;
t193 = t221 * t233;
t63 = -t130 * t192 - t131 * t193;
t299 = t147 * t193;
t66 = -t146 * t192 - t299;
t336 = m(6) * (t66 + t63);
t365 = t263 + t336 / 0.2e1;
t64 = -t136 * t192 - t137 * t193;
t335 = m(6) * (t66 + t64);
t364 = t263 + t335 / 0.2e1;
t77 = t165 * t175 - t174 * t166;
t103 = t184 * t195 - t194 * t185;
t361 = t369 * t246 + (t217 - t220) * t248;
t212 = Icges(6,4) * t282;
t172 = Icges(6,1) * t281 + t232 * Icges(6,5) - t212;
t271 = -Icges(6,2) * t281 + t172 - t212;
t170 = Icges(6,4) * t281 - Icges(6,2) * t282 + t232 * Icges(6,6);
t273 = t233 * t363 + t170;
t360 = -t246 * t271 - t248 * t273;
t171 = -Icges(6,5) * t233 + t220 * t232;
t272 = -t217 * t232 + t171;
t169 = -Icges(6,6) * t233 + t218 * t232;
t274 = t232 * t363 + t169;
t359 = -t246 * t272 - t248 * t274;
t358 = t232 ^ 2;
t357 = t233 ^ 2;
t355 = 0.4e1 * qJD(1);
t352 = 2 * qJD(3);
t343 = m(5) * t77;
t341 = m(5) * t91;
t340 = m(5) * t276;
t337 = m(6) * (t64 + t63);
t334 = m(6) * ((-t131 + t137) * t233 + (-t130 + t136) * t232) * t221;
t333 = m(6) * (t299 + (-t131 * t233 + (-t130 + t146) * t232) * t221);
t332 = m(6) * (t299 + (-t137 * t233 + (-t136 + t146) * t232) * t221);
t48 = t130 * t137 - t136 * t131;
t331 = m(6) * t48;
t329 = m(6) * t50;
t328 = m(6) * t308;
t326 = m(6) * t63;
t325 = m(6) * t64;
t324 = m(6) * t66;
t323 = -t232 / 0.2e1;
t322 = -t233 / 0.2e1;
t321 = t233 / 0.2e1;
t318 = m(4) * t103;
t316 = m(4) * t116;
t315 = m(4) * t275;
t284 = t232 * t248;
t153 = t171 * t284;
t154 = t172 * t284;
t155 = t169 * t282;
t216 = Icges(6,5) * t248 - Icges(6,6) * t246;
t286 = t216 * t232;
t167 = -Icges(6,3) * t233 + t286;
t297 = t170 * t246;
t254 = t172 * t248 - t297;
t296 = t171 * t248;
t298 = t169 * t246;
t85 = -t167 * t232 - t171 * t281 + t155;
t168 = Icges(6,5) * t281 - Icges(6,6) * t282 + Icges(6,3) * t232;
t86 = t168 * t232 + t254 * t233;
t17 = (t85 + t154 - t155 + (t167 - t297) * t232) * t232 + (-t153 - t86 + (t167 + t254) * t233 + (t296 + t298) * t232) * t233;
t83 = -t167 * t233 - t169 * t285 + t153;
t84 = t168 * t233 + t170 * t285 - t154;
t18 = (t155 - t84 + (t168 - t296) * t233) * t233 + (-t153 + t83 + (t168 + t298) * t232) * t232;
t46 = -t232 * t84 - t233 * t83;
t47 = -t232 * t86 - t233 * t85;
t5 = (-t18 / 0.2e1 - t47 / 0.2e1) * t233 + (t46 / 0.2e1 - t17 / 0.2e1) * t232;
t309 = t5 * qJD(5);
t260 = t337 / 0.2e1 + t263;
t256 = Icges(6,5) * t246 + Icges(6,6) * t248;
t253 = t232 * t17 / 0.2e1 + (-t216 * t233 - t361 * t232 - t274 * t246 + t272 * t248) * t322 + (t18 + t47) * t321 + (t361 * t233 + t273 * t246 - t271 * t248 - t286 + t46) * t323;
t252 = -t263 + (t321 + t322) * (t248 * t170 + t246 * t172);
t187 = t256 * t233;
t186 = t232 * t256;
t134 = -t192 * t232 - t193 * t233;
t62 = t263 + t324;
t59 = t263 + t325;
t53 = t263 + t326;
t42 = t332 / 0.2e1;
t40 = t333 / 0.2e1;
t38 = t334 / 0.2e1;
t29 = -t315 - t328 - t340;
t28 = t316 + t329 + t341;
t19 = t318 + t320 + t331 + t343;
t14 = -t332 / 0.2e1 + t364;
t13 = t42 + t364;
t12 = -t333 / 0.2e1 + t365;
t11 = t40 + t365;
t10 = -t334 / 0.2e1 + t260;
t9 = t38 + t260;
t8 = t42 - t335 / 0.2e1 + t252;
t7 = t40 - t336 / 0.2e1 + t252;
t6 = t38 - t337 / 0.2e1 + t252;
t2 = t266 + t267;
t1 = [qJD(2) * t19 + qJD(3) * t28 + qJD(5) * t53, t19 * qJD(1) + t2 * qJD(3) + t9 * qJD(5) + 0.2e1 * (t320 / 0.2e1 + t103 * t350 + t48 * t348 + t77 * t349) * qJD(2), t28 * qJD(1) + t2 * qJD(2) + t11 * qJD(5) + (t116 * t350 + t50 * t348 + t91 * t349) * t352, 0, t53 * qJD(1) + t9 * qJD(2) + t11 * qJD(3) + ((t130 * t233 - t131 * t232) * t367 + t253) * qJD(5); -t3 * qJD(3) + t10 * qJD(5) + (-t320 / 0.4e1 - t318 / 0.4e1 - t343 / 0.4e1 - t331 / 0.4e1) * t355, qJD(3) * t29 + qJD(5) * t59, -t370 + t29 * qJD(2) + t13 * qJD(5) + (-t275 * t350 - t276 * t349 - t308 * t348) * t352, 0, t10 * qJD(1) + t59 * qJD(2) + t13 * qJD(3) + ((t136 * t233 - t137 * t232) * t367 + t253) * qJD(5); t3 * qJD(2) + t12 * qJD(5) + (-t316 / 0.4e1 - t341 / 0.4e1 - t329 / 0.4e1) * t355, t370 + t14 * qJD(5) + 0.4e1 * (t328 / 0.4e1 + t340 / 0.4e1 + t315 / 0.4e1) * qJD(2), qJD(5) * t62, 0, t12 * qJD(1) + t14 * qJD(2) + t62 * qJD(3) + ((t146 * t233 - t147 * t232) * t367 + t253) * qJD(5); 0, 0, 0, 0, m(6) * t134 * qJD(5); (t252 - t326) * qJD(1) + t6 * qJD(2) + t7 * qJD(3) + t309, t6 * qJD(1) + (t252 - t325) * qJD(2) + t8 * qJD(3) + t309, t7 * qJD(1) + t8 * qJD(2) + (t252 - t324) * qJD(3) + t309, 0, (m(6) * ((t233 * t264 + t232 * (rSges(6,1) * t284 + t270)) * t134 + (t357 + t358) * t223 * t221) + (-t357 * t186 + (t360 * t232 + (t187 - t359) * t233) * t232) * t322 + (t358 * t187 + (t359 * t233 + (-t186 - t360) * t232) * t233) * t323) * qJD(5) + (qJD(1) + qJD(2) + qJD(3)) * t5;];
Cq = t1;
