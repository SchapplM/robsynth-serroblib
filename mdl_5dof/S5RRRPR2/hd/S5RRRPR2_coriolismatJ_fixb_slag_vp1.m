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
% m [6x1]
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:30:33
% EndTime: 2022-01-20 11:30:39
% DurationCPUTime: 2.84s
% Computational Cost: add. (27928->238), mult. (13707->315), div. (0->0), fcn. (11584->10), ass. (0->162)
t242 = qJ(1) + qJ(2);
t241 = qJ(3) + t242;
t236 = sin(t241);
t237 = cos(t241);
t199 = -rSges(4,1) * t236 - rSges(4,2) * t237;
t238 = sin(t242);
t313 = pkin(2) * t238;
t195 = t199 - t313;
t315 = sin(qJ(1)) * pkin(1);
t185 = t195 - t315;
t200 = t237 * rSges(4,1) - t236 * rSges(4,2);
t239 = cos(t242);
t312 = pkin(2) * t239;
t196 = t200 + t312;
t314 = cos(qJ(1)) * pkin(1);
t186 = t196 + t314;
t117 = -t200 * t185 + t186 * t199;
t275 = t200 * t195 - t196 * t199;
t235 = pkin(9) + t241;
t232 = sin(t235);
t233 = cos(t235);
t311 = pkin(3) * t236;
t183 = -rSges(5,1) * t232 - rSges(5,2) * t233 - t311;
t175 = t183 - t313;
t310 = pkin(3) * t237;
t184 = t233 * rSges(5,1) - t232 * rSges(5,2) + t310;
t176 = t184 + t312;
t278 = t184 * t175 - t176 * t183;
t245 = cos(qJ(5));
t307 = rSges(6,1) * t245;
t266 = pkin(4) + t307;
t243 = sin(qJ(5));
t286 = t232 * t243;
t270 = rSges(6,2) * t286 + t233 * rSges(6,3);
t146 = t233 * pkin(8) - t266 * t232 + t270 - t311;
t136 = t146 - t313;
t284 = t233 * t243;
t215 = rSges(6,2) * t284;
t147 = t310 - t215 + t266 * t233 + (rSges(6,3) + pkin(8)) * t232;
t137 = t147 + t312;
t308 = t147 * t136 - t137 * t146;
t351 = m(6) / 0.2e1;
t352 = m(5) / 0.2e1;
t353 = m(4) / 0.2e1;
t130 = t136 - t315;
t131 = t137 + t314;
t50 = -t147 * t130 + t131 * t146;
t166 = t175 - t315;
t167 = t176 + t314;
t91 = -t184 * t166 + t167 * t183;
t268 = (t308 + t50) * t351 + (t278 + t91) * t352 + (t275 + t117) * t353;
t269 = (-t308 + t50) * t351 + (-t278 + t91) * t352 + (-t275 + t117) * t353;
t3 = t268 - t269;
t368 = t3 * qJD(1);
t240 = Icges(6,4) * t245;
t218 = -Icges(6,2) * t243 + t240;
t219 = Icges(6,1) * t243 + t240;
t367 = t218 + t219;
t323 = m(3) * (t314 * (-rSges(3,1) * t238 - rSges(3,2) * t239) + (t239 * rSges(3,1) - t238 * rSges(3,2)) * t315);
t305 = Icges(6,4) * t243;
t217 = Icges(6,2) * t245 + t305;
t220 = Icges(6,1) * t245 - t305;
t263 = t367 * t245 / 0.2e1 + (-t217 / 0.2e1 + t220 / 0.2e1) * t243;
t221 = rSges(6,1) * t243 + rSges(6,2) * t245;
t193 = t221 * t232;
t116 = t130 * t193;
t194 = t221 * t233;
t63 = -t131 * t194 + t116;
t66 = t146 * t193 - t147 * t194;
t339 = m(6) * (t66 + t63);
t364 = t263 + t339 / 0.2e1;
t64 = t136 * t193 - t137 * t194;
t338 = m(6) * (t66 + t64);
t363 = t263 + t338 / 0.2e1;
t77 = -t176 * t166 + t167 * t175;
t103 = -t196 * t185 + t186 * t195;
t361 = t232 ^ 2;
t360 = t233 ^ 2;
t358 = 0.4e1 * qJD(1);
t355 = 2 * qJD(3);
t346 = m(5) * t77;
t344 = m(5) * t91;
t343 = m(5) * t278;
t340 = m(6) * (t64 + t63);
t337 = m(6) * (t116 + (-t136 * t232 + (-t131 + t137) * t233) * t221);
t336 = m(6) * (t116 + (-t146 * t232 + (-t131 + t147) * t233) * t221);
t335 = m(6) * ((-t137 + t147) * t233 + (t136 - t146) * t232) * t221;
t48 = -t137 * t130 + t131 * t136;
t334 = m(6) * t48;
t332 = m(6) * t50;
t331 = m(6) * t308;
t329 = m(6) * t63;
t328 = m(6) * t64;
t327 = m(6) * t66;
t326 = -t232 / 0.2e1;
t325 = t232 / 0.2e1;
t324 = -t233 / 0.2e1;
t321 = m(4) * t103;
t319 = m(4) * t117;
t318 = m(4) * t275;
t285 = t232 * t245;
t168 = Icges(6,5) * t285 - Icges(6,6) * t286 - Icges(6,3) * t233;
t171 = Icges(6,6) * t232 + t218 * t233;
t264 = t171 * t243 - t168;
t173 = Icges(6,5) * t232 + t220 * t233;
t153 = t173 * t285;
t216 = Icges(6,5) * t245 - Icges(6,6) * t243;
t288 = t216 * t233;
t169 = Icges(6,3) * t232 + t288;
t265 = t169 * t233 - t153;
t283 = t233 * t245;
t276 = t232 * t169 + t173 * t283;
t170 = Icges(6,4) * t285 - Icges(6,2) * t286 - Icges(6,6) * t233;
t213 = Icges(6,4) * t286;
t172 = Icges(6,1) * t285 - Icges(6,5) * t233 - t213;
t277 = -t232 * t168 - t172 * t283;
t85 = -t170 * t284 - t277;
t86 = -t171 * t284 + t276;
t17 = (t264 * t233 - t276 + t86) * t233 + (t264 * t232 + t265 + t85) * t232;
t298 = t170 * t243;
t84 = -t171 * t286 - t265;
t18 = (t84 - t153 + (t169 + t298) * t233 + t277) * t233 + t276 * t232;
t46 = t232 * t84 - t233 * (-(-t172 * t245 + t298) * t232 - t168 * t233);
t47 = t232 * t86 - t233 * t85;
t5 = (t47 / 0.2e1 - t18 / 0.2e1) * t233 + (t17 / 0.2e1 + t46 / 0.2e1) * t232;
t309 = t5 * qJD(5);
t274 = t219 * t232 + t170;
t273 = -t219 * t233 - t171;
t272 = -Icges(6,2) * t285 + t172 - t213;
t271 = -t217 * t233 + t173;
t261 = t340 / 0.2e1 + t263;
t257 = Icges(6,5) * t243 + Icges(6,6) * t245;
t254 = (-t193 * t233 + t194 * t232) * t221;
t248 = (-t217 + t220) * t245 - t367 * t243;
t253 = t233 * t18 / 0.2e1 + (t17 + t46) * t326 + (t216 * t232 + t248 * t233 + t273 * t243 + t271 * t245) * t325 + (t248 * t232 - t274 * t243 + t272 * t245 - t288 + t47) * t324;
t252 = -t263 + (t325 + t326) * (t245 * t170 + t243 * t172);
t250 = t272 * t243 + t274 * t245;
t249 = -t271 * t243 + t273 * t245;
t223 = -rSges(6,2) * t243 + t307;
t188 = t233 * t257;
t187 = t257 * t232;
t134 = -t193 * t232 - t194 * t233;
t62 = t263 + t327;
t59 = t263 + t328;
t53 = t263 + t329;
t42 = t335 / 0.2e1;
t40 = t336 / 0.2e1;
t38 = t337 / 0.2e1;
t29 = -t318 - t331 - t343;
t28 = t319 + t332 + t344;
t19 = t321 + t323 + t334 + t346;
t14 = -t335 / 0.2e1 + t363;
t13 = t42 + t363;
t12 = -t336 / 0.2e1 + t364;
t11 = t40 + t364;
t10 = -t337 / 0.2e1 + t261;
t9 = t38 + t261;
t8 = t42 - t338 / 0.2e1 + t252;
t7 = t40 - t339 / 0.2e1 + t252;
t6 = t38 - t340 / 0.2e1 + t252;
t2 = t268 + t269;
t1 = [qJD(2) * t19 + qJD(3) * t28 + qJD(5) * t53, t19 * qJD(1) + t2 * qJD(3) + t9 * qJD(5) + 0.2e1 * (t323 / 0.2e1 + t103 * t353 + t48 * t351 + t77 * t352) * qJD(2), t28 * qJD(1) + t2 * qJD(2) + t11 * qJD(5) + (t117 * t353 + t50 * t351 + t91 * t352) * t355, 0, t53 * qJD(1) + t9 * qJD(2) + t11 * qJD(3) + (((-t130 * t233 - t131 * t232) * t223 + t254) * m(6) + t253) * qJD(5); -t3 * qJD(3) + t10 * qJD(5) + (-t346 / 0.4e1 - t334 / 0.4e1 - t321 / 0.4e1 - t323 / 0.4e1) * t358, qJD(3) * t29 + qJD(5) * t59, -t368 + t29 * qJD(2) + t13 * qJD(5) + (-t275 * t353 - t278 * t352 - t308 * t351) * t355, 0, t10 * qJD(1) + t59 * qJD(2) + t13 * qJD(3) + (((-t136 * t233 - t137 * t232) * t223 + t254) * m(6) + t253) * qJD(5); t3 * qJD(2) + t12 * qJD(5) + (-t344 / 0.4e1 - t332 / 0.4e1 - t319 / 0.4e1) * t358, t368 + t14 * qJD(5) + 0.4e1 * (t331 / 0.4e1 + t343 / 0.4e1 + t318 / 0.4e1) * qJD(2), qJD(5) * t62, 0, t12 * qJD(1) + t14 * qJD(2) + t62 * qJD(3) + (((-t146 * t233 - t147 * t232) * t223 + t254) * m(6) + t253) * qJD(5); 0, 0, 0, 0, m(6) * t134 * qJD(5); (t252 - t329) * qJD(1) + t6 * qJD(2) + t7 * qJD(3) + t309, t6 * qJD(1) + (t252 - t328) * qJD(2) + t8 * qJD(3) + t309, t7 * qJD(1) + t8 * qJD(2) + (t252 - t327) * qJD(3) + t309, 0, (m(6) * ((t232 * (rSges(6,1) * t285 - t270) + t233 * (rSges(6,1) * t283 + t232 * rSges(6,3) - t215)) * t134 + (t360 + t361) * t223 * t221) + (-t361 * t188 + (t250 * t233 + (t187 + t249) * t232) * t233) * t325 + (-t360 * t187 + (t249 * t232 + (t188 + t250) * t233) * t232) * t324) * qJD(5) + (qJD(1) + qJD(2) + qJD(3)) * t5;];
Cq = t1;
