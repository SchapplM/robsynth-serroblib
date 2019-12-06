% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:19
% EndTime: 2019-12-05 16:17:26
% DurationCPUTime: 3.33s
% Computational Cost: add. (22570->226), mult. (17407->326), div. (0->0), fcn. (18339->8), ass. (0->150)
t225 = pkin(8) + qJ(2);
t224 = qJ(3) + t225;
t221 = cos(t224);
t217 = t221 * qJ(4);
t220 = sin(t224);
t227 = cos(pkin(9));
t241 = rSges(5,1) * t227 + pkin(3);
t226 = sin(pkin(9));
t269 = t220 * t226;
t174 = rSges(5,2) * t269 + t221 * rSges(5,3) - t241 * t220 + t217;
t291 = pkin(2) * sin(t225);
t164 = t174 - t291;
t268 = t221 * t226;
t175 = -rSges(5,2) * t268 + t241 * t221 + (rSges(5,3) + qJ(4)) * t220;
t290 = pkin(2) * cos(t225);
t165 = t175 + t290;
t103 = t164 * t221 + t165 * t220;
t115 = t174 * t221 + t175 * t220;
t315 = m(6) / 0.2e1;
t334 = m(5) / 0.2e1;
t229 = cos(qJ(5));
t228 = sin(qJ(5));
t264 = t227 * t228;
t198 = t220 * t229 - t221 * t264;
t263 = t227 * t229;
t199 = t220 * t228 + t221 * t263;
t238 = t199 * rSges(6,1) + t198 * rSges(6,2);
t321 = pkin(4) * t227 + pkin(3) + (rSges(6,3) + pkin(7)) * t226;
t128 = t220 * qJ(4) + t321 * t221 + t238;
t121 = t128 + t290;
t196 = t220 * t264 + t221 * t229;
t197 = t220 * t263 - t221 * t228;
t324 = -t197 * rSges(6,1) + t196 * rSges(6,2);
t328 = -t321 * t220 + t217 + t324;
t338 = t328 - t291;
t65 = t121 * t220 + t221 * t338;
t68 = t128 * t220 + t221 * t328;
t287 = (t68 + t65) * t315 + (t115 + t103) * t334;
t288 = (-t68 + t65) * t315 + (t103 - t115) * t334;
t6 = t288 - t287;
t346 = t6 * qJD(2);
t172 = -t196 * rSges(6,1) - t197 * rSges(6,2);
t173 = t198 * rSges(6,1) - t199 * rSges(6,2);
t56 = t128 * t173 - t172 * t328;
t345 = m(6) * t56;
t53 = t121 * t173 - t172 * t338;
t344 = m(6) * t53;
t47 = t121 * t328 - t128 * t338;
t156 = rSges(6,3) * t269 - t324;
t205 = -t227 * rSges(6,3) + (rSges(6,1) * t229 - rSges(6,2) * t228) * t226;
t342 = t227 * t156 + t205 * t269;
t277 = Icges(6,4) * t228;
t204 = -Icges(6,5) * t227 + (Icges(6,1) * t229 - t277) * t226;
t256 = t204 + (-Icges(6,2) * t229 - t277) * t226;
t339 = t228 * t256;
t186 = Icges(6,4) * t197;
t150 = -Icges(6,2) * t196 + Icges(6,6) * t269 + t186;
t185 = Icges(6,4) * t196;
t154 = -Icges(6,1) * t197 - Icges(6,5) * t269 + t185;
t337 = -t196 * t150 - t197 * t154;
t336 = t198 * t150 - t199 * t154;
t312 = m(5) * (-t175 * t164 + t165 * t174);
t296 = m(4) * (t290 * (-t220 * rSges(4,1) - t221 * rSges(4,2)) + (t221 * rSges(4,1) - t220 * rSges(4,2)) * t291);
t202 = -Icges(6,3) * t227 + (Icges(6,5) * t229 - Icges(6,6) * t228) * t226;
t276 = Icges(6,4) * t229;
t203 = -Icges(6,6) * t227 + (-Icges(6,2) * t228 + t276) * t226;
t331 = (-t196 * t203 + t197 * t204 + t202 * t269) * t227;
t147 = Icges(6,5) * t197 - Icges(6,6) * t196 + Icges(6,3) * t269;
t330 = t147 * t268;
t255 = qJD(2) + qJD(3);
t208 = (-Icges(6,5) * t228 - Icges(6,6) * t229) * t226;
t265 = t227 * t208;
t326 = -t265 / 0.2e1 - t226 * t339 / 0.2e1;
t59 = t330 + t336;
t325 = t59 - t330;
t149 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t268;
t278 = Icges(6,4) * t199;
t152 = Icges(6,2) * t198 + Icges(6,6) * t268 + t278;
t187 = Icges(6,4) * t198;
t155 = Icges(6,1) * t199 + Icges(6,5) * t268 + t187;
t322 = t149 * t269 - t196 * t152 + t197 * t155;
t320 = 0.4e1 * qJD(2);
t158 = rSges(6,3) * t268 + t238;
t131 = t227 * t158 + t205 * t268;
t310 = m(6) * ((t121 - t128) * t342 + (-t328 + t338) * t131);
t307 = m(6) * (t56 + t53);
t303 = m(6) * t47;
t301 = m(6) * t65;
t300 = m(6) * t68;
t299 = m(6) * (-t131 * t220 + t221 * t342);
t294 = m(5) * t103;
t293 = m(5) * t115;
t292 = m(6) * (-t220 * t172 - t221 * t173);
t109 = t292 / 0.2e1;
t69 = t299 / 0.2e1;
t22 = t109 + t69;
t289 = t22 * qJD(4);
t283 = -t131 * t173 - t172 * t342;
t281 = m(6) * qJD(5);
t57 = t147 * t269 + t337;
t280 = t220 * t57;
t275 = (t198 * t203 + t199 * t204 + t202 * t268) * t227;
t210 = (-Icges(6,1) * t228 - t276) * t226;
t257 = -t203 + t210;
t274 = (-t265 + (t229 * t257 - t339) * t226) * t227;
t266 = t226 * t229;
t262 = -Icges(6,1) * t196 - t150 - t186;
t261 = Icges(6,1) * t198 - t152 - t278;
t260 = -Icges(6,2) * t197 - t154 - t185;
t259 = -Icges(6,2) * t199 + t155 + t187;
t11 = t331 + ((-t322 + t325 - t336) * t221 - t280) * t226;
t60 = t149 * t268 + t198 * t152 + t199 * t155;
t12 = -t275 + ((-t337 + t57 + t60) * t221 + t325 * t220) * t226;
t31 = -t331 + (t221 * t322 + t280) * t226;
t32 = -t275 + (t220 * t59 + t221 * t60) * t226;
t2 = ((t11 / 0.2e1 + t31 / 0.2e1) * t221 + (-t32 / 0.2e1 + t12 / 0.2e1) * t220) * t226;
t23 = t69 - t292 / 0.2e1;
t254 = t23 * qJD(4) + t2 * qJD(5);
t248 = t269 / 0.2e1;
t246 = t268 / 0.2e1;
t243 = -t266 / 0.2e1;
t242 = t266 / 0.2e1;
t237 = t203 * t243 + t210 * t242 + t326;
t235 = t307 / 0.2e1 + t237;
t232 = t203 * t242 + t210 * t243 - t326;
t166 = -Icges(6,5) * t196 - Icges(6,6) * t197;
t54 = -t227 * t166 + (-t228 * t260 + t229 * t262) * t226;
t167 = Icges(6,5) * t198 - Icges(6,6) * t199;
t55 = -t227 * t167 + (-t228 * t259 + t229 * t261) * t226;
t84 = -t196 * t256 + t197 * t257 + t208 * t269;
t85 = t198 * t256 + t199 * t257 + t208 * t268;
t230 = -t12 * t269 / 0.2e1 - t274 - (t11 + t31) * t268 / 0.2e1 + (t55 + t85) * t246 + (t32 + t54 + t84) * t248;
t211 = (-rSges(6,1) * t228 - rSges(6,2) * t229) * t226;
t142 = -t227 * t173 - t211 * t268;
t141 = t227 * t172 + t211 * t269;
t107 = (t172 * t221 - t173 * t220) * t226;
t50 = t293 + t300;
t49 = t294 + t301;
t42 = t237 + t345;
t37 = t237 + t344;
t24 = t296 + t303 + t312;
t21 = t109 - t299 / 0.2e1;
t18 = t22 * qJD(5);
t17 = t21 * qJD(5);
t13 = t310 / 0.2e1;
t7 = t287 + t288;
t5 = -t310 / 0.2e1 + t235;
t4 = t13 + t235;
t3 = t13 - t307 / 0.2e1 + t232;
t1 = [0, 0, 0, 0, t107 * t281; 0, t24 * qJD(3) + t49 * qJD(4) + t37 * qJD(5), t24 * qJD(2) + t7 * qJD(4) + t4 * qJD(5) + 0.2e1 * (t47 * t315 + t312 / 0.2e1 + t296 / 0.2e1) * qJD(3), t49 * qJD(2) + t7 * qJD(3) + t18, t37 * qJD(2) + t4 * qJD(3) + (m(6) * (t142 * t121 + t141 * t338 + t283) + t230) * qJD(5) + t289; 0, -t6 * qJD(4) + t5 * qJD(5) + (-t303 / 0.4e1 - t312 / 0.4e1 - t296 / 0.4e1) * t320, t50 * qJD(4) + t42 * qJD(5), t50 * qJD(3) + t18 - t346, t5 * qJD(2) + t42 * qJD(3) + (m(6) * (t142 * t128 + t141 * t328 + t283) + t230) * qJD(5) + t289; 0, t6 * qJD(3) + t17 + (-t301 / 0.4e1 - t294 / 0.4e1) * t320, t346 + t17 + 0.4e1 * (-t300 / 0.4e1 - t293 / 0.4e1) * qJD(3), 0, (t141 * t220 - t142 * t221) * t281 + t255 * t21; 0, (t232 - t344) * qJD(2) + t3 * qJD(3) + t254, t3 * qJD(2) + (t232 - t345) * qJD(3) + t254, t255 * t23, (m(6) * ((t156 * t221 - t158 * t220) * t226 * t107 - t131 * t142 + t141 * t342) + ((t167 * t268 + t198 * t259 + t199 * t261) * t268 + (t166 * t268 + t198 * t260 + t199 * t262) * t269 - t85 * t227) * t246 + ((t167 * t269 - t196 * t259 + t197 * t261) * t268 + (t166 * t269 - t196 * t260 + t197 * t262) * t269 - t84 * t227) * t248 - t227 * (-t274 + (t220 * t54 + t221 * t55) * t226) / 0.2e1) * qJD(5) + t255 * t2;];
Cq = t1;
