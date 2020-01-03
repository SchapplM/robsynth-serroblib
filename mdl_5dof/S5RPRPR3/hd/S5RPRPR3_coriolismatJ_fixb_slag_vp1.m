% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:58
% EndTime: 2020-01-03 11:36:10
% DurationCPUTime: 3.59s
% Computational Cost: add. (22730->242), mult. (17575->339), div. (0->0), fcn. (18501->10), ass. (0->151)
t227 = qJ(1) + pkin(8);
t225 = qJ(3) + t227;
t222 = cos(t225);
t215 = t222 * qJ(4);
t221 = sin(t225);
t228 = sin(pkin(9));
t229 = cos(pkin(9));
t232 = cos(qJ(5));
t230 = sin(qJ(5));
t265 = t229 * t230;
t187 = -t221 * t265 - t222 * t232;
t264 = t229 * t232;
t188 = t221 * t264 - t222 * t230;
t241 = t188 * rSges(6,1) + t187 * rSges(6,2);
t123 = -t215 + (pkin(4) * t229 + pkin(3) + (rSges(6,3) + pkin(7)) * t228) * t221 + t241;
t242 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t227);
t114 = t123 + t242;
t344 = t114 - t123;
t255 = pkin(2) * cos(t227) + cos(qJ(1)) * pkin(1);
t189 = -t221 * t232 + t222 * t265;
t190 = t221 * t230 + t222 * t264;
t270 = t222 * t228;
t153 = t190 * rSges(6,1) - t189 * rSges(6,2) + rSges(6,3) * t270;
t256 = t222 * pkin(3) + t221 * qJ(4);
t269 = t222 * t229;
t334 = pkin(4) * t269 + pkin(7) * t270 + t153 + t256;
t337 = t334 + t255;
t110 = t337 * t221;
t167 = rSges(5,1) * t269 - rSges(5,2) * t270 + t221 * rSges(5,3) + t256;
t156 = t167 + t255;
t142 = t156 * t221;
t271 = t221 * t228;
t166 = -rSges(5,2) * t271 - t215 - t222 * rSges(5,3) + (rSges(5,1) * t229 + pkin(3)) * t221;
t155 = t166 + t242;
t157 = t167 * t221;
t316 = m(6) / 0.2e1;
t330 = m(5) / 0.2e1;
t338 = t221 * t334;
t290 = (t110 + t338 + (-t114 - t123) * t222) * t316 + (t142 + t157 + (-t155 - t166) * t222) * t330;
t291 = (-t344 * t222 + t110 - t338) * t316 + (t142 - t157 + (-t155 + t166) * t222) * t330;
t324 = t290 - t291;
t343 = t324 * qJD(1);
t164 = t187 * rSges(6,1) - t188 * rSges(6,2);
t165 = t189 * rSges(6,1) + t190 * rSges(6,2);
t51 = t114 * t164 - t165 * t337;
t342 = m(6) * t51;
t56 = t123 * t164 - t165 * t334;
t341 = m(6) * t56;
t45 = t114 * t334 - t123 * t337;
t144 = Icges(6,5) * t190 - Icges(6,6) * t189 + Icges(6,3) * t270;
t180 = Icges(6,4) * t190;
t148 = Icges(6,2) * t189 - Icges(6,6) * t270 - t180;
t179 = Icges(6,4) * t189;
t150 = Icges(6,1) * t190 + Icges(6,5) * t270 - t179;
t239 = t189 * t148 + t190 * t150;
t60 = t144 * t270 + t239;
t336 = t222 * t60;
t280 = Icges(6,4) * t230;
t195 = -Icges(6,5) * t229 + (Icges(6,1) * t232 - t280) * t228;
t257 = t195 + (-Icges(6,2) * t232 - t280) * t228;
t335 = t230 * t257;
t196 = -t229 * rSges(6,3) + (rSges(6,1) * t232 - rSges(6,2) * t230) * t228;
t333 = -t229 * t153 - t196 * t270;
t332 = -t187 * t148 + t188 * t150;
t313 = m(5) * (t155 * t167 - t166 * t156);
t298 = m(4) * (-t255 * (t221 * rSges(4,1) + t222 * rSges(4,2)) + t242 * (t222 * rSges(4,1) - t221 * rSges(4,2)));
t193 = -Icges(6,3) * t229 + (Icges(6,5) * t232 - Icges(6,6) * t230) * t228;
t279 = Icges(6,4) * t232;
t194 = -Icges(6,6) * t229 + (-Icges(6,2) * t230 + t279) * t228;
t327 = (-t189 * t194 + t190 * t195 + t193 * t270) * t229;
t254 = qJD(1) + qJD(3);
t201 = (-Icges(6,5) * t230 - Icges(6,6) * t232) * t228;
t266 = t229 * t201;
t323 = -t266 / 0.2e1 - t228 * t335 / 0.2e1;
t281 = Icges(6,4) * t188;
t146 = Icges(6,2) * t187 + Icges(6,6) * t271 + t281;
t178 = Icges(6,4) * t187;
t149 = Icges(6,1) * t188 + Icges(6,5) * t271 + t178;
t263 = -t189 * t146 + t190 * t149;
t320 = 0.4e1 * qJD(1);
t152 = rSges(6,3) * t271 + t241;
t125 = t229 * t152 + t196 * t271;
t311 = m(6) * (t344 * t333 + (-t334 + t337) * t125);
t308 = m(6) * (t56 + t51);
t304 = m(6) * t45;
t301 = m(6) * (-t114 * t222 + t110);
t300 = m(6) * (-t123 * t222 + t338);
t299 = m(6) * (t125 * t222 + t333 * t221);
t296 = m(5) * (-t155 * t222 + t142);
t295 = m(5) * (-t166 * t222 + t157);
t293 = m(6) * (-t221 * t164 + t222 * t165);
t109 = t293 / 0.2e1;
t69 = t299 / 0.2e1;
t22 = t109 + t69;
t292 = t22 * qJD(4);
t286 = -t125 * t164 - t165 * t333;
t284 = m(6) * qJD(5);
t278 = (t187 * t194 + t188 * t195 + t193 * t271) * t229;
t203 = (-Icges(6,1) * t230 - t279) * t228;
t258 = t194 - t203;
t277 = (-t266 + (-t232 * t258 - t335) * t228) * t229;
t143 = Icges(6,5) * t188 + Icges(6,6) * t187 + Icges(6,3) * t271;
t276 = t143 * t222;
t275 = t144 * t221;
t267 = t228 * t232;
t262 = -Icges(6,1) * t187 + t146 + t281;
t261 = -Icges(6,1) * t189 + t148 - t180;
t260 = -Icges(6,2) * t188 + t149 + t178;
t259 = Icges(6,2) * t190 - t150 + t179;
t57 = t143 * t271 + t187 * t146 + t188 * t149;
t58 = -t144 * t271 - t332;
t59 = -t143 * t270 - t263;
t11 = -t278 + ((-t239 + t57 + t60) * t221 + (t59 + (-t275 + t276) * t228 - t58 + t263) * t222) * t228;
t12 = -t327 + (t336 + (t263 + t58 + (t275 + t276) * t228 + t332) * t221) * t228;
t31 = -t278 + (t221 * t57 - t222 * t58) * t228;
t32 = t327 + (t221 * t59 - t336) * t228;
t2 = ((t31 / 0.2e1 - t11 / 0.2e1) * t222 + (t12 / 0.2e1 + t32 / 0.2e1) * t221) * t228;
t23 = t69 - t293 / 0.2e1;
t253 = t23 * qJD(4) + t2 * qJD(5);
t248 = t271 / 0.2e1;
t247 = -t270 / 0.2e1;
t244 = -t267 / 0.2e1;
t243 = t267 / 0.2e1;
t240 = t194 * t244 + t203 * t243 + t323;
t238 = t308 / 0.2e1 + t240;
t235 = t194 * t243 + t203 * t244 - t323;
t158 = Icges(6,5) * t187 - Icges(6,6) * t188;
t54 = -t229 * t158 + (-t230 * t260 - t232 * t262) * t228;
t159 = Icges(6,5) * t189 + Icges(6,6) * t190;
t55 = -t229 * t159 + (-t230 * t259 - t232 * t261) * t228;
t84 = t187 * t257 - t188 * t258 + t201 * t271;
t85 = t189 * t257 + t190 * t258 - t201 * t270;
t234 = t11 * t270 / 0.2e1 - t277 - (t12 + t32) * t271 / 0.2e1 + (t54 + t84) * t248 + (t31 + t55 + t85) * t247;
t204 = (-rSges(6,1) * t230 - rSges(6,2) * t232) * t228;
t138 = t229 * t165 - t204 * t270;
t137 = -t229 * t164 - t204 * t271;
t107 = (t164 * t222 + t165 * t221) * t228;
t50 = t295 + t300;
t49 = t296 + t301;
t42 = t240 + t341;
t37 = t240 + t342;
t24 = t298 + t304 + t313;
t21 = t109 - t299 / 0.2e1;
t18 = t22 * qJD(5);
t17 = t21 * qJD(5);
t13 = t311 / 0.2e1;
t8 = t290 + t291;
t5 = -t311 / 0.2e1 + t238;
t4 = t13 + t238;
t3 = t13 - t308 / 0.2e1 + t235;
t1 = [t24 * qJD(3) + t49 * qJD(4) + t37 * qJD(5), 0, t24 * qJD(1) + t8 * qJD(4) + t4 * qJD(5) + 0.2e1 * (t298 / 0.2e1 + t313 / 0.2e1 + t45 * t316) * qJD(3), t49 * qJD(1) + t8 * qJD(3) + t18, t37 * qJD(1) + t4 * qJD(3) + (m(6) * (t137 * t114 + t138 * t337 + t286) + t234) * qJD(5) + t292; 0, 0, 0, 0, t107 * t284; t324 * qJD(4) + t5 * qJD(5) + (-t298 / 0.4e1 - t313 / 0.4e1 - t304 / 0.4e1) * t320, 0, t50 * qJD(4) + t42 * qJD(5), t50 * qJD(3) + t18 + t343, t5 * qJD(1) + t42 * qJD(3) + (m(6) * (t137 * t123 + t138 * t334 + t286) + t234) * qJD(5) + t292; -t324 * qJD(3) + t17 + (-t296 / 0.4e1 - t301 / 0.4e1) * t320, 0, -t343 + t17 + 0.4e1 * (-t300 / 0.4e1 - t295 / 0.4e1) * qJD(3), 0, (-t137 * t221 - t138 * t222) * t284 + t254 * t21; (t235 - t342) * qJD(1) + t3 * qJD(3) + t253, 0, t3 * qJD(1) + (t235 - t341) * qJD(3) + t253, t254 * t23, (m(6) * ((t152 * t222 - t153 * t221) * t228 * t107 - t125 * t137 + t138 * t333) - t229 * (-t277 + (t221 * t54 - t222 * t55) * t228) / 0.2e1 + (-t84 * t229 + (t158 * t271 + t187 * t260 - t188 * t262) * t271 - (t159 * t271 + t187 * t259 - t188 * t261) * t270) * t248 + (-t85 * t229 + (-t158 * t270 + t189 * t260 + t190 * t262) * t271 - (-t159 * t270 + t189 * t259 + t190 * t261) * t270) * t247) * qJD(5) + t254 * t2;];
Cq = t1;
