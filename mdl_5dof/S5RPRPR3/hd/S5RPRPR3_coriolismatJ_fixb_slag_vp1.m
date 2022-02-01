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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:20:36
% EndTime: 2022-01-23 09:20:42
% DurationCPUTime: 3.38s
% Computational Cost: add. (22730->228), mult. (17575->328), div. (0->0), fcn. (18501->10), ass. (0->150)
t227 = qJ(1) + pkin(8);
t226 = qJ(3) + t227;
t223 = cos(t226);
t219 = t223 * qJ(4);
t222 = sin(t226);
t229 = cos(pkin(9));
t247 = rSges(5,1) * t229 + pkin(3);
t228 = sin(pkin(9));
t275 = t222 * t228;
t174 = rSges(5,2) * t275 + t223 * rSges(5,3) - t247 * t222 + t219;
t245 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t227);
t161 = t174 + t245;
t274 = t223 * t228;
t175 = -rSges(5,2) * t274 + t247 * t223 + (rSges(5,3) + qJ(4)) * t222;
t244 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t227);
t162 = t175 + t244;
t103 = t161 * t223 + t162 * t222;
t118 = t174 * t223 + t175 * t222;
t319 = m(6) / 0.2e1;
t338 = m(5) / 0.2e1;
t232 = cos(qJ(5));
t230 = sin(qJ(5));
t270 = t229 * t230;
t198 = t222 * t232 - t223 * t270;
t269 = t229 * t232;
t199 = t222 * t230 + t223 * t269;
t242 = t199 * rSges(6,1) + t198 * rSges(6,2);
t325 = pkin(4) * t229 + pkin(3) + (rSges(6,3) + pkin(7)) * t228;
t128 = t222 * qJ(4) + t325 * t223 + t242;
t116 = t128 + t244;
t196 = t222 * t270 + t223 * t232;
t197 = t222 * t269 - t223 * t230;
t328 = -t197 * rSges(6,1) + t196 * rSges(6,2);
t332 = -t325 * t222 + t219 + t328;
t342 = t332 + t245;
t65 = t116 * t222 + t223 * t342;
t68 = t128 * t222 + t223 * t332;
t293 = (t68 + t65) * t319 + (t118 + t103) * t338;
t294 = (-t68 + t65) * t319 + (t103 - t118) * t338;
t6 = t294 - t293;
t350 = t6 * qJD(1);
t172 = -t196 * rSges(6,1) - t197 * rSges(6,2);
t173 = t198 * rSges(6,1) - t199 * rSges(6,2);
t56 = t128 * t173 - t172 * t332;
t349 = m(6) * t56;
t51 = t116 * t173 - t172 * t342;
t348 = m(6) * t51;
t45 = t116 * t332 - t128 * t342;
t158 = rSges(6,3) * t275 - t328;
t205 = -t229 * rSges(6,3) + (rSges(6,1) * t232 - rSges(6,2) * t230) * t228;
t346 = t229 * t158 + t205 * t275;
t283 = Icges(6,4) * t230;
t204 = -Icges(6,5) * t229 + (Icges(6,1) * t232 - t283) * t228;
t262 = t204 + (-Icges(6,2) * t232 - t283) * t228;
t343 = t262 * t230;
t188 = Icges(6,4) * t197;
t152 = -Icges(6,2) * t196 + Icges(6,6) * t275 + t188;
t187 = Icges(6,4) * t196;
t156 = -Icges(6,1) * t197 - Icges(6,5) * t275 + t187;
t341 = -t196 * t152 - t197 * t156;
t340 = t198 * t152 - t199 * t156;
t316 = m(5) * (-t175 * t161 + t162 * t174);
t300 = m(4) * (t244 * (-t222 * rSges(4,1) - t223 * rSges(4,2)) - (t223 * rSges(4,1) - t222 * rSges(4,2)) * t245);
t202 = -Icges(6,3) * t229 + (Icges(6,5) * t232 - Icges(6,6) * t230) * t228;
t282 = Icges(6,4) * t232;
t203 = -Icges(6,6) * t229 + (-Icges(6,2) * t230 + t282) * t228;
t335 = (-t196 * t203 + t197 * t204 + t202 * t275) * t229;
t149 = Icges(6,5) * t197 - Icges(6,6) * t196 + Icges(6,3) * t275;
t334 = t149 * t274;
t261 = qJD(1) + qJD(3);
t210 = (-Icges(6,5) * t230 - Icges(6,6) * t232) * t228;
t271 = t229 * t210;
t330 = -t271 / 0.2e1 - t228 * t343 / 0.2e1;
t59 = t334 + t340;
t329 = t59 - t334;
t151 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t274;
t284 = Icges(6,4) * t199;
t154 = Icges(6,2) * t198 + Icges(6,6) * t274 + t284;
t189 = Icges(6,4) * t198;
t157 = Icges(6,1) * t199 + Icges(6,5) * t274 + t189;
t326 = t151 * t275 - t196 * t154 + t197 * t157;
t324 = 0.4e1 * qJD(1);
t160 = rSges(6,3) * t274 + t242;
t131 = t229 * t160 + t205 * t274;
t314 = m(6) * ((t116 - t128) * t346 + (-t332 + t342) * t131);
t311 = m(6) * (t56 + t51);
t308 = m(6) * t45;
t305 = m(6) * t65;
t304 = m(6) * t68;
t303 = m(6) * (-t131 * t222 + t223 * t346);
t298 = m(5) * t103;
t297 = m(5) * t118;
t296 = m(6) * (-t222 * t172 - t223 * t173);
t109 = t296 / 0.2e1;
t69 = t303 / 0.2e1;
t22 = t109 + t69;
t295 = t22 * qJD(4);
t289 = -t131 * t173 - t172 * t346;
t287 = m(6) * qJD(5);
t57 = t149 * t275 + t341;
t286 = t222 * t57;
t281 = (t198 * t203 + t199 * t204 + t202 * t274) * t229;
t212 = (-Icges(6,1) * t230 - t282) * t228;
t263 = -t203 + t212;
t280 = (-t271 + (t263 * t232 - t343) * t228) * t229;
t272 = t228 * t232;
t268 = -Icges(6,1) * t196 - t152 - t188;
t267 = Icges(6,1) * t198 - t154 - t284;
t266 = -Icges(6,2) * t197 - t156 - t187;
t265 = -Icges(6,2) * t199 + t157 + t189;
t11 = t335 + ((-t326 + t329 - t340) * t223 - t286) * t228;
t60 = t151 * t274 + t198 * t154 + t199 * t157;
t12 = -t281 + ((-t341 + t57 + t60) * t223 + t329 * t222) * t228;
t31 = -t335 + (t223 * t326 + t286) * t228;
t32 = -t281 + (t222 * t59 + t223 * t60) * t228;
t2 = ((t11 / 0.2e1 + t31 / 0.2e1) * t223 + (-t32 / 0.2e1 + t12 / 0.2e1) * t222) * t228;
t23 = t69 - t296 / 0.2e1;
t260 = t23 * qJD(4) + t2 * qJD(5);
t254 = t275 / 0.2e1;
t252 = t274 / 0.2e1;
t249 = -t272 / 0.2e1;
t248 = t272 / 0.2e1;
t241 = t203 * t249 + t212 * t248 + t330;
t239 = t311 / 0.2e1 + t241;
t236 = t203 * t248 + t212 * t249 - t330;
t166 = -Icges(6,5) * t196 - Icges(6,6) * t197;
t54 = -t229 * t166 + (-t266 * t230 + t268 * t232) * t228;
t167 = Icges(6,5) * t198 - Icges(6,6) * t199;
t55 = -t229 * t167 + (-t265 * t230 + t267 * t232) * t228;
t84 = -t262 * t196 + t263 * t197 + t210 * t275;
t85 = t262 * t198 + t263 * t199 + t210 * t274;
t234 = -t12 * t275 / 0.2e1 - t280 - (t11 + t31) * t274 / 0.2e1 + (t55 + t85) * t252 + (t32 + t54 + t84) * t254;
t213 = (-rSges(6,1) * t230 - rSges(6,2) * t232) * t228;
t142 = -t229 * t173 - t213 * t274;
t141 = t229 * t172 + t213 * t275;
t107 = (t172 * t223 - t173 * t222) * t228;
t50 = t297 + t304;
t49 = t298 + t305;
t42 = t241 + t349;
t37 = t241 + t348;
t24 = t300 + t308 + t316;
t21 = t109 - t303 / 0.2e1;
t18 = t22 * qJD(5);
t17 = t21 * qJD(5);
t13 = t314 / 0.2e1;
t7 = t293 + t294;
t5 = -t314 / 0.2e1 + t239;
t4 = t13 + t239;
t3 = t13 - t311 / 0.2e1 + t236;
t1 = [t24 * qJD(3) + t49 * qJD(4) + t37 * qJD(5), 0, t24 * qJD(1) + t7 * qJD(4) + t4 * qJD(5) + 0.2e1 * (t45 * t319 + t316 / 0.2e1 + t300 / 0.2e1) * qJD(3), qJD(1) * t49 + qJD(3) * t7 + t18, t37 * qJD(1) + t4 * qJD(3) + (m(6) * (t142 * t116 + t141 * t342 + t289) + t234) * qJD(5) + t295; 0, 0, 0, 0, t107 * t287; -t6 * qJD(4) + t5 * qJD(5) + (-t308 / 0.4e1 - t316 / 0.4e1 - t300 / 0.4e1) * t324, 0, qJD(4) * t50 + qJD(5) * t42, qJD(3) * t50 + t18 - t350, t5 * qJD(1) + t42 * qJD(3) + (m(6) * (t142 * t128 + t141 * t332 + t289) + t234) * qJD(5) + t295; t6 * qJD(3) + t17 + (-t305 / 0.4e1 - t298 / 0.4e1) * t324, 0, t350 + t17 + 0.4e1 * (-t304 / 0.4e1 - t297 / 0.4e1) * qJD(3), 0, (t141 * t222 - t142 * t223) * t287 + t261 * t21; (t236 - t348) * qJD(1) + t3 * qJD(3) + t260, 0, t3 * qJD(1) + (t236 - t349) * qJD(3) + t260, t261 * t23, (m(6) * ((t158 * t223 - t160 * t222) * t228 * t107 - t131 * t142 + t141 * t346) + ((t167 * t274 + t265 * t198 + t267 * t199) * t274 + (t166 * t274 + t266 * t198 + t268 * t199) * t275 - t85 * t229) * t252 + ((t167 * t275 - t265 * t196 + t267 * t197) * t274 + (t166 * t275 - t266 * t196 + t268 * t197) * t275 - t84 * t229) * t254 - t229 * (-t280 + (t222 * t54 + t223 * t55) * t228) / 0.2e1) * qJD(5) + t261 * t2;];
Cq = t1;
