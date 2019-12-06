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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:51:15
% EndTime: 2019-12-05 17:51:23
% DurationCPUTime: 3.74s
% Computational Cost: add. (22730->238), mult. (17575->341), div. (0->0), fcn. (18501->10), ass. (0->150)
t223 = qJ(1) + pkin(8);
t222 = qJ(3) + t223;
t219 = cos(t222);
t215 = t219 * qJ(4);
t218 = sin(t222);
t228 = cos(qJ(5));
t225 = cos(pkin(9));
t226 = sin(qJ(5));
t272 = t225 * t226;
t193 = t218 * t272 + t219 * t228;
t271 = t225 * t228;
t194 = t218 * t271 - t219 * t226;
t261 = t194 * rSges(6,1) - t193 * rSges(6,2);
t224 = sin(pkin(9));
t326 = pkin(4) * t225 + pkin(3) + (rSges(6,3) + pkin(7)) * t224;
t125 = t326 * t218 - t215 + t261;
t243 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t223);
t114 = t125 + t243;
t110 = t219 * t114;
t244 = rSges(5,1) * t225 + pkin(3);
t279 = t218 * t224;
t168 = -rSges(5,2) * t279 - t219 * rSges(5,3) + t244 * t218 - t215;
t156 = t168 + t243;
t143 = t219 * t156;
t276 = t219 * t224;
t169 = rSges(5,2) * t276 - t244 * t219 + (-rSges(5,3) - qJ(4)) * t218;
t242 = -cos(qJ(1)) * pkin(1) - pkin(2) * cos(t223);
t157 = t169 + t242;
t277 = t219 * t168;
t278 = t219 * t125;
t320 = m(6) / 0.2e1;
t195 = -t218 * t228 + t219 * t272;
t196 = t218 * t226 + t219 * t271;
t329 = -t196 * rSges(6,1) + t195 * rSges(6,2);
t333 = -t218 * qJ(4) - t326 * t219 + t329;
t337 = m(5) / 0.2e1;
t339 = t333 + t242;
t294 = (-t110 - t278 + (-t339 - t333) * t218) * t320 + (-t143 - t277 + (-t157 - t169) * t218) * t337;
t269 = t218 * t333 + t278;
t65 = -t218 * t339 - t110;
t295 = (t65 + t269) * t320 + (-t143 + t277 + (-t157 + t169) * t218) * t337;
t6 = t295 - t294;
t348 = t6 * qJD(1);
t166 = t193 * rSges(6,1) + t194 * rSges(6,2);
t167 = -t195 * rSges(6,1) - t196 * rSges(6,2);
t56 = -t125 * t166 - t167 * t333;
t347 = m(6) * t56;
t51 = -t114 * t166 - t167 * t339;
t346 = m(6) * t51;
t45 = -t114 * t333 + t125 * t339;
t154 = rSges(6,3) * t276 - t329;
t232 = t224 * (-t225 * rSges(6,3) + (rSges(6,1) * t228 - rSges(6,2) * t226) * t224);
t344 = t225 * t154 + t219 * t232;
t285 = Icges(6,4) * t226;
t201 = -Icges(6,5) * t225 + (Icges(6,1) * t228 - t285) * t224;
t258 = t201 + (-Icges(6,2) * t228 - t285) * t224;
t340 = t258 * t226;
t184 = Icges(6,4) * t196;
t148 = -Icges(6,2) * t195 + Icges(6,6) * t276 + t184;
t183 = Icges(6,4) * t195;
t152 = -Icges(6,1) * t196 - Icges(6,5) * t276 + t183;
t267 = t193 * t148 + t194 * t152;
t239 = -t195 * t148 - t196 * t152;
t317 = m(5) * (-t169 * t156 + t157 * t168);
t302 = m(4) * (t242 * (t218 * rSges(4,1) + t219 * rSges(4,2)) - (-t219 * rSges(4,1) + t218 * rSges(4,2)) * t243);
t199 = -Icges(6,3) * t225 + (Icges(6,5) * t228 - Icges(6,6) * t226) * t224;
t284 = Icges(6,4) * t228;
t200 = -Icges(6,6) * t225 + (-Icges(6,2) * t226 + t284) * t224;
t335 = (-t195 * t200 + t196 * t201 + t199 * t276) * t225;
t257 = qJD(1) + qJD(3);
t206 = (-Icges(6,5) * t226 - Icges(6,6) * t228) * t224;
t273 = t225 * t206;
t331 = -t273 / 0.2e1 - t224 * t340 / 0.2e1;
t144 = -Icges(6,5) * t194 + Icges(6,6) * t193 - Icges(6,3) * t279;
t286 = Icges(6,4) * t194;
t147 = Icges(6,2) * t193 - Icges(6,6) * t279 - t286;
t182 = Icges(6,4) * t193;
t150 = -Icges(6,1) * t194 - Icges(6,5) * t279 + t182;
t327 = t144 * t276 - t195 * t147 + t196 * t150;
t325 = 0.4e1 * qJD(1);
t153 = -rSges(6,3) * t279 - t261;
t127 = t225 * t153 - t218 * t232;
t315 = m(6) * ((-t114 + t125) * t344 + (-t333 + t339) * t127);
t312 = m(6) * (t56 + t51);
t309 = m(6) * t45;
t306 = m(6) * t65;
t305 = m(6) * t269;
t304 = m(6) * (-t127 * t219 - t218 * t344);
t300 = m(5) * (-t157 * t218 - t143);
t299 = m(5) * (-t169 * t218 - t277);
t297 = m(6) * (t218 * t166 - t219 * t167);
t109 = t297 / 0.2e1;
t69 = t304 / 0.2e1;
t22 = t109 + t69;
t296 = t22 * qJD(4);
t290 = -t127 * t166 - t167 * t344;
t288 = m(6) * qJD(5);
t283 = (t193 * t200 - t194 * t201 - t199 * t279) * t225;
t208 = (-Icges(6,1) * t226 - t284) * t224;
t259 = t200 - t208;
t282 = (-t273 + (-t259 * t228 - t340) * t224) * t225;
t274 = t224 * t228;
t266 = -t193 * t147 + t194 * t150;
t265 = -Icges(6,1) * t193 + t147 - t286;
t264 = Icges(6,1) * t195 + t148 + t184;
t263 = Icges(6,2) * t194 + t150 + t182;
t262 = -Icges(6,2) * t196 - t152 - t183;
t237 = t144 * t279 + t266;
t145 = Icges(6,5) * t196 - Icges(6,6) * t195 + Icges(6,3) * t276;
t60 = t145 * t276 + t239;
t11 = -t283 + (t267 * t219 + (t237 + t239 - t60) * t218) * t224;
t58 = -t145 * t279 + t267;
t12 = t335 + (-(-t267 - t327) * t218 + t237 * t219 + (-t239 - t266) * t219 - t218 * t58 + (-t145 * t218 ^ 2 + (-t144 * t218 - t145 * t219) * t219) * t224) * t224;
t31 = -t283 + (t218 * t237 + t219 * t58) * t224;
t32 = -t335 + (-t218 * t327 + t219 * t60) * t224;
t2 = ((-t31 / 0.2e1 + t11 / 0.2e1) * t219 + (-t12 / 0.2e1 - t32 / 0.2e1) * t218) * t224;
t23 = t69 - t297 / 0.2e1;
t256 = t23 * qJD(4) + t2 * qJD(5);
t251 = -t279 / 0.2e1;
t249 = t276 / 0.2e1;
t246 = -t274 / 0.2e1;
t245 = t274 / 0.2e1;
t240 = t200 * t246 + t208 * t245 + t331;
t236 = t312 / 0.2e1 + t240;
t233 = t200 * t245 + t208 * t246 - t331;
t160 = Icges(6,5) * t193 + Icges(6,6) * t194;
t54 = -t225 * t160 + (-t263 * t226 - t265 * t228) * t224;
t161 = -Icges(6,5) * t195 - Icges(6,6) * t196;
t55 = -t225 * t161 + (-t262 * t226 - t264 * t228) * t224;
t84 = t258 * t193 + t259 * t194 - t206 * t279;
t85 = -t258 * t195 - t259 * t196 + t206 * t276;
t231 = -t11 * t276 / 0.2e1 - t282 + (t54 + t84) * t251 + (t12 + t32) * t279 / 0.2e1 + (t31 + t55 + t85) * t249;
t209 = (-rSges(6,1) * t226 - rSges(6,2) * t228) * t224;
t140 = t225 * t167 + t209 * t276;
t139 = -t225 * t166 + t209 * t279;
t107 = (-t166 * t219 - t167 * t218) * t224;
t50 = t299 - t305;
t49 = t300 + t306;
t42 = t240 + t347;
t37 = t240 + t346;
t24 = t302 + t309 + t317;
t21 = t109 - t304 / 0.2e1;
t18 = t22 * qJD(5);
t17 = t21 * qJD(5);
t13 = t315 / 0.2e1;
t8 = t294 + t295;
t5 = -t315 / 0.2e1 + t236;
t4 = t13 + t236;
t3 = t13 - t312 / 0.2e1 + t233;
t1 = [t24 * qJD(3) + t49 * qJD(4) + t37 * qJD(5), 0, t24 * qJD(1) + t8 * qJD(4) + t4 * qJD(5) + 0.2e1 * (t302 / 0.2e1 + t317 / 0.2e1 + t45 * t320) * qJD(3), t49 * qJD(1) + t8 * qJD(3) + t18, t37 * qJD(1) + t4 * qJD(3) + (m(6) * (-t139 * t114 + t140 * t339 + t290) + t231) * qJD(5) + t296; 0, 0, 0, 0, t107 * t288; -t6 * qJD(4) + t5 * qJD(5) + (-t302 / 0.4e1 - t317 / 0.4e1 - t309 / 0.4e1) * t325, 0, t50 * qJD(4) + t42 * qJD(5), t50 * qJD(3) + t18 - t348, t5 * qJD(1) + t42 * qJD(3) + (m(6) * (-t139 * t125 + t140 * t333 + t290) + t231) * qJD(5) + t296; t6 * qJD(3) + t17 + (-t300 / 0.4e1 - t306 / 0.4e1) * t325, 0, t348 + t17 + 0.4e1 * (t305 / 0.4e1 - t299 / 0.4e1) * qJD(3), 0, (t139 * t218 + t140 * t219) * t288 + t257 * t21; (t233 - t346) * qJD(1) + t3 * qJD(3) + t256, 0, t3 * qJD(1) + (t233 - t347) * qJD(3) + t256, t257 * t23, (m(6) * ((-t153 * t219 - t154 * t218) * t224 * t107 - t127 * t139 + t140 * t344) - t225 * (-t282 + (-t218 * t54 + t219 * t55) * t224) / 0.2e1 + (-t84 * t225 - (-t160 * t279 + t193 * t263 + t194 * t265) * t279 + (-t161 * t279 + t193 * t262 + t194 * t264) * t276) * t251 + (-t85 * t225 - (t160 * t276 - t195 * t263 - t196 * t265) * t279 + (t161 * t276 - t195 * t262 - t196 * t264) * t276) * t249) * qJD(5) + t257 * t2;];
Cq = t1;
