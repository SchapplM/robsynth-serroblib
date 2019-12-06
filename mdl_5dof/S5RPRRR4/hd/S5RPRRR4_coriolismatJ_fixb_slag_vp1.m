% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:19
% EndTime: 2019-12-05 18:14:25
% DurationCPUTime: 2.34s
% Computational Cost: add. (24670->215), mult. (12027->283), div. (0->0), fcn. (10264->10), ass. (0->146)
t315 = m(6) / 0.2e1;
t316 = m(5) / 0.2e1;
t217 = qJ(1) + pkin(9);
t215 = qJ(3) + t217;
t212 = qJ(4) + t215;
t208 = sin(t212);
t209 = cos(t212);
t220 = cos(qJ(5));
t280 = rSges(6,1) * t220;
t239 = pkin(4) + t280;
t218 = sin(qJ(5));
t261 = t208 * t218;
t244 = -rSges(6,2) * t261 - t209 * rSges(6,3);
t135 = -t209 * pkin(8) + t239 * t208 + t244;
t210 = sin(t215);
t286 = pkin(3) * t210;
t129 = t135 + t286;
t235 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t217);
t119 = t129 + t235;
t258 = t209 * t218;
t194 = rSges(6,2) * t258;
t136 = t194 - t239 * t209 + (-rSges(6,3) - pkin(8)) * t208;
t211 = cos(t215);
t285 = pkin(3) * t211;
t130 = t136 - t285;
t234 = -cos(qJ(1)) * pkin(1) - pkin(2) * cos(t217);
t120 = t130 + t234;
t50 = -t136 * t119 + t120 * t135;
t53 = -t136 * t129 + t130 * t135;
t175 = rSges(5,1) * t208 + rSges(5,2) * t209;
t165 = t175 + t286;
t152 = t165 + t235;
t176 = -t209 * rSges(5,1) + t208 * rSges(5,2);
t166 = t176 - t285;
t153 = t166 + t234;
t84 = -t176 * t152 + t153 * t175;
t96 = -t176 * t165 + t166 * t175;
t282 = (-t53 + t50) * t315 + (-t96 + t84) * t316;
t283 = (t53 + t50) * t315 + (t96 + t84) * t316;
t4 = t283 - t282;
t338 = t4 * qJD(1);
t200 = rSges(6,1) * t218 + rSges(6,2) * t220;
t174 = t200 * t209;
t173 = t200 * t208;
t274 = t129 * t173;
t64 = t130 * t174 - t274;
t337 = m(6) * t64;
t289 = m(4) * (t234 * (rSges(4,1) * t210 + rSges(4,2) * t211) - (-t211 * rSges(4,1) + t210 * rSges(4,2)) * t235);
t278 = Icges(6,4) * t218;
t196 = Icges(6,2) * t220 + t278;
t199 = Icges(6,1) * t220 - t278;
t216 = Icges(6,4) * t220;
t328 = Icges(6,2) * t218 - t216;
t330 = Icges(6,1) * t218 + t216;
t236 = (-t328 / 0.2e1 + t330 / 0.2e1) * t220 + (-t196 / 0.2e1 + t199 / 0.2e1) * t218;
t329 = -t119 * t173 + t120 * t174;
t271 = t135 * t173;
t67 = t136 * t174 - t271;
t305 = m(6) * (t67 + t329);
t336 = t236 + t305 / 0.2e1;
t303 = m(6) * (t67 + t64);
t335 = t236 + t303 / 0.2e1;
t333 = m(6) * t67;
t202 = -rSges(6,2) * t218 + t280;
t332 = m(6) * t202;
t257 = t209 * t220;
t157 = Icges(6,4) * t257 - Icges(6,2) * t258 + Icges(6,6) * t208;
t192 = Icges(6,4) * t258;
t159 = Icges(6,1) * t257 + Icges(6,5) * t208 - t192;
t331 = (t157 * t218 - t159 * t220) * t209;
t77 = -t166 * t152 + t153 * t165;
t326 = (-t328 + t330) * t218 + (t196 - t199) * t220;
t245 = -Icges(6,2) * t257 + t159 - t192;
t247 = t209 * t330 + t157;
t325 = t245 * t218 + t247 * t220;
t191 = Icges(6,4) * t261;
t260 = t208 * t220;
t158 = -Icges(6,1) * t260 + Icges(6,5) * t209 + t191;
t246 = Icges(6,2) * t260 + t158 + t191;
t156 = Icges(6,6) * t209 + t328 * t208;
t248 = -t208 * t330 + t156;
t324 = -t246 * t218 - t248 * t220;
t323 = t208 ^ 2;
t322 = t209 ^ 2;
t321 = 0.4e1 * qJD(1);
t318 = 2 * qJD(4);
t312 = m(5) * t77;
t311 = m(5) * t84;
t309 = m(5) * t96;
t306 = m(6) * (t64 + t329);
t304 = m(6) * (t329 - t64);
t302 = m(6) * (t329 - t67);
t301 = m(6) * (-t274 + t271 + (t130 - t136) * t174);
t48 = -t130 * t119 + t120 * t129;
t300 = m(6) * t48;
t299 = m(6) * t50;
t297 = m(6) * t53;
t295 = m(6) * t329;
t292 = t208 / 0.2e1;
t291 = -t209 / 0.2e1;
t290 = t209 / 0.2e1;
t155 = Icges(6,5) * t257 - Icges(6,6) * t258 + Icges(6,3) * t208;
t237 = -t158 * t220 - t155;
t195 = Icges(6,5) * t220 - Icges(6,6) * t218;
t262 = t195 * t208;
t154 = Icges(6,3) * t209 - t262;
t249 = t208 * t154 + t158 * t257;
t250 = t209 * t154 + t156 * t261;
t268 = t156 * t218;
t74 = t209 * t155 + t157 * t261 - t159 * t260;
t75 = -t156 * t258 + t249;
t76 = t208 * t155 - t331;
t15 = (t250 + t76 + t331) * t209 + (-t75 + (t237 - t268) * t209 + t74 + t249) * t208;
t73 = -t158 * t260 + t250;
t16 = (t74 + (-t155 + t268) * t209 - t249) * t209 + (t237 * t208 + t250 - t73) * t208;
t46 = t208 * t74 + t209 * t73;
t47 = t208 * t76 + t209 * t75;
t2 = (t16 / 0.2e1 + t47 / 0.2e1) * t209 + (-t46 / 0.2e1 + t15 / 0.2e1) * t208;
t284 = t2 * qJD(5);
t232 = t306 / 0.2e1 + t236;
t228 = Icges(6,5) * t218 + Icges(6,6) * t220;
t167 = t208 * t228;
t224 = -t208 * t15 / 0.2e1 + (t16 + t47) * t291 + (t195 * t209 + t326 * t208 - t248 * t218 + t246 * t220) * t290 + (-t326 * t209 - t247 * t218 + t245 * t220 + t262 + t46) * t292;
t223 = -t236 + (t290 + t291) * (t220 * t157 + t218 * t159);
t168 = t228 * t209;
t125 = -t173 * t208 - t174 * t209;
t54 = t236 + t333;
t52 = t236 + t337;
t49 = t236 + t295;
t40 = t301 / 0.2e1;
t38 = t302 / 0.2e1;
t36 = t297 + t309;
t33 = t304 / 0.2e1;
t28 = t299 + t311;
t19 = t289 + t300 + t312;
t18 = -t301 / 0.2e1 + t335;
t17 = t40 + t335;
t12 = -t302 / 0.2e1 + t336;
t11 = t38 + t336;
t10 = -t304 / 0.2e1 + t232;
t9 = t33 + t232;
t8 = t40 - t303 / 0.2e1 + t223;
t7 = t38 - t305 / 0.2e1 + t223;
t6 = t33 - t306 / 0.2e1 + t223;
t3 = t282 + t283;
t1 = [qJD(3) * t19 + qJD(4) * t28 + qJD(5) * t49, 0, t19 * qJD(1) + t3 * qJD(4) + t9 * qJD(5) + 0.2e1 * (t289 / 0.2e1 + t48 * t315 + t77 * t316) * qJD(3), t28 * qJD(1) + t3 * qJD(3) + t11 * qJD(5) + (t50 * t315 + t84 * t316) * t318, t49 * qJD(1) + t9 * qJD(3) + t11 * qJD(4) + ((t119 * t209 + t120 * t208) * t332 + t224) * qJD(5); 0, 0, 0, 0, m(6) * t125 * qJD(5); t4 * qJD(4) + t10 * qJD(5) + (-t289 / 0.4e1 - t312 / 0.4e1 - t300 / 0.4e1) * t321, 0, qJD(4) * t36 + qJD(5) * t52, t338 + t36 * qJD(3) + t17 * qJD(5) + (t53 * t315 + t96 * t316) * t318, t10 * qJD(1) + t52 * qJD(3) + t17 * qJD(4) + ((t129 * t209 + t130 * t208) * t332 + t224) * qJD(5); -t4 * qJD(3) + t12 * qJD(5) + (-t311 / 0.4e1 - t299 / 0.4e1) * t321, 0, -t338 + t18 * qJD(5) + 0.4e1 * (-t309 / 0.4e1 - t297 / 0.4e1) * qJD(3), qJD(5) * t54, t12 * qJD(1) + t18 * qJD(3) + t54 * qJD(4) + ((t135 * t209 + t136 * t208) * t332 + t224) * qJD(5); (t223 - t295) * qJD(1) + t6 * qJD(3) + t7 * qJD(4) + t284, 0, t6 * qJD(1) + (t223 - t337) * qJD(3) + t8 * qJD(4) + t284, t7 * qJD(1) + t8 * qJD(3) + (t223 - t333) * qJD(4) + t284, (m(6) * ((t209 * (rSges(6,1) * t257 + t208 * rSges(6,3) - t194) - t208 * (-rSges(6,1) * t260 - t244)) * t125 + (t322 + t323) * t202 * t200) + (t322 * t167 + (t325 * t208 + (-t168 - t324) * t209) * t208) * t290 + (-t323 * t168 + (t324 * t209 + (t167 - t325) * t208) * t209) * t292) * qJD(5) + (qJD(1) + qJD(3) + qJD(4)) * t2;];
Cq = t1;
