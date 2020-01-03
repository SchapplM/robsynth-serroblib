% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:34
% DurationCPUTime: 2.75s
% Computational Cost: add. (12345->263), mult. (11029->374), div. (0->0), fcn. (10134->7), ass. (0->178)
t216 = sin(qJ(1));
t217 = cos(qJ(1));
t212 = pkin(7) + qJ(3);
t202 = qJ(4) + t212;
t196 = sin(t202);
t197 = cos(t202);
t168 = rSges(5,1) * t196 + rSges(5,2) * t197;
t200 = sin(t212);
t311 = pkin(3) * t200;
t220 = t168 + t311;
t341 = t220 * t217;
t342 = t220 * t216;
t337 = t216 * t342 + t217 * t341;
t201 = cos(t212);
t176 = rSges(4,1) * t200 + rSges(4,2) * t201;
t213 = t216 ^ 2;
t214 = t217 ^ 2;
t252 = t213 + t214;
t340 = t252 * t176;
t305 = -m(5) * t337 / 0.2e1 - m(4) * t340 / 0.2e1;
t332 = m(5) / 0.2e1;
t161 = t176 * t216;
t162 = t176 * t217;
t89 = t216 * t161 + t162 * t217;
t306 = t337 * t332 + m(4) * t89 / 0.2e1;
t23 = t306 - t305;
t351 = t23 * qJD(1);
t350 = t168 * t252;
t191 = Icges(5,4) * t197;
t165 = -Icges(5,2) * t196 + t191;
t166 = Icges(5,1) * t196 + t191;
t349 = t165 + t166;
t314 = t216 / 0.2e1;
t313 = -t217 / 0.2e1;
t347 = t217 / 0.2e1;
t291 = Icges(4,4) * t200;
t172 = Icges(4,2) * t201 + t291;
t175 = Icges(4,1) * t201 - t291;
t346 = (t175 / 0.2e1 - t172 / 0.2e1) * t200;
t303 = rSges(5,1) * t197;
t169 = -rSges(5,2) * t196 + t303;
t227 = Icges(5,5) * t196 + Icges(5,6) * t197;
t147 = t227 * t216;
t148 = t217 * t227;
t290 = Icges(5,4) * t196;
t167 = Icges(5,1) * t197 - t290;
t128 = Icges(5,5) * t216 + t167 * t217;
t164 = Icges(5,2) * t197 + t290;
t259 = -t164 * t217 + t128;
t278 = t196 * t216;
t181 = Icges(5,4) * t278;
t274 = t197 * t216;
t127 = Icges(5,1) * t274 - Icges(5,5) * t217 - t181;
t260 = -Icges(5,2) * t274 + t127 - t181;
t126 = Icges(5,6) * t216 + t165 * t217;
t261 = -t166 * t217 - t126;
t125 = Icges(5,4) * t274 - Icges(5,2) * t278 - Icges(5,6) * t217;
t262 = t166 * t216 + t125;
t334 = (-t259 * t216 + t260 * t217) * t196 + (t261 * t216 + t262 * t217) * t197;
t307 = (-t213 * t148 + (t216 * t147 + t334) * t217) * t314 + (-t214 * t147 + (t217 * t148 + t334) * t216) * t313;
t153 = t168 * t216;
t154 = t168 * t217;
t339 = t216 * t153 + t217 * t154;
t198 = cos(pkin(7)) * pkin(2) + pkin(1);
t310 = pkin(3) * t201;
t178 = t198 + t310;
t308 = -pkin(5) - qJ(2);
t251 = -pkin(6) + t308;
t193 = t216 * t251;
t199 = t216 * t308;
t249 = t217 * t308;
t254 = -t216 * t178 - t217 * t251;
t129 = rSges(5,1) * t274 - rSges(5,2) * t278 - t217 * rSges(5,3);
t277 = t196 * t217;
t243 = -rSges(5,2) * t277 + t216 * rSges(5,3);
t273 = t197 * t217;
t79 = t216 * t129 + t217 * (rSges(5,1) * t273 + t243);
t44 = -t216 * (t216 * t198 + t249 + t254) + (-t193 + t199 + (t178 - t198) * t217) * t217 + t79;
t6 = t307 + m(5) * (t337 * t169 - t339 * t44);
t345 = t6 * qJD(4);
t87 = -t129 + t254;
t88 = -t193 + (t178 + t303) * t217 + t243;
t344 = t88 * t216 + t217 * t87;
t304 = rSges(4,1) * t201;
t245 = t198 + t304;
t272 = t200 * t216;
t253 = rSges(4,2) * t272 + t217 * rSges(4,3);
t93 = -t216 * t245 - t249 + t253;
t271 = t200 * t217;
t244 = -rSges(4,2) * t271 + t216 * rSges(4,3);
t94 = t217 * t245 - t199 + t244;
t343 = t94 * t216 + t217 * t93;
t195 = Icges(4,4) * t201;
t173 = -Icges(4,2) * t200 + t195;
t174 = Icges(4,1) * t200 + t195;
t234 = t349 * t197 / 0.2e1 + (-t164 / 0.2e1 + t167 / 0.2e1) * t196;
t123 = Icges(5,5) * t274 - Icges(5,6) * t278 - Icges(5,3) * t217;
t163 = Icges(5,5) * t197 - Icges(5,6) * t196;
t282 = t163 * t217;
t124 = Icges(5,3) * t216 + t282;
t236 = t126 * t196 - t123;
t98 = t128 * t274;
t242 = t124 * t217 - t98;
t266 = t216 * t124 + t128 * t273;
t285 = t125 * t196;
t294 = -t216 * t123 - t127 * t273;
t55 = -t126 * t278 - t242;
t56 = -t125 * t277 - t294;
t57 = -t126 * t277 + t266;
t247 = ((t55 - t98 + (t124 + t285) * t217 + t294) * t217 + t266 * t216) * t313 + (t57 * t216 - t217 * t56) * t347 + ((t216 * t236 + t242 + t55 + t56) * t216 + (t216 * (-t127 * t197 + t285) - t266 + t57 + (t123 + t236) * t217) * t217) * t314;
t145 = Icges(4,5) * t216 + t175 * t217;
t255 = -t172 * t217 + t145;
t188 = Icges(4,4) * t272;
t270 = t201 * t216;
t144 = Icges(4,1) * t270 - Icges(4,5) * t217 - t188;
t256 = -Icges(4,2) * t270 + t144 - t188;
t143 = Icges(4,6) * t216 + t173 * t217;
t257 = -t174 * t217 - t143;
t142 = Icges(4,4) * t270 - Icges(4,2) * t272 - Icges(4,6) * t217;
t258 = t174 * t216 + t142;
t335 = (-t255 * t216 + t256 * t217) * t200 + (t257 * t216 + t258 * t217) * t201;
t333 = 0.4e1 * qJD(1);
t331 = m(3) * t252 * (rSges(3,3) + qJ(2));
t330 = m(4) * (t161 * t93 - t162 * t94);
t329 = m(4) * t343;
t222 = t344 * t169;
t324 = m(5) * (-t153 * t341 + t154 * t342 - t222);
t323 = m(5) * (-t222 + (t216 * t341 - t217 * t342) * t168);
t322 = m(5) * (-t341 * t88 + t342 * t87);
t321 = m(5) * (t153 * t87 - t154 * t88);
t320 = m(5) * t344;
t315 = -t216 / 0.2e1;
t283 = t142 * t200;
t269 = t201 * t217;
t223 = t339 * t332;
t232 = m(5) * t350;
t59 = t223 + t232 / 0.2e1;
t267 = t59 * qJD(1);
t140 = Icges(4,5) * t270 - Icges(4,6) * t272 - Icges(4,3) * t217;
t265 = -t216 * t140 - t144 * t269;
t229 = Icges(4,5) * t201 - Icges(4,6) * t200;
t141 = Icges(4,3) * t216 + t217 * t229;
t264 = t216 * t141 + t145 * t269;
t246 = -t169 - t310;
t103 = t145 * t270;
t237 = t141 * t217 - t103;
t235 = t143 * t200 - t140;
t228 = -Icges(4,5) * t200 - Icges(4,6) * t201;
t218 = (-t164 + t167) * t197 - t349 * t196;
t221 = -t247 + (t216 * t163 + t261 * t196 + t259 * t197 + t217 * t218) * t314 + (-t262 * t196 + t260 * t197 + t216 * t218 - t282) * t313;
t219 = -t234 + (t314 + t315) * (t125 * t197 + t127 * t196);
t177 = -rSges(4,2) * t200 + t304;
t156 = t228 * t217;
t155 = t228 * t216;
t114 = t246 * t217;
t112 = t246 * t216;
t73 = -t252 * t311 - t339;
t63 = -t143 * t271 + t264;
t62 = -t142 * t271 - t265;
t61 = -t143 * t272 - t237;
t58 = t223 - t232 / 0.2e1;
t41 = t63 * t216 - t217 * t62;
t40 = t61 * t216 - t217 * (-t216 * (-t144 * t201 + t283) - t140 * t217);
t30 = t323 / 0.2e1;
t29 = t234 + t321;
t28 = t324 / 0.2e1;
t27 = t320 + t329 + t331;
t24 = t305 + t306;
t19 = (t174 / 0.2e1 + t173 / 0.2e1) * t201 + t346 + t330 + t322 + t234;
t16 = (t61 - t103 + (t141 + t283) * t217 + t265) * t217 + t264 * t216;
t15 = (t217 * t235 - t264 + t63) * t217 + (t216 * t235 + t237 + t62) * t216;
t8 = m(5) * (t169 * t350 - t339 * t79) + t307;
t7 = t8 * qJD(4);
t4 = t28 - t323 / 0.2e1 + t247;
t3 = t30 - t324 / 0.2e1 + t247;
t2 = t28 + t30 + t221;
t1 = (t41 / 0.2e1 - t16 / 0.2e1) * t217 + (t15 / 0.2e1 + t40 / 0.2e1) * t216 + t247;
t5 = [t27 * qJD(2) + t19 * qJD(3) + t29 * qJD(4), qJD(1) * t27 + qJD(3) * t24 + qJD(4) * t58, t19 * qJD(1) + t24 * qJD(2) + t2 * qJD(4) + (t16 * t347 + (t257 * t200 + t255 * t201) * t314 + t221 + 0.2e1 * (t112 * t88 + t114 * t87) * t332 + m(4) * (-t343 * t177 + (-t161 * t217 + t162 * t216) * t176) + (t213 / 0.2e1 + t214 / 0.2e1) * t229 + (t15 + t40) * t315 + (-t258 * t200 + t256 * t201 + t41) * t313) * qJD(3), t29 * qJD(1) + t58 * qJD(2) + t2 * qJD(3) + ((-t222 + (-t153 * t217 + t154 * t216) * t168) * m(5) + t221) * qJD(4); t23 * qJD(3) + t59 * qJD(4) + (-t320 / 0.4e1 - t329 / 0.4e1 - t331 / 0.4e1) * t333, 0, t351 + m(5) * (-t112 * t217 + t114 * t216) * qJD(3), t267; (t219 - (t174 + t173) * t201 / 0.2e1 - t346) * qJD(1) - t23 * qJD(2) + t1 * qJD(3) + t4 * qJD(4) + (-t330 / 0.4e1 - t322 / 0.4e1) * t333, -t351, t1 * qJD(1) + (m(4) * (t177 * t340 - (t216 * (rSges(4,1) * t270 - t253) + t217 * (rSges(4,1) * t269 + t244)) * t89) + (t213 * t156 + (-t216 * t155 + t335) * t217) * t314 + (t155 * t214 + (-t156 * t217 + t335) * t216) * t313 + m(5) * (-t112 * t342 - t114 * t341 + t44 * t73) + t307) * qJD(3) + t345, t4 * qJD(1) + t6 * qJD(3) + t345; (t219 - t321) * qJD(1) - t59 * qJD(2) + t3 * qJD(3) + t247 * qJD(4), -t267, t3 * qJD(1) + ((t79 * t73 + (-t112 * t216 - t114 * t217) * t168) * m(5) + t307) * qJD(3) + t7, qJD(1) * t247 + qJD(3) * t8 + t7;];
Cq = t5;
