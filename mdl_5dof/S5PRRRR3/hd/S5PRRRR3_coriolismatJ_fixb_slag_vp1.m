% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:05
% EndTime: 2019-12-05 17:06:10
% DurationCPUTime: 2.31s
% Computational Cost: add. (24482->215), mult. (11831->290), div. (0->0), fcn. (10074->8), ass. (0->146)
t211 = pkin(9) + qJ(2);
t209 = qJ(3) + t211;
t206 = qJ(4) + t209;
t202 = sin(t206);
t203 = cos(t206);
t169 = -rSges(5,1) * t202 - rSges(5,2) * t203;
t204 = sin(t209);
t274 = pkin(3) * t204;
t157 = t169 - t274;
t170 = rSges(5,1) * t203 - t202 * rSges(5,2);
t205 = cos(t209);
t273 = pkin(3) * t205;
t158 = t170 + t273;
t242 = t170 * t157 - t158 * t169;
t213 = cos(qJ(5));
t268 = rSges(6,1) * t213;
t233 = pkin(4) + t268;
t212 = sin(qJ(5));
t253 = t202 * t212;
t235 = rSges(6,2) * t253 + t203 * rSges(6,3);
t131 = t203 * pkin(8) - t233 * t202 + t235;
t125 = t131 - t274;
t267 = t212 * rSges(6,2);
t187 = t203 * t267;
t132 = -t187 + t233 * t203 + (rSges(6,3) + pkin(8)) * t202;
t126 = t132 + t273;
t269 = t132 * t125 - t126 * t131;
t305 = m(6) / 0.2e1;
t306 = m(5) / 0.2e1;
t276 = pkin(2) * sin(t211);
t119 = t125 - t276;
t275 = pkin(2) * cos(t211);
t120 = t126 + t275;
t51 = -t132 * t119 + t120 * t131;
t155 = t157 - t276;
t156 = t158 + t275;
t91 = -t170 * t155 + t156 * t169;
t270 = (t269 + t51) * t305 + (t242 + t91) * t306;
t271 = (-t269 + t51) * t305 + (-t242 + t91) * t306;
t5 = t270 - t271;
t320 = t5 * qJD(2);
t210 = Icges(6,4) * t213;
t190 = -Icges(6,2) * t212 + t210;
t191 = Icges(6,1) * t212 + t210;
t319 = t190 + t191;
t279 = m(4) * (t275 * (-rSges(4,1) * t204 - rSges(4,2) * t205) + (rSges(4,1) * t205 - t204 * rSges(4,2)) * t276);
t265 = Icges(6,4) * t212;
t189 = Icges(6,2) * t213 + t265;
t192 = Icges(6,1) * t213 - t265;
t229 = t319 * t213 / 0.2e1 + (-t189 / 0.2e1 + t192 / 0.2e1) * t212;
t193 = t212 * rSges(6,1) + rSges(6,2) * t213;
t165 = t193 * t202;
t110 = t119 * t165;
t166 = t193 * t203;
t63 = -t120 * t166 + t110;
t67 = t131 * t165 - t132 * t166;
t295 = m(6) * (t67 + t63);
t316 = t229 + t295 / 0.2e1;
t64 = t125 * t165 - t126 * t166;
t294 = m(6) * (t67 + t64);
t315 = t229 + t294 / 0.2e1;
t79 = -t158 * t155 + t156 * t157;
t313 = t202 ^ 2;
t312 = t203 ^ 2;
t311 = 0.4e1 * qJD(2);
t308 = 2 * qJD(4);
t302 = m(5) * t79;
t300 = m(5) * t91;
t299 = m(5) * t242;
t296 = m(6) * (t64 + t63);
t293 = m(6) * (t110 + (-t125 * t202 + (-t120 + t126) * t203) * t193);
t292 = m(6) * (t110 + (-t131 * t202 + (-t120 + t132) * t203) * t193);
t291 = m(6) * ((-t126 + t132) * t203 + (t125 - t131) * t202) * t193;
t48 = -t126 * t119 + t120 * t125;
t290 = m(6) * t48;
t288 = m(6) * t51;
t287 = m(6) * t269;
t285 = m(6) * t63;
t284 = m(6) * t64;
t283 = m(6) * t67;
t282 = -t202 / 0.2e1;
t281 = t202 / 0.2e1;
t280 = -t203 / 0.2e1;
t252 = t202 * t213;
t146 = Icges(6,5) * t252 - Icges(6,6) * t253 - Icges(6,3) * t203;
t149 = Icges(6,6) * t202 + t190 * t203;
t247 = t212 * t149;
t230 = -t146 + t247;
t151 = Icges(6,5) * t202 + t192 * t203;
t137 = t151 * t252;
t188 = Icges(6,5) * t213 - Icges(6,6) * t212;
t250 = t203 * t188;
t147 = Icges(6,3) * t202 + t250;
t231 = t203 * t147 - t137;
t249 = t203 * t213;
t240 = t202 * t147 + t151 * t249;
t185 = Icges(6,4) * t253;
t150 = Icges(6,1) * t252 - Icges(6,5) * t203 - t185;
t241 = -t202 * t146 - t150 * t249;
t148 = Icges(6,4) * t252 - Icges(6,2) * t253 - Icges(6,6) * t203;
t248 = t212 * t148;
t75 = -t203 * t248 - t241;
t76 = -t203 * t247 + t240;
t15 = (t230 * t203 - t240 + t76) * t203 + (t230 * t202 + t231 + t75) * t202;
t74 = -t202 * t247 - t231;
t16 = (t74 - t137 + (t147 + t248) * t203 + t241) * t203 + t240 * t202;
t46 = t202 * t74 - t203 * (-(-t150 * t213 + t248) * t202 - t203 * t146);
t47 = t202 * t76 - t203 * t75;
t2 = (t47 / 0.2e1 - t16 / 0.2e1) * t203 + (t15 / 0.2e1 + t46 / 0.2e1) * t202;
t272 = t2 * qJD(5);
t239 = t191 * t202 + t148;
t238 = -t191 * t203 - t149;
t237 = -Icges(6,2) * t252 + t150 - t185;
t236 = -t189 * t203 + t151;
t228 = t296 / 0.2e1 + t229;
t224 = Icges(6,5) * t212 + Icges(6,6) * t213;
t220 = (-t165 * t203 + t166 * t202) * t193;
t215 = (-t189 + t192) * t213 - t319 * t212;
t219 = t203 * t16 / 0.2e1 + (t15 + t46) * t282 + (t202 * t188 + t215 * t203 + t238 * t212 + t236 * t213) * t281 + (t215 * t202 - t239 * t212 + t237 * t213 - t250 + t47) * t280;
t218 = -t229 + (t281 + t282) * (t213 * t148 + t212 * t150);
t217 = t237 * t212 + t239 * t213;
t216 = -t236 * t212 + t238 * t213;
t194 = -t267 + t268;
t160 = t203 * t224;
t159 = t224 * t202;
t121 = -t165 * t202 - t166 * t203;
t54 = t229 + t283;
t52 = t229 + t284;
t49 = t229 + t285;
t40 = t291 / 0.2e1;
t38 = t292 / 0.2e1;
t35 = t293 / 0.2e1;
t34 = -t287 - t299;
t29 = t288 + t300;
t19 = t279 + t290 + t302;
t18 = -t291 / 0.2e1 + t315;
t17 = t40 + t315;
t12 = -t292 / 0.2e1 + t316;
t11 = t38 + t316;
t10 = -t293 / 0.2e1 + t228;
t9 = t35 + t228;
t8 = t40 - t294 / 0.2e1 + t218;
t7 = t38 - t295 / 0.2e1 + t218;
t6 = t35 - t296 / 0.2e1 + t218;
t3 = t270 + t271;
t1 = [0, 0, 0, 0, m(6) * t121 * qJD(5); 0, qJD(3) * t19 + qJD(4) * t29 + qJD(5) * t49, t19 * qJD(2) + t3 * qJD(4) + t9 * qJD(5) + 0.2e1 * (t279 / 0.2e1 + t48 * t305 + t79 * t306) * qJD(3), t29 * qJD(2) + t3 * qJD(3) + t11 * qJD(5) + (t51 * t305 + t91 * t306) * t308, t49 * qJD(2) + t9 * qJD(3) + t11 * qJD(4) + (((-t119 * t203 - t120 * t202) * t194 + t220) * m(6) + t219) * qJD(5); 0, -t5 * qJD(4) + t10 * qJD(5) + (-t290 / 0.4e1 - t302 / 0.4e1 - t279 / 0.4e1) * t311, qJD(4) * t34 + qJD(5) * t52, -t320 + t34 * qJD(3) + t17 * qJD(5) + (-t242 * t306 - t269 * t305) * t308, t10 * qJD(2) + t52 * qJD(3) + t17 * qJD(4) + (((-t125 * t203 - t126 * t202) * t194 + t220) * m(6) + t219) * qJD(5); 0, t5 * qJD(3) + t12 * qJD(5) + (-t288 / 0.4e1 - t300 / 0.4e1) * t311, t320 + t18 * qJD(5) + 0.4e1 * (t299 / 0.4e1 + t287 / 0.4e1) * qJD(3), qJD(5) * t54, t12 * qJD(2) + t18 * qJD(3) + t54 * qJD(4) + (((-t131 * t203 - t132 * t202) * t194 + t220) * m(6) + t219) * qJD(5); 0, (t218 - t285) * qJD(2) + t6 * qJD(3) + t7 * qJD(4) + t272, t6 * qJD(2) + (t218 - t284) * qJD(3) + t8 * qJD(4) + t272, t7 * qJD(2) + t8 * qJD(3) + (t218 - t283) * qJD(4) + t272, (m(6) * ((t202 * (rSges(6,1) * t252 - t235) + t203 * (rSges(6,1) * t249 + t202 * rSges(6,3) - t187)) * t121 + (t312 + t313) * t194 * t193) + (-t313 * t160 + (t217 * t203 + (t159 + t216) * t202) * t203) * t281 + (-t312 * t159 + (t216 * t202 + (t160 + t217) * t203) * t202) * t280) * qJD(5) + (qJD(2) + qJD(3) + qJD(4)) * t2;];
Cq = t1;
