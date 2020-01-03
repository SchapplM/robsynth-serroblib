% Calculate vector of inverse dynamics joint torques for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR9_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR9_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:44
% DurationCPUTime: 5.47s
% Computational Cost: add. (9719->357), mult. (10522->461), div. (0->0), fcn. (10846->6), ass. (0->185)
t253 = pkin(8) + qJ(2);
t236 = sin(t253);
t237 = cos(t253);
t304 = sin(qJ(4));
t305 = cos(qJ(4));
t125 = -t236 * t305 + t237 * t304;
t172 = qJDD(2) - qJDD(4);
t173 = qJD(2) - qJD(4);
t170 = t237 * pkin(3);
t176 = qJD(2) ^ 2;
t163 = qJD(3) * t236;
t217 = qJD(2) * t236;
t260 = t237 * pkin(2) + t236 * qJ(3);
t218 = qJD(2) * t237;
t263 = qJ(3) * t218 + t163;
t195 = qJDD(2) * t260 - qJDD(3) * t237 + (-pkin(2) * t217 + t163 + t263) * qJD(2);
t230 = t236 * pkin(3);
t179 = qJDD(2) * t170 - t176 * t230 + t195;
t174 = sin(qJ(5));
t175 = cos(qJ(5));
t228 = rSges(6,1) * t174 + rSges(6,2) * t175;
t153 = -rSges(6,1) * t175 + rSges(6,2) * t174;
t141 = t153 * qJD(5);
t257 = qJD(5) * t141;
t124 = -t236 * t304 - t237 * t305;
t273 = t124 * t175;
t274 = t124 * t174;
t75 = -rSges(6,1) * t273 + rSges(6,2) * t274 + t125 * rSges(6,3);
t284 = -t124 * pkin(4) + pkin(7) * t125 + t75;
t255 = qJD(5) * t174;
t99 = t173 * t125;
t203 = t124 * t255 + t175 * t99;
t254 = qJD(5) * t175;
t204 = t124 * t254 - t174 * t99;
t98 = t173 * t124;
t30 = rSges(6,1) * t203 + t204 * rSges(6,2) + t98 * rSges(6,3);
t314 = t99 * pkin(4) + t98 * pkin(7) + t30;
t77 = qJD(5) * t98 + qJDD(5) * t125;
t9 = -t125 * t257 + t284 * t172 + t173 * t314 + t77 * t228 + t179;
t344 = g(2) - t9;
t166 = t237 * qJ(3);
t231 = t236 * pkin(2);
t131 = t231 - t166;
t164 = qJD(3) * t237;
t198 = -qJDD(2) * t131 + qJDD(3) * t236 + (-qJD(2) * t260 + 0.2e1 * t164) * qJD(2);
t178 = (-qJDD(2) * t236 - t176 * t237) * pkin(3) + t198;
t201 = t125 * t255 - t175 * t98;
t202 = t125 * t254 + t174 * t98;
t29 = rSges(6,1) * t201 + t202 * rSges(6,2) + t99 * rSges(6,3);
t324 = t98 * pkin(4) - t99 * pkin(7) - t29;
t270 = t125 * t175;
t271 = t125 * t174;
t74 = -rSges(6,1) * t270 + rSges(6,2) * t271 - t124 * rSges(6,3);
t39 = -t125 * pkin(4) - t124 * pkin(7) + t74;
t78 = qJD(5) * t99 - qJDD(5) * t124;
t8 = -t124 * t257 - t39 * t172 + t173 * t324 - t78 * t228 + t178;
t343 = t8 - g(1);
t266 = t125 * rSges(5,1) - t124 * rSges(5,2);
t50 = -t98 * rSges(5,1) - t99 * rSges(5,2);
t19 = t172 * t266 - t173 * t50 + t178;
t342 = t19 - g(1);
t341 = t173 * t39;
t282 = Icges(6,4) * t175;
t214 = -Icges(6,2) * t174 + t282;
t64 = -Icges(6,6) * t125 + t124 * t214;
t287 = t174 * t64;
t283 = Icges(6,4) * t174;
t216 = Icges(6,1) * t175 - t283;
t68 = -Icges(6,5) * t125 + t124 * t216;
t219 = t175 * t68 - t287;
t65 = Icges(6,6) * t124 + t125 * t214;
t288 = t174 * t65;
t69 = Icges(6,5) * t124 + t125 * t216;
t339 = -t175 * t69 + t288;
t212 = Icges(6,5) * t175 - Icges(6,6) * t174;
t61 = Icges(6,3) * t124 + t125 * t212;
t298 = -t124 * t61 - t270 * t69;
t60 = -Icges(6,3) * t125 + t124 * t212;
t338 = -(t60 - t288) * t125 + t298;
t256 = qJD(5) * t228;
t337 = t125 * t256 + t173 * t284;
t220 = -t174 * t68 - t175 * t64;
t222 = -t174 * t69 - t175 * t65;
t258 = qJD(5) * t125;
t259 = qJD(5) * t124;
t306 = -t173 / 0.2e1;
t336 = ((-t220 * t124 + t222 * t125) * qJD(5) + t220 * t259 - t222 * t258) * t306;
t100 = rSges(5,1) * t124 + rSges(5,2) * t125;
t51 = t99 * rSges(5,1) - t98 * rSges(5,2);
t20 = -t100 * t172 + t173 * t51 + t179;
t335 = t20 - g(2);
t334 = t125 * t60;
t333 = t125 * t61;
t329 = t61 + t287;
t328 = t173 * t266;
t211 = Icges(6,5) * t174 + Icges(6,6) * t175;
t83 = t211 * t124;
t82 = t211 * t125;
t323 = qJD(2) * t131 - t163 + t263;
t247 = t170 + t260;
t191 = t247 * qJD(2) - t164;
t213 = Icges(6,2) * t175 + t283;
t215 = Icges(6,1) * t174 + t282;
t210 = -t174 * t213 + t175 * t215;
t44 = t125 * t210 + t83;
t42 = t44 * t173;
t45 = t124 * t210 - t82;
t43 = t45 * t173;
t261 = t237 * rSges(4,1) + t236 * rSges(4,3);
t317 = t260 + t261;
t313 = t124 * t256 - t341;
t189 = -t231 - t230;
t311 = pkin(3) * t217 + t189 * qJD(2) + t323;
t292 = t213 * t125 - t69;
t294 = -t215 * t125 - t65;
t310 = t174 * t292 + t175 * t294;
t307 = -m(4) - m(5);
t297 = t124 * t60 + t270 * t68;
t296 = t273 * t69 - t333;
t295 = t273 * t68 - t334;
t293 = -t215 * t124 - t64;
t291 = t213 * t124 - t68;
t269 = t212 * t173;
t265 = -t213 + t216;
t264 = -t214 - t215;
t248 = m(2) + m(3) - t307;
t240 = t259 / 0.2e1;
t239 = -t258 / 0.2e1;
t238 = t258 / 0.2e1;
t235 = -t237 / 0.2e1;
t234 = t236 / 0.2e1;
t26 = Icges(6,4) * t203 + Icges(6,2) * t204 + Icges(6,6) * t98;
t28 = Icges(6,1) * t203 + Icges(6,4) * t204 + Icges(6,5) * t98;
t183 = qJD(5) * t220 + t174 * t26 - t175 * t28;
t25 = Icges(6,4) * t201 + Icges(6,2) * t202 + Icges(6,6) * t99;
t27 = Icges(6,1) * t201 + Icges(6,4) * t202 + Icges(6,5) * t99;
t184 = qJD(5) * t222 + t174 * t25 - t175 * t27;
t23 = Icges(6,5) * t201 + Icges(6,6) * t202 + Icges(6,3) * t99;
t24 = Icges(6,5) * t203 + Icges(6,6) * t204 + Icges(6,3) * t98;
t233 = -(-t124 * t23 + t125 * t184 - t339 * t98 - t61 * t99) * t124 + t125 * (-t124 * t24 + t125 * t183 + t219 * t98 - t60 * t99);
t232 = -t124 * (t124 * t184 + t125 * t23 + t339 * t99 - t61 * t98) + t125 * (t124 * t183 + t125 * t24 - t219 * t99 - t60 * t98);
t229 = t236 * rSges(4,1);
t15 = -t271 * t65 - t298;
t16 = -t271 * t64 + t297;
t227 = -t124 * t15 + t125 * t16;
t17 = -t274 * t65 + t296;
t18 = -t274 * t64 + t295;
t226 = -t124 * t17 + t125 * t18;
t225 = t124 * t30 + t125 * t29;
t192 = t163 + (-t230 - t131) * qJD(2);
t31 = t192 + t313;
t32 = t191 + t337;
t224 = -t124 * t31 - t125 * t32;
t223 = t124 * t75 + t125 * t74;
t107 = -t174 * t215 - t175 * t213;
t197 = t174 * t291 + t175 * t293;
t190 = (t174 * t264 + t175 * t265) * t173;
t188 = -t229 - t231;
t136 = rSges(3,1) * t237 - rSges(3,2) * t236;
t133 = rSges(3,1) * t236 + rSges(3,2) * t237;
t185 = t166 + t189;
t139 = t214 * qJD(5);
t140 = t216 * qJD(5);
t182 = qJD(5) * t107 - t139 * t174 + t140 * t175;
t11 = -t339 * qJD(5) - t174 * t27 - t175 * t25;
t12 = qJD(5) * t219 - t174 * t28 - t175 * t26;
t138 = t212 * qJD(5);
t13 = t124 * t138 + t125 * t182 + t210 * t98 - t211 * t99;
t14 = t124 * t182 - t125 * t138 - t210 * t99 - t211 * t98;
t6 = qJD(5) * t227 + t42;
t7 = qJD(5) * t226 + t43;
t177 = t173 * (-t140 * t174 + t213 * t255 + (-qJD(5) * t215 - t139) * t175) + t238 * t6 - (-t220 + t45) * t77 / 0.2e1 - (-t222 + t44) * t78 / 0.2e1 + (t12 + t14) * t239 + (-Icges(5,3) + t107) * t172 + (t11 + t13 + t7) * t240;
t168 = t237 * rSges(4,3);
t161 = rSges(4,3) * t218;
t132 = t229 - t168;
t89 = t228 * t124;
t88 = t228 * t125;
t59 = -t173 * t100 + t191;
t58 = t192 + t328;
t41 = qJDD(2) * t261 + qJD(2) * (-rSges(4,1) * t217 + t161) + t195;
t40 = -qJDD(2) * t132 - t176 * t261 + t198;
t33 = qJD(5) * t223 + qJD(1);
t10 = (-t124 * t65 - t125 * t64) * t174 + t296 + t297;
t5 = qJD(5) * t225 + t74 * t77 - t75 * t78 + qJDD(1);
t1 = [m(6) * t5 + t248 * qJDD(1) + (-m(6) - t248) * g(3); -m(3) * (-g(1) * t133 + g(2) * t136) - t177 + t336 + (t43 + ((t339 * t125 + t15 + t295) * t125 + (-t329 * t125 - t10 + t16) * t124) * qJD(5)) * t240 + (-t42 + ((t17 + (-t219 + t61) * t125) * t125 + (t329 * t124 + t18 - t295 - t334) * t124) * qJD(5)) * t239 + (m(3) * (t133 ^ 2 + t136 ^ 2) + Icges(3,3) + Icges(4,2)) * qJDD(2) + (-t344 * (t247 + t284) + (-t191 + t324) * t31 + (-t228 * t259 + t31 + t311 + t314 + t341) * t32 + t343 * (t185 - t39)) * m(6) + (t335 * (-t100 + t247) + (-t191 - t50) * t58 + (t311 + t51 + t58 - t328) * t59 + t342 * (t185 + t266)) * m(5) + ((t161 + (t132 + t188) * qJD(2) + t323) * (t317 * qJD(2) - t164) + (t41 - g(2)) * t317 + (t40 - g(1)) * (t166 + t168 + t188)) * m(4); (-m(6) + t307) * (g(1) * t236 - g(2) * t237) + 0.2e1 * (t234 * t8 + t235 * t9) * m(6) + 0.2e1 * (t19 * t234 + t20 * t235) * m(5) + 0.2e1 * (t234 * t40 + t235 * t41) * m(4); t177 + t336 + (t42 + ((t10 - t17) * t125 + (t219 * t124 - t18 + t338) * t124) * qJD(5)) * t239 + (-t43 + ((-t15 - t338) * t125 + (-t16 + (-t339 + t60) * t124 - t333) * t124) * qJD(5)) * t240 + (t343 * t39 + (t313 - t314) * t32 + (-t324 - t337) * t31 + t344 * t284) * m(6) + (t50 * t58 - t51 * t59 - (-t173 * t59 + t342) * t266 + (t173 * t58 + t335) * t100) * m(5); t98 * t7 / 0.2e1 + t125 * (qJD(5) * t232 + t14 * t173 + t17 * t78 + t172 * t45 + t18 * t77) / 0.2e1 + t77 * t226 / 0.2e1 + (t17 * t99 + t18 * t98 + t232) * t238 + t99 * t6 / 0.2e1 - t124 * (qJD(5) * t233 + t13 * t173 + t15 * t78 + t16 * t77 + t172 * t44) / 0.2e1 + t78 * t227 / 0.2e1 - (t15 * t99 + t16 * t98 + t233) * t259 / 0.2e1 + t172 * (t124 * t222 - t125 * t220) / 0.2e1 + t173 * (-t11 * t124 + t12 * t125 - t220 * t98 - t222 * t99) / 0.2e1 + ((t83 * t258 - t269) * t125 + (t190 + (-t310 * t124 + (-t82 + t197) * t125) * qJD(5)) * t124) * t239 + ((t82 * t259 + t269) * t124 + (t190 + (t197 * t125 + (-t310 - t83) * t124) * qJD(5)) * t125) * t240 + ((t174 * t265 - t175 * t264) * t173 + ((t124 * t292 - t125 * t291) * t175 + (-t124 * t294 + t125 * t293) * t174) * qJD(5)) * t306 + (t5 * t223 + t33 * (t74 * t98 - t75 * t99 + t225) + t224 * t141 - (-t8 * t124 - t125 * t9 + t31 * t99 - t32 * t98) * t228 - (-t31 * t88 + t32 * t89) * t173 - (t33 * (t124 * t89 + t125 * t88) + t224 * t153) * qJD(5) - g(1) * t89 - g(2) * t88 - g(3) * t153) * m(6);];
tau = t1;
