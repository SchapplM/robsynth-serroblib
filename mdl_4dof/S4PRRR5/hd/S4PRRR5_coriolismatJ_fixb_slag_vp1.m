% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:36
% EndTime: 2019-12-31 16:33:41
% DurationCPUTime: 4.50s
% Computational Cost: add. (23050->346), mult. (31522->527), div. (0->0), fcn. (34065->8), ass. (0->223)
t234 = sin(pkin(7));
t231 = t234 ^ 2;
t235 = cos(pkin(7));
t232 = t235 ^ 2;
t365 = t231 + t232;
t233 = qJ(2) + qJ(3);
t229 = sin(t233);
t230 = cos(t233);
t236 = sin(qJ(4));
t238 = cos(qJ(4));
t273 = rSges(5,1) * t238 - rSges(5,2) * t236;
t179 = -rSges(5,3) * t230 + t273 * t229;
t165 = t179 * t234;
t166 = t179 * t235;
t221 = rSges(4,1) * t229 + rSges(4,2) * t230;
t143 = t365 * t221;
t260 = -Icges(4,5) * t229 - Icges(4,6) * t230;
t237 = sin(qJ(2));
t239 = cos(qJ(2));
t364 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t239 * t237 + (-0.2e1 * t237 ^ 2 + 0.2e1 * t239 ^ 2) * Icges(3,4);
t340 = t234 / 0.2e1;
t339 = -t235 / 0.2e1;
t261 = -Icges(3,5) * t237 - Icges(3,6) * t239;
t215 = t261 * t234;
t216 = t261 * t235;
t306 = t230 * (-Icges(5,5) * t236 - Icges(5,6) * t238) * t229;
t317 = Icges(5,4) * t238;
t262 = -Icges(5,2) * t236 + t317;
t175 = -Icges(5,6) * t230 + t262 * t229;
t293 = -t175 + (-Icges(5,1) * t236 - t317) * t229;
t318 = Icges(5,4) * t236;
t267 = Icges(5,1) * t238 - t318;
t177 = -Icges(5,5) * t230 + t267 * t229;
t292 = t177 + (-Icges(5,2) * t238 - t318) * t229;
t222 = rSges(4,1) * t230 - rSges(4,2) * t229;
t128 = t365 * t222;
t101 = (-t128 + t222) * t143;
t300 = t235 * t238;
t303 = t234 * t236;
t211 = -t230 * t303 - t300;
t301 = t235 * t236;
t302 = t234 * t238;
t212 = t230 * t302 - t301;
t308 = t229 * t234;
t141 = rSges(5,1) * t212 + rSges(5,2) * t211 + rSges(5,3) * t308;
t213 = -t230 * t301 + t302;
t214 = t230 * t300 + t303;
t307 = t229 * t235;
t142 = rSges(5,1) * t214 + rSges(5,2) * t213 + rSges(5,3) * t307;
t224 = pkin(3) * t230 + pkin(6) * t229;
t105 = t234 * t141 + t235 * t142 + t224 * t365;
t223 = pkin(3) * t229 - pkin(6) * t230;
t111 = -t234 * t165 - t235 * t166 - t223 * t365;
t291 = -t179 - t223;
t144 = t291 * t234;
t180 = rSges(5,3) * t229 + t273 * t230;
t290 = -t180 - t224;
t145 = t290 * t234;
t146 = t291 * t235;
t147 = t290 * t235;
t54 = t105 * t111 + t144 * t145 + t146 * t147;
t357 = m(4) * t101 + m(5) * t54;
t335 = pkin(2) * t237;
t275 = t291 - t335;
t129 = t275 * t234;
t131 = t275 * t235;
t334 = pkin(2) * t239;
t294 = t365 * t334;
t89 = t105 + t294;
t50 = t89 * t111 + t129 * t145 + t131 * t147;
t113 = t128 + t294;
t281 = -t221 - t335;
t181 = t281 * t234;
t183 = t281 * t235;
t84 = -t113 * t143 + (-t181 * t234 - t183 * t235) * t222;
t356 = m(4) * t84 + m(5) * t50;
t148 = Icges(5,5) * t211 - Icges(5,6) * t212;
t199 = Icges(5,4) * t211;
t139 = Icges(5,1) * t212 + Icges(5,5) * t308 + t199;
t296 = -Icges(5,2) * t212 + t139 + t199;
t320 = Icges(5,4) * t212;
t137 = Icges(5,2) * t211 + Icges(5,6) * t308 + t320;
t298 = Icges(5,1) * t211 - t137 - t320;
t66 = t148 * t308 + t296 * t211 + t298 * t212;
t149 = Icges(5,5) * t213 - Icges(5,6) * t214;
t200 = Icges(5,4) * t213;
t140 = Icges(5,1) * t214 + Icges(5,5) * t307 + t200;
t295 = -Icges(5,2) * t214 + t140 + t200;
t319 = Icges(5,4) * t214;
t138 = Icges(5,2) * t213 + Icges(5,6) * t307 + t319;
t297 = Icges(5,1) * t213 - t138 - t319;
t67 = t149 * t308 + t295 * t211 + t297 * t212;
t43 = t234 * t67 - t235 * t66;
t68 = t148 * t307 + t296 * t213 + t298 * t214;
t69 = t149 * t307 + t295 * t213 + t297 * t214;
t44 = t234 * t69 - t235 * t68;
t331 = t43 * t339 + t44 * t340;
t259 = Icges(5,5) * t238 - Icges(5,6) * t236;
t173 = -Icges(5,3) * t230 + t259 * t229;
t354 = t365 * t260;
t351 = 2 * qJD(2);
t350 = m(4) / 0.2e1;
t349 = m(5) / 0.2e1;
t256 = t141 * t235 - t142 * t234;
t112 = t256 * t229;
t122 = t141 * t230 + t179 * t308;
t123 = -t142 * t230 - t179 * t307;
t284 = t112 * t111 + t122 * t147 + t123 * t145;
t103 = (t180 * t234 - t141) * t229;
t104 = (-t180 * t235 + t142) * t229;
t92 = t256 * t230 + (-t165 * t235 + t166 * t234) * t229;
t285 = t103 * t131 + t104 * t129 + t92 * t89;
t347 = m(5) * (t284 + t285);
t21 = t103 * t146 + t104 * t144 + t105 * t92 + t284;
t346 = m(5) * t21;
t210 = (-rSges(5,1) * t236 - rSges(5,2) * t238) * t229;
t154 = rSges(5,1) * t211 - rSges(5,2) * t212;
t155 = rSges(5,1) * t213 - rSges(5,2) * t214;
t124 = t154 * t234 + t155 * t235;
t76 = t89 * t124;
t85 = t105 * t124;
t345 = m(5) * (t76 + t85 + ((-t131 - t146) * t235 + (-t129 - t144) * t234) * t210);
t344 = m(5) * (t103 * t122 + t104 * t123 + t112 * t92);
t341 = -t230 / 0.2e1;
t338 = t235 / 0.2e1;
t77 = 0.2e1 * (t92 / 0.4e1 - t124 / 0.4e1) * m(5);
t299 = t77 * qJD(1);
t109 = t173 * t308 + t175 * t211 + t177 * t212;
t176 = Icges(5,6) * t229 + t262 * t230;
t178 = Icges(5,5) * t229 + t267 * t230;
t255 = -t175 * t236 + t177 * t238;
t252 = (Icges(5,3) * t229 + t259 * t230 - t255) * t230;
t311 = t173 * t230;
t136 = Icges(5,5) * t214 + Icges(5,6) * t213 + Icges(5,3) * t307;
t94 = t136 * t308 + t138 * t211 + t140 * t212;
t326 = t94 * t235;
t161 = t175 * t234;
t163 = t177 * t234;
t258 = -t137 * t236 + t139 * t238;
t254 = -t173 * t234 - t258;
t135 = Icges(5,5) * t212 + Icges(5,6) * t211 + Icges(5,3) * t308;
t313 = t135 * t230;
t242 = t254 * t229 + t313;
t62 = -t161 * t211 - t163 * t212 + t242 * t234;
t162 = t175 * t235;
t164 = t177 * t235;
t257 = -t138 * t236 + t140 * t238;
t253 = -t173 * t235 - t257;
t312 = t136 * t230;
t241 = t253 * t229 + t312;
t63 = -t162 * t211 - t164 * t212 + t241 * t234;
t93 = t135 * t308 + t137 * t211 + t139 * t212;
t15 = (-t211 * t176 - t212 * t178 + t326 + (t93 - t311) * t234) * t230 + (t63 * t235 + t109 + (t62 - t252) * t234) * t229;
t110 = t173 * t307 + t175 * t213 + t177 * t214;
t95 = t135 * t307 + t137 * t213 + t139 * t214;
t325 = t95 * t234;
t64 = -t161 * t213 - t163 * t214 + t242 * t235;
t65 = -t162 * t213 - t164 * t214 + t241 * t235;
t96 = t136 * t307 + t138 * t213 + t140 * t214;
t16 = (-t213 * t176 - t214 * t178 + t325 + (t96 - t311) * t235) * t230 + (t64 * t234 + t110 + (t65 - t252) * t235) * t229;
t117 = t255 * t229 - t311;
t100 = t257 * t229 - t312;
t99 = t258 * t229 - t313;
t272 = t100 * t235 + t99 * t234;
t79 = -t254 * t230 + (t161 * t236 - t163 * t238 + t135) * t229;
t80 = -t253 * t230 + (t162 * t236 - t164 * t238 + t136) * t229;
t22 = (t252 + t272) * t230 + (t80 * t235 + t79 * t234 - (-t176 * t236 + t178 * t238 + t173) * t230 + t117) * t229;
t48 = -t109 * t230 + (t234 * t93 + t326) * t229;
t49 = -t110 * t230 + (t235 * t96 + t325) * t229;
t53 = -t117 * t230 + t272 * t229;
t3 = t344 + (t49 * t338 + t48 * t340 - t22 / 0.2e1) * t230 + (t16 * t338 + t15 * t340 + t53 / 0.2e1) * t229;
t337 = t3 * qJD(4) + t299;
t78 = (t124 + t92) * t349;
t98 = -m(4) * t143 + m(5) * t111;
t330 = t98 * qJD(3) + t78 * qJD(4);
t305 = t230 * t234;
t304 = t230 * t235;
t287 = qJD(2) + qJD(3);
t286 = t345 / 0.2e1 + t331;
t283 = t308 / 0.2e1;
t282 = t307 / 0.2e1;
t280 = -t222 - t334;
t201 = t260 * t234;
t202 = t260 * t235;
t35 = t234 * t63 - t235 * t62;
t36 = t234 * t65 - t235 * t64;
t279 = (t36 + t231 * t202 + (-t234 * t201 + t354) * t235) * t340 + (t35 + t232 * t201 + (-t235 * t202 + t354) * t234) * t339;
t276 = t365 * t335;
t274 = t290 - t334;
t225 = rSges(3,1) * t237 + rSges(3,2) * t239;
t251 = t112 * t124 + (-t122 * t235 - t123 * t234) * t210;
t250 = t364 * t234 + t216;
t249 = -t364 * t235 + t215;
t248 = t15 * t339 + t16 * t340 + t36 * t282 + t35 * t283 + (t234 * t80 - t235 * t79) * t341 + (t234 * t94 - t235 * t93) * t305 / 0.2e1 + (t234 * t96 - t235 * t95) * t304 / 0.2e1 + t229 * (t100 * t234 - t235 * t99) / 0.2e1 - t331;
t28 = -(t292 * t211 + t293 * t212) * t230 + (t67 * t235 + (t66 - t306) * t234) * t229;
t29 = -(t292 * t213 + t293 * t214) * t230 + (t68 * t234 + (t69 - t306) * t235) * t229;
t82 = -t148 * t230 + (-t296 * t236 + t298 * t238) * t229;
t83 = -t149 * t230 + (-t295 * t236 + t297 * t238) * t229;
t243 = -t15 * t308 / 0.2e1 - t16 * t307 / 0.2e1 + t230 * t22 / 0.2e1 + t29 * t340 + t28 * t339 - t344 + t43 * t283 + t44 * t282 - t48 * t305 / 0.2e1 - t49 * t304 / 0.2e1 + (t234 * t83 - t235 * t82) * t341 - t229 * t53 / 0.2e1;
t184 = t280 * t235;
t182 = t280 * t234;
t156 = t365 * t225;
t132 = t274 * t235;
t130 = t274 * t234;
t127 = -t276 - t143;
t126 = -t155 * t230 - t210 * t307;
t125 = t154 * t230 + t210 * t308;
t118 = (t154 * t235 - t155 * t234) * t229;
t106 = -t276 + t111;
t74 = t77 * qJD(4);
t59 = t85 + (-t144 * t234 - t146 * t235) * t210;
t55 = t76 + (-t129 * t234 - t131 * t235) * t210;
t19 = t346 / 0.2e1;
t17 = t347 / 0.2e1;
t10 = m(5) * t59 + t331;
t9 = m(5) * t55 + t331;
t8 = t279 + t357;
t7 = t8 * qJD(3);
t6 = t279 + t356;
t5 = t19 - t347 / 0.2e1 + t286;
t4 = t17 - t346 / 0.2e1 + t286;
t1 = t17 + t19 - t345 / 0.2e1 + t248;
t2 = [0, (-m(3) * t156 / 0.2e1 + t127 * t350 + t106 * t349) * t351 + t330, qJD(2) * t98 + t330, m(5) * t118 * qJD(4) + t287 * t78; -t74, ((t232 * t215 + (t249 * t234 + (-t216 + t250) * t235) * t234) * t339 + m(4) * (t113 * t127 + t181 * t182 + t183 * t184) + m(5) * (t106 * t89 + t129 * t130 + t131 * t132) + (t231 * t216 + (t250 * t235 + (-t215 + t249) * t234) * t235) * t340 + t279 + m(3) * (-t156 + t225) * t365 * (rSges(3,1) * t239 - rSges(3,2) * t237)) * qJD(2) + t6 * qJD(3) + t9 * qJD(4), t6 * qJD(2) + t4 * qJD(4) + (t279 + 0.2e1 * (t84 + t101) * t350 + 0.2e1 * (t54 + t50) * t349 - t357) * qJD(3), -t299 + t9 * qJD(2) + t4 * qJD(3) + (t243 + m(5) * (t118 * t89 + t125 * t131 + t126 * t129 + t251)) * qJD(4); -t74, t7 + t5 * qJD(4) + ((t127 * t128 + (-t182 * t234 - t184 * t235) * t221 + t84) * t350 + (t105 * t106 + t130 * t144 + t132 * t146 + t50) * t349) * t351 + (t279 - t356) * qJD(2), qJD(2) * t8 + qJD(4) * t10 + t7, -t299 + t5 * qJD(2) + t10 * qJD(3) + (t243 + m(5) * (t105 * t118 + t125 * t146 + t126 * t144 + t251)) * qJD(4); t287 * t77, ((t106 * t112 + t122 * t132 + t123 * t130 + t285 - t55) * m(5) + t248) * qJD(2) + t1 * qJD(3) + t337, t1 * qJD(2) + ((t21 - t59) * m(5) + t248) * qJD(3) + t337, t287 * t3 + (m(5) * (t112 * t118 + t122 * t125 + t123 * t126) - t230 ^ 2 * t306 / 0.2e1 + (t29 * t338 + t28 * t340 + (t83 * t235 + t82 * t234 - (-t236 * t292 + t238 * t293) * t230) * t341) * t229) * qJD(4);];
Cq = t2;
