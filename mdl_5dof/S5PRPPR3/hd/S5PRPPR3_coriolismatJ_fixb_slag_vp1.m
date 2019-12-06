% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:19
% EndTime: 2019-12-05 15:26:33
% DurationCPUTime: 7.53s
% Computational Cost: add. (9898->357), mult. (13516->573), div. (0->0), fcn. (14296->8), ass. (0->226)
t194 = qJ(2) + pkin(8);
t190 = sin(t194);
t191 = cos(t194);
t350 = -Icges(4,5) * t190 + (Icges(5,5) - Icges(4,6)) * t191;
t195 = sin(pkin(7));
t192 = t195 ^ 2;
t196 = cos(pkin(7));
t193 = t196 ^ 2;
t259 = t192 + t193;
t348 = Icges(5,4) * t190 + t350;
t198 = sin(qJ(5));
t200 = cos(qJ(5));
t240 = rSges(6,1) * t198 + rSges(6,2) * t200;
t207 = -rSges(6,3) * t190 + t191 * t240;
t113 = t207 * t195;
t114 = t207 * t196;
t199 = sin(qJ(2));
t201 = cos(qJ(2));
t343 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t201 * t199 + (-0.2e1 * t199 ^ 2 + 0.2e1 * t201 ^ 2) * Icges(3,4);
t279 = t190 * t191;
t223 = -Icges(3,5) * t199 - Icges(3,6) * t201;
t168 = t223 * t195;
t169 = t223 * t196;
t342 = -Icges(5,2) + Icges(5,3);
t275 = t191 * t196;
t276 = t191 * t195;
t277 = t190 * t196;
t278 = t190 * t195;
t333 = 2 * Icges(5,6);
t336 = (-t196 * (-Icges(5,4) * t196 + t342 * t276 + t278 * t333) + (Icges(5,4) * t195 + t342 * t275 + t277 * t333) * t195) * t190 + t350 * t259;
t335 = t348 * t195 + t168;
t334 = t348 * t196 + t169;
t281 = (-Icges(6,5) * t200 + Icges(6,6) * t198) * t279;
t257 = m(5) / 0.4e1 + m(6) / 0.4e1;
t260 = t259 * t279;
t332 = t257 * (t260 - t279);
t331 = t259 * t190;
t288 = Icges(6,4) * t198;
t224 = Icges(6,2) * t200 + t288;
t205 = -Icges(6,6) * t190 + t191 * t224;
t266 = -t205 - (-Icges(6,1) * t200 + t288) * t191;
t287 = Icges(6,4) * t200;
t230 = Icges(6,1) * t198 + t287;
t206 = -Icges(6,5) * t190 + t191 * t230;
t265 = -t206 + (Icges(6,2) * t198 - t287) * t191;
t221 = Icges(6,5) * t198 + Icges(6,6) * t200;
t204 = -Icges(6,3) * t190 + t191 * t221;
t325 = 2 * qJD(2);
t324 = m(4) / 0.2e1;
t323 = m(5) / 0.2e1;
t322 = m(6) / 0.2e1;
t239 = rSges(5,2) * t190 + rSges(5,3) * t191;
t174 = pkin(3) * t190 - qJ(4) * t191;
t312 = pkin(2) * t199;
t254 = -t174 - t312;
t243 = t239 + t254;
t93 = t243 * t195;
t95 = t243 * t196;
t305 = t95 * t275 + t93 * t276;
t178 = -rSges(5,2) * t191 + rSges(5,3) * t190;
t177 = pkin(3) * t191 + qJ(4) * t190;
t311 = pkin(2) * t201;
t267 = t259 * t311;
t245 = t259 * t177 + t267;
t54 = t259 * t178 + t245;
t321 = m(5) * (t331 * t54 + t305);
t272 = t196 * t198;
t273 = t195 * t200;
t165 = t190 * t273 + t272;
t271 = t196 * t200;
t274 = t195 * t198;
t166 = -t190 * t274 + t271;
t92 = -rSges(6,1) * t166 + rSges(6,2) * t165 + rSges(6,3) * t276;
t67 = -t190 * t92 - t207 * t276;
t163 = t190 * t271 - t274;
t164 = t190 * t272 + t273;
t91 = rSges(6,1) * t164 + rSges(6,2) * t163 + rSges(6,3) * t275;
t68 = t190 * t91 + t207 * t275;
t307 = t67 * t275 + t68 * t276;
t237 = t195 * t91 - t196 * t92;
t45 = (t113 * t196 - t114 * t195) * t191 + t237 * t190;
t123 = rSges(6,3) * t191 + t190 * t240;
t52 = (t195 * t123 - t92) * t191;
t53 = (-t196 * t123 + t91) * t191;
t60 = t237 * t191;
t320 = m(6) * (-t191 * t45 + (t195 * t53 + t196 * t52 - t60) * t190 + t307);
t212 = -pkin(6) * t190 + t207 + t254;
t75 = t212 * t195;
t77 = t212 * t196;
t306 = t77 * t275 + t75 * t276;
t310 = pkin(6) * t191;
t43 = t195 * t92 + t196 * t91 + t259 * t310 + t245;
t319 = m(6) * (t331 * t43 + t306);
t318 = m(6) * (-t331 * t60 + t307);
t162 = (-rSges(6,1) * t200 + rSges(6,2) * t198) * t191;
t103 = rSges(6,1) * t163 - rSges(6,2) * t164;
t104 = rSges(6,1) * t165 + rSges(6,2) * t166;
t69 = t103 * t196 + t104 * t195;
t317 = m(6) * (-t162 * t331 - t191 * t69);
t316 = t190 / 0.2e1;
t315 = t195 / 0.2e1;
t314 = -t196 / 0.2e1;
t313 = t196 / 0.2e1;
t157 = Icges(6,4) * t163;
t89 = Icges(6,1) * t164 + Icges(6,5) * t275 + t157;
t304 = -Icges(6,2) * t164 + t157 + t89;
t303 = m(6) * qJD(5);
t85 = Icges(6,5) * t164 + Icges(6,6) * t163 + Icges(6,3) * t275;
t301 = t190 * t85;
t86 = -Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t276;
t300 = t190 * t86;
t289 = Icges(6,4) * t166;
t88 = Icges(6,2) * t165 + Icges(6,6) * t276 - t289;
t158 = Icges(6,4) * t165;
t90 = -Icges(6,1) * t166 + Icges(6,5) * t276 + t158;
t47 = t163 * t88 + t164 * t90 + t275 * t86;
t299 = t47 * t195;
t290 = Icges(6,4) * t164;
t87 = Icges(6,2) * t163 + Icges(6,6) * t275 + t290;
t48 = t165 * t87 - t166 * t89 + t276 * t85;
t298 = t48 * t196;
t297 = Icges(6,2) * t166 + t158 + t90;
t296 = Icges(6,1) * t163 - t290 - t87;
t295 = Icges(6,1) * t165 + t289 - t88;
t282 = t204 * t190;
t26 = t195 * t52 - t196 * t53;
t270 = t26 * qJD(3);
t36 = 0.2e1 * (t45 / 0.4e1 - t69 / 0.4e1) * m(6);
t269 = t36 * qJD(1);
t241 = 0.2e1 * t257 * t331;
t256 = t322 + t323;
t81 = -t190 * t256 + t241;
t268 = t81 * qJD(1);
t258 = qJD(5) * t191;
t255 = t26 * t322;
t176 = rSges(4,1) * t190 + rSges(4,2) * t191;
t253 = -t176 - t312;
t252 = -t177 - t311;
t179 = rSges(4,1) * t191 - rSges(4,2) * t190;
t251 = -t179 - t311;
t244 = t259 * t312;
t242 = -t178 + t252;
t185 = t199 * rSges(3,1) + rSges(3,2) * t201;
t236 = -t198 * t89 - t200 * t87;
t50 = t191 * t236 + t301;
t235 = -t198 * t90 - t200 * t88;
t51 = t191 * t235 + t300;
t238 = t51 * t195 + t50 * t196;
t218 = t198 * t206 + t200 * t205;
t217 = t204 * t195 - t235;
t216 = t204 * t196 - t236;
t97 = Icges(6,5) * t163 - Icges(6,6) * t164;
t32 = t163 * t304 + t164 * t296 + t275 * t97;
t98 = Icges(6,5) * t165 + Icges(6,6) * t166;
t33 = t163 * t297 + t164 * t295 + t275 * t98;
t17 = t195 * t32 - t196 * t33;
t34 = t165 * t304 - t166 * t296 + t276 * t97;
t35 = t165 * t297 - t166 * t295 + t276 * t98;
t18 = t195 * t34 - t196 * t35;
t215 = t17 * t315 + t18 * t314;
t79 = -t259 * t176 - t244;
t214 = t79 * t179;
t213 = -t123 + t252 - t310;
t211 = -t259 * t174 - t244;
t210 = (Icges(6,3) * t191 + t190 * t221 - t218) * t190;
t209 = t343 * t195 + t169;
t208 = -t343 * t196 + t168;
t203 = t191 * t216 - t301;
t202 = t191 * t217 - t300;
t126 = t251 * t196;
t125 = t251 * t195;
t121 = Icges(6,5) * t191 + t190 * t230;
t119 = Icges(6,6) * t191 + t190 * t224;
t112 = t206 * t196;
t111 = t206 * t195;
t110 = t205 * t196;
t109 = t205 * t195;
t105 = t259 * t185;
t96 = t242 * t196;
t94 = t242 * t195;
t80 = t241 + (m(5) + m(6)) * t316;
t78 = t213 * t196;
t76 = t213 * t195;
t72 = t103 * t190 - t162 * t275;
t71 = -t104 * t190 + t162 * t276;
t70 = 0.4e1 * t332;
t65 = (-t103 * t195 + t104 * t196) * t191;
t62 = t191 * t218 - t282;
t61 = t259 * t239 + t211;
t59 = -t165 * t205 + t166 * t206 - t204 * t276;
t58 = -t163 * t205 - t164 * t206 - t204 * t275;
t56 = t317 / 0.2e1;
t55 = -pkin(6) * t331 + t113 * t195 + t114 * t196 + t211;
t49 = t165 * t88 - t166 * t90 + t276 * t86;
t46 = t163 * t87 + t164 * t89 + t275 * t85;
t42 = t190 * t98 + (-t198 * t295 - t200 * t297) * t191;
t41 = t190 * t97 + (-t198 * t296 - t200 * t304) * t191;
t39 = (-t109 * t200 - t111 * t198 + t86) * t191 + t217 * t190;
t38 = (-t110 * t200 - t112 * t198 + t85) * t191 + t216 * t190;
t37 = (t45 + t69) * t322;
t31 = t109 * t163 + t111 * t164 + t196 * t202;
t30 = t110 * t163 + t112 * t164 + t196 * t203;
t29 = t109 * t165 - t111 * t166 + t195 * t202;
t28 = t110 * t165 - t112 * t166 + t195 * t203;
t25 = t318 / 0.2e1;
t24 = qJD(2) * t255;
t22 = t190 * t62 + t191 * t238;
t21 = t43 * t69 + (-t195 * t75 - t196 * t77) * t162;
t20 = t190 * t59 + (t195 * t49 + t298) * t191;
t19 = t190 * t58 + (t196 * t46 + t299) * t191;
t16 = t195 * t30 - t196 * t31;
t15 = t195 * t28 - t196 * t29;
t14 = t319 + t321;
t13 = (t165 * t265 + t166 * t266) * t190 + (t34 * t196 + (t35 + t281) * t195) * t191;
t12 = (t163 * t265 - t164 * t266) * t190 + (t33 * t195 + (t32 + t281) * t196) * t191;
t11 = -t45 * t60 + t52 * t67 + t53 * t68;
t9 = t320 / 0.2e1;
t8 = (t195 * t39 + t196 * t38 + t62) * t191 + (t210 + (-t119 * t200 - t121 * t198 - t204) * t191 - t238) * t190;
t7 = t25 + t9 - t317 / 0.2e1;
t6 = t56 + t25 - t320 / 0.2e1;
t5 = t56 + t9 - t318 / 0.2e1;
t4 = (t119 * t163 + t121 * t164 - t299 + (-t46 + t282) * t196) * t190 + (t31 * t195 + t58 + (t30 + t210) * t196) * t191;
t3 = (t119 * t165 - t121 * t166 - t298 + (-t49 + t282) * t195) * t190 + (t28 * t196 + t59 + (t29 + t210) * t195) * t191;
t2 = m(6) * t21 + t215;
t1 = m(6) * t11 + (t4 * t313 + t3 * t315 + t22 / 0.2e1) * t191 + (t19 * t314 - t195 * t20 / 0.2e1 + t8 / 0.2e1) * t190;
t10 = [0, t80 * qJD(4) + t37 * qJD(5) + (-m(3) * t105 / 0.2e1 + t79 * t324 + t61 * t323 + t55 * t322) * t325, 0, t80 * qJD(2), t37 * qJD(2) + t303 * t65; qJD(4) * t81 - qJD(5) * t36, t14 * qJD(4) + t2 * qJD(5) + (m(6) * (t43 * t55 + t75 * t76 + t77 * t78) + m(5) * (t54 * t61 + t93 * t94 + t95 * t96) + m(4) * (t267 * t79 + (t126 * t253 + t196 * t214) * t196 + (t125 * t253 + t195 * t214) * t195) + m(3) * (-t105 + t185) * t259 * (rSges(3,1) * t201 - t199 * rSges(3,2)) + (t16 + (t209 * t196 + (t208 - t335) * t195 + t336) * t196 + t334 * t192) * t315 + (t15 + (t208 * t195 + (t209 - t334) * t196 + t336) * t195 + t335 * t193) * t314) * qJD(2), -t26 * t303 / 0.2e1, t14 * qJD(2) + t6 * qJD(5) + t268 + (-0.4e1 * t332 + 0.2e1 * t256 * (-t191 * t331 + t260)) * qJD(4), -t269 + t2 * qJD(2) + t6 * qJD(4) - t270 * m(6) / 0.2e1 + (-t22 / 0.2e1 + (t17 / 0.2e1 - t4 / 0.2e1) * t196 + (t18 / 0.2e1 - t3 / 0.2e1) * t195) * t258 + (t12 * t315 + t13 * t314 + (-t8 / 0.2e1 + (-t42 / 0.2e1 + t19 / 0.2e1) * t196 + (t41 / 0.2e1 + t20 / 0.2e1) * t195) * t190 + (t43 * t65 - t60 * t69 + t71 * t77 + t72 * t75 + (-t195 * t68 - t196 * t67) * t162 - t11) * m(6)) * qJD(5); 0, ((t195 * t78 - t196 * t76) * t322 + (t195 * t96 - t196 * t94) * t323 + (-t125 * t196 + t126 * t195) * t324) * t325 + qJD(5) * t255, 0, 0, t24 + (t195 * t71 - t196 * t72) * t303; -t81 * qJD(2), -t268 + t70 * qJD(4) + t5 * qJD(5) + 0.4e1 * (-t319 / 0.4e1 - t321 / 0.4e1) * qJD(2) + ((-t191 * t55 + t306) * t322 + (-t191 * t61 + t305) * t323 + ((t195 * t76 + t196 * t78 + t43) * t322 + (t195 * t94 + t196 * t96 + t54) * t323) * t190) * t325, 0, t70 * qJD(2), t5 * qJD(2) + (-t191 * t65 + (t195 * t72 + t196 * t71) * t190) * t303; t36 * qJD(2), t269 + (t3 * t314 + t191 * (t195 * t50 - t196 * t51) / 0.2e1 + (t195 * t38 - t196 * t39) * t316 - (t195 * t46 - t196 * t47) * t277 / 0.2e1 + t16 * t275 / 0.2e1 - (t195 * t48 - t196 * t49) * t278 / 0.2e1 + t15 * t276 / 0.2e1 + t4 * t315 - t215) * qJD(2) + t7 * qJD(4) + t1 * qJD(5) + ((t43 * t45 + t52 * t77 + t53 * t75 - t55 * t60 + t67 * t78 + t68 * t76 - t21) * qJD(2) + t270 / 0.2e1) * m(6), t24, t7 * qJD(2), t1 * qJD(2) + (m(6) * (-t60 * t65 + t67 * t71 + t68 * t72) + t190 ^ 2 * t281 / 0.2e1) * qJD(5) + (t12 * t313 + t13 * t315 + (t41 * t196 + t42 * t195 + (t198 * t266 - t200 * t265) * t190) * t316) * t258;];
Cq = t10;
