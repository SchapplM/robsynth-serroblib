% Calculate joint inertia matrix for
% S6RRPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:16:48
% EndTime: 2019-03-09 12:16:57
% DurationCPUTime: 4.21s
% Computational Cost: add. (6790->432), mult. (17172->643), div. (0->0), fcn. (20663->8), ass. (0->209)
t313 = Icges(3,1) + Icges(4,1);
t311 = Icges(4,4) + Icges(3,5);
t185 = sin(qJ(2));
t312 = (Icges(3,4) - Icges(4,5)) * t185;
t310 = Icges(4,2) + Icges(3,3);
t188 = cos(qJ(2));
t309 = t311 * t188 + (-Icges(3,6) + Icges(4,6)) * t185;
t308 = t313 * t188 - t312;
t184 = sin(qJ(4));
t286 = cos(qJ(4));
t152 = t185 * t184 + t188 * t286;
t186 = sin(qJ(1));
t146 = t152 * t186;
t183 = sin(qJ(5));
t187 = cos(qJ(5));
t189 = cos(qJ(1));
t117 = t146 * t183 - t189 * t187;
t118 = t146 * t187 + t183 * t189;
t234 = t185 * t286;
t153 = -t188 * t184 + t234;
t145 = t153 * t186;
t300 = rSges(7,3) + qJ(6);
t303 = rSges(7,1) + pkin(5);
t282 = -t145 * rSges(7,2) + t117 * t300 + t118 * t303;
t108 = Icges(5,4) * t153 - Icges(5,2) * t152;
t307 = t108 / 0.2e1;
t306 = -t186 / 0.2e1;
t305 = t186 / 0.2e1;
t304 = -t189 / 0.2e1;
t296 = t189 / 0.2e1;
t256 = Icges(7,3) * t117;
t257 = Icges(7,6) * t145;
t262 = Icges(7,5) * t118;
t199 = t256 - t257 + t262;
t201 = Icges(7,4) * t118 - Icges(7,2) * t145 + Icges(7,6) * t117;
t267 = Icges(7,4) * t145;
t272 = Icges(7,1) * t118;
t204 = Icges(7,5) * t117 - t267 + t272;
t27 = t152 * t201 + (t183 * t199 + t187 * t204) * t153;
t200 = Icges(6,5) * t118 - Icges(6,6) * t117 - Icges(6,3) * t145;
t258 = Icges(6,6) * t145;
t260 = Icges(6,2) * t117;
t268 = Icges(6,4) * t118;
t202 = -t258 - t260 + t268;
t263 = Icges(6,5) * t145;
t273 = Icges(6,1) * t118;
t205 = -Icges(6,4) * t117 - t263 + t273;
t29 = t152 * t200 + (-t183 * t202 + t187 * t205) * t153;
t302 = t27 + t29;
t148 = t152 * t189;
t119 = t148 * t183 + t186 * t187;
t120 = t148 * t187 - t183 * t186;
t250 = t188 * t189;
t147 = t184 * t250 - t189 * t234;
t62 = Icges(7,5) * t120 + Icges(7,6) * t147 + Icges(7,3) * t119;
t64 = Icges(7,4) * t120 + Icges(7,2) * t147 + Icges(7,6) * t119;
t66 = Icges(7,1) * t120 + Icges(7,4) * t147 + Icges(7,5) * t119;
t28 = t152 * t64 + (t183 * t62 + t187 * t66) * t153;
t63 = Icges(6,5) * t120 - Icges(6,6) * t119 + Icges(6,3) * t147;
t65 = Icges(6,4) * t120 - Icges(6,2) * t119 + Icges(6,6) * t147;
t67 = Icges(6,1) * t120 - Icges(6,4) * t119 + Icges(6,5) * t147;
t30 = t152 * t63 + (-t183 * t65 + t187 * t67) * t153;
t301 = t28 + t30;
t254 = t153 * t183;
t83 = Icges(7,6) * t152 + (Icges(7,5) * t187 + Icges(7,3) * t183) * t153;
t84 = Icges(6,3) * t152 + (Icges(6,5) * t187 - Icges(6,6) * t183) * t153;
t85 = Icges(7,2) * t152 + (Icges(7,4) * t187 + Icges(7,6) * t183) * t153;
t87 = Icges(7,4) * t152 + (Icges(7,1) * t187 + Icges(7,5) * t183) * t153;
t88 = Icges(6,5) * t152 + (Icges(6,1) * t187 - Icges(6,4) * t183) * t153;
t299 = t83 * t254 + (t87 + t88) * t153 * t187 + (t84 + t85) * t152;
t295 = t310 * t186 + t309 * t189;
t294 = -t309 * t186 + t310 * t189;
t182 = t189 ^ 2;
t259 = Icges(5,6) * t189;
t261 = Icges(5,2) * t145;
t264 = Icges(5,5) * t189;
t269 = Icges(5,4) * t146;
t274 = Icges(5,1) * t146;
t291 = t145 ^ 2;
t193 = Icges(7,2) * t291 + (-0.2e1 * t267 + t272) * t118 + (t256 - 0.2e1 * t257 + 0.2e1 * t262) * t117;
t21 = t117 * t62 + t118 * t66 - t145 * t64;
t5 = t186 * t21 - t189 * t193;
t192 = Icges(6,3) * t291 + (-0.2e1 * t263 + t273) * t118 + (0.2e1 * t258 + t260 - 0.2e1 * t268) * t117;
t22 = -t117 * t65 + t118 * t67 - t145 * t63;
t6 = t186 * t22 - t189 * t192;
t92 = Icges(5,5) * t148 - Icges(5,6) * t147 - Icges(5,3) * t186;
t93 = Icges(5,4) * t148 - Icges(5,2) * t147 - Icges(5,6) * t186;
t94 = Icges(5,1) * t148 - Icges(5,4) * t147 - Icges(5,5) * t186;
t243 = -t6 - t5 - (t145 * t93 + t146 * t94 + t189 * t92) * t186 + (Icges(5,3) * t182 + (0.2e1 * t264 + t274) * t146 - (-0.2e1 * t259 - t261 - 0.2e1 * t269) * t145) * t189;
t203 = t259 + t261 + t269;
t206 = Icges(5,4) * t145 + t264 + t274;
t190 = t119 * t199 + t120 * t204 + t147 * t201;
t23 = t119 * t62 + t120 * t66 + t147 * t64;
t7 = t186 * t23 - t189 * t190;
t191 = -t119 * t202 + t120 * t205 + t147 * t200;
t24 = -t119 * t65 + t120 * t67 + t147 * t63;
t8 = t186 * t24 - t189 * t191;
t242 = t8 + t7 + (-t147 * t93 + t148 * t94 - t186 * t92) * t186 - (t148 * t206 - t147 * t203 - t186 * (Icges(5,5) * t146 + Icges(5,6) * t145 + Icges(5,3) * t189)) * t189;
t293 = -t242 * t186 - t243 * t189;
t230 = t148 * pkin(4) + pkin(9) * t147;
t231 = t146 * pkin(4) - t145 * pkin(9);
t279 = t186 * t231 + t189 * t230;
t281 = t147 * rSges(7,2) + t119 * t300 + t120 * t303;
t12 = -t282 * t186 - t281 * t189 - t279;
t32 = t117 * t83 + t118 * t87 - t145 * t85;
t1 = -t145 * t193 + t21 * t147 + t32 * t152;
t86 = Icges(6,6) * t152 + (Icges(6,4) * t187 - Icges(6,2) * t183) * t153;
t33 = -t117 * t86 + t118 * t88 - t145 * t84;
t2 = -t145 * t192 + t22 * t147 + t33 * t152;
t34 = t119 * t83 + t120 * t87 + t147 * t85;
t3 = -t145 * t190 + t23 * t147 + t34 * t152;
t35 = -t119 * t86 + t120 * t88 + t147 * t84;
t4 = -t145 * t191 + t24 * t147 + t35 * t152;
t292 = -(t6 / 0.2e1 + t5 / 0.2e1) * t145 + (t8 / 0.2e1 + t7 / 0.2e1) * t147 - (t296 * t302 + t301 * t306) * t152 + (t4 / 0.2e1 + t3 / 0.2e1) * t186 - (t2 / 0.2e1 + t1 / 0.2e1) * t189;
t181 = t186 ^ 2;
t290 = m(4) / 0.2e1;
t289 = m(6) / 0.2e1;
t288 = m(7) / 0.2e1;
t287 = -rSges(5,3) - pkin(8);
t285 = pkin(8) * t186;
t284 = pkin(8) * t189;
t283 = (-t86 * t254 + t299) * t152;
t277 = t189 * rSges(4,2);
t276 = t189 * rSges(3,3);
t275 = t152 * rSges(7,2) + (t183 * t300 + t187 * t303) * t153;
t270 = Icges(3,4) * t188;
t265 = Icges(4,5) * t188;
t255 = qJ(3) * t185;
t252 = t185 * t189;
t249 = t148 * rSges(5,1) - t147 * rSges(5,2);
t246 = pkin(2) * t250 + qJ(3) * t252;
t248 = t181 * (pkin(2) * t188 + t255) + t189 * t246;
t161 = pkin(2) * t185 - qJ(3) * t188;
t247 = -rSges(4,1) * t185 + rSges(4,3) * t188 - t161;
t245 = t189 * pkin(1) + t186 * pkin(7);
t244 = t181 + t182;
t71 = t120 * rSges(6,1) - t119 * rSges(6,2) + t147 * rSges(6,3);
t235 = rSges(4,1) * t250 + t186 * rSges(4,2) + rSges(4,3) * t252;
t233 = t311 * t185 / 0.2e1 + (-Icges(4,6) / 0.2e1 + Icges(3,6) / 0.2e1) * t188;
t232 = -pkin(3) * t185 - t161;
t172 = pkin(3) * t250;
t229 = t186 * (t186 * t188 * pkin(3) + t284) + t189 * (t172 - t285) + t248;
t228 = t245 + t246;
t90 = t152 * rSges(6,3) + (rSges(6,1) * t187 - rSges(6,2) * t183) * t153;
t227 = t232 - t90;
t114 = rSges(5,1) * t153 - rSges(5,2) * t152;
t226 = -t114 + t232;
t225 = rSges(3,1) * t188 - rSges(3,2) * t185;
t224 = -t146 * rSges(5,1) - t145 * rSges(5,2);
t81 = t226 * t186;
t82 = t226 * t189;
t223 = t186 * t81 + t189 * t82;
t58 = t186 * t224 - t189 * t249;
t222 = t172 + t228;
t219 = -Icges(3,2) * t185 + t270;
t216 = Icges(4,3) * t185 + t265;
t211 = t232 - t275;
t210 = rSges(3,1) * t250 - rSges(3,2) * t252 + t186 * rSges(3,3);
t209 = t30 / 0.2e1 + t35 / 0.2e1 + t34 / 0.2e1 + t28 / 0.2e1;
t208 = t33 / 0.2e1 + t32 / 0.2e1 + t27 / 0.2e1 + t29 / 0.2e1;
t178 = t189 * pkin(7);
t196 = t178 + (-t255 - pkin(1) + (-pkin(2) - pkin(3)) * t188) * t186;
t72 = t287 * t189 + t196 + t224;
t73 = t287 * t186 + t222 + t249;
t207 = m(5) * (t186 * t73 + t189 * t72);
t69 = t118 * rSges(6,1) - t117 * rSges(6,2) - t145 * rSges(6,3);
t31 = -t186 * t69 - t189 * t71 - t279;
t107 = Icges(5,5) * t153 - Icges(5,6) * t152;
t109 = Icges(5,1) * t153 - Icges(5,4) * t152;
t198 = t152 * t93 / 0.2e1 - t153 * t94 / 0.2e1 + t107 * t305 + t147 * t307 - t109 * t148 / 0.2e1 - t209;
t197 = t107 * t296 + t145 * t307 + t146 * t109 / 0.2e1 - t152 * t203 / 0.2e1 + t153 * t206 / 0.2e1 + t208;
t195 = t222 + t230 - t285;
t194 = t196 - t231 - t284;
t165 = rSges(2,1) * t189 - t186 * rSges(2,2);
t164 = -t186 * rSges(2,1) - rSges(2,2) * t189;
t163 = rSges(3,1) * t185 + rSges(3,2) * t188;
t124 = t247 * t189;
t123 = t247 * t186;
t122 = t210 + t245;
t121 = t276 + t178 + (-pkin(1) - t225) * t186;
t116 = pkin(4) * t153 + pkin(9) * t152;
t105 = t189 * t116;
t104 = t186 * t116;
t100 = t228 + t235;
t99 = t277 + t178 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t188 + (-rSges(4,3) - qJ(3)) * t185) * t186;
t91 = t189 * t210 + (t225 * t186 - t276) * t186;
t61 = t189 * t235 + (-t277 + (rSges(4,1) * t188 + rSges(4,3) * t185) * t186) * t186 + t248;
t60 = t189 * t90 + t105;
t59 = t186 * t90 + t104;
t57 = t227 * t189 - t105;
t56 = t227 * t186 - t104;
t53 = t275 * t189 + t105;
t52 = t275 * t186 + t104;
t49 = t211 * t189 - t105;
t48 = t211 * t186 - t104;
t47 = t195 + t71;
t46 = t194 - t69;
t45 = -t147 * t90 + t152 * t71;
t44 = -t145 * t90 - t152 * t69;
t43 = -t58 + t229;
t40 = t145 * t71 + t147 * t69;
t37 = t195 + t281;
t36 = t194 - t282;
t26 = -t275 * t147 + t281 * t152;
t25 = -t275 * t145 - t282 * t152;
t20 = -t31 + t229;
t15 = t281 * t145 + t282 * t147;
t11 = -t12 + t229;
t9 = [-t152 * t108 + Icges(2,3) + (-t183 * t86 + t109) * t153 + m(7) * (t36 ^ 2 + t37 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2) + m(5) * (t72 ^ 2 + t73 ^ 2) + m(4) * (t100 ^ 2 + t99 ^ 2) + m(3) * (t121 ^ 2 + t122 ^ 2) + m(2) * (t164 ^ 2 + t165 ^ 2) + ((Icges(3,2) + Icges(4,3)) * t188 + t312) * t188 + (t313 * t185 - t265 + t270) * t185 + t299; (t233 * t189 + (Icges(3,6) * t296 + Icges(4,6) * t304 + t216 * t305 + t219 * t306) * t188 + (t311 * t296 + t308 * t306) * t185 - t197) * t189 + (t233 * t186 + (Icges(3,6) * t305 + Icges(4,6) * t306 + t216 * t304 + t219 * t296) * t188 + (t308 * t296 + t311 * t305) * t185 - t198) * t186 + m(5) * (t72 * t82 + t73 * t81) + m(6) * (t46 * t57 + t47 * t56) + m(7) * (t36 * t49 + t37 * t48) + m(4) * (t100 * t123 + t124 * t99) + m(3) * (-t121 * t189 - t122 * t186) * t163; m(7) * (t11 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t20 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t43 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t123 ^ 2 + t124 ^ 2 + t61 ^ 2) + m(3) * (t244 * t163 ^ 2 + t91 ^ 2) + (t294 * t182 + t243) * t189 + (t295 * t181 + (t294 * t186 + t295 * t189) * t189 + t242) * t186; 0.2e1 * ((t186 * t37 + t189 * t36) * t288 + (t186 * t47 + t189 * t46) * t289 + t207 / 0.2e1 + (t100 * t186 + t189 * t99) * t290) * t185; m(7) * (-t188 * t11 + (t186 * t48 + t189 * t49) * t185) + m(6) * (-t188 * t20 + (t186 * t56 + t189 * t57) * t185) + m(5) * (t223 * t185 - t188 * t43) + m(4) * (-t188 * t61 + (t123 * t186 + t124 * t189) * t185); 0.2e1 * (t290 + m(5) / 0.2e1 + t289 + t288) * (t244 * t185 ^ 2 + t188 ^ 2); t197 * t189 + t198 * t186 + m(7) * (t36 * t53 + t37 * t52) + m(6) * (t46 * t60 + t47 * t59) + t114 * t207; m(7) * (t11 * t12 + t48 * t52 + t49 * t53) + m(6) * (t20 * t31 + t56 * t59 + t57 * t60) + m(5) * (t223 * t114 + t58 * t43) + t293; m(5) * (t244 * t185 * t114 - t58 * t188) + m(6) * (-t31 * t188 + (t186 * t59 + t189 * t60) * t185) + m(7) * (-t12 * t188 + (t186 * t52 + t189 * t53) * t185); m(7) * (t12 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t31 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t244 * t114 ^ 2 + t58 ^ 2) - t293; m(7) * (t25 * t36 + t26 * t37) + m(6) * (t44 * t46 + t45 * t47) + t209 * t147 - t208 * t145 + t283; m(7) * (t11 * t15 + t25 * t49 + t26 * t48) + m(6) * (t20 * t40 + t44 * t57 + t45 * t56) + t292; m(6) * (-t40 * t188 + (t186 * t45 + t189 * t44) * t185) + m(7) * (-t15 * t188 + (t186 * t26 + t189 * t25) * t185); m(7) * (t12 * t15 + t25 * t53 + t26 * t52) + m(6) * (t31 * t40 + t44 * t60 + t45 * t59) - t292; t283 * t152 + (t152 * t301 + t3 + t4) * t147 - (t152 * t302 + t1 + t2) * t145 + m(7) * (t15 ^ 2 + t25 ^ 2 + t26 ^ 2) + m(6) * (t40 ^ 2 + t44 ^ 2 + t45 ^ 2); m(7) * (t117 * t37 + t119 * t36); m(7) * (t11 * t254 + t117 * t48 + t119 * t49); m(7) * (-t188 * t254 + (t117 * t186 + t119 * t189) * t185); m(7) * (t117 * t52 + t119 * t53 + t12 * t254); m(7) * (t117 * t26 + t119 * t25 + t15 * t254); m(7) * (t153 ^ 2 * t183 ^ 2 + t117 ^ 2 + t119 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
