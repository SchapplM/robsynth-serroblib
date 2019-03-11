% Calculate joint inertia matrix for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:25:03
% EndTime: 2019-03-09 15:25:14
% DurationCPUTime: 5.15s
% Computational Cost: add. (6928->373), mult. (6339->544), div. (0->0), fcn. (6410->10), ass. (0->187)
t305 = Icges(5,4) + Icges(6,6);
t304 = Icges(5,1) + Icges(6,2);
t303 = -Icges(5,2) - Icges(6,3);
t188 = qJ(2) + qJ(3);
t175 = pkin(10) + t188;
t173 = cos(t175);
t302 = t305 * t173;
t172 = sin(t175);
t301 = t305 * t172;
t300 = Icges(6,4) - Icges(5,5);
t299 = Icges(6,5) - Icges(5,6);
t298 = t303 * t172 + t302;
t297 = t304 * t173 - t301;
t191 = sin(qJ(1));
t194 = cos(qJ(1));
t296 = t298 * t191 + t299 * t194;
t295 = -t299 * t191 + t298 * t194;
t294 = t297 * t191 + t300 * t194;
t293 = -t300 * t191 + t297 * t194;
t292 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t176 = sin(t188);
t177 = cos(t188);
t291 = Icges(4,5) * t177 - Icges(4,6) * t176 + t299 * t172 - t300 * t173;
t290 = t303 * t173 - t301;
t289 = t304 * t172 + t302;
t288 = -t191 * t291 + t194 * t292;
t287 = t191 * t292 + t194 * t291;
t261 = Icges(4,4) * t177;
t224 = -Icges(4,2) * t176 + t261;
t121 = Icges(4,6) * t191 + t224 * t194;
t262 = Icges(4,4) * t176;
t227 = Icges(4,1) * t177 - t262;
t123 = Icges(4,5) * t191 + t227 * t194;
t286 = t121 * t176 - t123 * t177 + t172 * t295 - t173 * t293;
t120 = -Icges(4,6) * t194 + t224 * t191;
t122 = -Icges(4,5) * t194 + t227 * t191;
t285 = t120 * t176 - t122 * t177 + t172 * t296 - t294 * t173;
t186 = t191 ^ 2;
t284 = t191 * pkin(7);
t253 = t173 * t194;
t192 = cos(qJ(6));
t249 = t194 * t192;
t189 = sin(qJ(6));
t252 = t191 * t189;
t132 = t172 * t249 - t252;
t250 = t194 * t189;
t251 = t191 * t192;
t133 = t172 * t250 + t251;
t71 = t133 * rSges(7,1) + t132 * rSges(7,2) + rSges(7,3) * t253;
t283 = t191 * pkin(5) + pkin(9) * t253 + t71;
t232 = rSges(4,1) * t177 - rSges(4,2) * t176;
t282 = Icges(4,5) * t176 + Icges(4,6) * t177 - t300 * t172 - t299 * t173;
t150 = Icges(4,2) * t177 + t262;
t151 = Icges(4,1) * t176 + t261;
t281 = -t150 * t176 + t151 * t177 + t290 * t172 + t289 * t173;
t187 = t194 ^ 2;
t134 = t172 * t251 + t250;
t135 = t172 * t252 - t249;
t254 = t173 * t191;
t65 = Icges(7,5) * t133 + Icges(7,6) * t132 + Icges(7,3) * t253;
t67 = Icges(7,4) * t133 + Icges(7,2) * t132 + Icges(7,6) * t253;
t69 = Icges(7,1) * t133 + Icges(7,4) * t132 + Icges(7,5) * t253;
t20 = t134 * t67 + t135 * t69 + t65 * t254;
t66 = Icges(7,5) * t135 + Icges(7,6) * t134 + Icges(7,3) * t254;
t68 = Icges(7,4) * t135 + Icges(7,2) * t134 + Icges(7,6) * t254;
t70 = Icges(7,1) * t135 + Icges(7,4) * t134 + Icges(7,5) * t254;
t21 = t134 * t68 + t135 * t70 + t66 * t254;
t9 = t20 * t191 - t21 * t194;
t280 = -t9 + t288 * t187 + (t286 * t191 + (-t285 + t287) * t194) * t191;
t279 = m(6) / 0.2e1;
t278 = m(7) / 0.2e1;
t195 = -pkin(8) - pkin(7);
t277 = t191 / 0.2e1;
t276 = -t194 / 0.2e1;
t190 = sin(qJ(2));
t275 = pkin(2) * t190;
t274 = pkin(3) * t176;
t193 = cos(qJ(2));
t174 = t193 * pkin(2) + pkin(1);
t153 = pkin(3) * t177 + t174;
t147 = t194 * t153;
t167 = t194 * t174;
t273 = t194 * (t147 - t167) + (t153 - t174) * t186;
t184 = t194 * pkin(7);
t272 = t191 * (t184 + (-pkin(1) + t174) * t191) + t194 * (-t194 * pkin(1) + t167 - t284);
t271 = rSges(3,1) * t193;
t269 = rSges(3,2) * t190;
t267 = t194 * rSges(3,3);
t28 = t172 * t65 + (-t189 * t69 - t192 * t67) * t173;
t266 = t28 * t191;
t29 = t172 * t66 + (-t189 * t70 - t192 * t68) * t173;
t265 = t29 * t194;
t264 = Icges(3,4) * t190;
t263 = Icges(3,4) * t193;
t256 = qJ(5) * t172;
t255 = t172 * t194;
t205 = t191 * rSges(4,3) + t194 * t232;
t64 = t191 * (-t194 * rSges(4,3) + t232 * t191) + t194 * t205;
t248 = pkin(4) * t253 + qJ(5) * t255;
t247 = t191 * rSges(3,3) + t194 * t271;
t244 = t186 + t187;
t243 = t279 + t278;
t18 = t132 * t67 + t133 * t69 + t65 * t253;
t19 = t132 * t68 + t133 * t70 + t66 * t253;
t8 = t18 * t191 - t19 * t194;
t242 = (t8 + t287 * t186 + ((-t286 + t288) * t191 + t285 * t194) * t194) * t191;
t152 = t176 * rSges(4,1) + t177 * rSges(4,2);
t241 = -t152 - t275;
t144 = t172 * pkin(4) - t173 * qJ(5);
t240 = -t144 - t274;
t146 = t172 * rSges(5,1) + t173 * rSges(5,2);
t239 = -t146 - t274;
t204 = rSges(5,1) * t253 - rSges(5,2) * t255 + t191 * rSges(5,3);
t231 = rSges(5,1) * t173 - rSges(5,2) * t172;
t32 = t191 * (-t194 * rSges(5,3) + t231 * t191) + t194 * t204 + t273;
t185 = -qJ(4) + t195;
t238 = -t191 * t185 + t147;
t237 = t186 * (pkin(4) * t173 + t256) + t194 * t248 + t273;
t85 = Icges(7,3) * t172 + (-Icges(7,5) * t189 - Icges(7,6) * t192) * t173;
t86 = Icges(7,6) * t172 + (-Icges(7,4) * t189 - Icges(7,2) * t192) * t173;
t87 = Icges(7,5) * t172 + (-Icges(7,1) * t189 - Icges(7,4) * t192) * t173;
t33 = t132 * t86 + t133 * t87 + t85 * t253;
t3 = t33 * t172 + (t18 * t194 + t19 * t191) * t173;
t34 = t134 * t86 + t135 * t87 + t85 * t254;
t4 = t34 * t172 + (t191 * t21 + t194 * t20) * t173;
t236 = t4 * t276 + t3 * t277 + t172 * (-t265 + t266) / 0.2e1 + t9 * t254 / 0.2e1 + t8 * t253 / 0.2e1;
t145 = -t172 * rSges(6,2) - t173 * rSges(6,3);
t235 = -t145 + t240;
t234 = -t274 - t275;
t233 = -t269 + t271;
t230 = -t135 * rSges(7,1) - t134 * rSges(7,2);
t229 = -t189 * t87 - t192 * t86;
t228 = Icges(3,1) * t193 - t264;
t225 = -Icges(3,2) * t190 + t263;
t221 = Icges(3,5) * t193 - Icges(3,6) * t190;
t203 = t191 * rSges(6,1) - rSges(6,2) * t253 + rSges(6,3) * t255;
t23 = t191 * (-t194 * rSges(6,1) + (-rSges(6,2) * t173 + rSges(6,3) * t172) * t191) + t194 * t203 + t237;
t202 = -t146 + t234;
t201 = t238 + t248;
t88 = t172 * rSges(7,3) + (-rSges(7,1) * t189 - rSges(7,2) * t192) * t173;
t200 = -pkin(9) * t172 + t240 - t88;
t199 = -t144 - t145 + t234;
t72 = rSges(7,3) * t254 - t230;
t12 = t237 + t283 * t194 + (-t194 * pkin(5) + pkin(9) * t254 + t72) * t191;
t198 = t280 * t194 + t242;
t197 = t200 - t275;
t196 = t266 / 0.2e1 - t265 / 0.2e1 + (t177 * t121 + t176 * t123 + t172 * t293 + t173 * t295 + t282 * t191 + t281 * t194 + t33) * t277 + (t177 * t120 + t176 * t122 + t294 * t172 + t173 * t296 + t281 * t191 - t282 * t194 + t34) * t276;
t166 = t194 * rSges(2,1) - t191 * rSges(2,2);
t165 = -t191 * rSges(2,1) - t194 * rSges(2,2);
t164 = t190 * rSges(3,1) + t193 * rSges(3,2);
t127 = Icges(3,3) * t191 + t221 * t194;
t126 = -Icges(3,3) * t194 + t221 * t191;
t117 = t241 * t194;
t116 = t241 * t191;
t96 = t284 + (pkin(1) - t269) * t194 + t247;
t95 = t267 + t184 + (-pkin(1) - t233) * t191;
t94 = t239 * t194;
t93 = t239 * t191;
t84 = t202 * t194;
t83 = t202 * t191;
t82 = t172 * t85;
t81 = -t191 * t195 + t167 + t205;
t80 = (rSges(4,3) - t195) * t194 + (-t174 - t232) * t191;
t77 = t194 * (-t194 * t269 + t247) + (t233 * t191 - t267) * t191;
t76 = t204 + t238;
t75 = (rSges(5,3) - t185) * t194 + (-t153 - t231) * t191;
t74 = t235 * t194;
t73 = t235 * t191;
t61 = t199 * t194;
t60 = t199 * t191;
t55 = t201 + t203;
t54 = (rSges(6,1) - t185) * t194 + (-t153 + (rSges(6,2) - pkin(4)) * t173 + (-rSges(6,3) - qJ(5)) * t172) * t191;
t49 = t200 * t194;
t48 = t200 * t191;
t43 = t197 * t194;
t42 = t197 * t191;
t41 = t172 * t71 - t88 * t253;
t40 = -t172 * t72 + t88 * t254;
t39 = t64 + t272;
t38 = t201 + t283;
t37 = (pkin(5) - t185) * t194 + (-t256 - t153 + (-rSges(7,3) - pkin(4) - pkin(9)) * t173) * t191 + t230;
t36 = (t229 * t173 + t82) * t172;
t35 = (-t191 * t71 + t194 * t72) * t173;
t22 = t32 + t272;
t13 = t23 + t272;
t11 = t12 + t272;
t1 = [t177 * t150 + t176 * t151 + t193 * (Icges(3,2) * t193 + t264) + t190 * (Icges(3,1) * t190 + t263) + Icges(2,3) + t82 + t289 * t172 + (t229 - t290) * t173 + m(7) * (t37 ^ 2 + t38 ^ 2) + m(6) * (t54 ^ 2 + t55 ^ 2) + m(5) * (t75 ^ 2 + t76 ^ 2) + m(4) * (t80 ^ 2 + t81 ^ 2) + m(3) * (t95 ^ 2 + t96 ^ 2) + m(2) * (t165 ^ 2 + t166 ^ 2); t196 + m(3) * (-t191 * t96 - t194 * t95) * t164 + m(4) * (t116 * t81 + t117 * t80) + m(5) * (t84 * t75 + t83 * t76) + m(6) * (t61 * t54 + t60 * t55) + m(7) * (t43 * t37 + t42 * t38) + (t187 / 0.2e1 + t186 / 0.2e1) * (Icges(3,5) * t190 + Icges(3,6) * t193) + (t193 * (Icges(3,6) * t191 + t225 * t194) + t190 * (Icges(3,5) * t191 + t228 * t194)) * t277 + (t193 * (-Icges(3,6) * t194 + t225 * t191) + t190 * (-Icges(3,5) * t194 + t228 * t191)) * t276; m(7) * (t11 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t13 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t22 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (t116 ^ 2 + t117 ^ 2 + t39 ^ 2) + t191 * t186 * t127 + m(3) * (t244 * t164 ^ 2 + t77 ^ 2) + t242 + (-t187 * t126 + (-t191 * t126 + t194 * t127) * t191 + t280) * t194; t196 + m(7) * (t49 * t37 + t48 * t38) + m(6) * (t74 * t54 + t73 * t55) + m(5) * (t94 * t75 + t93 * t76) + m(4) * (-t191 * t81 - t194 * t80) * t152; m(7) * (t12 * t11 + t48 * t42 + t49 * t43) + m(6) * (t23 * t13 + t73 * t60 + t74 * t61) + m(5) * (t32 * t22 + t93 * t83 + t94 * t84) + m(4) * (t64 * t39 + (-t116 * t191 - t117 * t194) * t152) + t198; m(7) * (t12 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(6) * (t23 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(5) * (t32 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(4) * (t244 * t152 ^ 2 + t64 ^ 2) + t198; m(7) * (t191 * t37 - t194 * t38) + m(6) * (t191 * t54 - t194 * t55) + m(5) * (t191 * t75 - t194 * t76); m(7) * (t191 * t43 - t194 * t42) + m(6) * (t191 * t61 - t194 * t60) + m(5) * (t191 * t84 - t194 * t83); m(7) * (t191 * t49 - t194 * t48) + m(6) * (t191 * t74 - t194 * t73) + m(5) * (t191 * t94 - t194 * t93); 0.2e1 * (m(5) / 0.2e1 + t243) * t244; 0.2e1 * ((t191 * t38 + t194 * t37) * t278 + (t191 * t55 + t194 * t54) * t279) * t172; m(7) * (-t173 * t11 + (t191 * t42 + t194 * t43) * t172) + m(6) * (-t173 * t13 + (t191 * t60 + t194 * t61) * t172); m(7) * (-t173 * t12 + (t191 * t48 + t194 * t49) * t172) + m(6) * (-t173 * t23 + (t191 * t73 + t194 * t74) * t172); 0; 0.2e1 * t243 * (t244 * t172 ^ 2 + t173 ^ 2); m(7) * (t40 * t37 + t41 * t38) + t36 + ((t28 / 0.2e1 + t33 / 0.2e1) * t194 + (t29 / 0.2e1 + t34 / 0.2e1) * t191) * t173; m(7) * (t35 * t11 + t40 * t43 + t41 * t42) + t236; m(7) * (t35 * t12 + t40 * t49 + t41 * t48) + t236; m(7) * (t40 * t191 - t41 * t194); m(7) * (-t35 * t173 + (t191 * t41 + t194 * t40) * t172); t172 * t36 + m(7) * (t35 ^ 2 + t40 ^ 2 + t41 ^ 2) + (t194 * t3 + t191 * t4 + t172 * (t191 * t29 + t194 * t28)) * t173;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
