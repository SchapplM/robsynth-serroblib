% Calculate joint inertia matrix for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:02:12
% EndTime: 2019-03-09 18:02:19
% DurationCPUTime: 2.81s
% Computational Cost: add. (11030->415), mult. (8109->605), div. (0->0), fcn. (8348->12), ass. (0->212)
t199 = sin(qJ(1));
t194 = t199 ^ 2;
t297 = t199 * pkin(7);
t196 = qJ(2) + qJ(3);
t183 = pkin(11) + t196;
t181 = qJ(5) + t183;
t175 = cos(t181);
t202 = cos(qJ(1));
t262 = t175 * t202;
t174 = sin(t181);
t263 = t174 * t202;
t200 = cos(qJ(6));
t259 = t199 * t200;
t197 = sin(qJ(6));
t260 = t197 * t202;
t134 = -t175 * t260 + t259;
t258 = t200 * t202;
t261 = t197 * t199;
t135 = t175 * t258 + t261;
t75 = t135 * rSges(7,1) + t134 * rSges(7,2) + rSges(7,3) * t263;
t296 = pkin(5) * t262 + pkin(10) * t263 + t75;
t178 = sin(t183);
t179 = cos(t183);
t184 = sin(t196);
t185 = cos(t196);
t295 = Icges(4,5) * t184 + Icges(5,5) * t178 + Icges(4,6) * t185 + Icges(5,6) * t179;
t268 = Icges(5,4) * t178;
t150 = Icges(5,2) * t179 + t268;
t267 = Icges(5,4) * t179;
t151 = Icges(5,1) * t178 + t267;
t270 = Icges(4,4) * t184;
t156 = Icges(4,2) * t185 + t270;
t269 = Icges(4,4) * t185;
t157 = Icges(4,1) * t184 + t269;
t294 = -t150 * t178 + t151 * t179 - t156 * t184 + t157 * t185;
t242 = rSges(5,1) * t179 - rSges(5,2) * t178;
t243 = rSges(4,1) * t185 - rSges(4,2) * t184;
t229 = Icges(5,5) * t179 - Icges(5,6) * t178;
t113 = -Icges(5,3) * t202 + t229 * t199;
t114 = Icges(5,3) * t199 + t229 * t202;
t230 = Icges(4,5) * t185 - Icges(4,6) * t184;
t123 = -Icges(4,3) * t202 + t230 * t199;
t124 = Icges(4,3) * t199 + t230 * t202;
t195 = t202 ^ 2;
t234 = -Icges(4,2) * t184 + t269;
t126 = Icges(4,6) * t199 + t234 * t202;
t238 = Icges(4,1) * t185 - t270;
t128 = Icges(4,5) * t199 + t238 * t202;
t222 = -t126 * t184 + t128 * t185;
t125 = -Icges(4,6) * t202 + t234 * t199;
t127 = -Icges(4,5) * t202 + t238 * t199;
t223 = t125 * t184 - t127 * t185;
t233 = -Icges(5,2) * t178 + t267;
t116 = Icges(5,6) * t199 + t233 * t202;
t237 = Icges(5,1) * t179 - t268;
t118 = Icges(5,5) * t199 + t237 * t202;
t224 = -t116 * t178 + t118 * t179;
t115 = -Icges(5,6) * t202 + t233 * t199;
t117 = -Icges(5,5) * t202 + t237 * t199;
t225 = t115 * t178 - t117 * t179;
t228 = Icges(6,5) * t175 - Icges(6,6) * t174;
t101 = -Icges(6,3) * t202 + t228 * t199;
t102 = Icges(6,3) * t199 + t228 * t202;
t265 = Icges(6,4) * t175;
t232 = -Icges(6,2) * t174 + t265;
t104 = Icges(6,6) * t199 + t232 * t202;
t266 = Icges(6,4) * t174;
t236 = Icges(6,1) * t175 - t266;
t106 = Icges(6,5) * t199 + t236 * t202;
t226 = -t104 * t174 + t106 * t175;
t103 = -Icges(6,6) * t202 + t232 * t199;
t105 = -Icges(6,5) * t202 + t236 * t199;
t227 = t103 * t174 - t105 * t175;
t132 = -t175 * t261 - t258;
t133 = t175 * t259 - t260;
t264 = t174 * t199;
t68 = Icges(7,5) * t133 + Icges(7,6) * t132 + Icges(7,3) * t264;
t70 = Icges(7,4) * t133 + Icges(7,2) * t132 + Icges(7,6) * t264;
t72 = Icges(7,1) * t133 + Icges(7,4) * t132 + Icges(7,5) * t264;
t17 = t132 * t70 + t133 * t72 + t68 * t264;
t69 = Icges(7,5) * t135 + Icges(7,6) * t134 + Icges(7,3) * t263;
t71 = Icges(7,4) * t135 + Icges(7,2) * t134 + Icges(7,6) * t263;
t73 = Icges(7,1) * t135 + Icges(7,4) * t134 + Icges(7,5) * t263;
t18 = t132 * t71 + t133 * t73 + t69 * t264;
t8 = -t17 * t202 + t18 * t199;
t288 = -t195 * t101 - (t226 * t199 + (-t102 + t227) * t202) * t199 - t8;
t293 = t288 + (-t113 - t123) * t195 + ((-t224 - t222) * t199 + (t114 - t225 + t124 - t223) * t202) * t199;
t203 = -pkin(8) - pkin(7);
t291 = t199 / 0.2e1;
t290 = -t202 / 0.2e1;
t19 = t134 * t70 + t135 * t72 + t68 * t263;
t20 = t134 * t71 + t135 * t73 + t69 * t263;
t9 = -t19 * t202 + t199 * t20;
t289 = (t194 * t102 + t9 + (t227 * t202 + (-t101 + t226) * t199) * t202) * t199;
t198 = sin(qJ(2));
t287 = pkin(2) * t198;
t286 = pkin(3) * t184;
t285 = pkin(5) * t175;
t201 = cos(qJ(2));
t182 = t201 * pkin(2) + pkin(1);
t159 = pkin(3) * t185 + t182;
t153 = t202 * t159;
t172 = t202 * t182;
t284 = t202 * (t153 - t172) + (t159 - t182) * t194;
t213 = rSges(6,1) * t262 - rSges(6,2) * t263 + t199 * rSges(6,3);
t241 = rSges(6,1) * t175 - rSges(6,2) * t174;
t57 = t199 * (-t202 * rSges(6,3) + t241 * t199) + t202 * t213;
t283 = rSges(3,1) * t201;
t280 = rSges(3,2) * t198;
t92 = -Icges(7,6) * t175 + (Icges(7,4) * t200 - Icges(7,2) * t197) * t174;
t277 = t197 * t92;
t276 = t202 * rSges(3,3);
t25 = -t175 * t68 + (-t197 * t70 + t200 * t72) * t174;
t275 = t25 * t202;
t26 = -t175 * t69 + (-t197 * t71 + t200 * t73) * t174;
t274 = t26 * t199;
t94 = -t175 * rSges(7,3) + (rSges(7,1) * t200 - rSges(7,2) * t197) * t174;
t273 = -pkin(5) * t174 + pkin(10) * t175 - t94;
t272 = Icges(3,4) * t198;
t271 = Icges(3,4) * t201;
t192 = t202 * pkin(7);
t257 = t199 * (t192 + (-pkin(1) + t182) * t199) + t202 * (-t202 * pkin(1) + t172 - t297);
t215 = t199 * rSges(4,3) + t243 * t202;
t76 = t199 * (-t202 * rSges(4,3) + t243 * t199) + t202 * t215;
t255 = t199 * rSges(3,3) + t202 * t283;
t252 = t194 + t195;
t193 = -qJ(4) + t203;
t251 = t199 * (t194 * t114 + (t225 * t202 + (-t113 + t224) * t199) * t202) + t199 * (t194 * t124 + (t223 * t202 + (-t123 + t222) * t199) * t202) + t289;
t158 = rSges(4,1) * t184 + rSges(4,2) * t185;
t250 = -t158 - t287;
t249 = -rSges(5,1) * t178 - rSges(5,2) * t179 - t286;
t136 = pkin(4) * t179 + t159;
t129 = t202 * t136;
t248 = t202 * (t129 - t153) + t284 + (t136 - t159) * t194;
t214 = t199 * rSges(5,3) + t242 * t202;
t35 = t199 * (-t202 * rSges(5,3) + t242 * t199) + t202 * t214 + t284;
t186 = -pkin(9) + t193;
t247 = -t186 * t199 + t129;
t240 = -t133 * rSges(7,1) - t132 * rSges(7,2);
t74 = rSges(7,3) * t264 - t240;
t30 = t199 * t74 + t194 * (pkin(10) * t174 + t285) + t296 * t202;
t91 = -Icges(7,3) * t175 + (Icges(7,5) * t200 - Icges(7,6) * t197) * t174;
t93 = -Icges(7,5) * t175 + (Icges(7,1) * t200 - Icges(7,4) * t197) * t174;
t33 = t132 * t92 + t133 * t93 + t91 * t264;
t3 = -t33 * t175 + (t17 * t199 + t18 * t202) * t174;
t34 = t134 * t92 + t135 * t93 + t91 * t263;
t4 = -t34 * t175 + (t19 * t199 + t20 * t202) * t174;
t246 = t3 * t290 + t4 * t291 - t175 * (t274 - t275) / 0.2e1 + t8 * t264 / 0.2e1 + t9 * t263 / 0.2e1;
t245 = -pkin(4) * t178 - t286;
t244 = -t280 + t283;
t14 = t248 + t57;
t239 = Icges(3,1) * t201 - t272;
t235 = -Icges(3,2) * t198 + t271;
t231 = Icges(3,5) * t201 - Icges(3,6) * t198;
t145 = Icges(6,2) * t175 + t266;
t146 = Icges(6,1) * t174 + t265;
t219 = -t145 * t174 + t146 * t175;
t216 = t288 * t202 + t289;
t212 = t249 - t287;
t147 = rSges(6,1) * t174 + rSges(6,2) * t175;
t211 = -t147 + t245;
t210 = t245 + t273;
t144 = Icges(6,5) * t174 + Icges(6,6) * t175;
t209 = -t275 / 0.2e1 + t274 / 0.2e1 + (t104 * t175 + t106 * t174 + t199 * t144 + t219 * t202 + t34) * t291 + (t103 * t175 + t105 * t174 - t202 * t144 + t219 * t199 + t33) * t290;
t208 = t245 - t287;
t12 = t30 + t248;
t207 = -t147 + t208;
t206 = t293 * t202 + t251;
t205 = t208 + t273;
t204 = t209 + (t116 * t179 + t118 * t178 + t126 * t185 + t128 * t184 + t295 * t199 + t294 * t202) * t291 + (t115 * t179 + t117 * t178 + t125 * t185 + t127 * t184 + t294 * t199 - t295 * t202) * t290;
t171 = rSges(2,1) * t202 - rSges(2,2) * t199;
t170 = -rSges(2,1) * t199 - rSges(2,2) * t202;
t169 = rSges(3,1) * t198 + rSges(3,2) * t201;
t138 = Icges(3,3) * t199 + t231 * t202;
t137 = -Icges(3,3) * t202 + t231 * t199;
t122 = t250 * t202;
t121 = t250 * t199;
t108 = t297 + (pkin(1) - t280) * t202 + t255;
t107 = t276 + t192 + (-pkin(1) - t244) * t199;
t100 = t249 * t202;
t99 = t249 * t199;
t90 = t212 * t202;
t89 = t212 * t199;
t88 = -t199 * t203 + t172 + t215;
t87 = (rSges(4,3) - t203) * t202 + (-t182 - t243) * t199;
t84 = t211 * t202;
t83 = t211 * t199;
t82 = t174 * t200 * t93;
t81 = t202 * (-t202 * t280 + t255) + (t244 * t199 - t276) * t199;
t80 = -t199 * t193 + t153 + t214;
t79 = (rSges(5,3) - t193) * t202 + (-t159 - t242) * t199;
t78 = t207 * t202;
t77 = t207 * t199;
t63 = t273 * t202;
t62 = t273 * t199;
t61 = t213 + t247;
t60 = (rSges(6,3) - t186) * t202 + (-t136 - t241) * t199;
t50 = t210 * t202;
t49 = t210 * t199;
t46 = t205 * t202;
t45 = t205 * t199;
t42 = t76 + t257;
t41 = -t175 * t75 - t94 * t263;
t40 = t175 * t74 + t94 * t264;
t39 = t247 + t296;
t38 = -t202 * t186 + (-t285 - t136 + (-rSges(7,3) - pkin(10)) * t174) * t199 + t240;
t37 = -t174 * t277 - t175 * t91 + t82;
t36 = (-t199 * t75 + t202 * t74) * t174;
t27 = t35 + t257;
t13 = t14 + t257;
t10 = t12 + t257;
t1 = [t179 * t150 + t178 * t151 + t185 * t156 + t184 * t157 + t201 * (Icges(3,2) * t201 + t272) + t198 * (Icges(3,1) * t198 + t271) + Icges(2,3) + t82 + (-t91 + t145) * t175 + (t146 - t277) * t174 + m(7) * (t38 ^ 2 + t39 ^ 2) + m(6) * (t60 ^ 2 + t61 ^ 2) + m(5) * (t79 ^ 2 + t80 ^ 2) + m(4) * (t87 ^ 2 + t88 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2) + m(2) * (t170 ^ 2 + t171 ^ 2); (t201 * (-Icges(3,6) * t202 + t235 * t199) + t198 * (-Icges(3,5) * t202 + t239 * t199)) * t290 + (t201 * (Icges(3,6) * t199 + t235 * t202) + t198 * (Icges(3,5) * t199 + t239 * t202)) * t291 + m(3) * (-t107 * t202 - t108 * t199) * t169 + t204 + (t195 / 0.2e1 + t194 / 0.2e1) * (Icges(3,5) * t198 + Icges(3,6) * t201) + m(7) * (t38 * t46 + t39 * t45) + m(6) * (t60 * t78 + t61 * t77) + m(5) * (t79 * t90 + t80 * t89) + m(4) * (t121 * t88 + t122 * t87); m(7) * (t10 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t13 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t27 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(4) * (t121 ^ 2 + t122 ^ 2 + t42 ^ 2) + m(3) * (t252 * t169 ^ 2 + t81 ^ 2) + t199 * t194 * t138 + t251 + (-t195 * t137 + (-t199 * t137 + t202 * t138) * t199 + t293) * t202; m(7) * (t38 * t50 + t39 * t49) + m(6) * (t60 * t84 + t61 * t83) + m(5) * (t100 * t79 + t80 * t99) + m(4) * (-t199 * t88 - t202 * t87) * t158 + t204; m(7) * (t10 * t12 + t45 * t49 + t46 * t50) + m(6) * (t13 * t14 + t77 * t83 + t78 * t84) + m(5) * (t100 * t90 + t27 * t35 + t89 * t99) + m(4) * (t76 * t42 + (-t121 * t199 - t122 * t202) * t158) + t206; m(7) * (t12 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t14 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(5) * (t100 ^ 2 + t35 ^ 2 + t99 ^ 2) + m(4) * (t252 * t158 ^ 2 + t76 ^ 2) + t206; m(7) * (t199 * t38 - t202 * t39) + m(6) * (t199 * t60 - t202 * t61) + m(5) * (t199 * t79 - t202 * t80); m(7) * (t199 * t46 - t202 * t45) + m(6) * (t199 * t78 - t202 * t77) + m(5) * (t199 * t90 - t202 * t89); m(7) * (t199 * t50 - t202 * t49) + m(6) * (t199 * t84 - t202 * t83) + m(5) * (t100 * t199 - t202 * t99); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t252; m(7) * (t38 * t63 + t39 * t62) + m(6) * (-t199 * t61 - t202 * t60) * t147 + t209; m(7) * (t10 * t30 + t45 * t62 + t46 * t63) + m(6) * (t57 * t13 + (-t199 * t77 - t202 * t78) * t147) + t216; m(7) * (t12 * t30 + t49 * t62 + t50 * t63) + m(6) * (t57 * t14 + (-t199 * t83 - t202 * t84) * t147) + t216; m(7) * (t199 * t63 - t202 * t62); m(6) * (t252 * t147 ^ 2 + t57 ^ 2) + m(7) * (t30 ^ 2 + t62 ^ 2 + t63 ^ 2) + t216; m(7) * (t38 * t40 + t39 * t41) - t37 * t175 + ((t34 / 0.2e1 + t26 / 0.2e1) * t202 + (t33 / 0.2e1 + t25 / 0.2e1) * t199) * t174; m(7) * (t10 * t36 + t40 * t46 + t41 * t45) + t246; m(7) * (t12 * t36 + t40 * t50 + t41 * t49) + t246; m(7) * (t199 * t40 - t202 * t41); m(7) * (t30 * t36 + t40 * t63 + t41 * t62) + t246; t175 ^ 2 * t37 + m(7) * (t36 ^ 2 + t40 ^ 2 + t41 ^ 2) + (t202 * t4 + t199 * t3 - t175 * (t199 * t25 + t202 * t26)) * t174;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
