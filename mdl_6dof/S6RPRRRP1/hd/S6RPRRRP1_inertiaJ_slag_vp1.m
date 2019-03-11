% Calculate joint inertia matrix for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:32
% EndTime: 2019-03-09 05:55:39
% DurationCPUTime: 3.32s
% Computational Cost: add. (9505->362), mult. (8043->512), div. (0->0), fcn. (8678->10), ass. (0->183)
t181 = qJ(1) + pkin(10);
t176 = sin(t181);
t177 = cos(t181);
t186 = cos(qJ(5));
t182 = qJ(3) + qJ(4);
t179 = cos(t182);
t183 = sin(qJ(5));
t229 = t179 * t183;
t142 = t176 * t229 + t177 * t186;
t228 = t179 * t186;
t143 = t176 * t228 - t177 * t183;
t178 = sin(t182);
t235 = t176 * t178;
t78 = Icges(7,5) * t143 + Icges(7,6) * t235 + Icges(7,3) * t142;
t84 = Icges(6,4) * t143 - Icges(6,2) * t142 + Icges(6,6) * t235;
t278 = t78 - t84;
t144 = -t176 * t186 + t177 * t229;
t145 = t176 * t183 + t177 * t228;
t234 = t177 * t178;
t79 = Icges(7,5) * t145 + Icges(7,6) * t234 + Icges(7,3) * t144;
t85 = Icges(6,4) * t145 - Icges(6,2) * t144 + Icges(6,6) * t234;
t277 = t79 - t85;
t80 = Icges(6,5) * t143 - Icges(6,6) * t142 + Icges(6,3) * t235;
t82 = Icges(7,4) * t143 + Icges(7,2) * t235 + Icges(7,6) * t142;
t276 = t80 + t82;
t81 = Icges(6,5) * t145 - Icges(6,6) * t144 + Icges(6,3) * t234;
t83 = Icges(7,4) * t145 + Icges(7,2) * t234 + Icges(7,6) * t144;
t275 = t81 + t83;
t86 = Icges(7,1) * t143 + Icges(7,4) * t235 + Icges(7,5) * t142;
t88 = Icges(6,1) * t143 - Icges(6,4) * t142 + Icges(6,5) * t235;
t274 = t86 + t88;
t87 = Icges(7,1) * t145 + Icges(7,4) * t234 + Icges(7,5) * t144;
t89 = Icges(6,1) * t145 - Icges(6,4) * t144 + Icges(6,5) * t234;
t273 = t87 + t89;
t260 = rSges(7,3) + qJ(6);
t263 = rSges(7,1) + pkin(5);
t272 = -t260 * t142 - t263 * t143;
t271 = t142 * t278 + t274 * t143 + t276 * t235;
t270 = t142 * t277 + t143 * t273 + t235 * t275;
t269 = t144 * t278 + t274 * t145 + t276 * t234;
t268 = t144 * t277 + t145 * t273 + t234 * t275;
t129 = -Icges(7,6) * t179 + (Icges(7,5) * t186 + Icges(7,3) * t183) * t178;
t131 = -Icges(7,2) * t179 + (Icges(7,4) * t186 + Icges(7,6) * t183) * t178;
t133 = -Icges(7,4) * t179 + (Icges(7,1) * t186 + Icges(7,5) * t183) * t178;
t54 = t129 * t142 + t131 * t235 + t133 * t143;
t130 = -Icges(6,3) * t179 + (Icges(6,5) * t186 - Icges(6,6) * t183) * t178;
t132 = -Icges(6,6) * t179 + (Icges(6,4) * t186 - Icges(6,2) * t183) * t178;
t134 = -Icges(6,5) * t179 + (Icges(6,1) * t186 - Icges(6,4) * t183) * t178;
t55 = t130 * t235 - t132 * t142 + t134 * t143;
t267 = -t54 - t55;
t56 = t129 * t144 + t131 * t234 + t133 * t145;
t57 = t130 * t234 - t132 * t144 + t134 * t145;
t266 = -t56 - t57;
t265 = t267 * t179 + (t271 * t176 + t270 * t177) * t178;
t264 = t266 * t179 + (t269 * t176 + t268 * t177) * t178;
t262 = t270 * t176 - t271 * t177;
t261 = t268 * t176 - t269 * t177;
t241 = rSges(7,2) * t235 - t272;
t259 = rSges(7,2) * t234 + t260 * t144 + t263 * t145;
t258 = -t130 - t131;
t231 = t178 * t183;
t257 = t129 * t231 + (t133 + t134) * t178 * t186;
t201 = Icges(5,5) * t179 - Icges(5,6) * t178;
t113 = -Icges(5,3) * t177 + t176 * t201;
t114 = Icges(5,3) * t176 + t177 * t201;
t175 = t177 ^ 2;
t236 = Icges(5,4) * t179;
t203 = -Icges(5,2) * t178 + t236;
t116 = Icges(5,6) * t176 + t177 * t203;
t237 = Icges(5,4) * t178;
t205 = Icges(5,1) * t179 - t237;
t118 = Icges(5,5) * t176 + t177 * t205;
t199 = -t116 * t178 + t118 * t179;
t115 = -Icges(5,6) * t177 + t176 * t203;
t117 = -Icges(5,5) * t177 + t176 * t205;
t200 = t115 * t178 - t117 * t179;
t256 = -t175 * t113 - (t199 * t176 + (-t114 + t200) * t177) * t176 - t262;
t174 = t176 ^ 2;
t255 = t176 / 0.2e1;
t254 = -t177 / 0.2e1;
t184 = sin(qJ(3));
t252 = pkin(3) * t184;
t251 = pkin(4) * t179;
t185 = sin(qJ(1));
t250 = t185 * pkin(1);
t249 = -t132 * t231 + t179 * t258 + t257;
t187 = cos(qJ(3));
t248 = rSges(4,1) * t187;
t247 = rSges(4,2) * t184;
t246 = t177 * rSges(4,3);
t39 = -t179 * t82 + (t183 * t78 + t186 * t86) * t178;
t245 = t39 * t177;
t40 = -t179 * t83 + (t183 * t79 + t186 * t87) * t178;
t244 = t40 * t176;
t41 = -t179 * t80 + (-t183 * t84 + t186 * t88) * t178;
t243 = t41 * t177;
t42 = -t179 * t81 + (-t183 * t85 + t186 * t89) * t178;
t242 = t42 * t176;
t239 = Icges(4,4) * t184;
t238 = Icges(4,4) * t187;
t233 = t177 * t179;
t189 = -pkin(8) - pkin(7);
t232 = t177 * t189;
t173 = pkin(3) * t187 + pkin(2);
t155 = t177 * t173;
t172 = t177 * pkin(7);
t227 = t176 * (t232 + t172 + (-pkin(2) + t173) * t176) + t177 * (-t177 * pkin(2) + t155 + (-pkin(7) - t189) * t176);
t194 = rSges(5,1) * t233 - rSges(5,2) * t234 + t176 * rSges(5,3);
t209 = rSges(5,1) * t179 - rSges(5,2) * t178;
t72 = t176 * (-t177 * rSges(5,3) + t176 * t209) + t177 * t194;
t222 = pkin(4) * t233 + pkin(9) * t234;
t225 = t174 * (pkin(9) * t178 + t251) + t177 * t222;
t224 = -t179 * rSges(7,2) + (t183 * t260 + t186 * t263) * t178;
t136 = -t179 * rSges(6,3) + (rSges(6,1) * t186 - rSges(6,2) * t183) * t178;
t154 = pkin(4) * t178 - pkin(9) * t179;
t223 = -t136 - t154;
t221 = t176 * rSges(4,3) + t177 * t248;
t220 = t174 + t175;
t219 = (t174 * t114 + (t200 * t177 + (-t113 + t199) * t176) * t177 + t261) * t176;
t218 = -t154 - t224;
t95 = t145 * rSges(6,1) - t144 * rSges(6,2) + rSges(6,3) * t234;
t153 = rSges(5,1) * t178 + rSges(5,2) * t179;
t215 = -t153 - t252;
t214 = -t154 - t252;
t213 = -t173 - t251;
t208 = -t143 * rSges(6,1) + t142 * rSges(6,2);
t93 = rSges(6,3) * t235 - t208;
t43 = t176 * t93 + t177 * t95 + t225;
t212 = -t136 + t214;
t188 = cos(qJ(1));
t180 = t188 * pkin(1);
t211 = -t176 * t189 + t155 + t180;
t210 = -t247 + t248;
t207 = -t232 - t250;
t206 = Icges(4,1) * t187 - t239;
t204 = -Icges(4,2) * t184 + t238;
t202 = Icges(4,5) * t187 - Icges(4,6) * t184;
t151 = Icges(5,2) * t179 + t237;
t152 = Icges(5,1) * t178 + t236;
t196 = -t151 * t178 + t152 * t179;
t195 = t214 - t224;
t22 = t176 * t241 + t177 * t259 + t225;
t193 = t211 + t222;
t192 = t177 * t256 + t219;
t191 = t264 * t255 + t265 * t254 - (t244 - t245 + t242 - t243) * t179 / 0.2e1 + t262 * t235 / 0.2e1 + t261 * t234 / 0.2e1;
t150 = Icges(5,5) * t178 + Icges(5,6) * t179;
t190 = -t245 / 0.2e1 + t244 / 0.2e1 - t243 / 0.2e1 + t242 / 0.2e1 + (t116 * t179 + t118 * t178 + t176 * t150 + t177 * t196 - t266) * t255 + (t115 * t179 + t117 * t178 - t177 * t150 + t176 * t196 - t267) * t254;
t167 = rSges(2,1) * t188 - rSges(2,2) * t185;
t166 = -rSges(2,1) * t185 - rSges(2,2) * t188;
t165 = rSges(4,1) * t184 + rSges(4,2) * t187;
t148 = rSges(3,1) * t177 - rSges(3,2) * t176 + t180;
t147 = -rSges(3,1) * t176 - rSges(3,2) * t177 - t250;
t124 = Icges(4,3) * t176 + t177 * t202;
t123 = -Icges(4,3) * t177 + t176 * t202;
t120 = t215 * t177;
t119 = t215 * t176;
t105 = t176 * pkin(7) + t180 + (pkin(2) - t247) * t177 + t221;
t104 = t246 - t250 + t172 + (-pkin(2) - t210) * t176;
t101 = t223 * t177;
t100 = t223 * t176;
t97 = t194 + t211;
t96 = -t250 + (rSges(5,3) - t189) * t177 + (-t173 - t209) * t176;
t91 = t212 * t177;
t90 = t212 * t176;
t77 = t177 * (-t177 * t247 + t221) + (t176 * t210 - t246) * t176;
t71 = t218 * t177;
t70 = t218 * t176;
t65 = t195 * t177;
t64 = t195 * t176;
t61 = -t136 * t234 - t179 * t95;
t60 = t136 * t235 + t179 * t93;
t59 = t193 + t95;
t58 = ((-rSges(6,3) - pkin(9)) * t178 + t213) * t176 + t207 + t208;
t53 = (-t176 * t95 + t177 * t93) * t178;
t48 = t72 + t227;
t47 = t193 + t259;
t46 = ((-rSges(7,2) - pkin(9)) * t178 + t213) * t176 + t207 + t272;
t45 = -t179 * t259 - t224 * t234;
t44 = t179 * t241 + t224 * t235;
t24 = (-t176 * t259 + t177 * t241) * t178;
t23 = t43 + t227;
t21 = t22 + t227;
t1 = [t187 * (Icges(4,2) * t187 + t239) + t184 * (Icges(4,1) * t184 + t238) + Icges(2,3) + Icges(3,3) + (-t132 * t183 + t152) * t178 + (t151 + t258) * t179 + m(7) * (t46 ^ 2 + t47 ^ 2) + m(6) * (t58 ^ 2 + t59 ^ 2) + m(5) * (t96 ^ 2 + t97 ^ 2) + m(4) * (t104 ^ 2 + t105 ^ 2) + m(3) * (t147 ^ 2 + t148 ^ 2) + m(2) * (t166 ^ 2 + t167 ^ 2) + t257; 0; m(3) + m(4) + m(5) + m(6) + m(7); (t174 / 0.2e1 + t175 / 0.2e1) * (Icges(4,5) * t184 + Icges(4,6) * t187) + t190 + m(4) * (-t104 * t177 - t105 * t176) * t165 + m(7) * (t46 * t65 + t47 * t64) + m(6) * (t58 * t91 + t59 * t90) + m(5) * (t119 * t97 + t120 * t96) + (t187 * (Icges(4,6) * t176 + t177 * t204) + t184 * (Icges(4,5) * t176 + t177 * t206)) * t255 + (t187 * (-Icges(4,6) * t177 + t176 * t204) + t184 * (-Icges(4,5) * t177 + t176 * t206)) * t254; m(4) * t77 + m(5) * t48 + m(6) * t23 + m(7) * t21; m(7) * (t21 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(6) * (t23 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(5) * (t119 ^ 2 + t120 ^ 2 + t48 ^ 2) + m(4) * (t165 ^ 2 * t220 + t77 ^ 2) + t176 * t174 * t124 + t219 + (-t175 * t123 + (-t176 * t123 + t177 * t124) * t176 + t256) * t177; m(7) * (t46 * t71 + t47 * t70) + m(6) * (t100 * t59 + t101 * t58) + m(5) * (-t176 * t97 - t177 * t96) * t153 + t190; m(5) * t72 + m(6) * t43 + m(7) * t22; m(7) * (t21 * t22 + t64 * t70 + t65 * t71) + m(6) * (t100 * t90 + t101 * t91 + t23 * t43) + m(5) * (t72 * t48 + (-t119 * t176 - t120 * t177) * t153) + t192; m(7) * (t22 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (t100 ^ 2 + t101 ^ 2 + t43 ^ 2) + m(5) * (t153 ^ 2 * t220 + t72 ^ 2) + t192; -t249 * t179 + m(7) * (t44 * t46 + t45 * t47) + m(6) * (t58 * t60 + t59 * t61) + ((t42 / 0.2e1 + t40 / 0.2e1 + t57 / 0.2e1 + t56 / 0.2e1) * t177 + (t54 / 0.2e1 + t41 / 0.2e1 + t39 / 0.2e1 + t55 / 0.2e1) * t176) * t178; m(6) * t53 + m(7) * t24; m(7) * (t21 * t24 + t44 * t65 + t45 * t64) + m(6) * (t23 * t53 + t60 * t91 + t61 * t90) + t191; m(7) * (t22 * t24 + t44 * t71 + t45 * t70) + m(6) * (t100 * t61 + t101 * t60 + t43 * t53) + t191; m(7) * (t24 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t53 ^ 2 + t60 ^ 2 + t61 ^ 2) + t249 * t179 ^ 2 + (((-t40 - t42) * t179 + t264) * t177 + ((-t39 - t41) * t179 + t265) * t176) * t178; m(7) * (t142 * t47 + t144 * t46); m(7) * t231; m(7) * (t142 * t64 + t144 * t65 + t21 * t231); m(7) * (t142 * t70 + t144 * t71 + t22 * t231); m(7) * (t142 * t45 + t144 * t44 + t231 * t24); m(7) * (t178 ^ 2 * t183 ^ 2 + t142 ^ 2 + t144 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
