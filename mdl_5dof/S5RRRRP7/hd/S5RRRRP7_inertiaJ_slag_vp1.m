% Calculate joint inertia matrix for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP7_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP7_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:08
% EndTime: 2019-12-31 21:56:17
% DurationCPUTime: 2.99s
% Computational Cost: add. (5730->339), mult. (7641->492), div. (0->0), fcn. (8302->8), ass. (0->178)
t176 = qJ(2) + qJ(3);
t170 = cos(t176);
t180 = cos(qJ(4));
t182 = cos(qJ(1));
t222 = t182 * t180;
t177 = sin(qJ(4));
t179 = sin(qJ(1));
t225 = t179 * t177;
t143 = t170 * t225 + t222;
t223 = t182 * t177;
t224 = t179 * t180;
t144 = t170 * t224 - t223;
t169 = sin(t176);
t230 = t169 * t179;
t79 = Icges(6,5) * t144 + Icges(6,6) * t230 + Icges(6,3) * t143;
t85 = Icges(5,4) * t144 - Icges(5,2) * t143 + Icges(5,6) * t230;
t273 = t79 - t85;
t145 = t170 * t223 - t224;
t146 = t170 * t222 + t225;
t228 = t169 * t182;
t80 = Icges(6,5) * t146 + Icges(6,6) * t228 + Icges(6,3) * t145;
t86 = Icges(5,4) * t146 - Icges(5,2) * t145 + Icges(5,6) * t228;
t272 = t80 - t86;
t81 = Icges(5,5) * t144 - Icges(5,6) * t143 + Icges(5,3) * t230;
t83 = Icges(6,4) * t144 + Icges(6,2) * t230 + Icges(6,6) * t143;
t271 = t81 + t83;
t82 = Icges(5,5) * t146 - Icges(5,6) * t145 + Icges(5,3) * t228;
t84 = Icges(6,4) * t146 + Icges(6,2) * t228 + Icges(6,6) * t145;
t270 = t82 + t84;
t87 = Icges(6,1) * t144 + Icges(6,4) * t230 + Icges(6,5) * t143;
t89 = Icges(5,1) * t144 - Icges(5,4) * t143 + Icges(5,5) * t230;
t269 = t87 + t89;
t88 = Icges(6,1) * t146 + Icges(6,4) * t228 + Icges(6,5) * t145;
t90 = Icges(5,1) * t146 - Icges(5,4) * t145 + Icges(5,5) * t228;
t268 = t88 + t90;
t255 = rSges(6,3) + qJ(5);
t258 = rSges(6,1) + pkin(4);
t267 = -t255 * t143 - t258 * t144;
t266 = t273 * t143 + t269 * t144 + t271 * t230;
t265 = t272 * t143 + t268 * t144 + t270 * t230;
t264 = t273 * t145 + t269 * t146 + t271 * t228;
t263 = t272 * t145 + t268 * t146 + t270 * t228;
t113 = -Icges(6,6) * t170 + (Icges(6,5) * t180 + Icges(6,3) * t177) * t169;
t115 = -Icges(6,2) * t170 + (Icges(6,4) * t180 + Icges(6,6) * t177) * t169;
t117 = -Icges(6,4) * t170 + (Icges(6,1) * t180 + Icges(6,5) * t177) * t169;
t50 = t143 * t113 + t115 * t230 + t144 * t117;
t114 = -Icges(5,3) * t170 + (Icges(5,5) * t180 - Icges(5,6) * t177) * t169;
t116 = -Icges(5,6) * t170 + (Icges(5,4) * t180 - Icges(5,2) * t177) * t169;
t118 = -Icges(5,5) * t170 + (Icges(5,1) * t180 - Icges(5,4) * t177) * t169;
t51 = t114 * t230 - t143 * t116 + t144 * t118;
t262 = -t50 - t51;
t52 = t145 * t113 + t115 * t228 + t146 * t117;
t53 = t114 * t228 - t145 * t116 + t146 * t118;
t261 = -t52 - t53;
t260 = t262 * t170 + (t266 * t179 + t265 * t182) * t169;
t259 = t261 * t170 + (t264 * t179 + t263 * t182) * t169;
t257 = t265 * t179 - t266 * t182;
t256 = t263 * t179 - t264 * t182;
t237 = rSges(6,2) * t230 - t267;
t254 = rSges(6,2) * t228 + t255 * t145 + t258 * t146;
t253 = -t114 - t115;
t231 = t169 * t177;
t252 = t113 * t231 + (t117 + t118) * t169 * t180;
t195 = Icges(4,5) * t170 - Icges(4,6) * t169;
t123 = -Icges(4,3) * t182 + t195 * t179;
t124 = Icges(4,3) * t179 + t195 * t182;
t175 = t182 ^ 2;
t232 = Icges(4,4) * t170;
t197 = -Icges(4,2) * t169 + t232;
t126 = Icges(4,6) * t179 + t197 * t182;
t233 = Icges(4,4) * t169;
t199 = Icges(4,1) * t170 - t233;
t128 = Icges(4,5) * t179 + t199 * t182;
t193 = -t126 * t169 + t128 * t170;
t125 = -Icges(4,6) * t182 + t197 * t179;
t127 = -Icges(4,5) * t182 + t199 * t179;
t194 = t125 * t169 - t127 * t170;
t251 = -t175 * t123 - (t193 * t179 + (-t124 + t194) * t182) * t179 - t257;
t174 = t179 ^ 2;
t249 = t179 / 0.2e1;
t248 = -t182 / 0.2e1;
t178 = sin(qJ(2));
t247 = pkin(2) * t178;
t246 = pkin(3) * t170;
t226 = t177 * t116;
t245 = -t169 * t226 + t253 * t170 + t252;
t181 = cos(qJ(2));
t244 = rSges(3,1) * t181;
t243 = rSges(3,2) * t178;
t242 = t182 * rSges(3,3);
t39 = -t170 * t83 + (t177 * t79 + t180 * t87) * t169;
t241 = t39 * t182;
t40 = -t170 * t84 + (t177 * t80 + t180 * t88) * t169;
t240 = t40 * t179;
t41 = -t170 * t81 + (-t177 * t85 + t180 * t89) * t169;
t239 = t41 * t182;
t42 = -t170 * t82 + (-t177 * t86 + t180 * t90) * t169;
t238 = t42 * t179;
t235 = Icges(3,4) * t178;
t234 = Icges(3,4) * t181;
t227 = t170 * t182;
t183 = -pkin(7) - pkin(6);
t221 = t182 * t183;
t168 = t181 * pkin(2) + pkin(1);
t160 = t182 * t168;
t173 = t182 * pkin(6);
t219 = t179 * (t221 + t173 + (-pkin(1) + t168) * t179) + t182 * (-t182 * pkin(1) + t160 + (-pkin(6) - t183) * t179);
t188 = rSges(4,1) * t227 - rSges(4,2) * t228 + t179 * rSges(4,3);
t202 = rSges(4,1) * t170 - rSges(4,2) * t169;
t72 = t179 * (-t182 * rSges(4,3) + t202 * t179) + t182 * t188;
t218 = -t170 * rSges(6,2) + (t255 * t177 + t258 * t180) * t169;
t120 = -t170 * rSges(5,3) + (rSges(5,1) * t180 - rSges(5,2) * t177) * t169;
t152 = t169 * pkin(3) - t170 * pkin(8);
t217 = -t120 - t152;
t215 = pkin(3) * t227 + pkin(8) * t228;
t216 = t174 * (pkin(8) * t169 + t246) + t182 * t215;
t214 = t179 * rSges(3,3) + t182 * t244;
t213 = t174 + t175;
t212 = (t174 * t124 + (t194 * t182 + (-t123 + t193) * t179) * t182 + t256) * t179;
t211 = -t152 - t218;
t94 = t146 * rSges(5,1) - t145 * rSges(5,2) + rSges(5,3) * t228;
t151 = t169 * rSges(4,1) + t170 * rSges(4,2);
t208 = -t151 - t247;
t207 = -t152 - t247;
t206 = -t168 - t246;
t205 = -t179 * t183 + t160;
t201 = -t144 * rSges(5,1) + t143 * rSges(5,2);
t92 = rSges(5,3) * t230 - t201;
t45 = t179 * t92 + t182 * t94 + t216;
t204 = -t120 + t207;
t203 = -t243 + t244;
t200 = Icges(3,1) * t181 - t235;
t198 = -Icges(3,2) * t178 + t234;
t196 = Icges(3,5) * t181 - Icges(3,6) * t178;
t149 = Icges(4,2) * t170 + t233;
t150 = Icges(4,1) * t169 + t232;
t190 = -t149 * t169 + t150 * t170;
t189 = t207 - t218;
t22 = t237 * t179 + t254 * t182 + t216;
t187 = t205 + t215;
t186 = t251 * t182 + t212;
t185 = -(t240 - t241 + t238 - t239) * t170 / 0.2e1 + t259 * t249 + t260 * t248 + t257 * t230 / 0.2e1 + t256 * t228 / 0.2e1;
t148 = Icges(4,5) * t169 + Icges(4,6) * t170;
t184 = t238 / 0.2e1 + t240 / 0.2e1 - t239 / 0.2e1 - t241 / 0.2e1 + (t170 * t126 + t169 * t128 + t179 * t148 + t190 * t182 - t261) * t249 + (t170 * t125 + t169 * t127 - t182 * t148 + t190 * t179 - t262) * t248;
t159 = t182 * rSges(2,1) - t179 * rSges(2,2);
t158 = -t179 * rSges(2,1) - t182 * rSges(2,2);
t157 = t178 * rSges(3,1) + t181 * rSges(3,2);
t132 = Icges(3,3) * t179 + t196 * t182;
t131 = -Icges(3,3) * t182 + t196 * t179;
t122 = t208 * t182;
t121 = t208 * t179;
t108 = t179 * pkin(6) + (pkin(1) - t243) * t182 + t214;
t107 = t242 + t173 + (-pkin(1) - t203) * t179;
t103 = t188 + t205;
t102 = (rSges(4,3) - t183) * t182 + (-t168 - t202) * t179;
t97 = t217 * t182;
t96 = t217 * t179;
t95 = t182 * (-t182 * t243 + t214) + (t203 * t179 - t242) * t179;
t74 = t204 * t182;
t73 = t204 * t179;
t67 = t211 * t182;
t66 = t211 * t179;
t65 = t189 * t182;
t64 = t189 * t179;
t63 = t187 + t94;
t62 = -t221 + ((-rSges(5,3) - pkin(8)) * t169 + t206) * t179 + t201;
t61 = -t120 * t228 - t170 * t94;
t60 = t120 * t230 + t170 * t92;
t57 = t72 + t219;
t56 = (-t179 * t94 + t182 * t92) * t169;
t55 = t187 + t254;
t54 = -t221 + ((-rSges(6,2) - pkin(8)) * t169 + t206) * t179 + t267;
t44 = -t170 * t254 - t218 * t228;
t43 = t237 * t170 + t218 * t230;
t24 = t45 + t219;
t23 = (-t179 * t254 + t237 * t182) * t169;
t21 = t22 + t219;
t1 = [t181 * (Icges(3,2) * t181 + t235) + t178 * (Icges(3,1) * t178 + t234) + Icges(2,3) + (t150 - t226) * t169 + (t149 + t253) * t170 + m(6) * (t54 ^ 2 + t55 ^ 2) + m(5) * (t62 ^ 2 + t63 ^ 2) + m(4) * (t102 ^ 2 + t103 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2) + m(2) * (t158 ^ 2 + t159 ^ 2) + t252; (t174 / 0.2e1 + t175 / 0.2e1) * (Icges(3,5) * t178 + Icges(3,6) * t181) + t184 + m(6) * (t65 * t54 + t64 * t55) + m(5) * (t74 * t62 + t73 * t63) + m(4) * (t122 * t102 + t121 * t103) + m(3) * (-t107 * t182 - t108 * t179) * t157 + (t181 * (Icges(3,6) * t179 + t198 * t182) + t178 * (Icges(3,5) * t179 + t200 * t182)) * t249 + (t181 * (-Icges(3,6) * t182 + t198 * t179) + t178 * (-Icges(3,5) * t182 + t200 * t179)) * t248; m(6) * (t21 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(5) * (t24 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(4) * (t121 ^ 2 + t122 ^ 2 + t57 ^ 2) + m(3) * (t213 * t157 ^ 2 + t95 ^ 2) + t179 * t174 * t132 + t212 + (-t175 * t131 + (-t179 * t131 + t182 * t132) * t179 + t251) * t182; t184 + m(6) * (t67 * t54 + t66 * t55) + m(5) * (t97 * t62 + t96 * t63) + m(4) * (-t102 * t182 - t103 * t179) * t151; m(6) * (t22 * t21 + t66 * t64 + t67 * t65) + m(5) * (t45 * t24 + t96 * t73 + t97 * t74) + m(4) * (t72 * t57 + (-t121 * t179 - t122 * t182) * t151) + t186; m(6) * (t22 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t45 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(4) * (t213 * t151 ^ 2 + t72 ^ 2) + t186; -t245 * t170 + m(6) * (t43 * t54 + t44 * t55) + m(5) * (t60 * t62 + t61 * t63) + ((t42 / 0.2e1 + t40 / 0.2e1 + t53 / 0.2e1 + t52 / 0.2e1) * t182 + (t50 / 0.2e1 + t41 / 0.2e1 + t39 / 0.2e1 + t51 / 0.2e1) * t179) * t169; m(6) * (t23 * t21 + t43 * t65 + t44 * t64) + m(5) * (t56 * t24 + t60 * t74 + t61 * t73) + t185; m(6) * (t23 * t22 + t43 * t67 + t44 * t66) + m(5) * (t56 * t45 + t60 * t97 + t61 * t96) + t185; m(6) * (t23 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t56 ^ 2 + t60 ^ 2 + t61 ^ 2) + t245 * t170 ^ 2 + (t259 * t182 + t260 * t179 + ((-t40 - t42) * t182 + (-t39 - t41) * t179) * t170) * t169; m(6) * (t143 * t55 + t145 * t54); m(6) * (t143 * t64 + t145 * t65 + t21 * t231); m(6) * (t143 * t66 + t145 * t67 + t22 * t231); m(6) * (t143 * t44 + t145 * t43 + t23 * t231); m(6) * (t169 ^ 2 * t177 ^ 2 + t143 ^ 2 + t145 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
