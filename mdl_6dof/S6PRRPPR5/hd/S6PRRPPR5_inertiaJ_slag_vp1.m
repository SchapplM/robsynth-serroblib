% Calculate joint inertia matrix for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:13
% EndTime: 2019-03-08 21:18:22
% DurationCPUTime: 5.09s
% Computational Cost: add. (13765->540), mult. (30858->784), div. (0->0), fcn. (38980->12), ass. (0->242)
t268 = m(6) + m(7);
t269 = m(5) + t268;
t215 = sin(pkin(10));
t218 = cos(pkin(10));
t222 = cos(qJ(2));
t219 = cos(pkin(6));
t221 = sin(qJ(2));
t252 = t219 * t221;
t200 = t215 * t222 + t218 * t252;
t216 = sin(pkin(6));
t262 = sin(qJ(3));
t233 = t216 * t262;
t263 = cos(qJ(3));
t191 = t200 * t263 - t218 * t233;
t202 = -t215 * t252 + t218 * t222;
t193 = t202 * t263 + t215 * t233;
t234 = t216 * t263;
t204 = t219 * t262 + t221 * t234;
t192 = t202 * t262 - t215 * t234;
t251 = t219 * t222;
t201 = t215 * t251 + t218 * t221;
t213 = pkin(11) + qJ(6);
t211 = sin(t213);
t212 = cos(t213);
t150 = t192 * t212 - t201 * t211;
t151 = t192 * t211 + t201 * t212;
t190 = t200 * t262 + t218 * t234;
t199 = t215 * t221 - t218 * t251;
t148 = t190 * t212 - t199 * t211;
t149 = t190 * t211 + t199 * t212;
t90 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t191;
t92 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t191;
t94 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t191;
t38 = t150 * t92 + t151 * t94 + t193 * t90;
t91 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t193;
t93 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t193;
t95 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t193;
t39 = t150 * t93 + t151 * t95 + t193 * t91;
t203 = -t219 * t263 + t221 * t233;
t253 = t216 * t222;
t186 = t203 * t212 + t211 * t253;
t187 = t203 * t211 - t212 * t253;
t111 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t204;
t112 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t204;
t113 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t204;
t53 = t111 * t193 + t112 * t150 + t113 * t151;
t2 = t191 * t38 + t193 * t39 + t204 * t53;
t267 = t2 / 0.2e1;
t266 = t191 / 0.2e1;
t265 = t193 / 0.2e1;
t264 = t204 / 0.2e1;
t217 = cos(pkin(11));
t261 = pkin(5) * t217;
t214 = sin(pkin(11));
t258 = t190 * t214;
t96 = rSges(7,1) * t149 + rSges(7,2) * t148 + rSges(7,3) * t191;
t260 = pkin(5) * t258 + pkin(9) * t191 + t261 * t199 + t96;
t257 = t192 * t214;
t97 = rSges(7,1) * t151 + rSges(7,2) * t150 + rSges(7,3) * t193;
t259 = pkin(5) * t257 + pkin(9) * t193 + t261 * t201 + t97;
t256 = t203 * t214;
t255 = t215 * t216;
t254 = t216 * t218;
t114 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t204;
t249 = pkin(5) * t256 + pkin(9) * t204 - t261 * t253 + t114;
t131 = rSges(5,1) * t201 - rSges(5,2) * t193 + rSges(5,3) * t192;
t147 = pkin(3) * t193 + qJ(4) * t192;
t248 = -t131 - t147;
t146 = pkin(3) * t191 + qJ(4) * t190;
t135 = t201 * t146;
t161 = t199 * pkin(4) + t191 * qJ(5);
t247 = t201 * t161 + t135;
t185 = pkin(3) * t204 + qJ(4) * t203;
t246 = t146 * t253 + t199 * t185;
t184 = pkin(2) * t202 + pkin(8) * t201;
t182 = t219 * t184;
t245 = t219 * t147 + t182;
t183 = pkin(2) * t200 + pkin(8) * t199;
t244 = -t146 - t183;
t162 = t201 * pkin(4) + t193 * qJ(5);
t243 = -t147 - t162;
t177 = -rSges(5,1) * t253 - t204 * rSges(5,2) + t203 * rSges(5,3);
t242 = -t177 - t185;
t241 = t183 * t255 + t184 * t254;
t194 = -pkin(4) * t253 + t204 * qJ(5);
t240 = -t185 - t194;
t157 = t192 * t217 - t201 * t214;
t158 = t201 * t217 + t257;
t107 = rSges(6,1) * t158 + rSges(6,2) * t157 + rSges(6,3) * t193;
t238 = -t107 + t243;
t188 = t203 * t217 + t214 * t253;
t189 = -t217 * t253 + t256;
t134 = rSges(6,1) * t189 + rSges(6,2) * t188 + rSges(6,3) * t204;
t237 = -t134 + t240;
t236 = t219 * t162 + t245;
t235 = -t161 + t244;
t178 = t204 * rSges(4,1) - t203 * rSges(4,2) - rSges(4,3) * t253;
t205 = (pkin(2) * t221 - pkin(8) * t222) * t216;
t232 = (-t178 - t205) * t216;
t230 = t243 - t259;
t229 = t240 - t249;
t228 = t146 * t255 + t147 * t254 + t241;
t227 = t161 * t253 + t199 * t194 + t246;
t226 = (-t205 + t242) * t216;
t225 = (-t205 + t237) * t216;
t224 = t161 * t255 + t162 * t254 + t228;
t223 = (-t205 + t229) * t216;
t198 = t219 * rSges(3,3) + (rSges(3,1) * t221 + rSges(3,2) * t222) * t216;
t197 = Icges(3,5) * t219 + (Icges(3,1) * t221 + Icges(3,4) * t222) * t216;
t196 = Icges(3,6) * t219 + (Icges(3,4) * t221 + Icges(3,2) * t222) * t216;
t195 = Icges(3,3) * t219 + (Icges(3,5) * t221 + Icges(3,6) * t222) * t216;
t176 = Icges(4,1) * t204 - Icges(4,4) * t203 - Icges(4,5) * t253;
t175 = Icges(4,4) * t204 - Icges(4,2) * t203 - Icges(4,6) * t253;
t174 = Icges(4,5) * t204 - Icges(4,6) * t203 - Icges(4,3) * t253;
t173 = -Icges(5,1) * t253 - Icges(5,4) * t204 + Icges(5,5) * t203;
t172 = -Icges(5,4) * t253 - Icges(5,2) * t204 + Icges(5,6) * t203;
t171 = -Icges(5,5) * t253 - Icges(5,6) * t204 + Icges(5,3) * t203;
t170 = rSges(3,1) * t202 - rSges(3,2) * t201 + rSges(3,3) * t255;
t169 = rSges(3,1) * t200 - rSges(3,2) * t199 - rSges(3,3) * t254;
t168 = Icges(3,1) * t202 - Icges(3,4) * t201 + Icges(3,5) * t255;
t167 = Icges(3,1) * t200 - Icges(3,4) * t199 - Icges(3,5) * t254;
t166 = Icges(3,4) * t202 - Icges(3,2) * t201 + Icges(3,6) * t255;
t165 = Icges(3,4) * t200 - Icges(3,2) * t199 - Icges(3,6) * t254;
t164 = Icges(3,5) * t202 - Icges(3,6) * t201 + Icges(3,3) * t255;
t163 = Icges(3,5) * t200 - Icges(3,6) * t199 - Icges(3,3) * t254;
t156 = t199 * t217 + t258;
t155 = t190 * t217 - t199 * t214;
t140 = -t169 * t219 - t198 * t254;
t139 = t170 * t219 - t198 * t255;
t133 = rSges(4,1) * t193 - rSges(4,2) * t192 + rSges(4,3) * t201;
t132 = rSges(4,1) * t191 - rSges(4,2) * t190 + rSges(4,3) * t199;
t130 = rSges(5,1) * t199 - rSges(5,2) * t191 + rSges(5,3) * t190;
t129 = Icges(6,1) * t189 + Icges(6,4) * t188 + Icges(6,5) * t204;
t128 = Icges(6,4) * t189 + Icges(6,2) * t188 + Icges(6,6) * t204;
t127 = Icges(6,5) * t189 + Icges(6,6) * t188 + Icges(6,3) * t204;
t126 = Icges(4,1) * t193 - Icges(4,4) * t192 + Icges(4,5) * t201;
t125 = Icges(4,1) * t191 - Icges(4,4) * t190 + Icges(4,5) * t199;
t124 = Icges(5,1) * t201 - Icges(5,4) * t193 + Icges(5,5) * t192;
t123 = Icges(5,1) * t199 - Icges(5,4) * t191 + Icges(5,5) * t190;
t122 = Icges(4,4) * t193 - Icges(4,2) * t192 + Icges(4,6) * t201;
t121 = Icges(4,4) * t191 - Icges(4,2) * t190 + Icges(4,6) * t199;
t120 = Icges(5,4) * t201 - Icges(5,2) * t193 + Icges(5,6) * t192;
t119 = Icges(5,4) * t199 - Icges(5,2) * t191 + Icges(5,6) * t190;
t118 = Icges(4,5) * t193 - Icges(4,6) * t192 + Icges(4,3) * t201;
t117 = Icges(4,5) * t191 - Icges(4,6) * t190 + Icges(4,3) * t199;
t116 = Icges(5,5) * t201 - Icges(5,6) * t193 + Icges(5,3) * t192;
t115 = Icges(5,5) * t199 - Icges(5,6) * t191 + Icges(5,3) * t190;
t110 = (t169 * t215 + t170 * t218) * t216;
t106 = rSges(6,1) * t156 + rSges(6,2) * t155 + rSges(6,3) * t191;
t105 = Icges(6,1) * t158 + Icges(6,4) * t157 + Icges(6,5) * t193;
t104 = Icges(6,1) * t156 + Icges(6,4) * t155 + Icges(6,5) * t191;
t103 = Icges(6,4) * t158 + Icges(6,2) * t157 + Icges(6,6) * t193;
t102 = Icges(6,4) * t156 + Icges(6,2) * t155 + Icges(6,6) * t191;
t101 = Icges(6,5) * t158 + Icges(6,6) * t157 + Icges(6,3) * t193;
t100 = Icges(6,5) * t156 + Icges(6,6) * t155 + Icges(6,3) * t191;
t99 = -t133 * t253 - t201 * t178;
t98 = t132 * t253 + t199 * t178;
t88 = -t174 * t253 - t203 * t175 + t204 * t176;
t87 = t203 * t171 - t204 * t172 - t173 * t253;
t86 = t132 * t201 - t133 * t199;
t85 = (-t132 - t183) * t219 + t218 * t232;
t84 = t219 * t133 + t215 * t232 + t182;
t83 = t174 * t201 - t175 * t192 + t176 * t193;
t82 = t174 * t199 - t175 * t190 + t176 * t191;
t81 = t171 * t192 - t172 * t193 + t173 * t201;
t80 = t171 * t190 - t172 * t191 + t173 * t199;
t79 = (t132 * t215 + t133 * t218) * t216 + t241;
t78 = t242 * t201 + t248 * t253;
t77 = t130 * t253 + t199 * t177 + t246;
t76 = -t118 * t253 - t203 * t122 + t204 * t126;
t75 = -t117 * t253 - t203 * t121 + t204 * t125;
t74 = t203 * t116 - t204 * t120 - t124 * t253;
t73 = t203 * t115 - t204 * t119 - t123 * t253;
t72 = -t114 * t193 + t204 * t97;
t71 = t114 * t191 - t204 * t96;
t70 = (-t130 + t244) * t219 + t218 * t226;
t69 = t219 * t131 + t215 * t226 + t245;
t68 = t127 * t204 + t128 * t188 + t129 * t189;
t67 = t118 * t201 - t122 * t192 + t126 * t193;
t66 = t117 * t201 - t121 * t192 + t125 * t193;
t65 = t118 * t199 - t122 * t190 + t126 * t191;
t64 = t117 * t199 - t121 * t190 + t125 * t191;
t63 = t116 * t192 - t120 * t193 + t124 * t201;
t62 = t115 * t192 - t119 * t193 + t123 * t201;
t61 = t116 * t190 - t120 * t191 + t124 * t199;
t60 = t115 * t190 - t119 * t191 + t123 * t199;
t59 = t111 * t204 + t112 * t186 + t113 * t187;
t58 = t201 * t130 + t199 * t248 + t135;
t57 = -t191 * t97 + t193 * t96;
t56 = t127 * t193 + t128 * t157 + t129 * t158;
t55 = t127 * t191 + t128 * t155 + t129 * t156;
t54 = (t130 * t215 + t131 * t218) * t216 + t228;
t52 = t111 * t191 + t112 * t148 + t113 * t149;
t51 = t201 * t237 + t238 * t253;
t50 = t106 * t253 + t199 * t134 + t227;
t49 = (-t106 + t235) * t219 + t218 * t225;
t48 = t219 * t107 + t215 * t225 + t236;
t47 = t101 * t204 + t103 * t188 + t105 * t189;
t46 = t100 * t204 + t102 * t188 + t104 * t189;
t45 = t186 * t93 + t187 * t95 + t204 * t91;
t44 = t186 * t92 + t187 * t94 + t204 * t90;
t43 = t101 * t193 + t103 * t157 + t105 * t158;
t42 = t100 * t193 + t102 * t157 + t104 * t158;
t41 = t101 * t191 + t103 * t155 + t105 * t156;
t40 = t100 * t191 + t102 * t155 + t104 * t156;
t37 = t148 * t93 + t149 * t95 + t191 * t91;
t36 = t148 * t92 + t149 * t94 + t191 * t90;
t35 = t201 * t106 + t199 * t238 + t247;
t34 = (t106 * t215 + t107 * t218) * t216 + t224;
t33 = t201 * t229 + t230 * t253;
t32 = t249 * t199 + t260 * t253 + t227;
t31 = (t235 - t260) * t219 + t218 * t223;
t30 = t215 * t223 + t259 * t219 + t236;
t29 = t88 * t219 + (t215 * t76 - t218 * t75) * t216;
t28 = t87 * t219 + (t215 * t74 - t218 * t73) * t216;
t27 = t75 * t199 + t76 * t201 - t88 * t253;
t26 = t73 * t199 + t74 * t201 - t87 * t253;
t25 = t83 * t219 + (t215 * t67 - t218 * t66) * t216;
t24 = t82 * t219 + (t215 * t65 - t218 * t64) * t216;
t23 = t81 * t219 + (t215 * t63 - t218 * t62) * t216;
t22 = t80 * t219 + (t215 * t61 - t218 * t60) * t216;
t21 = t66 * t199 + t67 * t201 - t83 * t253;
t20 = t64 * t199 + t65 * t201 - t82 * t253;
t19 = t62 * t199 + t63 * t201 - t81 * t253;
t18 = t60 * t199 + t61 * t201 - t80 * t253;
t17 = t230 * t199 + t260 * t201 + t247;
t16 = (t260 * t215 + t259 * t218) * t216 + t224;
t15 = t68 * t219 + (t215 * t47 - t218 * t46) * t216;
t14 = t46 * t199 + t47 * t201 - t68 * t253;
t13 = t59 * t219 + (t215 * t45 - t218 * t44) * t216;
t12 = t44 * t199 + t45 * t201 - t59 * t253;
t11 = t191 * t44 + t193 * t45 + t204 * t59;
t10 = t56 * t219 + (t215 * t43 - t218 * t42) * t216;
t9 = t55 * t219 + (t215 * t41 - t218 * t40) * t216;
t8 = t42 * t199 + t43 * t201 - t56 * t253;
t7 = t40 * t199 + t41 * t201 - t55 * t253;
t6 = t53 * t219 + (t215 * t39 - t218 * t38) * t216;
t5 = t52 * t219 + (t215 * t37 - t218 * t36) * t216;
t4 = t38 * t199 + t39 * t201 - t53 * t253;
t3 = t36 * t199 + t37 * t201 - t52 * t253;
t1 = t191 * t36 + t193 * t37 + t204 * t52;
t89 = [m(2) + m(3) + m(4) + t269; m(3) * t110 + m(4) * t79 + m(5) * t54 + m(6) * t34 + m(7) * t16; m(7) * (t16 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t34 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t54 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(4) * (t79 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(3) * (t110 ^ 2 + t139 ^ 2 + t140 ^ 2) + (t6 + t25 + t23 + t10 + (t164 * t255 - t166 * t201 + t168 * t202) * t255) * t255 + (-t5 - t24 - t22 - t9 + (-t163 * t254 - t165 * t199 + t167 * t200) * t254 + (-t163 * t255 + t164 * t254 + t165 * t201 + t166 * t199 - t167 * t202 - t168 * t200) * t255) * t254 + (-(-t195 * t254 - t199 * t196 + t200 * t197) * t254 + (t195 * t255 - t201 * t196 + t202 * t197) * t255 + t13 + t29 + t28 + t15 + ((t166 * t222 + t168 * t221) * t215 - (t165 * t222 + t167 * t221) * t218) * t216 ^ 2 + ((-t163 * t218 + t164 * t215 + t196 * t222 + t197 * t221) * t216 + t219 * t195) * t219) * t219; m(4) * t86 + m(5) * t58 + m(6) * t35 + m(7) * t17; (t26 / 0.2e1 + t27 / 0.2e1 + t14 / 0.2e1 + t12 / 0.2e1) * t219 + (t25 / 0.2e1 + t23 / 0.2e1 + t10 / 0.2e1 + t6 / 0.2e1) * t201 + (t22 / 0.2e1 + t24 / 0.2e1 + t9 / 0.2e1 + t5 / 0.2e1) * t199 + m(7) * (t16 * t17 + t30 * t33 + t31 * t32) + m(6) * (t34 * t35 + t48 * t51 + t49 * t50) + m(5) * (t54 * t58 + t69 * t78 + t70 * t77) + m(4) * (t79 * t86 + t84 * t99 + t85 * t98) + ((-t13 / 0.2e1 - t15 / 0.2e1 - t28 / 0.2e1 - t29 / 0.2e1) * t222 + (-t3 / 0.2e1 - t7 / 0.2e1 - t18 / 0.2e1 - t20 / 0.2e1) * t218 + (t4 / 0.2e1 + t8 / 0.2e1 + t19 / 0.2e1 + t21 / 0.2e1) * t215) * t216; (-t12 - t14 - t26 - t27) * t253 + (t19 + t8 + t21 + t4) * t201 + (t18 + t7 + t20 + t3) * t199 + m(5) * (t58 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t86 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(7) * (t17 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t35 ^ 2 + t50 ^ 2 + t51 ^ 2); t203 * t269; m(7) * (t16 * t203 + t190 * t30 + t192 * t31) + m(6) * (t190 * t48 + t192 * t49 + t203 * t34) + m(5) * (t190 * t69 + t192 * t70 + t203 * t54); m(7) * (t17 * t203 + t190 * t33 + t192 * t32) + m(6) * (t190 * t51 + t192 * t50 + t203 * t35) + m(5) * (t190 * t78 + t192 * t77 + t203 * t58); (t190 ^ 2 + t192 ^ 2 + t203 ^ 2) * t269; t204 * t268; m(7) * (t16 * t204 + t191 * t30 + t193 * t31) + m(6) * (t191 * t48 + t193 * t49 + t204 * t34); m(7) * (t17 * t204 + t191 * t33 + t193 * t32) + m(6) * (t191 * t51 + t193 * t50 + t204 * t35); (t190 * t191 + t192 * t193 + t203 * t204) * t268; (t191 ^ 2 + t193 ^ 2 + t204 ^ 2) * t268; m(7) * t57; m(7) * (t16 * t57 + t30 * t72 + t31 * t71) + t6 * t265 + t219 * t11 / 0.2e1 + t5 * t266 + t13 * t264 + (t215 * t267 - t218 * t1 / 0.2e1) * t216; m(7) * (t17 * t57 + t32 * t71 + t33 * t72) + t3 * t266 + t4 * t265 + t201 * t267 + t199 * t1 / 0.2e1 + t12 * t264 - t11 * t253 / 0.2e1; m(7) * (t190 * t72 + t192 * t71 + t203 * t57); m(7) * (t191 * t72 + t193 * t71 + t204 * t57); t193 * t2 + t191 * t1 + t204 * t11 + m(7) * (t57 ^ 2 + t71 ^ 2 + t72 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t89(1) t89(2) t89(4) t89(7) t89(11) t89(16); t89(2) t89(3) t89(5) t89(8) t89(12) t89(17); t89(4) t89(5) t89(6) t89(9) t89(13) t89(18); t89(7) t89(8) t89(9) t89(10) t89(14) t89(19); t89(11) t89(12) t89(13) t89(14) t89(15) t89(20); t89(16) t89(17) t89(18) t89(19) t89(20) t89(21);];
Mq  = res;
