% Calculate joint inertia matrix for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:14:51
% EndTime: 2019-12-05 17:15:13
% DurationCPUTime: 4.16s
% Computational Cost: add. (17976->421), mult. (31028->626), div. (0->0), fcn. (39272->12), ass. (0->213)
t203 = sin(pkin(10));
t205 = cos(pkin(10));
t209 = sin(qJ(2));
t206 = cos(pkin(5));
t212 = cos(qJ(2));
t241 = t206 * t212;
t191 = t203 * t209 - t205 * t241;
t193 = t203 * t241 + t205 * t209;
t204 = sin(pkin(5));
t244 = t204 * t212;
t242 = t206 * t209;
t192 = t203 * t212 + t205 * t242;
t240 = qJ(3) + qJ(4);
t202 = sin(t240);
t222 = cos(t240);
t248 = t204 * t205;
t176 = t192 * t222 - t202 * t248;
t207 = sin(qJ(5));
t210 = cos(qJ(5));
t148 = -t176 * t207 + t191 * t210;
t149 = t176 * t210 + t191 * t207;
t217 = t204 * t222;
t175 = t192 * t202 + t205 * t217;
t100 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t175;
t102 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t175;
t98 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t175;
t46 = t148 * t100 + t149 * t102 + t175 * t98;
t194 = -t203 * t242 + t205 * t212;
t249 = t203 * t204;
t178 = t194 * t222 + t202 * t249;
t150 = -t178 * t207 + t193 * t210;
t151 = t178 * t210 + t193 * t207;
t177 = t194 * t202 - t203 * t217;
t101 = Icges(6,4) * t151 + Icges(6,2) * t150 + Icges(6,6) * t177;
t103 = Icges(6,1) * t151 + Icges(6,4) * t150 + Icges(6,5) * t177;
t99 = Icges(6,5) * t151 + Icges(6,6) * t150 + Icges(6,3) * t177;
t47 = t148 * t101 + t149 * t103 + t175 * t99;
t190 = t206 * t202 + t209 * t217;
t179 = -t190 * t207 - t210 * t244;
t180 = t190 * t210 - t207 * t244;
t246 = t204 * t209;
t189 = t202 * t246 - t206 * t222;
t117 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t189;
t118 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t189;
t119 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t189;
t58 = t175 * t117 + t148 * t118 + t149 * t119;
t11 = t46 * t191 + t47 * t193 - t58 * t244;
t121 = Icges(5,5) * t176 - Icges(5,6) * t175 + Icges(5,3) * t191;
t123 = Icges(5,4) * t176 - Icges(5,2) * t175 + Icges(5,6) * t191;
t125 = Icges(5,1) * t176 - Icges(5,4) * t175 + Icges(5,5) * t191;
t67 = t191 * t121 - t175 * t123 + t176 * t125;
t122 = Icges(5,5) * t178 - Icges(5,6) * t177 + Icges(5,3) * t193;
t124 = Icges(5,4) * t178 - Icges(5,2) * t177 + Icges(5,6) * t193;
t126 = Icges(5,1) * t178 - Icges(5,4) * t177 + Icges(5,5) * t193;
t68 = t191 * t122 - t175 * t124 + t176 * t126;
t152 = Icges(5,5) * t190 - Icges(5,6) * t189 - Icges(5,3) * t244;
t153 = Icges(5,4) * t190 - Icges(5,2) * t189 - Icges(5,6) * t244;
t154 = Icges(5,1) * t190 - Icges(5,4) * t189 - Icges(5,5) * t244;
t84 = t191 * t152 - t175 * t153 + t176 * t154;
t268 = t67 * t191 + t68 * t193 - t84 * t244 + t11;
t48 = t150 * t100 + t151 * t102 + t177 * t98;
t49 = t150 * t101 + t151 * t103 + t177 * t99;
t59 = t177 * t117 + t150 * t118 + t151 * t119;
t12 = t48 * t191 + t49 * t193 - t59 * t244;
t69 = t193 * t121 - t177 * t123 + t178 * t125;
t70 = t193 * t122 - t177 * t124 + t178 * t126;
t85 = t193 * t152 - t177 * t153 + t178 * t154;
t267 = t69 * t191 + t70 * t193 - t85 * t244 + t12;
t15 = t58 * t206 + (t203 * t47 - t205 * t46) * t204;
t266 = t15 + t84 * t206 + (t203 * t68 - t205 * t67) * t204;
t16 = t59 * t206 + (t203 * t49 - t205 * t48) * t204;
t265 = t16 + t85 * t206 + (t203 * t70 - t205 * t69) * t204;
t54 = t179 * t100 + t180 * t102 + t189 * t98;
t55 = t179 * t101 + t180 * t103 + t189 * t99;
t64 = t189 * t117 + t179 * t118 + t180 * t119;
t21 = t54 * t191 + t55 * t193 - t64 * t244;
t77 = -t121 * t244 - t189 * t123 + t190 * t125;
t78 = -t122 * t244 - t189 * t124 + t190 * t126;
t89 = -t152 * t244 - t189 * t153 + t190 * t154;
t264 = t77 * t191 + t78 * t193 - t89 * t244 + t21;
t23 = t64 * t206 + (t203 * t55 - t205 * t54) * t204;
t263 = t23 + t89 * t206 + (t203 * t78 - t205 * t77) * t204;
t104 = t149 * rSges(6,1) + t148 * rSges(6,2) + t175 * rSges(6,3);
t239 = t176 * pkin(4) + t175 * pkin(9) + t104;
t120 = t180 * rSges(6,1) + t179 * rSges(6,2) + t189 * rSges(6,3);
t262 = t190 * pkin(4) + t189 * pkin(9) + t120;
t261 = t175 / 0.2e1;
t260 = t177 / 0.2e1;
t259 = t189 / 0.2e1;
t258 = t191 / 0.2e1;
t257 = t193 / 0.2e1;
t256 = t203 / 0.2e1;
t255 = -t205 / 0.2e1;
t254 = t206 / 0.2e1;
t211 = cos(qJ(3));
t253 = t211 * pkin(3);
t250 = t239 * t193;
t208 = sin(qJ(3));
t247 = t204 * t208;
t245 = t204 * t211;
t243 = t206 * t208;
t105 = t151 * rSges(6,1) + t150 * rSges(6,2) + t177 * rSges(6,3);
t238 = t178 * pkin(4) + t177 * pkin(9) + t105;
t229 = t205 * t247;
t129 = -pkin(3) * t229 + pkin(8) * t191 + t253 * t192;
t169 = pkin(3) * t243 + (-pkin(8) * t212 + t253 * t209) * t204;
t237 = t129 * t244 + t191 * t169;
t230 = t203 * t247;
t130 = pkin(3) * t230 + pkin(8) * t193 + t253 * t194;
t174 = t194 * pkin(2) + t193 * pkin(7);
t172 = t206 * t174;
t236 = t206 * t130 + t172;
t128 = t178 * rSges(5,1) - t177 * rSges(5,2) + t193 * rSges(5,3);
t234 = -t128 - t130;
t173 = t192 * pkin(2) + t191 * pkin(7);
t233 = -t129 - t173;
t127 = t176 * rSges(5,1) - t175 * rSges(5,2) + t191 * rSges(5,3);
t155 = t190 * rSges(5,1) - t189 * rSges(5,2) - rSges(5,3) * t244;
t95 = t127 * t244 + t191 * t155;
t232 = -t155 - t169;
t231 = t173 * t249 + t174 * t248;
t228 = -t130 - t238;
t227 = -t169 - t262;
t224 = -t244 / 0.2e1;
t223 = t268 * t191 + t267 * t193;
t195 = t206 * t211 - t208 * t246;
t196 = t209 * t245 + t243;
t167 = t196 * rSges(4,1) + t195 * rSges(4,2) - rSges(4,3) * t244;
t197 = (pkin(2) * t209 - pkin(7) * t212) * t204;
t221 = (-t167 - t197) * t204;
t61 = t262 * t191 + t239 * t244;
t220 = t129 * t249 + t130 * t248 + t231;
t219 = (-t197 + t232) * t204;
t18 = t54 * t175 + t55 * t177 + t64 * t189;
t3 = t46 * t175 + t47 * t177 + t58 * t189;
t4 = t48 * t175 + t49 * t177 + t59 * t189;
t218 = t11 * t261 + t12 * t260 + t18 * t224 + t21 * t259 + t4 * t257 + t3 * t258;
t216 = (-t197 + t227) * t204;
t215 = -t244 * t264 + t223;
t214 = t266 * t258 + t265 * t257 + t264 * t254 + t267 * t249 / 0.2e1 - t268 * t248 / 0.2e1 + t263 * t224;
t188 = t206 * rSges(3,3) + (rSges(3,1) * t209 + rSges(3,2) * t212) * t204;
t187 = Icges(3,5) * t206 + (Icges(3,1) * t209 + Icges(3,4) * t212) * t204;
t186 = Icges(3,6) * t206 + (Icges(3,4) * t209 + Icges(3,2) * t212) * t204;
t185 = Icges(3,3) * t206 + (Icges(3,5) * t209 + Icges(3,6) * t212) * t204;
t184 = t194 * t211 + t230;
t183 = -t194 * t208 + t203 * t245;
t182 = t192 * t211 - t229;
t181 = -t192 * t208 - t205 * t245;
t166 = Icges(4,1) * t196 + Icges(4,4) * t195 - Icges(4,5) * t244;
t165 = Icges(4,4) * t196 + Icges(4,2) * t195 - Icges(4,6) * t244;
t164 = Icges(4,5) * t196 + Icges(4,6) * t195 - Icges(4,3) * t244;
t163 = t194 * rSges(3,1) - t193 * rSges(3,2) + rSges(3,3) * t249;
t162 = t192 * rSges(3,1) - t191 * rSges(3,2) - rSges(3,3) * t248;
t161 = Icges(3,1) * t194 - Icges(3,4) * t193 + Icges(3,5) * t249;
t160 = Icges(3,1) * t192 - Icges(3,4) * t191 - Icges(3,5) * t248;
t159 = Icges(3,4) * t194 - Icges(3,2) * t193 + Icges(3,6) * t249;
t158 = Icges(3,4) * t192 - Icges(3,2) * t191 - Icges(3,6) * t248;
t157 = Icges(3,5) * t194 - Icges(3,6) * t193 + Icges(3,3) * t249;
t156 = Icges(3,5) * t192 - Icges(3,6) * t191 - Icges(3,3) * t248;
t142 = -t206 * t162 - t188 * t248;
t141 = t206 * t163 - t188 * t249;
t139 = t184 * rSges(4,1) + t183 * rSges(4,2) + t193 * rSges(4,3);
t138 = t182 * rSges(4,1) + t181 * rSges(4,2) + t191 * rSges(4,3);
t137 = Icges(4,1) * t184 + Icges(4,4) * t183 + Icges(4,5) * t193;
t136 = Icges(4,1) * t182 + Icges(4,4) * t181 + Icges(4,5) * t191;
t135 = Icges(4,4) * t184 + Icges(4,2) * t183 + Icges(4,6) * t193;
t134 = Icges(4,4) * t182 + Icges(4,2) * t181 + Icges(4,6) * t191;
t133 = Icges(4,5) * t184 + Icges(4,6) * t183 + Icges(4,3) * t193;
t132 = Icges(4,5) * t182 + Icges(4,6) * t181 + Icges(4,3) * t191;
t111 = (t162 * t203 + t163 * t205) * t204;
t110 = t193 * t129;
t109 = t193 * t127;
t107 = -t139 * t244 - t193 * t167;
t106 = t138 * t244 + t191 * t167;
t96 = -t128 * t244 - t193 * t155;
t93 = -t164 * t244 + t195 * t165 + t196 * t166;
t92 = t193 * t138 - t191 * t139;
t91 = (-t138 - t173) * t206 + t205 * t221;
t90 = t206 * t139 + t203 * t221 + t172;
t88 = t193 * t164 + t183 * t165 + t184 * t166;
t87 = t191 * t164 + t181 * t165 + t182 * t166;
t86 = -t191 * t128 + t109;
t83 = (t138 * t203 + t139 * t205) * t204 + t231;
t82 = -t133 * t244 + t195 * t135 + t196 * t137;
t81 = -t132 * t244 + t195 * t134 + t196 * t136;
t80 = t189 * t105 - t177 * t120;
t79 = -t189 * t104 + t175 * t120;
t76 = t232 * t193 + t234 * t244;
t75 = t95 + t237;
t74 = t193 * t133 + t183 * t135 + t184 * t137;
t73 = t193 * t132 + t183 * t134 + t184 * t136;
t72 = t191 * t133 + t181 * t135 + t182 * t137;
t71 = t191 * t132 + t181 * t134 + t182 * t136;
t66 = (-t127 + t233) * t206 + t205 * t219;
t65 = t206 * t128 + t203 * t219 + t236;
t63 = t177 * t104 - t175 * t105;
t62 = -t193 * t262 - t238 * t244;
t60 = t234 * t191 + t109 + t110;
t57 = (t127 * t203 + t128 * t205) * t204 + t220;
t56 = -t238 * t191 + t250;
t53 = t227 * t193 + t228 * t244;
t52 = t61 + t237;
t51 = (t233 - t239) * t206 + t205 * t216;
t50 = t203 * t216 + t238 * t206 + t236;
t45 = t228 * t191 + t110 + t250;
t44 = (t239 * t203 + t238 * t205) * t204 + t220;
t43 = t93 * t206 + (t203 * t82 - t205 * t81) * t204;
t42 = t81 * t191 + t82 * t193 - t93 * t244;
t39 = t88 * t206 + (t203 * t74 - t205 * t73) * t204;
t38 = t87 * t206 + (t203 * t72 - t205 * t71) * t204;
t35 = t73 * t191 + t74 * t193 - t88 * t244;
t34 = t71 * t191 + t72 * t193 - t87 * t244;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t111 + m(4) * t83 + m(5) * t57 + m(6) * t44; m(6) * (t44 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t57 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(4) * (t83 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(3) * (t111 ^ 2 + t141 ^ 2 + t142 ^ 2) + (t39 + (t157 * t249 - t193 * t159 + t194 * t161) * t249 + t265) * t249 + (-t38 + (-t156 * t248 - t191 * t158 + t192 * t160) * t248 + (-t156 * t249 + t157 * t248 + t193 * t158 + t191 * t159 - t194 * t160 - t192 * t161) * t249 - t266) * t248 + (-(-t185 * t248 - t191 * t186 + t192 * t187) * t248 + (t185 * t249 - t193 * t186 + t194 * t187) * t249 + t43 + ((t159 * t212 + t161 * t209) * t203 - (t158 * t212 + t160 * t209) * t205) * t204 ^ 2 + ((-t156 * t205 + t157 * t203 + t186 * t212 + t187 * t209) * t204 + t206 * t185) * t206 + t263) * t206; m(4) * t92 + m(5) * t60 + m(6) * t45; t214 + m(6) * (t45 * t44 + t53 * t50 + t52 * t51) + m(5) * (t60 * t57 + t76 * t65 + t75 * t66) + m(4) * (t106 * t91 + t107 * t90 + t92 * t83) + (-t212 * t43 / 0.2e1 + t35 * t256 + t34 * t255) * t204 + t38 * t258 + t39 * t257 + t42 * t254; t191 * t34 + t193 * t35 + (-t42 - t264) * t244 + m(6) * (t45 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(5) * (t60 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(4) * (t106 ^ 2 + t107 ^ 2 + t92 ^ 2) + t223; m(5) * t86 + m(6) * t56; t214 + m(6) * (t56 * t44 + t62 * t50 + t61 * t51) + m(5) * (t86 * t57 + t96 * t65 + t95 * t66); m(6) * (t56 * t45 + t61 * t52 + t62 * t53) + m(5) * (t86 * t60 + t95 * t75 + t96 * t76) + t215; m(6) * (t56 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t86 ^ 2 + t95 ^ 2 + t96 ^ 2) + t215; m(6) * t63; t15 * t261 + t16 * t260 + m(6) * (t63 * t44 + t80 * t50 + t79 * t51) + t23 * t259 + t18 * t254 + (t3 * t255 + t4 * t256) * t204; m(6) * (t63 * t45 + t79 * t52 + t80 * t53) + t218; m(6) * (t63 * t56 + t79 * t61 + t80 * t62) + t218; m(6) * (t63 ^ 2 + t79 ^ 2 + t80 ^ 2) + t177 * t4 + t175 * t3 + t189 * t18;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
