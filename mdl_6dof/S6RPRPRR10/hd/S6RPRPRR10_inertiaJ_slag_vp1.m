% Calculate joint inertia matrix for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:40
% EndTime: 2019-03-09 04:07:45
% DurationCPUTime: 2.55s
% Computational Cost: add. (6918->416), mult. (7765->602), div. (0->0), fcn. (8241->10), ass. (0->204)
t188 = sin(qJ(1));
t267 = -t188 / 0.2e1;
t190 = cos(qJ(1));
t253 = t190 / 0.2e1;
t266 = rSges(5,3) + qJ(4);
t185 = cos(pkin(10));
t169 = t185 * pkin(4) + pkin(3);
t186 = -pkin(8) - qJ(4);
t189 = cos(qJ(3));
t228 = t188 * t189;
t187 = sin(qJ(3));
t236 = t187 * t188;
t184 = sin(pkin(10));
t237 = t184 * t190;
t218 = pkin(4) * t237 + t169 * t236 + t186 * t228;
t180 = pkin(10) + qJ(5);
t171 = cos(t180);
t170 = sin(t180);
t232 = t188 * t170;
t121 = t171 * t190 - t187 * t232;
t231 = t188 * t171;
t122 = t170 * t190 + t187 * t231;
t77 = t122 * rSges(6,1) + t121 * rSges(6,2) - rSges(6,3) * t228;
t265 = -t77 - t218;
t235 = t187 * t190;
t166 = pkin(3) * t235;
t175 = t190 * qJ(2);
t229 = t188 * t185;
t140 = t184 * t235 + t229;
t230 = t188 * t184;
t141 = -t185 * t235 + t230;
t206 = -t141 * rSges(5,1) - t140 * rSges(5,2);
t213 = t266 * t189;
t259 = -pkin(1) - pkin(7);
t58 = t188 * t259 - t190 * t213 + t166 + t175 + t206;
t223 = t190 * pkin(1) + t188 * qJ(2);
t216 = t190 * pkin(7) + t223;
t138 = t185 * t190 - t187 * t230;
t139 = t187 * t229 + t237;
t260 = -t139 * rSges(5,1) - t138 * rSges(5,2) - pkin(3) * t236;
t59 = -t188 * t213 + t216 - t260;
t264 = m(5) * (t188 * t58 - t190 * t59);
t227 = t189 * t190;
t224 = t169 * t235 + t186 * t227;
t251 = pkin(4) * t184;
t123 = t170 * t235 + t231;
t124 = -t171 * t235 + t232;
t78 = t124 * rSges(6,1) + t123 * rSges(6,2) + rSges(6,3) * t227;
t50 = t175 + (-t251 + t259) * t188 - t78 + t224;
t51 = t216 - t265;
t263 = m(6) * (t188 * t50 - t190 * t51);
t149 = pkin(5) * t170 + t251;
t179 = -pkin(9) + t186;
t172 = qJ(6) + t180;
t167 = sin(t172);
t168 = cos(t172);
t233 = t188 * t168;
t113 = t167 * t235 + t233;
t234 = t188 * t167;
t114 = -t168 * t235 + t234;
t205 = -t114 * rSges(7,1) - t113 * rSges(7,2);
t145 = pkin(5) * t171 + t169;
t238 = t145 * t187;
t46 = t175 + (t238 + (-rSges(7,3) + t179) * t189) * t190 + (-t149 + t259) * t188 + t205;
t219 = t145 * t236 + t190 * t149 + t179 * t228;
t111 = t168 * t190 - t187 * t234;
t112 = t167 * t190 + t187 * t233;
t69 = t112 * rSges(7,1) + t111 * rSges(7,2) - rSges(7,3) * t228;
t47 = t69 + t216 + t219;
t262 = m(7) * (t188 * t46 - t190 * t47);
t261 = (rSges(4,1) * t187 + rSges(4,2) * t189) * t190;
t181 = t188 ^ 2;
t183 = t190 ^ 2;
t63 = Icges(7,5) * t112 + Icges(7,6) * t111 - Icges(7,3) * t228;
t65 = Icges(7,4) * t112 + Icges(7,2) * t111 - Icges(7,6) * t228;
t67 = Icges(7,1) * t112 + Icges(7,4) * t111 - Icges(7,5) * t228;
t31 = t187 * t63 + (-t167 * t65 + t168 * t67) * t189;
t64 = Icges(7,5) * t114 + Icges(7,6) * t113 + Icges(7,3) * t227;
t66 = Icges(7,4) * t114 + Icges(7,2) * t113 + Icges(7,6) * t227;
t68 = Icges(7,1) * t114 + Icges(7,4) * t113 + Icges(7,5) * t227;
t32 = t187 * t64 + (-t167 * t66 + t168 * t68) * t189;
t101 = Icges(7,6) * t187 + (Icges(7,4) * t168 - Icges(7,2) * t167) * t189;
t240 = t101 * t167;
t100 = Icges(7,3) * t187 + (Icges(7,5) * t168 - Icges(7,6) * t167) * t189;
t102 = Icges(7,5) * t187 + (Icges(7,1) * t168 - Icges(7,4) * t167) * t189;
t248 = t189 * t168 * t102 + t187 * t100;
t48 = (-t189 * t240 + t248) * t187;
t21 = t113 * t65 + t114 * t67 + t227 * t63;
t22 = t113 * t66 + t114 * t68 + t227 * t64;
t38 = t100 * t227 + t113 * t101 + t114 * t102;
t5 = t38 * t187 + (-t188 * t21 + t190 * t22) * t189;
t258 = t5 * t227 + t187 * (t48 + (-t188 * t31 + t190 * t32) * t189);
t118 = Icges(5,6) * t187 + (Icges(5,4) * t185 - Icges(5,2) * t184) * t189;
t257 = t118 / 0.2e1;
t119 = Icges(5,5) * t187 + (Icges(5,1) * t185 - Icges(5,4) * t184) * t189;
t256 = t119 / 0.2e1;
t255 = t187 / 0.2e1;
t254 = t188 / 0.2e1;
t157 = rSges(4,1) * t189 - rSges(4,2) * t187;
t252 = m(4) * t157;
t56 = -t218 + t219;
t250 = -t56 - t69;
t70 = rSges(7,3) * t227 - t205;
t249 = (-t179 * t189 - t238) * t190 + (t149 - t251) * t188 + t224 + t70;
t105 = Icges(6,3) * t187 + (Icges(6,5) * t171 - Icges(6,6) * t170) * t189;
t107 = Icges(6,5) * t187 + (Icges(6,1) * t171 - Icges(6,4) * t170) * t189;
t247 = t189 * t171 * t107 + t187 * t105;
t103 = t187 * rSges(7,3) + (rSges(7,1) * t168 - rSges(7,2) * t167) * t189;
t52 = t103 * t228 + t187 * t69;
t88 = (t145 - t169) * t189 + (-t179 + t186) * t187;
t246 = t103 + t88;
t209 = qJ(4) * t227 - t166;
t134 = t190 * t209;
t245 = t134 + t190 * (pkin(4) * t230 - t209 - t224);
t104 = (-pkin(3) + t169) * t189 + (-qJ(4) - t186) * t187;
t156 = t189 * pkin(3) + t187 * qJ(4);
t144 = t188 * t156;
t243 = t188 * t104 + t144;
t242 = Icges(4,4) * t187;
t241 = Icges(4,4) * t189;
t106 = Icges(6,6) * t187 + (Icges(6,4) * t171 - Icges(6,2) * t170) * t189;
t239 = t106 * t170;
t226 = -t104 - t156;
t222 = t181 + t183;
t71 = Icges(6,5) * t122 + Icges(6,6) * t121 - Icges(6,3) * t228;
t73 = Icges(6,4) * t122 + Icges(6,2) * t121 - Icges(6,6) * t228;
t75 = Icges(6,1) * t122 + Icges(6,4) * t121 - Icges(6,5) * t228;
t33 = t187 * t71 + (-t170 * t73 + t171 * t75) * t189;
t41 = -t105 * t228 + t106 * t121 + t107 * t122;
t221 = -t41 / 0.2e1 - t33 / 0.2e1;
t72 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t227;
t74 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t227;
t76 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t227;
t34 = t187 * t72 + (-t170 * t74 + t171 * t76) * t189;
t42 = t105 * t227 + t123 * t106 + t124 * t107;
t220 = t42 / 0.2e1 + t34 / 0.2e1;
t217 = rSges(4,1) * t236 + rSges(4,2) * t228 + t190 * rSges(4,3);
t215 = -t228 / 0.2e1;
t214 = t227 / 0.2e1;
t212 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t19 = t111 * t65 + t112 * t67 - t228 * t63;
t20 = t111 * t66 + t112 * t68 - t228 * t64;
t11 = t20 * t188 + t19 * t190;
t12 = t22 * t188 + t190 * t21;
t37 = -t100 * t228 + t101 * t111 + t102 * t112;
t4 = t37 * t187 + (-t188 * t19 + t190 * t20) * t189;
t211 = t11 * t215 + t12 * t214 + t4 * t253 + t5 * t254 + (t32 * t188 + t31 * t190) * t255;
t210 = t48 + (t31 + t37) * t215 + (t32 + t38) * t214;
t208 = -t228 * t4 + t258;
t23 = t187 * t56 + t228 * t88 + t52;
t96 = t103 * t227;
t24 = -t187 * t249 + t227 * t88 + t96;
t204 = t24 * t188 - t190 * t23;
t44 = t188 * t246 + t243;
t45 = (t226 - t246) * t190;
t203 = t44 * t188 - t190 * t45;
t53 = -t187 * t70 + t96;
t200 = t53 * t188 - t190 * t52;
t110 = t187 * rSges(6,3) + (rSges(6,1) * t171 - rSges(6,2) * t170) * t189;
t54 = t110 * t228 + t187 * t77;
t55 = t110 * t227 - t187 * t78;
t199 = t55 * t188 - t190 * t54;
t60 = t110 * t188 + t243;
t61 = (-t110 + t226) * t190;
t197 = t60 * t188 - t190 * t61;
t120 = t187 * rSges(5,3) + (rSges(5,1) * t185 - rSges(5,2) * t184) * t189;
t89 = t120 * t188 + t144;
t90 = (-t120 - t156) * t190;
t196 = t89 * t188 - t190 * t90;
t195 = Icges(4,1) * t187 + t241;
t194 = Icges(4,2) * t189 + t242;
t193 = Icges(4,5) * t187 + Icges(4,6) * t189;
t158 = rSges(2,1) * t190 - t188 * rSges(2,2);
t155 = -t188 * rSges(2,1) - rSges(2,2) * t190;
t151 = Icges(4,5) * t189 - Icges(4,6) * t187;
t133 = -rSges(3,2) * t190 + t188 * rSges(3,3) + t223;
t132 = t190 * rSges(3,3) + t175 + (rSges(3,2) - pkin(1)) * t188;
t126 = Icges(4,3) * t188 - t190 * t193;
t125 = Icges(4,3) * t190 + t188 * t193;
t94 = t216 + t217;
t93 = t175 + t261 + (-rSges(4,3) + t259) * t188;
t86 = Icges(5,1) * t141 + Icges(5,4) * t140 + Icges(5,5) * t227;
t85 = Icges(5,1) * t139 + Icges(5,4) * t138 - Icges(5,5) * t228;
t84 = Icges(5,4) * t141 + Icges(5,2) * t140 + Icges(5,6) * t227;
t83 = Icges(5,4) * t139 + Icges(5,2) * t138 - Icges(5,6) * t228;
t82 = Icges(5,5) * t141 + Icges(5,6) * t140 + Icges(5,3) * t227;
t81 = Icges(5,5) * t139 + Icges(5,6) * t138 - Icges(5,3) * t228;
t80 = -t188 * t217 + (t188 * rSges(4,3) - t261) * t190;
t49 = (-t189 * t239 + t247) * t187;
t43 = (-t188 * t78 - t190 * t77) * t189;
t40 = (-t188 * t70 - t190 * t69) * t189;
t39 = t134 + t190 * (rSges(5,3) * t227 - t206) + (t228 * t266 + t260) * t188;
t28 = t123 * t74 + t124 * t76 + t227 * t72;
t27 = t123 * t73 + t124 * t75 + t227 * t71;
t26 = t121 * t74 + t122 * t76 - t228 * t72;
t25 = t121 * t73 + t122 * t75 - t228 * t71;
t18 = t188 * t265 + t190 * t78 + t245;
t17 = (-t188 * t249 + t190 * t250) * t189;
t15 = t249 * t190 + (-t218 + t250) * t188 + t245;
t14 = t28 * t188 + t190 * t27;
t13 = t26 * t188 + t190 * t25;
t7 = t42 * t187 + (-t188 * t27 + t190 * t28) * t189;
t6 = t41 * t187 + (-t188 * t25 + t190 * t26) * t189;
t1 = [Icges(3,1) + Icges(2,3) + (-t241 + (Icges(5,5) * t185 - Icges(5,6) * t184) * t189 + (Icges(4,2) + Icges(5,3)) * t187) * t187 + (Icges(4,1) * t189 - t118 * t184 + t119 * t185 - t239 - t240 - t242) * t189 + m(7) * (t46 ^ 2 + t47 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t93 ^ 2 + t94 ^ 2) + m(3) * (t132 ^ 2 + t133 ^ 2) + m(2) * (t155 ^ 2 + t158 ^ 2) + t247 + t248; t262 + t263 + t264 + m(4) * (t188 * t93 - t190 * t94) + m(3) * (t188 * t132 - t133 * t190); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t212) * t222; m(7) * (t44 * t46 + t45 * t47) + m(6) * (t50 * t60 + t51 * t61) + m(5) * (t58 * t89 + t59 * t90) + (t37 / 0.2e1 + t31 / 0.2e1 + t138 * t257 + t139 * t256 - t94 * t252 + t151 * t253 + (-Icges(4,6) * t190 / 0.2e1 + t194 * t267 + t81 / 0.2e1) * t187 - t221) * t190 + (t38 / 0.2e1 + t32 / 0.2e1 + t140 * t257 + t141 * t256 + t93 * t252 + t151 * t254 + (Icges(4,6) * t267 + t194 * t253 + t82 / 0.2e1) * t187 + t220) * t188 + ((Icges(4,5) * t188 - t184 * t84 + t185 * t86 - t190 * t195) * t254 + (Icges(4,5) * t190 - t184 * t83 + t185 * t85 + t188 * t195) * t253) * t189; m(5) * t196 + m(6) * t197 + m(7) * t203 + t222 * t252; m(7) * (t15 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t18 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t39 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(4) * (t157 ^ 2 * t222 + t80 ^ 2) + (t183 * t125 + t11 + t13 + (t138 * t83 + t139 * t85 - t228 * t81) * t190) * t190 + (t12 + t14 + t181 * t126 + (t140 * t84 + t141 * t86 + t227 * t82) * t188 + (t188 * t125 + t190 * t126 + t138 * t84 + t139 * t86 + t140 * t83 + t141 * t85 + t227 * t81 - t228 * t82) * t190) * t188; 0.2e1 * (-t262 / 0.2e1 - t263 / 0.2e1 - t264 / 0.2e1) * t189; -0.2e1 * t212 * t222 * t189; m(7) * (t187 * t15 - t189 * t203) + m(6) * (t187 * t18 - t189 * t197) + m(5) * (t187 * t39 - t189 * t196); 0.2e1 * t212 * (t189 ^ 2 * t222 + t187 ^ 2); t49 + m(7) * (t23 * t47 + t24 * t46) + m(6) * (t50 * t55 + t51 * t54) + (t188 * t221 + t190 * t220) * t189 + t210; m(6) * t199 + m(7) * t204; t6 * t253 + t7 * t254 + (t34 * t188 + t33 * t190) * t255 + (t13 * t267 + t14 * t253) * t189 + m(7) * (t15 * t17 + t23 * t45 + t24 * t44) + m(6) * (t18 * t43 + t54 * t61 + t55 * t60) + t211; m(6) * (t43 * t187 - t189 * t199) + m(7) * (t17 * t187 - t189 * t204); t187 * t49 + m(7) * (t17 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(6) * (t43 ^ 2 + t54 ^ 2 + t55 ^ 2) + ((t187 * t34 + t7) * t190 + (-t187 * t33 - t4 - t6) * t188) * t189 + t258; m(7) * (t46 * t53 + t47 * t52) + t210; m(7) * t200; m(7) * (t15 * t40 + t44 * t53 + t45 * t52) + t211; m(7) * (t40 * t187 - t189 * t200); m(7) * (t17 * t40 + t23 * t52 + t24 * t53) + t208; m(7) * (t40 ^ 2 + t52 ^ 2 + t53 ^ 2) + t208;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
