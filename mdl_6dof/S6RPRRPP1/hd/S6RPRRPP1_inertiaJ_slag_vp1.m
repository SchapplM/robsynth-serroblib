% Calculate joint inertia matrix for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:55
% EndTime: 2019-03-09 04:28:02
% DurationCPUTime: 2.72s
% Computational Cost: add. (7879->389), mult. (7506->560), div. (0->0), fcn. (8051->10), ass. (0->196)
t180 = sin(qJ(3));
t257 = Icges(4,5) * t180;
t256 = t257 / 0.2e1;
t175 = qJ(4) + pkin(10);
t170 = sin(t175);
t172 = cos(t175);
t176 = qJ(1) + pkin(9);
t173 = cos(t176);
t171 = sin(t176);
t183 = cos(qJ(3));
t223 = t171 * t183;
t118 = t170 * t223 + t172 * t173;
t119 = -t170 * t173 + t172 * t223;
t253 = rSges(7,3) + qJ(6);
t254 = rSges(7,1) + pkin(5);
t255 = -t253 * t118 - t254 * t119;
t123 = -Icges(6,3) * t183 + (Icges(6,5) * t172 - Icges(6,6) * t170) * t180;
t124 = -Icges(7,2) * t183 + (Icges(7,4) * t172 + Icges(7,6) * t170) * t180;
t179 = sin(qJ(4));
t182 = cos(qJ(4));
t137 = -Icges(5,3) * t183 + (Icges(5,5) * t182 - Icges(5,6) * t179) * t180;
t252 = -t123 - t124 - t137;
t224 = t171 * t180;
t63 = Icges(7,5) * t119 + Icges(7,6) * t224 + Icges(7,3) * t118;
t67 = Icges(7,4) * t119 + Icges(7,2) * t224 + Icges(7,6) * t118;
t71 = Icges(7,1) * t119 + Icges(7,4) * t224 + Icges(7,5) * t118;
t19 = t118 * t63 + t119 * t71 + t224 * t67;
t219 = t173 * t183;
t120 = t170 * t219 - t171 * t172;
t121 = t170 * t171 + t172 * t219;
t220 = t173 * t180;
t64 = Icges(7,5) * t121 + Icges(7,6) * t220 + Icges(7,3) * t120;
t68 = Icges(7,4) * t121 + Icges(7,2) * t220 + Icges(7,6) * t120;
t72 = Icges(7,1) * t121 + Icges(7,4) * t220 + Icges(7,5) * t120;
t20 = t118 * t64 + t119 * t72 + t224 * t68;
t65 = Icges(6,5) * t119 - Icges(6,6) * t118 + Icges(6,3) * t224;
t69 = Icges(6,4) * t119 - Icges(6,2) * t118 + Icges(6,6) * t224;
t73 = Icges(6,1) * t119 - Icges(6,4) * t118 + Icges(6,5) * t224;
t21 = -t118 * t69 + t119 * t73 + t224 * t65;
t66 = Icges(6,5) * t121 - Icges(6,6) * t120 + Icges(6,3) * t220;
t70 = Icges(6,4) * t121 - Icges(6,2) * t120 + Icges(6,6) * t220;
t74 = Icges(6,1) * t121 - Icges(6,4) * t120 + Icges(6,5) * t220;
t22 = -t118 * t70 + t119 * t74 + t224 * t66;
t217 = t179 * t183;
t133 = -t171 * t217 - t173 * t182;
t216 = t182 * t183;
t221 = t173 * t179;
t134 = t171 * t216 - t221;
t83 = Icges(5,5) * t134 + Icges(5,6) * t133 + Icges(5,3) * t224;
t85 = Icges(5,4) * t134 + Icges(5,2) * t133 + Icges(5,6) * t224;
t87 = Icges(5,1) * t134 + Icges(5,4) * t133 + Icges(5,5) * t224;
t27 = t133 * t85 + t134 * t87 + t224 * t83;
t135 = t171 * t182 - t173 * t217;
t225 = t171 * t179;
t136 = t173 * t216 + t225;
t84 = Icges(5,5) * t136 + Icges(5,6) * t135 + Icges(5,3) * t220;
t86 = Icges(5,4) * t136 + Icges(5,2) * t135 + Icges(5,6) * t220;
t88 = Icges(5,1) * t136 + Icges(5,4) * t135 + Icges(5,5) * t220;
t28 = t133 * t86 + t134 * t88 + t224 * t84;
t122 = -Icges(7,6) * t183 + (Icges(7,5) * t172 + Icges(7,3) * t170) * t180;
t126 = -Icges(7,4) * t183 + (Icges(7,1) * t172 + Icges(7,5) * t170) * t180;
t42 = t118 * t122 + t119 * t126 + t124 * t224;
t125 = -Icges(6,6) * t183 + (Icges(6,4) * t172 - Icges(6,2) * t170) * t180;
t127 = -Icges(6,5) * t183 + (Icges(6,1) * t172 - Icges(6,4) * t170) * t180;
t43 = -t118 * t125 + t119 * t127 + t123 * t224;
t138 = -Icges(5,6) * t183 + (Icges(5,4) * t182 - Icges(5,2) * t179) * t180;
t139 = -Icges(5,5) * t183 + (Icges(5,1) * t182 - Icges(5,4) * t179) * t180;
t49 = t133 * t138 + t134 * t139 + t137 * t224;
t251 = (-t42 - t43 - t49) * t183 + ((t20 + t22 + t28) * t173 + (t19 + t21 + t27) * t171) * t180;
t23 = t120 * t63 + t121 * t71 + t220 * t67;
t24 = t120 * t64 + t121 * t72 + t220 * t68;
t25 = -t120 * t69 + t121 * t73 + t220 * t65;
t26 = -t120 * t70 + t121 * t74 + t220 * t66;
t29 = t135 * t85 + t136 * t87 + t220 * t83;
t30 = t135 * t86 + t136 * t88 + t220 * t84;
t44 = t120 * t122 + t121 * t126 + t124 * t220;
t45 = -t120 * t125 + t121 * t127 + t123 * t220;
t50 = t135 * t138 + t136 * t139 + t137 * t220;
t250 = (-t44 - t45 - t50) * t183 + ((t24 + t26 + t30) * t173 + (t23 + t25 + t29) * t171) * t180;
t31 = -t183 * t67 + (t170 * t63 + t172 * t71) * t180;
t33 = -t183 * t65 + (-t170 * t69 + t172 * t73) * t180;
t37 = -t183 * t83 + (-t179 * t85 + t182 * t87) * t180;
t249 = -t31 - t33 - t37;
t32 = -t183 * t68 + (t170 * t64 + t172 * t72) * t180;
t34 = -t183 * t66 + (-t170 * t70 + t172 * t74) * t180;
t38 = -t183 * t84 + (-t179 * t86 + t182 * t88) * t180;
t248 = t32 + t34 + t38;
t226 = t170 * t180;
t247 = t122 * t226 + (t182 * t139 + (t126 + t127) * t172) * t180;
t168 = t171 ^ 2;
t169 = t173 ^ 2;
t246 = t183 ^ 2;
t245 = m(6) / 0.2e1;
t244 = m(7) / 0.2e1;
t243 = -m(6) - m(7);
t242 = t171 / 0.2e1;
t240 = -t183 / 0.2e1;
t151 = rSges(4,1) * t180 + rSges(4,2) * t183;
t239 = m(4) * t151;
t238 = pkin(3) * t183;
t237 = pkin(8) * t180;
t181 = sin(qJ(1));
t236 = t181 * pkin(1);
t167 = t182 * pkin(4) + pkin(3);
t235 = -pkin(3) + t167;
t234 = rSges(7,2) * t224 - t255;
t233 = rSges(7,2) * t220 + t253 * t120 + t254 * t121;
t78 = t121 * rSges(6,1) - t120 * rSges(6,2) + rSges(6,3) * t220;
t178 = -qJ(5) - pkin(8);
t218 = t178 * t180;
t188 = pkin(4) * t225 + t167 * t219 - t173 * t218;
t209 = pkin(3) * t219 + pkin(8) * t220;
t90 = t188 - t209;
t232 = -t78 - t90;
t117 = (pkin(8) + t178) * t183 + t235 * t180;
t210 = -pkin(4) * t221 - t171 * t218;
t89 = (t183 * t235 - t237) * t171 + t210;
t231 = t117 * t224 + t183 * t89;
t230 = t173 * rSges(4,3);
t228 = Icges(4,4) * t183;
t227 = t138 * t179;
t214 = t168 * (t237 + t238) + t173 * t209;
t129 = -t183 * rSges(6,3) + (rSges(6,1) * t172 - rSges(6,2) * t170) * t180;
t213 = -t117 - t129;
t212 = -t183 * rSges(7,2) + (t253 * t170 + t254 * t172) * t180;
t141 = -t183 * rSges(5,3) + (rSges(5,1) * t182 - rSges(5,2) * t179) * t180;
t158 = t180 * pkin(3) - t183 * pkin(8);
t211 = -t141 - t158;
t208 = t168 + t169;
t207 = -t125 * t226 - t180 * t227 + t252 * t183 + t247;
t206 = -t90 - t233;
t205 = -t117 - t212;
t204 = -t158 + t213;
t92 = t136 * rSges(5,1) + t135 * rSges(5,2) + rSges(5,3) * t220;
t184 = cos(qJ(1));
t174 = t184 * pkin(1);
t203 = t173 * pkin(2) + t171 * pkin(7) + t174;
t202 = t173 * pkin(7) - t236;
t201 = -t167 * t183 - pkin(2);
t200 = t171 * t89 + t173 * t90 + t214;
t199 = -t158 + t205;
t198 = rSges(4,1) * t183 - rSges(4,2) * t180;
t197 = -t134 * rSges(5,1) - t133 * rSges(5,2);
t196 = -t119 * rSges(6,1) + t118 * rSges(6,2);
t194 = -Icges(4,2) * t180 + t228;
t193 = Icges(4,5) * t183 - Icges(4,6) * t180;
t190 = t202 - t210;
t189 = rSges(4,1) * t219 - rSges(4,2) * t220 + t171 * rSges(4,3);
t187 = t37 / 0.2e1 + t33 / 0.2e1 + t31 / 0.2e1 + t49 / 0.2e1 + t43 / 0.2e1 + t42 / 0.2e1;
t186 = t50 / 0.2e1 + t45 / 0.2e1 + t44 / 0.2e1 + t38 / 0.2e1 + t34 / 0.2e1 + t32 / 0.2e1;
t185 = t188 + t203;
t177 = t180 ^ 2;
t153 = rSges(2,1) * t184 - t181 * rSges(2,2);
t152 = -t181 * rSges(2,1) - rSges(2,2) * t184;
t148 = Icges(4,6) * t183 + t257;
t143 = rSges(3,1) * t173 - rSges(3,2) * t171 + t174;
t142 = -rSges(3,1) * t171 - rSges(3,2) * t173 - t236;
t104 = Icges(4,3) * t171 + t173 * t193;
t103 = -Icges(4,3) * t173 + t171 * t193;
t98 = t211 * t173;
t97 = t211 * t171;
t96 = t189 + t203;
t95 = t230 + (-pkin(2) - t198) * t171 + t202;
t91 = rSges(5,3) * t224 - t197;
t79 = t89 * t220;
t76 = rSges(6,3) * t224 - t196;
t62 = t173 * t189 + (t171 * t198 - t230) * t171;
t61 = t204 * t173;
t60 = t204 * t171;
t58 = -t141 * t220 - t183 * t92;
t57 = t141 * t224 + t183 * t91;
t56 = t203 + t92 + t209;
t55 = (-t238 - pkin(2) + (-rSges(5,3) - pkin(8)) * t180) * t171 + t197 + t202;
t52 = t199 * t173;
t51 = t199 * t171;
t48 = t185 + t78;
t47 = (-rSges(6,3) * t180 + t201) * t171 + t190 + t196;
t46 = (-t171 * t92 + t173 * t91) * t180;
t41 = t185 + t233;
t40 = (-rSges(7,2) * t180 + t201) * t171 + t190 + t255;
t39 = t171 * t91 + t173 * t92 + t214;
t36 = t183 * t232 + t213 * t220;
t35 = t129 * t224 + t183 * t76 + t231;
t18 = t79 + (t171 * t232 + t173 * t76) * t180;
t17 = t183 * t206 + t205 * t220;
t16 = t183 * t234 + t212 * t224 + t231;
t15 = t171 * t76 + t173 * t78 + t200;
t14 = t79 + (t171 * t206 + t173 * t234) * t180;
t13 = t171 * t234 + t173 * t233 + t200;
t12 = t171 * t30 - t173 * t29;
t11 = t171 * t28 - t173 * t27;
t10 = t171 * t26 - t173 * t25;
t9 = t171 * t24 - t173 * t23;
t8 = t171 * t22 - t173 * t21;
t7 = t171 * t20 - t173 * t19;
t1 = [Icges(2,3) + Icges(3,3) + (Icges(4,1) * t180 - t125 * t170 - t227 + t228) * t180 + (Icges(4,4) * t180 + Icges(4,2) * t183 + t252) * t183 + m(6) * (t47 ^ 2 + t48 ^ 2) + m(7) * (t40 ^ 2 + t41 ^ 2) + m(5) * (t55 ^ 2 + t56 ^ 2) + m(4) * (t95 ^ 2 + t96 ^ 2) + m(3) * (t142 ^ 2 + t143 ^ 2) + m(2) * (t152 ^ 2 + t153 ^ 2) + t247; 0; m(3) + m(4) + m(5) - t243; m(6) * (t47 * t61 + t48 * t60) + m(7) * (t40 * t52 + t41 * t51) + m(5) * (t55 * t98 + t56 * t97) + (t171 * t194 * t240 - t95 * t239 - t187 + (-Icges(4,6) * t240 + t256 + t148 / 0.2e1) * t173) * t173 + (-t96 * t239 + t183 * (Icges(4,6) * t171 + t173 * t194) / 0.2e1 + t171 * t256 + t148 * t242 + t186) * t171; m(4) * t62 + m(5) * t39 + m(6) * t15 + m(7) * t13; m(6) * (t15 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(7) * (t13 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t39 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(4) * (t151 ^ 2 * t208 + t62 ^ 2) + (-t169 * t103 - t11 - t7 - t8) * t173 + (t168 * t104 + t10 + t12 + t9 + (-t171 * t103 + t173 * t104) * t173) * t171; -t207 * t183 + m(6) * (t35 * t47 + t36 * t48) + m(7) * (t16 * t40 + t17 * t41) + m(5) * (t55 * t57 + t56 * t58) + (t171 * t187 + t173 * t186) * t180; m(5) * t46 + m(6) * t18 + m(7) * t14; m(6) * (t15 * t18 + t35 * t61 + t36 * t60) + m(7) * (t13 * t14 + t16 * t52 + t17 * t51) + m(5) * (t39 * t46 + t57 * t98 + t58 * t97) + ((t9 / 0.2e1 + t12 / 0.2e1 + t10 / 0.2e1) * t173 + (t11 / 0.2e1 + t8 / 0.2e1 + t7 / 0.2e1) * t171) * t180 + t250 * t242 - t251 * t173 / 0.2e1 + (t248 * t171 + t249 * t173) * t240; m(7) * (t14 ^ 2 + t16 ^ 2 + t17 ^ 2) + m(6) * (t18 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t46 ^ 2 + t57 ^ 2 + t58 ^ 2) + t207 * t246 + ((-t248 * t183 + t250) * t173 + (t249 * t183 + t251) * t171) * t180; 0.2e1 * ((t171 * t48 + t173 * t47) * t245 + (t171 * t41 + t173 * t40) * t244) * t180; t243 * t183; m(6) * (-t183 * t15 + (t171 * t60 + t173 * t61) * t180) + m(7) * (-t183 * t13 + (t171 * t51 + t173 * t52) * t180); m(7) * (-t183 * t14 + (t16 * t173 + t17 * t171) * t180) + m(6) * (-t183 * t18 + (t171 * t36 + t173 * t35) * t180); 0.2e1 * (t245 + t244) * (t177 * t208 + t246); m(7) * (t118 * t41 + t120 * t40); m(7) * t226; m(7) * (t118 * t51 + t120 * t52 + t13 * t226); m(7) * (t118 * t17 + t120 * t16 + t14 * t226); m(7) * (t118 * t171 + t120 * t173 - t170 * t183) * t180; m(7) * (t170 ^ 2 * t177 + t118 ^ 2 + t120 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
