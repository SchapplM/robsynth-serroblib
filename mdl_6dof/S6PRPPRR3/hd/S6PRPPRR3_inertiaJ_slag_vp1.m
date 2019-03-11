% Calculate joint inertia matrix for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPPRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPPRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:30
% EndTime: 2019-03-08 19:21:39
% DurationCPUTime: 4.74s
% Computational Cost: add. (13144->431), mult. (34057->649), div. (0->0), fcn. (44451->12), ass. (0->189)
t230 = Icges(3,1) + Icges(4,1);
t228 = Icges(3,4) - Icges(4,5);
t227 = Icges(4,4) + Icges(3,5);
t229 = Icges(3,2) + Icges(4,3);
t226 = Icges(3,6) - Icges(4,6);
t225 = Icges(3,3) + Icges(4,2);
t224 = Icges(5,3) + t225;
t170 = sin(pkin(10));
t172 = cos(pkin(10));
t176 = sin(qJ(2));
t173 = cos(pkin(6));
t178 = cos(qJ(2));
t196 = t173 * t178;
t160 = t170 * t176 - t172 * t196;
t197 = t173 * t176;
t161 = t170 * t178 + t172 * t197;
t171 = sin(pkin(6));
t200 = t171 * t172;
t223 = t229 * t160 - t228 * t161 + t226 * t200;
t162 = t170 * t196 + t172 * t176;
t163 = -t170 * t197 + t172 * t178;
t201 = t170 * t171;
t222 = t229 * t162 - t228 * t163 - t226 * t201;
t221 = -t228 * t160 + t230 * t161 - t227 * t200;
t220 = -t228 * t162 + t230 * t163 + t227 * t201;
t184 = m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1;
t219 = 0.2e1 * t184;
t218 = t225 * t173 + (t176 * t227 + t178 * t226) * t171;
t217 = -t226 * t173 + (-t228 * t176 - t229 * t178) * t171;
t216 = t227 * t173 + (t230 * t176 + t228 * t178) * t171;
t169 = sin(pkin(11));
t202 = cos(pkin(11));
t135 = -t162 * t202 + t163 * t169;
t136 = t162 * t169 + t163 * t202;
t215 = -Icges(5,5) * t136 + Icges(5,6) * t135 - t162 * t226 + t163 * t227 + t201 * t224;
t133 = -t160 * t202 + t161 * t169;
t134 = t160 * t169 + t161 * t202;
t214 = Icges(5,5) * t134 - Icges(5,6) * t133 + t160 * t226 - t161 * t227 + t200 * t224;
t157 = (t169 * t176 + t178 * t202) * t171;
t158 = (-t169 * t178 + t176 * t202) * t171;
t109 = Icges(5,5) * t158 - Icges(5,6) * t157 - Icges(5,3) * t173;
t213 = -t109 + t218;
t212 = t173 ^ 2;
t175 = sin(qJ(5));
t206 = cos(qJ(5));
t186 = t171 * t206;
t105 = t134 * t175 - t172 * t186;
t107 = t136 * t175 + t170 * t186;
t143 = t158 * t175 + t173 * t206;
t199 = t171 * t175;
t106 = t134 * t206 + t172 * t199;
t174 = sin(qJ(6));
t177 = cos(qJ(6));
t78 = -t106 * t174 + t133 * t177;
t79 = t106 * t177 + t133 * t174;
t50 = Icges(7,5) * t79 + Icges(7,6) * t78 + Icges(7,3) * t105;
t52 = Icges(7,4) * t79 + Icges(7,2) * t78 + Icges(7,6) * t105;
t54 = Icges(7,1) * t79 + Icges(7,4) * t78 + Icges(7,5) * t105;
t108 = t136 * t206 - t170 * t199;
t80 = -t108 * t174 + t135 * t177;
t81 = t108 * t177 + t135 * t174;
t19 = t107 * t50 + t52 * t80 + t54 * t81;
t51 = Icges(7,5) * t81 + Icges(7,6) * t80 + Icges(7,3) * t107;
t53 = Icges(7,4) * t81 + Icges(7,2) * t80 + Icges(7,6) * t107;
t55 = Icges(7,1) * t81 + Icges(7,4) * t80 + Icges(7,5) * t107;
t20 = t107 * t51 + t53 * t80 + t55 * t81;
t144 = t158 * t206 - t173 * t175;
t102 = -t144 * t174 + t157 * t177;
t103 = t144 * t177 + t157 * t174;
t69 = Icges(7,5) * t103 + Icges(7,6) * t102 + Icges(7,3) * t143;
t70 = Icges(7,4) * t103 + Icges(7,2) * t102 + Icges(7,6) * t143;
t71 = Icges(7,1) * t103 + Icges(7,4) * t102 + Icges(7,5) * t143;
t29 = t107 * t69 + t70 * t80 + t71 * t81;
t2 = t105 * t19 + t107 * t20 + t143 * t29;
t211 = t2 / 0.2e1;
t21 = t102 * t52 + t103 * t54 + t143 * t50;
t22 = t102 * t53 + t103 * t55 + t143 * t51;
t36 = t102 * t70 + t103 * t71 + t143 * t69;
t7 = t105 * t21 + t107 * t22 + t143 * t36;
t210 = t7 / 0.2e1;
t209 = t105 / 0.2e1;
t208 = t107 / 0.2e1;
t207 = t143 / 0.2e1;
t58 = rSges(7,1) * t79 + rSges(7,2) * t78 + rSges(7,3) * t105;
t205 = pkin(5) * t106 + pkin(9) * t105 + t58;
t59 = rSges(7,1) * t81 + rSges(7,2) * t80 + rSges(7,3) * t107;
t204 = pkin(5) * t108 + pkin(9) * t107 + t59;
t72 = rSges(7,1) * t103 + rSges(7,2) * t102 + rSges(7,3) * t143;
t203 = pkin(5) * t144 + pkin(9) * t143 + t72;
t198 = t171 * t178;
t138 = pkin(2) * t161 + qJ(3) * t160;
t139 = pkin(2) * t163 + qJ(3) * t162;
t195 = t138 * t201 + t139 * t200;
t137 = t173 * t139;
t146 = pkin(3) * t163 - qJ(4) * t201;
t194 = t173 * t146 + t137;
t145 = pkin(3) * t161 + qJ(4) * t200;
t193 = -t138 - t145;
t164 = (pkin(2) * t176 - qJ(3) * t178) * t171;
t192 = -pkin(3) * t171 * t176 + qJ(4) * t173 - t164;
t191 = m(5) + m(6) + m(7);
t94 = pkin(4) * t136 + pkin(8) * t135;
t190 = t173 * t94 + t194;
t93 = pkin(4) * t134 + pkin(8) * t133;
t189 = -t93 + t193;
t188 = -m(4) - t191;
t187 = -pkin(4) * t158 - pkin(8) * t157 + t192;
t185 = (-t173 * rSges(4,2) - (rSges(4,1) * t176 - rSges(4,3) * t178) * t171 - t164) * t171;
t183 = t145 * t201 + t146 * t200 + t195;
t182 = (-rSges(5,1) * t158 + rSges(5,2) * t157 + rSges(5,3) * t173 + t192) * t171;
t98 = rSges(6,1) * t144 - rSges(6,2) * t143 + rSges(6,3) * t157;
t181 = (-t98 + t187) * t171;
t180 = t94 * t200 + t93 * t201 + t183;
t179 = (t187 - t203) * t171;
t168 = t171 ^ 2;
t154 = t173 * rSges(3,3) + (rSges(3,1) * t176 + rSges(3,2) * t178) * t171;
t128 = rSges(3,1) * t163 - rSges(3,2) * t162 + rSges(3,3) * t201;
t127 = rSges(4,1) * t163 + rSges(4,2) * t201 + rSges(4,3) * t162;
t126 = rSges(3,1) * t161 - rSges(3,2) * t160 - rSges(3,3) * t200;
t125 = rSges(4,1) * t161 - rSges(4,2) * t200 + rSges(4,3) * t160;
t111 = Icges(5,1) * t158 - Icges(5,4) * t157 - Icges(5,5) * t173;
t110 = Icges(5,4) * t158 - Icges(5,2) * t157 - Icges(5,6) * t173;
t100 = -t126 * t173 - t154 * t200;
t99 = t128 * t173 - t154 * t201;
t97 = Icges(6,1) * t144 - Icges(6,4) * t143 + Icges(6,5) * t157;
t96 = Icges(6,4) * t144 - Icges(6,2) * t143 + Icges(6,6) * t157;
t95 = Icges(6,5) * t144 - Icges(6,6) * t143 + Icges(6,3) * t157;
t89 = rSges(5,1) * t136 - rSges(5,2) * t135 - rSges(5,3) * t201;
t88 = rSges(5,1) * t134 - rSges(5,2) * t133 + rSges(5,3) * t200;
t87 = Icges(5,1) * t136 - Icges(5,4) * t135 - Icges(5,5) * t201;
t86 = Icges(5,1) * t134 - Icges(5,4) * t133 + Icges(5,5) * t200;
t85 = Icges(5,4) * t136 - Icges(5,2) * t135 - Icges(5,6) * t201;
t84 = Icges(5,4) * t134 - Icges(5,2) * t133 + Icges(5,6) * t200;
t77 = (t126 * t170 + t128 * t172) * t171;
t74 = (-t125 - t138) * t173 + t172 * t185;
t73 = t127 * t173 + t170 * t185 + t137;
t68 = rSges(6,1) * t108 - rSges(6,2) * t107 + rSges(6,3) * t135;
t67 = rSges(6,1) * t106 - rSges(6,2) * t105 + rSges(6,3) * t133;
t66 = Icges(6,1) * t108 - Icges(6,4) * t107 + Icges(6,5) * t135;
t65 = Icges(6,1) * t106 - Icges(6,4) * t105 + Icges(6,5) * t133;
t64 = Icges(6,4) * t108 - Icges(6,2) * t107 + Icges(6,6) * t135;
t63 = Icges(6,4) * t106 - Icges(6,2) * t105 + Icges(6,6) * t133;
t62 = Icges(6,5) * t108 - Icges(6,6) * t107 + Icges(6,3) * t135;
t61 = Icges(6,5) * t106 - Icges(6,6) * t105 + Icges(6,3) * t133;
t60 = (t125 * t170 + t127 * t172) * t171 + t195;
t57 = (-t88 + t193) * t173 + t172 * t182;
t56 = t170 * t182 + t173 * t89 + t194;
t49 = -t135 * t98 + t157 * t68;
t48 = t133 * t98 - t157 * t67;
t47 = -t143 * t96 + t144 * t97 + t157 * t95;
t46 = (t170 * t88 + t172 * t89) * t171 + t183;
t45 = -t133 * t68 + t135 * t67;
t44 = -t107 * t96 + t108 * t97 + t135 * t95;
t43 = -t105 * t96 + t106 * t97 + t133 * t95;
t42 = -t107 * t72 + t143 * t59;
t41 = t105 * t72 - t143 * t58;
t40 = (-t67 + t189) * t173 + t172 * t181;
t39 = t170 * t181 + t173 * t68 + t190;
t38 = -t143 * t64 + t144 * t66 + t157 * t62;
t37 = -t143 * t63 + t144 * t65 + t157 * t61;
t35 = -t107 * t64 + t108 * t66 + t135 * t62;
t34 = -t107 * t63 + t108 * t65 + t135 * t61;
t33 = -t105 * t64 + t106 * t66 + t133 * t62;
t32 = -t105 * t63 + t106 * t65 + t133 * t61;
t31 = -t105 * t59 + t107 * t58;
t30 = (t170 * t67 + t172 * t68) * t171 + t180;
t28 = t105 * t69 + t70 * t78 + t71 * t79;
t27 = -t135 * t203 + t157 * t204;
t26 = t133 * t203 - t157 * t205;
t25 = (t189 - t205) * t173 + t172 * t179;
t24 = t170 * t179 + t173 * t204 + t190;
t23 = -t133 * t204 + t135 * t205;
t18 = t105 * t51 + t53 * t78 + t55 * t79;
t17 = t105 * t50 + t52 * t78 + t54 * t79;
t16 = (t170 * t205 + t172 * t204) * t171 + t180;
t15 = t173 * t47 + (t170 * t38 - t172 * t37) * t171;
t14 = t133 * t37 + t135 * t38 + t157 * t47;
t13 = t173 * t44 + (t170 * t35 - t172 * t34) * t171;
t12 = t173 * t43 + (t170 * t33 - t172 * t32) * t171;
t11 = t133 * t34 + t135 * t35 + t157 * t44;
t10 = t133 * t32 + t135 * t33 + t157 * t43;
t9 = t173 * t36 + (t170 * t22 - t172 * t21) * t171;
t8 = t133 * t21 + t135 * t22 + t157 * t36;
t6 = t173 * t29 + (t170 * t20 - t172 * t19) * t171;
t5 = t173 * t28 + (-t17 * t172 + t170 * t18) * t171;
t4 = t133 * t19 + t135 * t20 + t157 * t29;
t3 = t133 * t17 + t135 * t18 + t157 * t28;
t1 = t105 * t17 + t107 * t18 + t143 * t28;
t75 = [m(2) + m(3) - t188; m(3) * t77 + m(4) * t60 + m(5) * t46 + m(6) * t30 + m(7) * t16; m(7) * (t16 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t30 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t46 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(4) * (t60 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(3) * (t100 ^ 2 + t77 ^ 2 + t99 ^ 2) + (t9 + t15 + (-t109 * t173 - t110 * t157 + t111 * t158) * t173 + t218 * t212 + ((t176 * t216 - t178 * t217) * t173 + (t157 * t84 - t158 * t86 + t214 * t173 + (-t221 * t176 + t223 * t178) * t171) * t172 + (-t157 * t85 + t158 * t87 + t215 * t173 + (t220 * t176 - t222 * t178) * t171) * t170) * t171) * t173 + (t6 + t13 + (-t135 * t85 + t136 * t87 + t222 * t162 + t220 * t163 + t215 * t201) * t201 + (-t110 * t135 + t111 * t136 + t162 * t217 + t163 * t216 + t201 * t213) * t173) * t201 + (-t5 - t12 + (-t133 * t84 + t134 * t86 + t223 * t160 + t221 * t161 + t214 * t200) * t200 + (t110 * t133 - t111 * t134 - t160 * t217 - t161 * t216 + t200 * t213) * t173 + (t133 * t85 - t134 * t87 + t135 * t84 - t136 * t86 - t222 * t160 - t220 * t161 - t223 * t162 - t221 * t163 + t215 * t200 + t214 * t201) * t201) * t200; t188 * t198; m(7) * (-t16 * t198 + t160 * t24 + t162 * t25) + m(6) * (t160 * t39 + t162 * t40 - t198 * t30) + m(5) * (t160 * t56 + t162 * t57 - t198 * t46) + m(4) * (t160 * t73 + t162 * t74 - t198 * t60); 0.2e1 * (m(4) / 0.2e1 + t184) * (t168 * t178 ^ 2 + t160 ^ 2 + t162 ^ 2); -t191 * t173; m(7) * (-t16 * t173 + (-t170 * t25 + t172 * t24) * t171) + m(6) * (-t173 * t30 + (-t170 * t40 + t172 * t39) * t171) + m(5) * (-t173 * t46 + (-t170 * t57 + t172 * t56) * t171); (t160 * t172 - t162 * t170 + t196) * t171 * t219; (t212 + (t170 ^ 2 + t172 ^ 2) * t168) * t219; m(6) * t45 + m(7) * t23; (t8 / 0.2e1 + t14 / 0.2e1) * t173 + (t9 / 0.2e1 + t15 / 0.2e1) * t157 + (t6 / 0.2e1 + t13 / 0.2e1) * t135 + (t5 / 0.2e1 + t12 / 0.2e1) * t133 + m(7) * (t16 * t23 + t24 * t27 + t25 * t26) + m(6) * (t30 * t45 + t39 * t49 + t40 * t48) + ((-t3 / 0.2e1 - t10 / 0.2e1) * t172 + (t4 / 0.2e1 + t11 / 0.2e1) * t170) * t171; m(6) * (t49 * t160 + t48 * t162 - t198 * t45) + m(7) * (t27 * t160 + t26 * t162 - t198 * t23); m(6) * (-t173 * t45 + (-t170 * t48 + t172 * t49) * t171) + m(7) * (-t173 * t23 + (-t170 * t26 + t172 * t27) * t171); (t8 + t14) * t157 + (t4 + t11) * t135 + (t3 + t10) * t133 + m(7) * (t23 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t45 ^ 2 + t48 ^ 2 + t49 ^ 2); m(7) * t31; m(7) * (t16 * t31 + t24 * t42 + t25 * t41) + t9 * t207 + t5 * t209 + t173 * t210 + t6 * t208 + (t170 * t211 - t172 * t1 / 0.2e1) * t171; m(7) * (t42 * t160 + t41 * t162 - t198 * t31); m(7) * (-t173 * t31 + (-t170 * t41 + t172 * t42) * t171); t157 * t210 + t133 * t1 / 0.2e1 + t3 * t209 + t4 * t208 + m(7) * (t23 * t31 + t26 * t41 + t27 * t42) + t135 * t211 + t8 * t207; t107 * t2 + t105 * t1 + t143 * t7 + m(7) * (t31 ^ 2 + t41 ^ 2 + t42 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t75(1) t75(2) t75(4) t75(7) t75(11) t75(16); t75(2) t75(3) t75(5) t75(8) t75(12) t75(17); t75(4) t75(5) t75(6) t75(9) t75(13) t75(18); t75(7) t75(8) t75(9) t75(10) t75(14) t75(19); t75(11) t75(12) t75(13) t75(14) t75(15) t75(20); t75(16) t75(17) t75(18) t75(19) t75(20) t75(21);];
Mq  = res;
