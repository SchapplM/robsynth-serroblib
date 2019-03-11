% Calculate joint inertia matrix for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:25
% EndTime: 2019-03-09 02:35:30
% DurationCPUTime: 2.26s
% Computational Cost: add. (6355->340), mult. (6495->496), div. (0->0), fcn. (6924->10), ass. (0->170)
t151 = pkin(10) + qJ(4);
t144 = cos(t151);
t221 = Icges(5,5) * t144;
t143 = sin(t151);
t220 = Icges(5,6) * t143;
t219 = t221 / 0.2e1 - t220 / 0.2e1;
t161 = cos(qJ(1));
t218 = (rSges(5,1) * t143 + rSges(5,2) * t144) * t161;
t159 = sin(qJ(1));
t152 = t159 ^ 2;
t153 = t161 ^ 2;
t195 = t144 * t161;
t154 = qJ(5) + qJ(6);
t145 = sin(t154);
t146 = cos(t154);
t189 = t161 * t146;
t194 = t159 * t145;
t105 = -t143 * t194 + t189;
t190 = t161 * t145;
t193 = t159 * t146;
t106 = t143 * t193 + t190;
t196 = t144 * t159;
t56 = Icges(7,5) * t106 + Icges(7,6) * t105 - Icges(7,3) * t196;
t58 = Icges(7,4) * t106 + Icges(7,2) * t105 - Icges(7,6) * t196;
t60 = Icges(7,1) * t106 + Icges(7,4) * t105 - Icges(7,5) * t196;
t26 = t143 * t56 + (-t145 * t58 + t146 * t60) * t144;
t107 = t143 * t190 + t193;
t108 = -t143 * t189 + t194;
t57 = Icges(7,5) * t108 + Icges(7,6) * t107 + Icges(7,3) * t195;
t59 = Icges(7,4) * t108 + Icges(7,2) * t107 + Icges(7,6) * t195;
t61 = Icges(7,1) * t108 + Icges(7,4) * t107 + Icges(7,5) * t195;
t27 = t143 * t57 + (-t145 * t59 + t146 * t61) * t144;
t89 = Icges(7,6) * t143 + (Icges(7,4) * t146 - Icges(7,2) * t145) * t144;
t203 = t145 * t89;
t88 = Icges(7,3) * t143 + (Icges(7,5) * t146 - Icges(7,6) * t145) * t144;
t90 = Icges(7,5) * t143 + (Icges(7,1) * t146 - Icges(7,4) * t145) * t144;
t206 = t144 * t146 * t90 + t143 * t88;
t42 = (-t144 * t203 + t206) * t143;
t20 = t107 * t58 + t108 * t60 + t56 * t195;
t21 = t107 * t59 + t108 * t61 + t57 * t195;
t38 = t107 * t89 + t108 * t90 + t88 * t195;
t5 = t38 * t143 + (-t159 * t20 + t161 * t21) * t144;
t217 = t5 * t195 + t143 * (t42 + (-t159 * t26 + t161 * t27) * t144);
t215 = t143 / 0.2e1;
t213 = t159 / 0.2e1;
t212 = t161 / 0.2e1;
t123 = t144 * rSges(5,1) - t143 * rSges(5,2);
t211 = m(5) * t123;
t137 = t152 + t153;
t128 = m(5) * t137;
t155 = sin(pkin(10));
t210 = pkin(3) * t155;
t162 = -pkin(9) - pkin(8);
t209 = -pkin(8) - t162;
t62 = t106 * rSges(7,1) + t105 * rSges(7,2) - rSges(7,3) * t196;
t197 = t143 * t159;
t135 = pkin(4) * t197;
t111 = -pkin(8) * t196 + t135;
t160 = cos(qJ(5));
t142 = t160 * pkin(5) + pkin(4);
t158 = sin(qJ(5));
t188 = t161 * t158;
t182 = pkin(5) * t188 + t142 * t197 + t162 * t196;
t65 = -t111 + t182;
t208 = -t62 - t65;
t136 = t161 * t143 * pkin(4);
t192 = t159 * t158;
t198 = t142 * t143;
t169 = -t108 * rSges(7,1) - t107 * rSges(7,2);
t63 = rSges(7,3) * t195 - t169;
t207 = t63 + pkin(5) * t192 + t136 + (t209 * t144 - t198) * t161;
t91 = t143 * rSges(7,3) + (rSges(7,1) * t146 - rSges(7,2) * t145) * t144;
t47 = t143 * t62 + t91 * t196;
t92 = Icges(6,3) * t143 + (Icges(6,5) * t160 - Icges(6,6) * t158) * t144;
t94 = Icges(6,5) * t143 + (Icges(6,1) * t160 - Icges(6,4) * t158) * t144;
t205 = t144 * t160 * t94 + t143 * t92;
t86 = (-pkin(4) + t142) * t144 + t209 * t143;
t204 = t86 + t91;
t93 = Icges(6,6) * t143 + (Icges(6,4) * t160 - Icges(6,2) * t158) * t144;
t202 = t158 * t93;
t201 = rSges(4,3) + qJ(3);
t191 = t159 * t160;
t187 = t161 * t160;
t114 = -t143 * t192 + t187;
t115 = t143 * t191 + t188;
t186 = t115 * rSges(6,1) + t114 * rSges(6,2);
t185 = t161 * pkin(1) + t159 * qJ(2);
t67 = Icges(6,5) * t115 + Icges(6,6) * t114 - Icges(6,3) * t196;
t69 = Icges(6,4) * t115 + Icges(6,2) * t114 - Icges(6,6) * t196;
t71 = Icges(6,1) * t115 + Icges(6,4) * t114 - Icges(6,5) * t196;
t32 = t143 * t67 + (-t158 * t69 + t160 * t71) * t144;
t40 = t114 * t93 + t115 * t94 - t92 * t196;
t184 = -t32 / 0.2e1 - t40 / 0.2e1;
t116 = t143 * t188 + t191;
t117 = -t143 * t187 + t192;
t68 = Icges(6,5) * t117 + Icges(6,6) * t116 + Icges(6,3) * t195;
t70 = Icges(6,4) * t117 + Icges(6,2) * t116 + Icges(6,6) * t195;
t72 = Icges(6,1) * t117 + Icges(6,4) * t116 + Icges(6,5) * t195;
t33 = t143 * t68 + (-t158 * t70 + t160 * t72) * t144;
t41 = t116 * t93 + t117 * t94 + t92 * t195;
t183 = t33 / 0.2e1 + t41 / 0.2e1;
t181 = rSges(5,1) * t197 + rSges(5,2) * t196 + t161 * rSges(5,3);
t148 = t161 * qJ(2);
t157 = -pkin(7) - qJ(3);
t180 = t159 * t157 + t161 * t210 + t148;
t179 = t144 * (-rSges(6,3) - pkin(8));
t178 = -t196 / 0.2e1;
t177 = t195 / 0.2e1;
t18 = t105 * t58 + t106 * t60 - t56 * t196;
t19 = t105 * t59 + t106 * t61 - t57 * t196;
t11 = t19 * t159 + t18 * t161;
t12 = t21 * t159 + t20 * t161;
t37 = t105 * t89 + t106 * t90 - t88 * t196;
t4 = t37 * t143 + (-t159 * t18 + t161 * t19) * t144;
t176 = t11 * t178 + t12 * t177 + t4 * t212 + t5 * t213 + (t27 * t159 + t26 * t161) * t215;
t175 = t128 + (m(4) + m(6) + m(7)) * t137;
t174 = t42 + (t26 + t37) * t178 + (t27 + t38) * t177;
t173 = -t4 * t196 + t217;
t156 = cos(pkin(10));
t172 = rSges(4,1) * t155 + rSges(4,2) * t156;
t170 = -t117 * rSges(6,1) - t116 * rSges(6,2);
t164 = Icges(5,5) * t143 + Icges(5,6) * t144;
t163 = -t161 * t157 + t159 * t210 + t185;
t132 = t161 * rSges(2,1) - t159 * rSges(2,2);
t131 = -t159 * rSges(2,1) - t161 * rSges(2,2);
t125 = t144 * pkin(4) + t143 * pkin(8);
t120 = -t220 + t221;
t118 = t159 * t125;
t113 = -t161 * rSges(3,2) + t159 * rSges(3,3) + t185;
t112 = t161 * rSges(3,3) + t148 + (rSges(3,2) - pkin(1)) * t159;
t104 = t161 * (pkin(8) * t195 - t136);
t97 = Icges(5,3) * t159 - t164 * t161;
t96 = Icges(5,3) * t161 + t164 * t159;
t95 = t143 * rSges(6,3) + (rSges(6,1) * t160 - rSges(6,2) * t158) * t144;
t85 = t172 * t159 + t201 * t161 + t185;
t84 = t148 + t172 * t161 + (-pkin(1) - t201) * t159;
t81 = t91 * t195;
t78 = t163 + t181;
t77 = t218 + (-rSges(5,3) - pkin(1)) * t159 + t180;
t76 = (-t125 - t95) * t161;
t75 = t159 * t95 + t118;
t74 = rSges(6,3) * t195 - t170;
t73 = -rSges(6,3) * t196 + t186;
t64 = -t159 * t181 + (t159 * rSges(5,3) - t218) * t161;
t54 = (-t125 - t204) * t161;
t53 = t204 * t159 + t118;
t52 = t159 * t179 + t135 + t163 + t186;
t51 = -t159 * pkin(1) + t161 * t179 + t136 + t170 + t180;
t50 = -t143 * t74 + t95 * t195;
t49 = t143 * t73 + t95 * t196;
t48 = -t143 * t63 + t81;
t46 = (-t144 * t202 + t205) * t143;
t45 = t163 + t62 + t182;
t44 = (-pkin(5) * t158 - pkin(1)) * t159 + (t198 + (-rSges(7,3) + t162) * t144) * t161 + t169 + t180;
t43 = (-t159 * t74 - t161 * t73) * t144;
t39 = (-t159 * t63 - t161 * t62) * t144;
t36 = t161 * t74 + t104 + (-t111 - t73) * t159;
t31 = t116 * t70 + t117 * t72 + t68 * t195;
t30 = t116 * t69 + t117 * t71 + t67 * t195;
t29 = t114 * t70 + t115 * t72 - t68 * t196;
t28 = t114 * t69 + t115 * t71 - t67 * t196;
t25 = -t207 * t143 + t86 * t195 + t81;
t24 = t143 * t65 + t86 * t196 + t47;
t17 = (-t207 * t159 + t208 * t161) * t144;
t16 = t104 + t207 * t161 + (-t111 + t208) * t159;
t15 = t31 * t159 + t30 * t161;
t14 = t29 * t159 + t28 * t161;
t8 = t41 * t143 + (-t159 * t30 + t161 * t31) * t144;
t7 = t40 * t143 + (-t159 * t28 + t161 * t29) * t144;
t1 = [Icges(4,1) * t156 ^ 2 + Icges(3,1) + Icges(2,3) + (-0.2e1 * Icges(4,4) * t156 + Icges(4,2) * t155) * t155 + (Icges(5,1) * t144 - t202 - t203) * t144 + m(7) * (t44 ^ 2 + t45 ^ 2) + m(6) * (t51 ^ 2 + t52 ^ 2) + m(5) * (t77 ^ 2 + t78 ^ 2) + m(4) * (t84 ^ 2 + t85 ^ 2) + m(3) * (t112 ^ 2 + t113 ^ 2) + m(2) * (t131 ^ 2 + t132 ^ 2) + t205 + t206 + (-0.2e1 * Icges(5,4) * t144 + Icges(5,2) * t143) * t143; m(7) * (t159 * t44 - t161 * t45) + m(6) * (t159 * t51 - t161 * t52) + m(5) * (t159 * t77 - t161 * t78) + m(4) * (t159 * t84 - t161 * t85) + m(3) * (t159 * t112 - t161 * t113); m(3) * t137 + t175; m(7) * (t159 * t45 + t161 * t44) + m(6) * (t159 * t52 + t161 * t51) + m(5) * (t159 * t78 + t161 * t77) + m(4) * (t159 * t85 + t161 * t84); 0; t175; m(7) * (t53 * t44 + t54 * t45) + m(6) * (t75 * t51 + t76 * t52) + (t37 / 0.2e1 + t26 / 0.2e1 - t78 * t211 + t120 * t212 - t184 + t219 * t161) * t161 + (t38 / 0.2e1 + t27 / 0.2e1 + t77 * t211 + t120 * t213 + t183 + t219 * t159) * t159; m(6) * (t75 * t159 - t76 * t161) + m(7) * (t53 * t159 - t54 * t161) + t123 * t128; m(6) * (t76 * t159 + t75 * t161) + m(7) * (t54 * t159 + t53 * t161); m(7) * (t16 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(6) * (t36 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t137 * t123 ^ 2 + t64 ^ 2) + (t152 * t97 + t12 + t15) * t159 + (t153 * t96 + t11 + t14 + (t159 * t96 + t161 * t97) * t159) * t161; t46 + m(7) * (t24 * t45 + t25 * t44) + m(6) * (t49 * t52 + t50 * t51) + (t184 * t159 + t183 * t161) * t144 + t174; m(6) * (t50 * t159 - t49 * t161) + m(7) * (t25 * t159 - t24 * t161); m(6) * (t49 * t159 + t50 * t161) + m(7) * (t24 * t159 + t25 * t161); (t33 * t159 + t32 * t161) * t215 + t8 * t213 + t7 * t212 + (t15 * t212 - t159 * t14 / 0.2e1) * t144 + m(7) * (t17 * t16 + t24 * t54 + t25 * t53) + m(6) * (t43 * t36 + t49 * t76 + t50 * t75) + t176; t143 * t46 + m(7) * (t17 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t43 ^ 2 + t49 ^ 2 + t50 ^ 2) + ((t143 * t33 + t8) * t161 + (-t143 * t32 - t4 - t7) * t159) * t144 + t217; m(7) * (t48 * t44 + t47 * t45) + t174; m(7) * (t48 * t159 - t47 * t161); m(7) * (t47 * t159 + t48 * t161); m(7) * (t39 * t16 + t47 * t54 + t48 * t53) + t176; m(7) * (t39 * t17 + t47 * t24 + t48 * t25) + t173; m(7) * (t39 ^ 2 + t47 ^ 2 + t48 ^ 2) + t173;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
