% Calculate joint inertia matrix for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR14_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR14_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:56
% EndTime: 2019-12-31 19:17:07
% DurationCPUTime: 3.75s
% Computational Cost: add. (21813->403), mult. (60471->586), div. (0->0), fcn. (79686->14), ass. (0->193)
t184 = sin(qJ(1));
t180 = cos(pkin(5));
t220 = cos(pkin(11));
t206 = t184 * t220;
t178 = sin(pkin(11));
t186 = cos(qJ(1));
t214 = t186 * t178;
t195 = t180 * t206 + t214;
t179 = sin(pkin(5));
t221 = cos(pkin(6));
t208 = t179 * t221;
t219 = sin(pkin(6));
t157 = t184 * t208 + t195 * t219;
t205 = t186 * t220;
t215 = t184 * t178;
t196 = -t180 * t205 + t215;
t156 = -t186 * t208 + t196 * t219;
t237 = m(3) / 0.2e1;
t236 = m(4) / 0.2e1;
t235 = m(5) / 0.2e1;
t234 = m(6) / 0.2e1;
t165 = t180 * t214 + t206;
t183 = sin(qJ(3));
t193 = t196 * t221;
t207 = t179 * t219;
t227 = cos(qJ(3));
t148 = t165 * t227 + (-t186 * t207 - t193) * t183;
t182 = sin(qJ(4));
t226 = cos(qJ(4));
t129 = t148 * t182 - t156 * t226;
t166 = -t180 * t215 + t205;
t191 = t195 * t221;
t150 = t166 * t227 + (t184 * t207 - t191) * t183;
t131 = t150 * t182 - t157 * t226;
t200 = t221 * t220;
t155 = t180 * t219 * t183 + (t178 * t227 + t183 * t200) * t179;
t164 = t180 * t221 - t207 * t220;
t145 = t155 * t182 - t164 * t226;
t130 = t148 * t226 + t156 * t182;
t202 = t227 * t219;
t199 = t179 * t202;
t147 = t165 * t183 + t186 * t199 + t193 * t227;
t181 = sin(qJ(5));
t185 = cos(qJ(5));
t98 = -t130 * t181 + t147 * t185;
t99 = t130 * t185 + t147 * t181;
t60 = Icges(6,5) * t99 + Icges(6,6) * t98 + Icges(6,3) * t129;
t62 = Icges(6,4) * t99 + Icges(6,2) * t98 + Icges(6,6) * t129;
t64 = Icges(6,1) * t99 + Icges(6,4) * t98 + Icges(6,5) * t129;
t16 = t129 * t60 + t62 * t98 + t64 * t99;
t132 = t150 * t226 + t157 * t182;
t149 = t166 * t183 - t184 * t199 + t191 * t227;
t100 = -t132 * t181 + t149 * t185;
t101 = t132 * t185 + t149 * t181;
t61 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t131;
t63 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t131;
t65 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t131;
t17 = t129 * t61 + t63 * t98 + t65 * t99;
t146 = t155 * t226 + t164 * t182;
t218 = t178 * t179;
t154 = -t179 * t200 * t227 - t180 * t202 + t183 * t218;
t123 = -t146 * t181 + t154 * t185;
t124 = t146 * t185 + t154 * t181;
t77 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t145;
t78 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t145;
t79 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t145;
t26 = t129 * t77 + t78 * t98 + t79 * t99;
t1 = t129 * t16 + t131 * t17 + t145 * t26;
t233 = t1 / 0.2e1;
t18 = t100 * t62 + t101 * t64 + t131 * t60;
t19 = t100 * t63 + t101 * t65 + t131 * t61;
t27 = t100 * t78 + t101 * t79 + t131 * t77;
t2 = t129 * t18 + t131 * t19 + t145 * t27;
t232 = t2 / 0.2e1;
t21 = t123 * t62 + t124 * t64 + t145 * t60;
t22 = t123 * t63 + t124 * t65 + t145 * t61;
t33 = t123 * t78 + t124 * t79 + t145 * t77;
t30 = t33 * t145;
t7 = t21 * t129 + t22 * t131 + t30;
t231 = t7 / 0.2e1;
t230 = t129 / 0.2e1;
t229 = t131 / 0.2e1;
t228 = t145 / 0.2e1;
t225 = pkin(4) * t130;
t201 = -rSges(6,1) * t99 - rSges(6,2) * t98;
t66 = rSges(6,3) * t129 - t201;
t224 = pkin(10) * t129 + t225 + t66;
t67 = rSges(6,1) * t101 + rSges(6,2) * t100 + rSges(6,3) * t131;
t223 = pkin(4) * t132 + pkin(10) * t131 + t67;
t80 = rSges(6,1) * t124 + rSges(6,2) * t123 + rSges(6,3) * t145;
t222 = pkin(4) * t146 + pkin(10) * t145 + t80;
t217 = t179 * t184;
t216 = t179 * t186;
t213 = t186 * pkin(1) + qJ(2) * t217;
t102 = Icges(5,5) * t146 - Icges(5,6) * t145 + Icges(5,3) * t154;
t103 = Icges(5,4) * t146 - Icges(5,2) * t145 + Icges(5,6) * t154;
t104 = Icges(5,1) * t146 - Icges(5,4) * t145 + Icges(5,5) * t154;
t53 = t102 * t154 - t103 * t145 + t104 * t146;
t212 = t21 / 0.2e1 + t26 / 0.2e1;
t211 = t27 / 0.2e1 + t22 / 0.2e1;
t133 = Icges(4,5) * t155 - Icges(4,6) * t154 + Icges(4,3) * t164;
t134 = Icges(4,4) * t155 - Icges(4,2) * t154 + Icges(4,6) * t164;
t135 = Icges(4,1) * t155 - Icges(4,4) * t154 + Icges(4,5) * t164;
t210 = t133 * t164 - t134 * t154 + t135 * t155;
t88 = rSges(5,1) * t132 - rSges(5,2) * t131 + t149 * rSges(5,3);
t113 = t150 * rSges(4,1) - t149 * rSges(4,2) + t157 * rSges(4,3);
t209 = -t184 * pkin(1) + qJ(2) * t216;
t119 = pkin(3) * t148 + t147 * pkin(9);
t120 = t150 * pkin(3) + pkin(9) * t149;
t81 = Icges(5,5) * t130 - Icges(5,6) * t129 + Icges(5,3) * t147;
t83 = Icges(5,4) * t130 - Icges(5,2) * t129 + Icges(5,6) * t147;
t85 = Icges(5,1) * t130 - Icges(5,4) * t129 + Icges(5,5) * t147;
t39 = -t145 * t83 + t146 * t85 + t154 * t81;
t46 = t102 * t147 - t103 * t129 + t104 * t130;
t198 = t46 / 0.2e1 + t39 / 0.2e1 + t212;
t82 = Icges(5,5) * t132 - Icges(5,6) * t131 + Icges(5,3) * t149;
t84 = Icges(5,4) * t132 - Icges(5,2) * t131 + Icges(5,6) * t149;
t86 = Icges(5,1) * t132 - Icges(5,4) * t131 + Icges(5,5) * t149;
t40 = -t145 * t84 + t146 * t86 + t154 * t82;
t47 = t102 * t149 - t103 * t131 + t104 * t132;
t197 = t40 / 0.2e1 + t47 / 0.2e1 + t211;
t112 = rSges(4,1) * t148 - rSges(4,2) * t147 + rSges(4,3) * t156;
t87 = rSges(5,1) * t130 - rSges(5,2) * t129 + rSges(5,3) * t147;
t194 = -pkin(2) * t165 - pkin(8) * t156 + t209;
t189 = -t119 + t194;
t188 = t166 * pkin(2) + pkin(8) * t157 + t213;
t187 = t120 + t188;
t172 = rSges(2,1) * t186 - rSges(2,2) * t184;
t171 = -rSges(2,1) * t184 - rSges(2,2) * t186;
t144 = t166 * rSges(3,1) - rSges(3,2) * t195 + rSges(3,3) * t217 + t213;
t143 = -t165 * rSges(3,1) + rSges(3,2) * t196 + rSges(3,3) * t216 + t209;
t137 = pkin(3) * t155 + pkin(9) * t154;
t136 = rSges(4,1) * t155 - rSges(4,2) * t154 + rSges(4,3) * t164;
t122 = t156 * t137;
t115 = t164 * t120;
t114 = t157 * t119;
t111 = Icges(4,1) * t150 - Icges(4,4) * t149 + Icges(4,5) * t157;
t110 = Icges(4,1) * t148 - Icges(4,4) * t147 + Icges(4,5) * t156;
t109 = Icges(4,4) * t150 - Icges(4,2) * t149 + Icges(4,6) * t157;
t108 = Icges(4,4) * t148 - Icges(4,2) * t147 + Icges(4,6) * t156;
t107 = Icges(4,5) * t150 - Icges(4,6) * t149 + Icges(4,3) * t157;
t106 = Icges(4,5) * t148 - Icges(4,6) * t147 + Icges(4,3) * t156;
t105 = rSges(5,1) * t146 - rSges(5,2) * t145 + rSges(5,3) * t154;
t92 = t188 + t113;
t91 = -t112 + t194;
t76 = t113 * t164 - t136 * t157;
t75 = -t112 * t164 + t136 * t156;
t71 = t112 * t157 - t113 * t156;
t70 = t210 * t164;
t69 = t133 * t157 - t134 * t149 + t135 * t150;
t68 = t133 * t156 - t134 * t147 + t135 * t148;
t59 = t187 + t88;
t58 = t189 - t87;
t57 = -t105 * t149 + t154 * t88;
t56 = t105 * t147 - t154 * t87;
t55 = t107 * t164 - t109 * t154 + t111 * t155;
t54 = t106 * t164 - t108 * t154 + t110 * t155;
t52 = -t147 * t88 + t149 * t87;
t51 = t53 * t164;
t50 = t53 * t154;
t49 = t164 * t88 + t115 + (-t105 - t137) * t157;
t48 = t105 * t156 + t122 + (-t119 - t87) * t164;
t45 = t187 + t223;
t44 = -t225 + (-rSges(6,3) - pkin(10)) * t129 + t189 + t201;
t43 = -t131 * t80 + t145 * t67;
t42 = t129 * t80 - t145 * t66;
t41 = t157 * t87 + t114 + (-t120 - t88) * t156;
t38 = -t131 * t84 + t132 * t86 + t149 * t82;
t37 = -t131 * t83 + t132 * t85 + t149 * t81;
t36 = -t129 * t84 + t130 * t86 + t147 * t82;
t35 = -t129 * t83 + t130 * t85 + t147 * t81;
t34 = -t129 * t67 + t131 * t66;
t32 = t33 * t164;
t31 = t33 * t154;
t29 = -t149 * t222 + t154 * t223;
t28 = t147 * t222 - t154 * t224;
t25 = t115 + t223 * t164 + (-t137 - t222) * t157;
t24 = t122 + t222 * t156 + (-t119 - t224) * t164;
t23 = -t147 * t223 + t149 * t224;
t20 = t114 + t224 * t157 + (-t120 - t223) * t156;
t15 = t39 * t156 + t40 * t157 + t51;
t14 = t39 * t147 + t40 * t149 + t50;
t13 = t156 * t37 + t157 * t38 + t164 * t47;
t12 = t156 * t35 + t157 * t36 + t164 * t46;
t11 = t147 * t37 + t149 * t38 + t154 * t47;
t10 = t147 * t35 + t149 * t36 + t154 * t46;
t9 = t21 * t156 + t22 * t157 + t32;
t8 = t21 * t147 + t22 * t149 + t31;
t6 = t156 * t18 + t157 * t19 + t164 * t27;
t5 = t156 * t16 + t157 * t17 + t164 * t26;
t4 = t147 * t18 + t149 * t19 + t154 * t27;
t3 = t147 * t16 + t149 * t17 + t154 * t26;
t72 = [m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t58 ^ 2 + t59 ^ 2) + m(4) * (t91 ^ 2 + t92 ^ 2) + m(3) * (t143 ^ 2 + t144 ^ 2) + m(2) * (t171 ^ 2 + t172 ^ 2) + Icges(2,3) + (Icges(3,5) * t180 + (Icges(3,1) * t178 + Icges(3,4) * t220) * t179) * t218 + t179 * t220 * (Icges(3,6) * t180 + (Icges(3,4) * t178 + Icges(3,2) * t220) * t179) + t180 * (Icges(3,3) * t180 + (Icges(3,5) * t178 + Icges(3,6) * t220) * t179) + t210 + t53 + t33; 0.2e1 * ((t184 * t44 - t186 * t45) * t234 + (t184 * t58 - t186 * t59) * t235 + (t184 * t91 - t186 * t92) * t236 + (t143 * t184 - t144 * t186) * t237) * t179; 0.2e1 * (t237 + t236 + t235 + t234) * (t180 ^ 2 + (t184 ^ 2 + t186 ^ 2) * t179 ^ 2); t32 + t51 + t70 + m(6) * (t24 * t44 + t25 * t45) + m(5) * (t48 * t58 + t49 * t59) + m(4) * (t75 * t91 + t76 * t92) + (t55 / 0.2e1 + t69 / 0.2e1 + t197) * t157 + (t68 / 0.2e1 + t54 / 0.2e1 + t198) * t156; m(4) * (t71 * t180 + (t184 * t75 - t186 * t76) * t179) + m(5) * (t41 * t180 + (t184 * t48 - t186 * t49) * t179) + m(6) * (t20 * t180 + (t184 * t24 - t186 * t25) * t179); (t15 + t70 + t9) * t164 + (t13 + t6 + (t107 * t157 - t109 * t149 + t111 * t150) * t157 + (t55 + t69) * t164) * t157 + (t5 + t12 + (t106 * t156 - t108 * t147 + t110 * t148) * t156 + (t54 + t68) * t164 + (t106 * t157 + t107 * t156 - t108 * t149 - t109 * t147 + t110 * t150 + t111 * t148) * t157) * t156 + m(6) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t41 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(4) * (t71 ^ 2 + t75 ^ 2 + t76 ^ 2); t31 + t50 + m(6) * (t28 * t44 + t29 * t45) + m(5) * (t56 * t58 + t57 * t59) + t197 * t149 + t198 * t147; m(5) * (t52 * t180 + (t184 * t56 - t186 * t57) * t179) + m(6) * (t23 * t180 + (t184 * t28 - t186 * t29) * t179); (t8 / 0.2e1 + t14 / 0.2e1) * t164 + (t4 / 0.2e1 + t11 / 0.2e1) * t157 + (t3 / 0.2e1 + t10 / 0.2e1) * t156 + (t9 / 0.2e1 + t15 / 0.2e1) * t154 + (t6 / 0.2e1 + t13 / 0.2e1) * t149 + (t5 / 0.2e1 + t12 / 0.2e1) * t147 + m(6) * (t20 * t23 + t24 * t28 + t25 * t29) + m(5) * (t41 * t52 + t48 * t56 + t49 * t57); (t8 + t14) * t154 + (t4 + t11) * t149 + (t3 + t10) * t147 + m(6) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(5) * (t52 ^ 2 + t56 ^ 2 + t57 ^ 2); m(6) * (t42 * t44 + t43 * t45) + t30 + t211 * t131 + t212 * t129; m(6) * (t34 * t180 + (t184 * t42 - t186 * t43) * t179); t5 * t230 + t9 * t228 + t157 * t232 + t156 * t233 + t6 * t229 + m(6) * (t20 * t34 + t24 * t42 + t25 * t43) + t164 * t231; m(6) * (t23 * t34 + t28 * t42 + t29 * t43) + t3 * t230 + t4 * t229 + t149 * t232 + t147 * t233 + t154 * t231 + t8 * t228; m(6) * (t34 ^ 2 + t42 ^ 2 + t43 ^ 2) + t131 * t2 + t129 * t1 + t145 * t7;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t72(1), t72(2), t72(4), t72(7), t72(11); t72(2), t72(3), t72(5), t72(8), t72(12); t72(4), t72(5), t72(6), t72(9), t72(13); t72(7), t72(8), t72(9), t72(10), t72(14); t72(11), t72(12), t72(13), t72(14), t72(15);];
Mq = res;
