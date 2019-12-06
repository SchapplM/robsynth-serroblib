% Calculate joint inertia matrix for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:37:48
% EndTime: 2019-12-05 18:37:53
% DurationCPUTime: 1.94s
% Computational Cost: add. (3662->265), mult. (2982->378), div. (0->0), fcn. (2730->10), ass. (0->141)
t214 = Icges(4,3) + Icges(5,3);
t132 = qJ(2) + qJ(3);
t119 = pkin(9) + t132;
t114 = sin(t119);
t115 = cos(t119);
t120 = sin(t132);
t121 = cos(t132);
t213 = Icges(4,5) * t121 + Icges(5,5) * t115 - Icges(4,6) * t120 - Icges(5,6) * t114;
t134 = sin(qJ(1));
t136 = cos(qJ(1));
t212 = -t134 * t213 + t136 * t214;
t211 = t134 * t214 + t136 * t213;
t183 = Icges(5,4) * t115;
t151 = -Icges(5,2) * t114 + t183;
t61 = -Icges(5,6) * t136 + t134 * t151;
t184 = Icges(5,4) * t114;
t155 = Icges(5,1) * t115 - t184;
t63 = -Icges(5,5) * t136 + t134 * t155;
t185 = Icges(4,4) * t121;
t152 = -Icges(4,2) * t120 + t185;
t71 = -Icges(4,6) * t136 + t134 * t152;
t186 = Icges(4,4) * t120;
t156 = Icges(4,1) * t121 - t186;
t73 = -Icges(4,5) * t136 + t134 * t156;
t210 = t114 * t61 - t115 * t63 + t120 * t71 - t121 * t73;
t62 = Icges(5,6) * t134 + t136 * t151;
t64 = Icges(5,5) * t134 + t136 * t155;
t72 = Icges(4,6) * t134 + t136 * t152;
t74 = Icges(4,5) * t134 + t136 * t156;
t209 = t114 * t62 - t115 * t64 + t120 * t72 - t121 * t74;
t130 = t134 ^ 2;
t208 = t134 * pkin(6);
t207 = Icges(4,5) * t120 + Icges(5,5) * t114 + Icges(4,6) * t121 + Icges(5,6) * t115;
t89 = Icges(5,2) * t115 + t184;
t90 = Icges(5,1) * t114 + t183;
t95 = Icges(4,2) * t121 + t186;
t96 = Icges(4,1) * t120 + t185;
t206 = -t114 * t89 + t115 * t90 - t120 * t95 + t121 * t96;
t117 = qJ(5) + t119;
t110 = sin(t117);
t111 = cos(t117);
t169 = rSges(6,1) * t111 - rSges(6,2) * t110;
t170 = rSges(5,1) * t115 - rSges(5,2) * t114;
t171 = rSges(4,1) * t121 - rSges(4,2) * t120;
t205 = (t211 * t130 + ((-t209 + t212) * t134 + t210 * t136) * t136) * t134;
t131 = t136 ^ 2;
t204 = t212 * t131 + (t209 * t134 + (-t210 + t211) * t136) * t134;
t137 = -pkin(7) - pkin(6);
t203 = t134 / 0.2e1;
t202 = -t136 / 0.2e1;
t133 = sin(qJ(2));
t201 = pkin(2) * t133;
t200 = pkin(3) * t120;
t135 = cos(qJ(2));
t118 = pkin(2) * t135 + pkin(1);
t108 = t136 * t118;
t98 = pkin(3) * t121 + t118;
t92 = t136 * t98;
t199 = t136 * (-t108 + t92) + (-t118 + t98) * t130;
t145 = t134 * rSges(6,3) + t136 * t169;
t22 = t134 * (-t136 * rSges(6,3) + t134 * t169) + t136 * t145;
t128 = t136 * pkin(6);
t198 = t134 * (t128 + (-pkin(1) + t118) * t134) + t136 * (-t136 * pkin(1) + t108 - t208);
t144 = t134 * rSges(4,3) + t136 * t171;
t29 = t134 * (-t136 * rSges(4,3) + t134 * t171) + t136 * t144;
t197 = rSges(3,1) * t135;
t193 = rSges(3,2) * t133;
t189 = t136 * rSges(3,3);
t188 = Icges(3,4) * t133;
t187 = Icges(3,4) * t135;
t182 = Icges(6,4) * t110;
t181 = Icges(6,4) * t111;
t180 = rSges(3,3) * t134 + t136 * t197;
t177 = t130 + t131;
t129 = -qJ(4) + t137;
t150 = -Icges(6,2) * t110 + t181;
t52 = Icges(6,6) * t134 + t136 * t150;
t154 = Icges(6,1) * t111 - t182;
t54 = Icges(6,5) * t134 + t136 * t154;
t167 = -t110 * t52 + t111 * t54;
t51 = -Icges(6,6) * t136 + t134 * t150;
t53 = -Icges(6,5) * t136 + t134 * t154;
t168 = t110 * t51 - t111 * t53;
t146 = Icges(6,5) * t111 - Icges(6,6) * t110;
t49 = -Icges(6,3) * t136 + t134 * t146;
t50 = Icges(6,3) * t134 + t136 * t146;
t3 = t134 * (t130 * t50 + (t168 * t136 + (t167 - t49) * t134) * t136);
t4 = t131 * t49 + (t167 * t134 + (t168 - t50) * t136) * t134;
t176 = -t136 * t4 + t3;
t97 = rSges(4,1) * t120 + rSges(4,2) * t121;
t175 = -t97 - t201;
t174 = -rSges(5,1) * t114 - rSges(5,2) * t115 - t200;
t85 = Icges(6,2) * t111 + t182;
t86 = Icges(6,1) * t110 + t181;
t166 = -t110 * t85 + t111 * t86;
t84 = Icges(6,5) * t110 + Icges(6,6) * t111;
t173 = (t110 * t54 + t111 * t52 + t134 * t84 + t136 * t166) * t203 + (t110 * t53 + t111 * t51 + t134 * t166 - t136 * t84) * t202;
t143 = t134 * rSges(5,3) + t136 * t170;
t10 = t134 * (-t136 * rSges(5,3) + t134 * t170) + t136 * t143 + t199;
t172 = -t193 + t197;
t76 = pkin(4) * t115 + t98;
t75 = t136 * t76;
t2 = t136 * (t75 - t92) + t22 + t199 + (t76 - t98) * t130;
t157 = Icges(3,1) * t135 - t188;
t153 = -Icges(3,2) * t133 + t187;
t149 = Icges(3,5) * t135 - Icges(3,6) * t133;
t142 = t174 - t201;
t87 = rSges(6,1) * t110 + rSges(6,2) * t111;
t141 = -pkin(4) * t114 - t200 - t87;
t140 = t3 + (-t4 + t204) * t136 + t205;
t139 = t141 - t201;
t138 = t173 + (t114 * t64 + t115 * t62 + t120 * t74 + t121 * t72 + t134 * t207 + t136 * t206) * t203 + (t114 * t63 + t115 * t61 + t120 * t73 + t121 * t71 + t134 * t206 - t136 * t207) * t202;
t122 = -pkin(8) + t129;
t107 = rSges(2,1) * t136 - rSges(2,2) * t134;
t106 = -rSges(2,1) * t134 - rSges(2,2) * t136;
t105 = rSges(3,1) * t133 + rSges(3,2) * t135;
t78 = Icges(3,3) * t134 + t136 * t149;
t77 = -Icges(3,3) * t136 + t134 * t149;
t68 = t175 * t136;
t67 = t175 * t134;
t56 = t208 + (pkin(1) - t193) * t136 + t180;
t55 = t189 + t128 + (-pkin(1) - t172) * t134;
t48 = t174 * t136;
t47 = t174 * t134;
t42 = t142 * t136;
t41 = t142 * t134;
t40 = -t134 * t137 + t108 + t144;
t39 = (rSges(4,3) - t137) * t136 + (-t118 - t171) * t134;
t36 = t141 * t136;
t35 = t141 * t134;
t34 = t136 * (-t136 * t193 + t180) + (t134 * t172 - t189) * t134;
t33 = -t129 * t134 + t143 + t92;
t32 = (rSges(5,3) - t129) * t136 + (-t170 - t98) * t134;
t31 = t139 * t136;
t30 = t139 * t134;
t26 = -t122 * t134 + t145 + t75;
t25 = (rSges(6,3) - t122) * t136 + (-t169 - t76) * t134;
t11 = t29 + t198;
t7 = t10 + t198;
t1 = t2 + t198;
t5 = [t135 * (Icges(3,2) * t135 + t188) + t133 * (Icges(3,1) * t133 + t187) + t110 * t86 + t111 * t85 + t114 * t90 + t115 * t89 + t120 * t96 + t121 * t95 + Icges(2,3) + m(5) * (t32 ^ 2 + t33 ^ 2) + m(6) * (t25 ^ 2 + t26 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(3) * (t55 ^ 2 + t56 ^ 2) + m(2) * (t106 ^ 2 + t107 ^ 2); (t131 / 0.2e1 + t130 / 0.2e1) * (Icges(3,5) * t133 + Icges(3,6) * t135) + m(5) * (t32 * t42 + t33 * t41) + m(6) * (t25 * t31 + t26 * t30) + m(4) * (t39 * t68 + t40 * t67) + m(3) * (-t134 * t56 - t136 * t55) * t105 + (t133 * (-Icges(3,5) * t136 + t134 * t157) + t135 * (-Icges(3,6) * t136 + t134 * t153)) * t202 + (t133 * (Icges(3,5) * t134 + t136 * t157) + t135 * (Icges(3,6) * t134 + t136 * t153)) * t203 + t138; m(6) * (t1 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t41 ^ 2 + t42 ^ 2 + t7 ^ 2) + m(4) * (t11 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(3) * (t105 ^ 2 * t177 + t34 ^ 2) + t134 * t130 * t78 + t176 + (-t131 * t77 + (-t134 * t77 + t136 * t78) * t134 + t204) * t136 + t205; m(4) * (-t134 * t40 - t136 * t39) * t97 + m(5) * (t32 * t48 + t33 * t47) + m(6) * (t25 * t36 + t26 * t35) + t138; m(6) * (t1 * t2 + t30 * t35 + t31 * t36) + m(5) * (t10 * t7 + t41 * t47 + t42 * t48) + m(4) * (t29 * t11 + (-t134 * t67 - t136 * t68) * t97) + t140; m(6) * (t2 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t10 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(4) * (t177 * t97 ^ 2 + t29 ^ 2) + t140; m(5) * (t134 * t32 - t136 * t33) + m(6) * (t134 * t25 - t136 * t26); m(6) * (t134 * t31 - t136 * t30) + m(5) * (t134 * t42 - t136 * t41); m(6) * (t134 * t36 - t136 * t35) + m(5) * (t134 * t48 - t136 * t47); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t177; m(6) * (-t134 * t26 - t136 * t25) * t87 + t173; m(6) * (t22 * t1 + (-t134 * t30 - t136 * t31) * t87) + t176; m(6) * (t22 * t2 + (-t134 * t35 - t136 * t36) * t87) + t176; 0; m(6) * (t177 * t87 ^ 2 + t22 ^ 2) + t176;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
