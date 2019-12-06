% Calculate joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:31
% EndTime: 2019-12-05 18:05:41
% DurationCPUTime: 2.82s
% Computational Cost: add. (4343->321), mult. (6023->453), div. (0->0), fcn. (6406->8), ass. (0->154)
t143 = qJ(3) + qJ(4);
t135 = sin(t143);
t136 = cos(t143);
t144 = sin(pkin(8));
t145 = cos(pkin(8));
t102 = -Icges(6,5) * t145 + (Icges(6,1) * t136 - Icges(6,4) * t135) * t144;
t103 = -Icges(5,5) * t145 + (Icges(5,1) * t136 - Icges(5,4) * t135) * t144;
t213 = t102 + t103;
t149 = cos(qJ(1));
t147 = sin(qJ(1));
t179 = t145 * t147;
t114 = t135 * t179 + t136 * t149;
t115 = t135 * t149 - t136 * t179;
t174 = t147 * t144;
t56 = Icges(6,5) * t115 + Icges(6,6) * t114 - Icges(6,3) * t174;
t58 = Icges(5,5) * t115 + Icges(5,6) * t114 - Icges(5,3) * t174;
t212 = -t56 - t58;
t178 = t145 * t149;
t116 = -t135 * t178 + t136 * t147;
t117 = t135 * t147 + t136 * t178;
t170 = t149 * t144;
t57 = Icges(6,5) * t117 + Icges(6,6) * t116 + Icges(6,3) * t170;
t59 = Icges(5,5) * t117 + Icges(5,6) * t116 + Icges(5,3) * t170;
t211 = t57 + t59;
t60 = Icges(6,4) * t115 + Icges(6,2) * t114 - Icges(6,6) * t174;
t62 = Icges(5,4) * t115 + Icges(5,2) * t114 - Icges(5,6) * t174;
t210 = -t60 - t62;
t61 = Icges(6,4) * t117 + Icges(6,2) * t116 + Icges(6,6) * t170;
t63 = Icges(5,4) * t117 + Icges(5,2) * t116 + Icges(5,6) * t170;
t209 = t61 + t63;
t64 = Icges(6,1) * t115 + Icges(6,4) * t114 - Icges(6,5) * t174;
t66 = Icges(5,1) * t115 + Icges(5,4) * t114 - Icges(5,5) * t174;
t208 = t64 + t66;
t65 = Icges(6,1) * t117 + Icges(6,4) * t116 + Icges(6,5) * t170;
t67 = Icges(5,1) * t117 + Icges(5,4) * t116 + Icges(5,5) * t170;
t207 = t65 + t67;
t98 = -Icges(6,3) * t145 + (Icges(6,5) * t136 - Icges(6,6) * t135) * t144;
t99 = -Icges(5,3) * t145 + (Icges(5,5) * t136 - Icges(5,6) * t135) * t144;
t206 = -t98 - t99;
t150 = -pkin(7) - pkin(6);
t139 = -qJ(5) + t150;
t205 = -rSges(6,3) + t139;
t201 = t213 * t136 * t144;
t100 = -Icges(6,6) * t145 + (Icges(6,4) * t136 - Icges(6,2) * t135) * t144;
t101 = -Icges(5,6) * t145 + (Icges(5,4) * t136 - Icges(5,2) * t135) * t144;
t203 = -t100 - t101;
t202 = t203 * t135;
t197 = t144 * t202 + t206 * t145 + t201;
t204 = t197 * t145;
t200 = t212 * t145 + (t210 * t135 + t208 * t136) * t144;
t199 = t211 * t145 + (t209 * t135 - t207 * t136) * t144;
t198 = t203 * t116 - t213 * t117 + t206 * t170;
t146 = sin(qJ(3));
t190 = pkin(3) * t146;
t126 = pkin(4) * t135 + t190;
t131 = t150 * t170;
t156 = -rSges(6,1) * t117 - rSges(6,2) * t116;
t148 = cos(qJ(3));
t134 = t148 * pkin(3) + pkin(2);
t125 = pkin(4) * t136 + t134;
t165 = t125 - t134;
t196 = t131 + (t126 - t190) * t147 + (-t139 * t144 + t165 * t145) * t149 + rSges(6,3) * t170 - t156;
t195 = (-t150 + t205) * t145 + (rSges(6,1) * t136 - rSges(6,2) * t135 + t165) * t144;
t194 = t115 * rSges(6,1) + t114 * rSges(6,2) + t149 * t126 + t139 * t174;
t188 = pkin(2) - t134;
t193 = pkin(6) * t144 + t188 * t145;
t141 = t147 ^ 2;
t142 = t149 ^ 2;
t192 = ((t209 * t116 + t207 * t117 + t211 * t170) * t170 + t198 * t145 + (t210 * t116 - t208 * t117 + t212 * t170) * t174) * t170;
t191 = t204 + (t200 * t147 + t199 * t149) * t144;
t175 = t146 * t149;
t164 = pkin(3) * t175 + t150 * t174;
t186 = rSges(6,3) * t174 + t165 * t179 + t164 - t194;
t105 = -rSges(5,3) * t145 + (rSges(5,1) * t136 - rSges(5,2) * t135) * t144;
t157 = -rSges(5,1) * t117 - rSges(5,2) * t116;
t71 = rSges(5,3) * t170 - t157;
t45 = t105 * t170 + t145 * t71;
t168 = t115 * rSges(5,1) + t114 * rSges(5,2);
t69 = -rSges(5,3) * t174 + t168;
t81 = t193 * t147 + t164;
t184 = -t69 - t81;
t176 = t146 * t147;
t82 = pkin(3) * t176 - t193 * t149 - t131;
t97 = (pkin(6) + t150) * t145 - t188 * t144;
t183 = t145 * t82 + t97 * t170;
t182 = t195 * t174;
t111 = -Icges(4,6) * t145 + (Icges(4,4) * t148 - Icges(4,2) * t146) * t144;
t177 = t146 * t111;
t173 = t147 * t148;
t172 = t147 * t149;
t171 = t148 * t149;
t120 = t145 * t176 + t171;
t121 = -t145 * t173 + t175;
t167 = t121 * rSges(4,1) + t120 * rSges(4,2);
t163 = t142 + t141;
t162 = -t81 + t186;
t159 = -t125 * t145 - pkin(1);
t13 = t196 * t145 + t195 * t170;
t122 = -t145 * t175 + t173;
t123 = t145 * t171 + t176;
t158 = -rSges(4,1) * t123 - rSges(4,2) * t122;
t155 = -rSges(3,1) * t145 + rSges(3,2) * t144 - pkin(1);
t154 = -rSges(5,3) * t144 - t134 * t145 - pkin(1);
t153 = -pkin(2) * t145 - pkin(1) + (-rSges(4,3) - pkin(6)) * t144;
t30 = t100 * t114 + t102 * t115 - t98 * t174;
t31 = t101 * t114 + t103 * t115 - t99 * t174;
t152 = -(t30 + t31 + t200) * t174 / 0.2e1 + (-t198 - t199) * t170 / 0.2e1;
t3 = -t30 * t145 - (t114 * t60 + t115 * t64 - t56 * t174) * t174 + (t114 * t61 + t115 * t65 - t57 * t174) * t170;
t4 = -t31 * t145 - (t114 * t62 + t115 * t66 - t58 * t174) * t174 + (t114 * t63 + t115 * t67 - t59 * t174) * t170;
t151 = t191 * t145 + (-t3 - t4) * t174 + t192;
t137 = t149 * qJ(2);
t129 = -rSges(2,1) * t149 + rSges(2,2) * t147;
t128 = -rSges(2,1) * t147 - rSges(2,2) * t149;
t113 = -rSges(4,3) * t145 + (rSges(4,1) * t148 - rSges(4,2) * t146) * t144;
t112 = -Icges(4,5) * t145 + (Icges(4,1) * t148 - Icges(4,4) * t146) * t144;
t110 = -Icges(4,3) * t145 + (Icges(4,5) * t148 - Icges(4,6) * t146) * t144;
t96 = t144 * t148 * t112;
t95 = (-rSges(3,3) - qJ(2)) * t147 + t155 * t149;
t94 = rSges(3,3) * t149 + t155 * t147 + t137;
t91 = t105 * t174;
t88 = t97 * t174;
t84 = rSges(4,3) * t170 - t158;
t83 = -rSges(4,3) * t174 + t167;
t80 = Icges(4,1) * t123 + Icges(4,4) * t122 + Icges(4,5) * t170;
t79 = Icges(4,1) * t121 + Icges(4,4) * t120 - Icges(4,5) * t174;
t78 = Icges(4,4) * t123 + Icges(4,2) * t122 + Icges(4,6) * t170;
t77 = Icges(4,4) * t121 + Icges(4,2) * t120 - Icges(4,6) * t174;
t76 = Icges(4,5) * t123 + Icges(4,6) * t122 + Icges(4,3) * t170;
t75 = Icges(4,5) * t121 + Icges(4,6) * t120 - Icges(4,3) * t174;
t53 = -qJ(2) * t147 + t149 * t153 + t158;
t52 = t147 * t153 + t137 + t167;
t49 = t113 * t170 + t145 * t84;
t48 = t113 * t174 - t145 * t83;
t46 = -t145 * t110 - t144 * t177 + t96;
t44 = -t145 * t69 + t91;
t43 = t131 + (-qJ(2) - t190) * t147 + t154 * t149 + t157;
t42 = t154 * t147 + t137 + t164 + t168;
t39 = (-t147 * t84 - t149 * t83) * t144;
t38 = t110 * t170 + t111 * t122 + t112 * t123;
t37 = -t110 * t174 + t111 * t120 + t112 * t121;
t36 = (-qJ(2) - t126) * t147 + (t205 * t144 + t159) * t149 + t156;
t35 = t137 + (-rSges(6,3) * t144 + t159) * t147 + t194;
t34 = (-t147 * t71 - t149 * t69) * t144;
t25 = -t145 * t76 + (-t146 * t78 + t148 * t80) * t144;
t24 = -t145 * t75 + (-t146 * t77 + t148 * t79) * t144;
t23 = t183 + t45;
t22 = t184 * t145 + t88 + t91;
t12 = t186 * t145 + t182;
t11 = (t184 * t149 + (-t71 - t82) * t147) * t144;
t10 = (-t147 * t196 + t186 * t149) * t144;
t9 = t13 + t183;
t8 = t162 * t145 + t182 + t88;
t7 = (t162 * t149 + (-t82 - t196) * t147) * t144;
t1 = [Icges(2,3) + t96 + (Icges(3,2) * t145 - t110 + t206) * t145 + (Icges(3,1) * t144 + 0.2e1 * Icges(3,4) * t145 - t177 + t202) * t144 + m(6) * (t35 ^ 2 + t36 ^ 2) + m(2) * (t128 ^ 2 + t129 ^ 2) + m(3) * (t94 ^ 2 + t95 ^ 2) + m(4) * (t52 ^ 2 + t53 ^ 2) + m(5) * (t42 ^ 2 + t43 ^ 2) + t201; m(6) * (t147 * t35 + t149 * t36) + m(3) * (t147 * t94 + t149 * t95) + m(4) * (t147 * t52 + t149 * t53) + m(5) * (t147 * t42 + t149 * t43); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t163; (-t46 - t197) * t145 + m(6) * (t35 * t8 + t36 * t9) + m(4) * (t48 * t52 + t49 * t53) + m(5) * (t22 * t42 + t23 * t43) + ((t25 / 0.2e1 + t38 / 0.2e1) * t149 + (-t24 / 0.2e1 - t37 / 0.2e1) * t147) * t144 + t152; m(4) * (t147 * t48 + t149 * t49) + m(5) * (t147 * t22 + t149 * t23) + m(6) * (t147 * t8 + t149 * t9); (t46 * t145 + t191) * t145 + (-t147 * t4 - t147 * t3 + (-t147 * (-(t120 * t77 + t121 * t79) * t147 + (t120 * t78 + t121 * t80) * t149) + t149 * (-(t122 * t77 + t123 * t79) * t147 + (t122 * t78 + t123 * t80) * t149) + (-t147 * (t141 * t75 - t76 * t172) + t149 * (t142 * t76 - t75 * t172)) * t144) * t144 + ((-t25 - t38) * t149 + (t24 + t37) * t147) * t145) * t144 + m(6) * (t7 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t11 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(4) * (t39 ^ 2 + t48 ^ 2 + t49 ^ 2) + t192; -t204 + m(6) * (t12 * t35 + t13 * t36) + m(5) * (t42 * t44 + t43 * t45) + t152; m(5) * (t147 * t44 + t149 * t45) + m(6) * (t12 * t147 + t13 * t149); m(6) * (t10 * t7 + t12 * t8 + t13 * t9) + m(5) * (t11 * t34 + t22 * t44 + t23 * t45) + t151; m(5) * (t34 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t10 ^ 2 + t12 ^ 2 + t13 ^ 2) + t151; m(6) * (-t147 * t36 + t149 * t35) * t144; 0; m(6) * (-t145 * t7 + (-t147 * t9 + t149 * t8) * t144); m(6) * (-t10 * t145 + (t12 * t149 - t13 * t147) * t144); m(6) * (t163 * t144 ^ 2 + t145 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
