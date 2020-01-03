% Calculate joint inertia matrix for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR13_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR13_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR13_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:33
% EndTime: 2019-12-31 19:14:39
% DurationCPUTime: 1.95s
% Computational Cost: add. (3774->307), mult. (6001->451), div. (0->0), fcn. (6476->8), ass. (0->156)
t148 = cos(qJ(3));
t203 = Icges(4,5) * t148;
t145 = sin(qJ(3));
t202 = Icges(4,6) * t145;
t201 = t203 / 0.2e1 - t202 / 0.2e1;
t149 = cos(qJ(1));
t200 = (rSges(4,1) * t145 + rSges(4,2) * t148) * t149;
t146 = sin(qJ(1));
t141 = t146 ^ 2;
t142 = t149 ^ 2;
t199 = -pkin(1) - pkin(6);
t173 = t148 * t149;
t143 = qJ(4) + qJ(5);
t134 = sin(t143);
t135 = cos(t143);
t178 = t145 * t146;
t102 = -t134 * t178 + t135 * t149;
t103 = t134 * t149 + t135 * t178;
t175 = t146 * t148;
t56 = Icges(6,5) * t103 + Icges(6,6) * t102 - Icges(6,3) * t175;
t58 = Icges(6,4) * t103 + Icges(6,2) * t102 - Icges(6,6) * t175;
t60 = Icges(6,1) * t103 + Icges(6,4) * t102 - Icges(6,5) * t175;
t24 = t145 * t56 + (-t134 * t58 + t135 * t60) * t148;
t177 = t145 * t149;
t104 = t134 * t177 + t135 * t146;
t105 = t134 * t146 - t135 * t177;
t57 = Icges(6,5) * t105 + Icges(6,6) * t104 + Icges(6,3) * t173;
t59 = Icges(6,4) * t105 + Icges(6,2) * t104 + Icges(6,6) * t173;
t61 = Icges(6,1) * t105 + Icges(6,4) * t104 + Icges(6,5) * t173;
t25 = t145 * t57 + (-t134 * t59 + t135 * t61) * t148;
t86 = Icges(6,6) * t145 + (Icges(6,4) * t135 - Icges(6,2) * t134) * t148;
t185 = t134 * t86;
t85 = Icges(6,3) * t145 + (Icges(6,5) * t135 - Icges(6,6) * t134) * t148;
t87 = Icges(6,5) * t145 + (Icges(6,1) * t135 - Icges(6,4) * t134) * t148;
t188 = t148 * t135 * t87 + t145 * t85;
t43 = (-t148 * t185 + t188) * t145;
t20 = t104 * t58 + t105 * t60 + t56 * t173;
t21 = t104 * t59 + t105 * t61 + t57 * t173;
t38 = t104 * t86 + t105 * t87 + t85 * t173;
t5 = t145 * t38 + (-t146 * t20 + t149 * t21) * t148;
t198 = t5 * t173 + t145 * (t43 + (-t146 * t24 + t149 * t25) * t148);
t196 = t145 / 0.2e1;
t195 = t146 / 0.2e1;
t193 = t149 / 0.2e1;
t123 = rSges(4,1) * t148 - rSges(4,2) * t145;
t192 = m(4) * t123;
t150 = -pkin(8) - pkin(7);
t191 = -pkin(7) - t150;
t62 = t103 * rSges(6,1) + t102 * rSges(6,2) - rSges(6,3) * t175;
t130 = pkin(3) * t178;
t115 = -pkin(7) * t175 + t130;
t147 = cos(qJ(4));
t133 = pkin(4) * t147 + pkin(3);
t144 = sin(qJ(4));
t179 = t144 * t149;
t167 = pkin(4) * t179 + t133 * t178 + t150 * t175;
t73 = -t115 + t167;
t190 = -t62 - t73;
t132 = pkin(3) * t177;
t180 = t144 * t146;
t181 = t133 * t145;
t156 = -t105 * rSges(6,1) - t104 * rSges(6,2);
t63 = rSges(6,3) * t173 - t156;
t189 = t63 + pkin(4) * t180 + t132 + (t191 * t148 - t181) * t149;
t88 = rSges(6,3) * t145 + (rSges(6,1) * t135 - rSges(6,2) * t134) * t148;
t46 = t145 * t62 + t88 * t175;
t92 = Icges(5,3) * t145 + (Icges(5,5) * t147 - Icges(5,6) * t144) * t148;
t98 = Icges(5,5) * t145 + (Icges(5,1) * t147 - Icges(5,4) * t144) * t148;
t187 = t148 * t147 * t98 + t145 * t92;
t84 = (-pkin(3) + t133) * t148 + t191 * t145;
t186 = t84 + t88;
t95 = Icges(5,6) * t145 + (Icges(5,4) * t147 - Icges(5,2) * t144) * t148;
t184 = t144 * t95;
t176 = t146 * t147;
t174 = t147 * t149;
t111 = -t144 * t178 + t174;
t112 = t145 * t176 + t179;
t172 = t112 * rSges(5,1) + t111 * rSges(5,2);
t171 = t149 * pkin(1) + t146 * qJ(2);
t170 = t141 + t142;
t65 = Icges(5,5) * t112 + Icges(5,6) * t111 - Icges(5,3) * t175;
t67 = Icges(5,4) * t112 + Icges(5,2) * t111 - Icges(5,6) * t175;
t69 = Icges(5,1) * t112 + Icges(5,4) * t111 - Icges(5,5) * t175;
t32 = t145 * t65 + (-t144 * t67 + t147 * t69) * t148;
t40 = t111 * t95 + t112 * t98 - t92 * t175;
t169 = -t32 / 0.2e1 - t40 / 0.2e1;
t113 = t144 * t177 + t176;
t114 = -t145 * t174 + t180;
t66 = Icges(5,5) * t114 + Icges(5,6) * t113 + Icges(5,3) * t173;
t68 = Icges(5,4) * t114 + Icges(5,2) * t113 + Icges(5,6) * t173;
t70 = Icges(5,1) * t114 + Icges(5,4) * t113 + Icges(5,5) * t173;
t33 = t145 * t66 + (-t144 * t68 + t147 * t70) * t148;
t41 = t113 * t95 + t114 * t98 + t92 * t173;
t168 = t41 / 0.2e1 + t33 / 0.2e1;
t166 = rSges(4,1) * t178 + rSges(4,2) * t175 + t149 * rSges(4,3);
t165 = t149 * pkin(6) + t171;
t164 = (-rSges(5,3) - pkin(7)) * t148;
t163 = -t175 / 0.2e1;
t162 = t173 / 0.2e1;
t18 = t102 * t58 + t103 * t60 - t56 * t175;
t19 = t102 * t59 + t103 * t61 - t57 * t175;
t11 = t146 * t19 + t149 * t18;
t12 = t146 * t21 + t149 * t20;
t37 = t102 * t86 + t103 * t87 - t85 * t175;
t4 = t145 * t37 + (-t146 * t18 + t149 * t19) * t148;
t161 = t11 * t163 + t12 * t162 + t4 * t193 + t5 * t195 + (t25 * t146 + t24 * t149) * t196;
t160 = t43 + (t24 + t37) * t163 + (t25 + t38) * t162;
t159 = -t4 * t175 + t198;
t157 = -rSges(5,1) * t114 - rSges(5,2) * t113;
t151 = Icges(4,5) * t145 + Icges(4,6) * t148;
t137 = t149 * qJ(2);
t125 = t148 * pkin(3) + t145 * pkin(7);
t124 = rSges(2,1) * t149 - rSges(2,2) * t146;
t122 = -rSges(2,1) * t146 - rSges(2,2) * t149;
t119 = -t202 + t203;
t116 = t146 * t125;
t108 = t149 * (pkin(7) * t173 - t132);
t107 = -rSges(3,2) * t149 + rSges(3,3) * t146 + t171;
t106 = rSges(3,3) * t149 + t137 + (rSges(3,2) - pkin(1)) * t146;
t101 = rSges(5,3) * t145 + (rSges(5,1) * t147 - rSges(5,2) * t144) * t148;
t94 = Icges(4,3) * t146 - t151 * t149;
t93 = Icges(4,3) * t149 + t151 * t146;
t81 = t88 * t173;
t78 = t165 + t166;
t77 = t137 + t200 + (-rSges(4,3) + t199) * t146;
t76 = (-t101 - t125) * t149;
t75 = t101 * t146 + t116;
t72 = rSges(5,3) * t173 - t157;
t71 = -rSges(5,3) * t175 + t172;
t64 = -t146 * t166 + (t146 * rSges(4,3) - t200) * t149;
t54 = (-t125 - t186) * t149;
t53 = t186 * t146 + t116;
t52 = t146 * t164 + t130 + t165 + t172;
t51 = t199 * t146 + t149 * t164 + t132 + t137 + t157;
t50 = t101 * t173 - t145 * t72;
t49 = t101 * t175 + t145 * t71;
t48 = (-t148 * t184 + t187) * t145;
t47 = -t145 * t63 + t81;
t45 = t62 + t165 + t167;
t44 = t137 + (t181 + (-rSges(6,3) + t150) * t148) * t149 + (-pkin(4) * t144 + t199) * t146 + t156;
t42 = (-t146 * t72 - t149 * t71) * t148;
t39 = (-t146 * t63 - t149 * t62) * t148;
t34 = t149 * t72 + t108 + (-t115 - t71) * t146;
t31 = -t189 * t145 + t84 * t173 + t81;
t30 = t145 * t73 + t84 * t175 + t46;
t29 = t113 * t68 + t114 * t70 + t66 * t173;
t28 = t113 * t67 + t114 * t69 + t65 * t173;
t27 = t111 * t68 + t112 * t70 - t66 * t175;
t26 = t111 * t67 + t112 * t69 - t65 * t175;
t17 = (-t189 * t146 + t190 * t149) * t148;
t16 = t108 + t189 * t149 + (-t115 + t190) * t146;
t15 = t146 * t29 + t149 * t28;
t14 = t146 * t27 + t149 * t26;
t8 = t145 * t41 + (-t146 * t28 + t149 * t29) * t148;
t7 = t145 * t40 + (-t146 * t26 + t149 * t27) * t148;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t148 - t184 - t185) * t148 + m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t51 ^ 2 + t52 ^ 2) + m(4) * (t77 ^ 2 + t78 ^ 2) + m(3) * (t106 ^ 2 + t107 ^ 2) + m(2) * (t122 ^ 2 + t124 ^ 2) + t187 + t188 + (-0.2e1 * Icges(4,4) * t148 + Icges(4,2) * t145) * t145; m(6) * (t146 * t44 - t149 * t45) + m(5) * (t146 * t51 - t149 * t52) + m(4) * (t146 * t77 - t149 * t78) + m(3) * (t106 * t146 - t107 * t149); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t170; m(6) * (t44 * t53 + t45 * t54) + m(5) * (t51 * t75 + t52 * t76) + (t37 / 0.2e1 + t24 / 0.2e1 - t78 * t192 + t119 * t193 - t169 + t201 * t149) * t149 + (t38 / 0.2e1 + t25 / 0.2e1 + t77 * t192 + t119 * t195 + t168 + t201 * t146) * t146; m(5) * (t146 * t75 - t149 * t76) + m(6) * (t146 * t53 - t149 * t54) + t170 * t192; m(6) * (t16 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t34 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(4) * (t170 * t123 ^ 2 + t64 ^ 2) + (t141 * t94 + t12 + t15) * t146 + (t142 * t93 + t11 + t14 + (t146 * t93 + t149 * t94) * t146) * t149; t48 + m(6) * (t30 * t45 + t31 * t44) + m(5) * (t49 * t52 + t50 * t51) + (t169 * t146 + t168 * t149) * t148 + t160; m(5) * (t146 * t50 - t149 * t49) + m(6) * (t146 * t31 - t149 * t30); t7 * t193 + t8 * t195 + (t33 * t146 + t32 * t149) * t196 + (t15 * t193 - t146 * t14 / 0.2e1) * t148 + m(6) * (t16 * t17 + t30 * t54 + t31 * t53) + m(5) * (t34 * t42 + t49 * t76 + t50 * t75) + t161; t145 * t48 + m(6) * (t17 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(5) * (t42 ^ 2 + t49 ^ 2 + t50 ^ 2) + ((t145 * t33 + t8) * t149 + (-t145 * t32 - t4 - t7) * t146) * t148 + t198; m(6) * (t44 * t47 + t45 * t46) + t160; m(6) * (t146 * t47 - t149 * t46); m(6) * (t16 * t39 + t46 * t54 + t47 * t53) + t161; m(6) * (t17 * t39 + t30 * t46 + t31 * t47) + t159; m(6) * (t39 ^ 2 + t46 ^ 2 + t47 ^ 2) + t159;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
