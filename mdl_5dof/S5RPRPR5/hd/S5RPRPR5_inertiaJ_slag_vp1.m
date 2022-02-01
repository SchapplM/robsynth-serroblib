% Calculate joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:54
% EndTime: 2022-01-23 09:24:56
% DurationCPUTime: 1.69s
% Computational Cost: add. (4066->326), mult. (4668->461), div. (0->0), fcn. (4929->10), ass. (0->158)
t139 = qJ(3) + pkin(9);
t130 = sin(t139);
t131 = cos(t139);
t142 = sin(pkin(8));
t143 = cos(pkin(8));
t84 = -Icges(5,3) * t143 + (Icges(5,5) * t131 - Icges(5,6) * t130) * t142;
t145 = sin(qJ(3));
t147 = cos(qJ(3));
t98 = -Icges(4,3) * t143 + (Icges(4,5) * t147 - Icges(4,6) * t145) * t142;
t196 = -t84 - t98;
t85 = -Icges(5,6) * t143 + (Icges(5,4) * t131 - Icges(5,2) * t130) * t142;
t99 = -Icges(4,6) * t143 + (Icges(4,4) * t147 - Icges(4,2) * t145) * t142;
t195 = -t130 * t85 - t145 * t99;
t100 = -Icges(4,5) * t143 + (Icges(4,1) * t147 - Icges(4,4) * t145) * t142;
t175 = t142 * t147;
t86 = -Icges(5,5) * t143 + (Icges(5,1) * t131 - Icges(5,4) * t130) * t142;
t194 = t142 * t131 * t86 + t100 * t175;
t189 = pkin(3) * t147;
t160 = -pkin(2) - t189;
t117 = pkin(4) * t131 - t160;
t190 = pkin(3) * t145;
t118 = pkin(4) * t130 + t190;
t144 = qJ(4) + pkin(6);
t137 = -pkin(7) - t144;
t146 = sin(qJ(1));
t133 = t146 * qJ(2);
t148 = cos(qJ(1));
t164 = t148 * pkin(1) + t133;
t173 = t143 * t148;
t174 = t142 * t148;
t132 = qJ(5) + t139;
t127 = sin(t132);
t128 = cos(t132);
t170 = t146 * t128;
t94 = -t127 * t173 + t170;
t171 = t146 * t127;
t95 = t128 * t173 + t171;
t53 = t95 * rSges(6,1) + t94 * rSges(6,2) + rSges(6,3) * t174;
t26 = t117 * t173 + t146 * t118 - t137 * t174 + t164 + t53;
t141 = t148 ^ 2;
t193 = m(5) / 0.2e1;
t192 = m(6) / 0.2e1;
t176 = t142 * t146;
t79 = -Icges(6,3) * t143 + (Icges(6,5) * t128 - Icges(6,6) * t127) * t142;
t80 = -Icges(6,6) * t143 + (Icges(6,4) * t128 - Icges(6,2) * t127) * t142;
t81 = -Icges(6,5) * t143 + (Icges(6,1) * t128 - Icges(6,4) * t127) * t142;
t92 = -t128 * t148 - t143 * t171;
t93 = -t127 * t148 + t143 * t170;
t20 = t79 * t176 + t80 * t92 + t81 * t93;
t21 = t79 * t174 + t94 * t80 + t95 * t81;
t46 = Icges(6,5) * t93 + Icges(6,6) * t92 + Icges(6,3) * t176;
t47 = Icges(6,5) * t95 + Icges(6,6) * t94 + Icges(6,3) * t174;
t48 = Icges(6,4) * t93 + Icges(6,2) * t92 + Icges(6,6) * t176;
t49 = Icges(6,4) * t95 + Icges(6,2) * t94 + Icges(6,6) * t174;
t50 = Icges(6,1) * t93 + Icges(6,4) * t92 + Icges(6,5) * t176;
t51 = Icges(6,1) * t95 + Icges(6,4) * t94 + Icges(6,5) * t174;
t191 = ((t47 * t176 + t49 * t92 + t51 * t93) * t174 + (t46 * t176 + t48 * t92 + t50 * t93) * t176 - t20 * t143) * t176 + ((t47 * t174 + t94 * t49 + t95 * t51) * t174 + (t46 * t174 + t94 * t48 + t95 * t50) * t176 - t21 * t143) * t174;
t188 = pkin(2) * t143 + pkin(1);
t187 = t195 * t142 + t196 * t143 + t194;
t107 = pkin(3) * t175 + (pkin(6) - t144) * t143;
t106 = t142 * t144 + t143 * t189 + t188;
t119 = pkin(6) * t142 + t188;
t134 = t148 * qJ(2);
t129 = qJ(2) + t190;
t177 = t129 * t148;
t61 = -t177 + t134 + (t106 - t119) * t146;
t186 = t107 * t176 + t143 * t61;
t155 = -t119 * t148 - t133;
t178 = t106 * t148 + t129 * t146;
t62 = t155 + t178;
t168 = t146 * t131;
t104 = -t130 * t173 + t168;
t169 = t146 * t130;
t105 = t131 * t173 + t169;
t64 = t105 * rSges(5,1) + t104 * rSges(5,2) + rSges(5,3) * t174;
t185 = -t62 - t64;
t102 = -t131 * t148 - t143 * t169;
t103 = -t130 * t148 + t143 * t168;
t55 = Icges(5,5) * t103 + Icges(5,6) * t102 + Icges(5,3) * t176;
t165 = t147 * t148;
t167 = t146 * t145;
t112 = -t143 * t167 - t165;
t166 = t146 * t147;
t172 = t145 * t148;
t113 = t143 * t166 - t172;
t65 = Icges(4,5) * t113 + Icges(4,6) * t112 + Icges(4,3) * t176;
t184 = t65 + t55;
t56 = Icges(5,5) * t105 + Icges(5,6) * t104 + Icges(5,3) * t174;
t114 = -t143 * t172 + t166;
t115 = t143 * t165 + t167;
t66 = Icges(4,5) * t115 + Icges(4,6) * t114 + Icges(4,3) * t174;
t183 = t66 + t56;
t153 = -t93 * rSges(6,1) - t92 * rSges(6,2);
t52 = rSges(6,3) * t176 - t153;
t82 = -rSges(6,3) * t143 + (rSges(6,1) * t128 - rSges(6,2) * t127) * t142;
t32 = t143 * t52 + t82 * t176;
t182 = t127 * t80;
t74 = t142 * t128 * t81;
t27 = -t142 * t182 - t143 * t79 + t74;
t179 = t27 * t143;
t163 = t146 ^ 2 + t141;
t162 = t193 + t192;
t161 = t178 - t62 - t26;
t72 = t115 * rSges(4,1) + t114 * rSges(4,2) + rSges(4,3) * t174;
t10 = -t143 * t46 + (-t127 * t48 + t128 * t50) * t142;
t11 = -t143 * t47 + (-t127 * t49 + t128 * t51) * t142;
t157 = (t10 + t20) * t176 / 0.2e1 + (t11 + t21) * t174 / 0.2e1;
t156 = -t117 * t143 - pkin(1);
t3 = -t179 + (t10 * t146 + t11 * t148) * t142;
t154 = -t143 * t3 + t191;
t152 = rSges(3,1) * t143 - rSges(3,2) * t142;
t151 = -rSges(4,1) * t113 - rSges(4,2) * t112;
t150 = -t103 * rSges(5,1) - t102 * rSges(5,2);
t121 = rSges(2,1) * t148 - t146 * rSges(2,2);
t120 = -t146 * rSges(2,1) - rSges(2,2) * t148;
t101 = -rSges(4,3) * t143 + (rSges(4,1) * t147 - rSges(4,2) * t145) * t142;
t88 = -rSges(5,3) * t143 + (rSges(5,1) * t131 - rSges(5,2) * t130) * t142;
t78 = t146 * rSges(3,3) + t152 * t148 + t164;
t77 = rSges(3,3) * t148 + t134 + (-pkin(1) - t152) * t146;
t73 = (t137 + t144) * t143 + (t117 + t160) * t142;
t71 = rSges(4,3) * t176 - t151;
t70 = Icges(4,1) * t115 + Icges(4,4) * t114 + Icges(4,5) * t174;
t69 = Icges(4,1) * t113 + Icges(4,4) * t112 + Icges(4,5) * t176;
t68 = Icges(4,4) * t115 + Icges(4,2) * t114 + Icges(4,6) * t174;
t67 = Icges(4,4) * t113 + Icges(4,2) * t112 + Icges(4,6) * t176;
t63 = rSges(5,3) * t176 - t150;
t60 = Icges(5,1) * t105 + Icges(5,4) * t104 + Icges(5,5) * t174;
t59 = Icges(5,1) * t103 + Icges(5,4) * t102 + Icges(5,5) * t176;
t58 = Icges(5,4) * t105 + Icges(5,2) * t104 + Icges(5,6) * t174;
t57 = Icges(5,4) * t103 + Icges(5,2) * t102 + Icges(5,6) * t176;
t45 = t61 * t174;
t43 = t134 + (-rSges(4,3) * t142 - t119) * t146 + t151;
t42 = -t155 + t72;
t41 = t52 * t174;
t40 = -t101 * t174 - t143 * t72;
t39 = t101 * t176 + t143 * t71;
t37 = t64 + t178;
t36 = t177 + (-rSges(5,3) * t142 - t106) * t146 + t150;
t34 = -t134 + (-t118 + t129) * t148 + (-t137 * t142 - t106 - t156) * t146;
t33 = -t143 * t53 - t82 * t174;
t30 = (-t146 * t72 + t148 * t71) * t142;
t29 = t115 * t100 + t114 * t99 + t98 * t174;
t28 = t100 * t113 + t112 * t99 + t98 * t176;
t25 = t118 * t148 + t134 + ((-rSges(6,3) + t137) * t142 + t156) * t146 + t153;
t24 = t104 * t85 + t105 * t86 + t84 * t174;
t23 = t102 * t85 + t103 * t86 + t84 * t176;
t22 = -t53 * t176 + t41;
t17 = -t143 * t66 + (-t145 * t68 + t147 * t70) * t142;
t16 = -t143 * t65 + (-t145 * t67 + t147 * t69) * t142;
t15 = t185 * t143 + (-t107 - t88) * t174;
t14 = t143 * t63 + t88 * t176 + t186;
t13 = -t143 * t56 + (-t130 * t58 + t131 * t60) * t142;
t12 = -t143 * t55 + (-t130 * t57 + t131 * t59) * t142;
t7 = t45 + (t185 * t146 + t148 * t63) * t142;
t6 = t161 * t143 + (-t107 - t73 - t82) * t174;
t5 = t143 * t34 + t73 * t176 + t186 + t32;
t4 = t41 + t45 + (t161 * t146 + t148 * t34) * t142;
t1 = [Icges(2,3) + t74 + (Icges(3,2) * t143 + t196 - t79) * t143 + (Icges(3,1) * t142 + 0.2e1 * Icges(3,4) * t143 - t182 + t195) * t142 + m(6) * (t25 ^ 2 + t26 ^ 2) + m(5) * (t36 ^ 2 + t37 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2) + m(3) * (t77 ^ 2 + t78 ^ 2) + m(2) * (t120 ^ 2 + t121 ^ 2) + t194; m(6) * (t146 * t25 - t148 * t26) + m(5) * (t146 * t36 - t148 * t37) + m(4) * (t146 * t43 - t148 * t42) + m(3) * (t146 * t77 - t148 * t78); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + t162) * t163; (-t27 - t187) * t143 + m(6) * (t25 * t5 + t26 * t6) + m(5) * (t14 * t36 + t15 * t37) + m(4) * (t39 * t43 + t40 * t42) + ((t29 / 0.2e1 + t24 / 0.2e1 + t17 / 0.2e1 + t13 / 0.2e1) * t148 + (t16 / 0.2e1 + t12 / 0.2e1 + t28 / 0.2e1 + t23 / 0.2e1) * t146) * t142 + t157; m(4) * (t39 * t146 - t148 * t40) + m(5) * (t14 * t146 - t148 * t15) + m(6) * (t5 * t146 - t148 * t6); m(6) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t14 ^ 2 + t15 ^ 2 + t7 ^ 2) + m(4) * (t30 ^ 2 + t39 ^ 2 + t40 ^ 2) + (t187 * t143 - t3) * t143 + (((t104 * t58 + t105 * t60 + t114 * t68 + t115 * t70 + t183 * t174) * t141 + ((t102 * t57 + t103 * t59 + t112 * t67 + t113 * t69 + t184 * t176) * t146 + (t102 * t58 + t103 * t60 + t104 * t57 + t105 * t59 + t112 * t68 + t113 * t70 + t114 * t67 + t115 * t69 + (t183 * t146 + t184 * t148) * t142) * t148) * t146) * t142 + ((-t13 - t17 - t24 - t29) * t148 + (-t12 - t16 - t23 - t28) * t146) * t143) * t142 + t191; 0.2e1 * ((t146 * t26 + t148 * t25) * t192 + (t146 * t37 + t148 * t36) * t193) * t142; 0; m(6) * (-t143 * t4 + (t146 * t6 + t148 * t5) * t142) + m(5) * (-t143 * t7 + (t14 * t148 + t146 * t15) * t142); 0.2e1 * t162 * (t163 * t142 ^ 2 + t143 ^ 2); m(6) * (t25 * t32 + t26 * t33) - t179 + t157; m(6) * (t32 * t146 - t148 * t33); m(6) * (t22 * t4 + t32 * t5 + t33 * t6) + t154; m(6) * (-t22 * t143 + (t146 * t33 + t148 * t32) * t142); m(6) * (t22 ^ 2 + t32 ^ 2 + t33 ^ 2) + t154;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
