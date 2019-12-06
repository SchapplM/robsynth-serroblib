% Calculate joint inertia matrix for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:34
% EndTime: 2019-12-05 17:09:40
% DurationCPUTime: 1.73s
% Computational Cost: add. (7752->262), mult. (8136->410), div. (0->0), fcn. (8918->10), ass. (0->147)
t141 = sin(pkin(9));
t137 = t141 ^ 2;
t142 = cos(pkin(9));
t138 = t142 ^ 2;
t205 = t137 + t138;
t140 = qJ(2) + qJ(3);
t134 = sin(t140);
t136 = cos(t140);
t145 = cos(qJ(4));
t195 = t145 * pkin(4);
t150 = pkin(8) * t134 + t195 * t136;
t143 = sin(qJ(4));
t183 = t141 * t143;
t139 = qJ(4) + qJ(5);
t133 = sin(t139);
t181 = t142 * t133;
t135 = cos(t139);
t184 = t141 * t135;
t114 = -t136 * t181 + t184;
t180 = t142 * t135;
t185 = t141 * t133;
t115 = t136 * t180 + t185;
t186 = t134 * t142;
t67 = t115 * rSges(6,1) + t114 * rSges(6,2) + rSges(6,3) * t186;
t204 = pkin(4) * t183 + t150 * t142 + t67;
t156 = Icges(4,5) * t136 - Icges(4,6) * t134;
t104 = -Icges(4,3) * t142 + t156 * t141;
t105 = Icges(4,3) * t141 + t156 * t142;
t112 = -t136 * t185 - t180;
t113 = t136 * t184 - t181;
t187 = t134 * t141;
t60 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t187;
t62 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t187;
t64 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t187;
t27 = t112 * t62 + t113 * t64 + t60 * t187;
t61 = Icges(6,5) * t115 + Icges(6,6) * t114 + Icges(6,3) * t186;
t63 = Icges(6,4) * t115 + Icges(6,2) * t114 + Icges(6,6) * t186;
t65 = Icges(6,1) * t115 + Icges(6,4) * t114 + Icges(6,5) * t186;
t28 = t112 * t63 + t113 * t65 + t61 * t187;
t15 = t28 * t141 - t27 * t142;
t158 = Icges(4,4) * t136 - Icges(4,2) * t134;
t160 = Icges(4,1) * t136 - Icges(4,4) * t134;
t154 = -(Icges(4,6) * t141 + t158 * t142) * t134 + (Icges(4,5) * t141 + t160 * t142) * t136;
t155 = (-Icges(4,6) * t142 + t158 * t141) * t134 - (-Icges(4,5) * t142 + t160 * t141) * t136;
t178 = t142 * t145;
t122 = -t136 * t183 - t178;
t179 = t142 * t143;
t182 = t141 * t145;
t123 = t136 * t182 - t179;
t77 = Icges(5,5) * t123 + Icges(5,6) * t122 + Icges(5,3) * t187;
t79 = Icges(5,4) * t123 + Icges(5,2) * t122 + Icges(5,6) * t187;
t81 = Icges(5,1) * t123 + Icges(5,4) * t122 + Icges(5,5) * t187;
t38 = t122 * t79 + t123 * t81 + t77 * t187;
t124 = -t136 * t179 + t182;
t125 = t136 * t178 + t183;
t78 = Icges(5,5) * t125 + Icges(5,6) * t124 + Icges(5,3) * t186;
t80 = Icges(5,4) * t125 + Icges(5,2) * t124 + Icges(5,6) * t186;
t82 = Icges(5,1) * t125 + Icges(5,4) * t124 + Icges(5,5) * t186;
t39 = t122 * t80 + t123 * t82 + t78 * t187;
t21 = t39 * t141 - t38 * t142;
t202 = -t15 - t21 - t138 * t104 - (t154 * t141 + (-t105 + t155) * t142) * t141;
t201 = t136 ^ 2;
t200 = -t136 / 0.2e1;
t199 = t141 / 0.2e1;
t198 = -t142 / 0.2e1;
t144 = sin(qJ(2));
t197 = pkin(2) * t144;
t66 = t113 * rSges(6,1) + t112 * rSges(6,2) + rSges(6,3) * t187;
t95 = -t136 * rSges(6,3) + (rSges(6,1) * t135 - rSges(6,2) * t133) * t134;
t48 = t136 * t66 + t95 * t187;
t89 = -pkin(8) * t136 + t195 * t134;
t191 = -t89 - t95;
t146 = cos(qJ(2));
t190 = t205 * t146 * pkin(2);
t68 = t205 * (rSges(4,1) * t136 - rSges(4,2) * t134);
t90 = -Icges(6,3) * t136 + (Icges(6,5) * t135 - Icges(6,6) * t133) * t134;
t189 = t136 * t90;
t98 = -Icges(5,3) * t136 + (Icges(5,5) * t145 - Icges(5,6) * t143) * t134;
t188 = t136 * t98;
t101 = -t136 * rSges(5,3) + (rSges(5,1) * t145 - rSges(5,2) * t143) * t134;
t128 = t134 * pkin(3) - t136 * pkin(7);
t177 = -t101 - t128;
t176 = t205 * (pkin(3) * t136 + pkin(7) * t134);
t29 = t114 * t62 + t115 * t64 + t60 * t186;
t30 = t114 * t63 + t115 * t65 + t61 * t186;
t16 = t30 * t141 - t29 * t142;
t40 = t124 * t79 + t125 * t81 + t77 * t186;
t41 = t124 * t80 + t125 * t82 + t78 * t186;
t22 = t41 * t141 - t40 * t142;
t174 = (t137 * t105 + t16 + t22 + (t155 * t142 + (-t104 + t154) * t141) * t142) * t141;
t173 = -t128 + t191;
t172 = t187 / 0.2e1;
t171 = t186 / 0.2e1;
t127 = t134 * rSges(4,1) + t136 * rSges(4,2);
t170 = -t127 - t197;
t169 = -t128 - t197;
t83 = t123 * rSges(5,1) + t122 * rSges(5,2) + rSges(5,3) * t187;
t84 = t125 * rSges(5,1) + t124 * rSges(5,2) + rSges(5,3) * t186;
t44 = t141 * t83 + t142 * t84 + t176;
t33 = -t136 * t60 + (-t133 * t62 + t135 * t64) * t134;
t34 = -t136 * t61 + (-t133 * t63 + t135 * t65) * t134;
t91 = -Icges(6,6) * t136 + (Icges(6,4) * t135 - Icges(6,2) * t133) * t134;
t92 = -Icges(6,5) * t136 + (Icges(6,1) * t135 - Icges(6,4) * t133) * t134;
t5 = -(t112 * t91 + t113 * t92) * t136 + (t28 * t142 + (t27 - t189) * t141) * t134;
t6 = -(t114 * t91 + t115 * t92) * t136 + (t29 * t141 + (t30 - t189) * t142) * t134;
t168 = -t136 * (t201 * t90 + (t34 * t142 + t33 * t141 - (-t133 * t91 + t135 * t92) * t136) * t134) + t6 * t186 + t5 * t187;
t167 = t15 * t172 + t16 * t171 + t5 * t198 + t6 * t199 + (t34 * t141 - t33 * t142) * t200;
t166 = -t101 + t169;
t162 = t169 + t191;
t157 = Icges(3,5) * t146 - Icges(3,6) * t144;
t75 = -pkin(4) * t179 + t150 * t141;
t25 = t176 + t204 * t142 + (t66 + t75) * t141;
t151 = t202 * t142 + t174;
t100 = -Icges(5,5) * t136 + (Icges(5,1) * t145 - Icges(5,4) * t143) * t134;
t99 = -Icges(5,6) * t136 + (Icges(5,4) * t145 - Icges(5,2) * t143) * t134;
t10 = -(t125 * t100 + t124 * t99) * t136 + (t40 * t141 + (t41 - t188) * t142) * t134;
t42 = -t136 * t77 + (-t143 * t79 + t145 * t81) * t134;
t43 = -t136 * t78 + (-t143 * t80 + t145 * t82) * t134;
t9 = -(t123 * t100 + t122 * t99) * t136 + (t39 * t142 + (t38 - t188) * t141) * t134;
t149 = t10 * t199 + t22 * t171 + t21 * t172 + t9 * t198 + t167 + (t43 * t141 - t42 * t142) * t200;
t130 = t144 * rSges(3,1) + t146 * rSges(3,2);
t117 = Icges(3,3) * t141 + t157 * t142;
t116 = -Icges(3,3) * t142 + t157 * t141;
t103 = t170 * t142;
t102 = t170 * t141;
t87 = t177 * t142;
t86 = t177 * t141;
t85 = t205 * (rSges(3,1) * t146 - rSges(3,2) * t144);
t72 = t166 * t142;
t71 = t166 * t141;
t56 = t66 * t186;
t55 = t173 * t142;
t54 = t173 * t141;
t53 = -t101 * t186 - t136 * t84;
t52 = t101 * t187 + t136 * t83;
t51 = t162 * t142;
t50 = t162 * t141;
t49 = -t136 * t67 - t95 * t186;
t47 = t68 + t190;
t46 = (-t141 * t84 + t142 * t83) * t134;
t45 = -t67 * t187 + t56;
t37 = -t136 * t204 + t191 * t186;
t36 = t136 * t75 + t89 * t187 + t48;
t31 = t44 + t190;
t26 = t56 + (-t141 * t204 + t142 * t75) * t134;
t24 = t25 + t190;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t85 + m(4) * t47 + m(5) * t31 + m(6) * t24; m(6) * (t24 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t31 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t102 ^ 2 + t103 ^ 2 + t47 ^ 2) + m(3) * (t130 ^ 2 * t205 + t85 ^ 2) + t141 * t137 * t117 + t174 + (-t138 * t116 + (-t141 * t116 + t142 * t117) * t141 + t202) * t142; m(4) * t68 + m(5) * t44 + m(6) * t25; m(6) * (t25 * t24 + t54 * t50 + t55 * t51) + m(5) * (t44 * t31 + t86 * t71 + t87 * t72) + m(4) * (t68 * t47 + (-t102 * t141 - t103 * t142) * t127) + t151; m(6) * (t25 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t44 ^ 2 + t86 ^ 2 + t87 ^ 2) + m(4) * (t127 ^ 2 * t205 + t68 ^ 2) + t151; m(5) * t46 + m(6) * t26; m(6) * (t26 * t24 + t36 * t51 + t37 * t50) + m(5) * (t46 * t31 + t52 * t72 + t53 * t71) + t149; m(6) * (t26 * t25 + t36 * t55 + t37 * t54) + m(5) * (t46 * t44 + t52 * t87 + t53 * t86) + t149; t9 * t187 + t10 * t186 + m(6) * (t26 ^ 2 + t36 ^ 2 + t37 ^ 2) - t136 * (t201 * t98 + (t43 * t142 + t42 * t141 - (t100 * t145 - t143 * t99) * t136) * t134) + m(5) * (t46 ^ 2 + t52 ^ 2 + t53 ^ 2) + t168; m(6) * t45; m(6) * (t45 * t24 + t48 * t51 + t49 * t50) + t167; m(6) * (t45 * t25 + t48 * t55 + t49 * t54) + t167; m(6) * (t45 * t26 + t48 * t36 + t49 * t37) + t168; m(6) * (t45 ^ 2 + t48 ^ 2 + t49 ^ 2) + t168;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
