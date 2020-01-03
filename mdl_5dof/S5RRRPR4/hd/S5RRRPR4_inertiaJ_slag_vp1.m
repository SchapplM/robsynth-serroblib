% Calculate joint inertia matrix for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:00
% EndTime: 2019-12-31 21:11:04
% DurationCPUTime: 1.37s
% Computational Cost: add. (2680->194), mult. (3172->285), div. (0->0), fcn. (3281->8), ass. (0->101)
t125 = cos(qJ(3));
t199 = t125 ^ 2;
t194 = Icges(5,4) + Icges(4,5);
t192 = Icges(4,6) - Icges(5,6);
t120 = qJ(1) + qJ(2);
t116 = sin(t120);
t114 = t116 ^ 2;
t117 = cos(t120);
t115 = t117 ^ 2;
t154 = t114 + t115;
t193 = Icges(5,2) + Icges(4,3);
t122 = sin(qJ(3));
t191 = -t192 * t122 + t194 * t125;
t186 = t193 * t116 + t191 * t117;
t185 = -t191 * t116 + t193 * t117;
t121 = sin(qJ(5));
t124 = cos(qJ(5));
t86 = -t125 * t121 + t122 * t124;
t73 = t86 * t116;
t172 = Icges(6,2) * t73;
t134 = t122 * t121 + t125 * t124;
t74 = t134 * t116;
t174 = Icges(6,4) * t74;
t129 = Icges(6,6) * t117 + t172 + t174;
t130 = Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t117;
t75 = t86 * t117;
t76 = t134 * t117;
t40 = Icges(6,4) * t76 + Icges(6,2) * t75 - Icges(6,6) * t116;
t41 = Icges(6,1) * t76 + Icges(6,4) * t75 - Icges(6,5) * t116;
t53 = Icges(6,5) * t86 - Icges(6,6) * t134;
t54 = Icges(6,4) * t86 - Icges(6,2) * t134;
t55 = Icges(6,1) * t86 - Icges(6,4) * t134;
t152 = -(-t116 * t53 - t134 * t40 + t86 * t41 + t75 * t54 + t76 * t55) * t116 / 0.2e1 + (t117 * t53 - t129 * t134 + t86 * t130 + t73 * t54 + t74 * t55) * t117 / 0.2e1;
t184 = 0.2e1 * t122;
t183 = m(5) / 0.2e1;
t96 = t122 * rSges(4,1) + t125 * rSges(4,2);
t182 = m(4) * t96;
t179 = -rSges(6,3) - pkin(8);
t123 = sin(qJ(1));
t178 = t123 * pkin(1);
t160 = qJ(4) * t122;
t157 = t117 * t125;
t158 = t117 * t122;
t166 = pkin(3) * t157 + qJ(4) * t158;
t177 = t114 * (pkin(3) * t125 + t160) + t117 * t166;
t176 = t76 * rSges(6,1) + t75 * rSges(6,2);
t94 = t122 * pkin(3) - t125 * qJ(4);
t175 = -t122 * rSges(5,1) + t125 * rSges(5,3) - t94;
t173 = Icges(6,5) * t74;
t171 = Icges(6,6) * t73;
t161 = Icges(6,3) * t117;
t159 = t116 * t125;
t156 = t116 * t122 * rSges(4,2) + t117 * rSges(4,3);
t155 = t117 * pkin(2) + t116 * pkin(7);
t153 = rSges(5,1) * t157 + t116 * rSges(5,2) + rSges(5,3) * t158;
t81 = t117 * rSges(3,1) - t116 * rSges(3,2);
t151 = t155 + t166;
t56 = t86 * rSges(6,1) - rSges(6,2) * t134;
t150 = -pkin(4) * t122 - t56 - t94;
t149 = -t74 * rSges(6,1) - t73 * rSges(6,2);
t39 = Icges(6,5) * t76 + Icges(6,6) * t75 - Icges(6,3) * t116;
t148 = t117 * ((t117 * t39 + t73 * t40 + t74 * t41) * t116 - (Icges(6,1) * t74 ^ 2 + (t172 + 0.2e1 * t174) * t73 + (t161 + 0.2e1 * t171 + 0.2e1 * t173) * t117) * t117) - t116 * ((-t116 * t39 + t75 * t40 + t76 * t41) * t116 - (t76 * t130 + t75 * t129 - t116 * (t161 + t171 + t173)) * t117);
t80 = -t116 * rSges(3,1) - t117 * rSges(3,2);
t31 = t150 * t116;
t32 = t150 * t117;
t147 = t116 * t31 + t117 * t32;
t133 = rSges(4,1) * t157 - rSges(4,2) * t158 + t116 * rSges(4,3);
t112 = t117 * pkin(7);
t22 = t112 + t179 * t117 + (-t160 - pkin(2) + (-pkin(3) - pkin(4)) * t125) * t116 + t149;
t20 = t22 - t178;
t126 = cos(qJ(1));
t118 = t126 * pkin(1);
t104 = pkin(4) * t157;
t23 = t179 * t116 + t104 + t151 + t176;
t21 = t118 + t23;
t132 = m(6) * (t116 * t21 + t117 * t20);
t131 = m(6) * (t116 * t23 + t117 * t22);
t38 = t151 + t153;
t49 = t133 + t155;
t128 = -t134 * t54 + t86 * t55 + Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t199 + ((Icges(4,1) + Icges(5,1)) * t122 + (2 * Icges(4,4) - 2 * Icges(5,5)) * t125) * t122;
t48 = t112 + (-rSges(4,1) * t125 - pkin(2)) * t116 + t156;
t127 = -t152 + t154 * (t194 * t122 + t192 * t125);
t109 = t117 * rSges(5,2);
t37 = t109 + t112 + (-pkin(2) + (-rSges(5,1) - pkin(3)) * t125 + (-rSges(5,3) - qJ(4)) * t122) * t116;
t98 = t126 * rSges(2,1) - t123 * rSges(2,2);
t97 = -t123 * rSges(2,1) - t126 * rSges(2,2);
t79 = t118 + t81;
t78 = t80 - t178;
t51 = t175 * t117;
t50 = t175 * t116;
t47 = t118 + t49;
t46 = t48 - t178;
t43 = -t116 * rSges(6,3) + t176;
t42 = t117 * rSges(6,3) - t149;
t30 = t118 + t38;
t29 = t37 - t178;
t28 = t116 * (rSges(4,1) * t159 - t156) + t117 * t133;
t19 = t117 * t153 + (-t109 + (rSges(5,1) * t125 + rSges(5,3) * t122) * t116) * t116 + t177;
t18 = -t116 * t42 - t117 * t43;
t5 = (t104 + t43) * t117 + (pkin(4) * t159 + t42) * t116 + t177;
t1 = [Icges(2,3) + m(6) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2) + m(4) * (t46 ^ 2 + t47 ^ 2) + m(3) * (t78 ^ 2 + t79 ^ 2) + m(2) * (t97 ^ 2 + t98 ^ 2) + t128; m(6) * (t22 * t20 + t23 * t21) + m(5) * (t37 * t29 + t38 * t30) + m(4) * (t48 * t46 + t49 * t47) + m(3) * (t80 * t78 + t81 * t79) + t128; m(6) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t37 ^ 2 + t38 ^ 2) + m(4) * (t48 ^ 2 + t49 ^ 2) + m(3) * (t80 ^ 2 + t81 ^ 2) + t128; t127 + m(6) * (t32 * t20 + t31 * t21) + m(5) * (t51 * t29 + t50 * t30) + (-t116 * t47 - t117 * t46) * t182; t127 + m(6) * (t32 * t22 + t31 * t23) + m(5) * (t51 * t37 + t50 * t38) + (-t116 * t49 - t117 * t48) * t182; m(6) * (t31 ^ 2 + t32 ^ 2 + t5 ^ 2) + m(5) * (t19 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(4) * (t154 * t96 ^ 2 + t28 ^ 2) - t148 + t186 * t116 * t114 + (t185 * t115 + (t185 * t116 + t186 * t117) * t116) * t117; (t132 / 0.2e1 + (t116 * t30 + t117 * t29) * t183) * t184; (t131 / 0.2e1 + (t116 * t38 + t117 * t37) * t183) * t184; m(6) * (t147 * t122 - t125 * t5) + m(5) * (-t125 * t19 + (t116 * t50 + t117 * t51) * t122); 0.2e1 * (t183 + m(6) / 0.2e1) * (t154 * t122 ^ 2 + t199); t56 * t132 + t152; t56 * t131 + t152; m(6) * (t147 * t56 + t18 * t5) + t148; m(6) * (t154 * t56 * t122 - t18 * t125); m(6) * (t154 * t56 ^ 2 + t18 ^ 2) - t148;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
