% Calculate joint inertia matrix for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP6_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP6_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:10
% EndTime: 2019-12-31 17:18:14
% DurationCPUTime: 1.40s
% Computational Cost: add. (1782->250), mult. (4120->370), div. (0->0), fcn. (4437->6), ass. (0->128)
t114 = -qJ(4) - pkin(6);
t176 = rSges(5,3) - t114;
t116 = sin(qJ(2));
t175 = Icges(3,5) * t116;
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t119 = cos(qJ(2));
t71 = -Icges(5,6) * t119 + (Icges(5,4) * t118 - Icges(5,2) * t115) * t116;
t72 = -Icges(4,6) * t119 + (Icges(4,4) * t118 - Icges(4,2) * t115) * t116;
t174 = (t71 + t72) * t115;
t173 = t175 / 0.2e1;
t67 = -Icges(5,3) * t119 + (Icges(5,5) * t118 - Icges(5,6) * t115) * t116;
t68 = -Icges(4,3) * t119 + (Icges(4,5) * t118 - Icges(4,6) * t115) * t116;
t172 = t67 + t68;
t117 = sin(qJ(1));
t145 = t116 * t117;
t120 = cos(qJ(1));
t138 = t120 * t118;
t141 = t117 * t119;
t87 = -t115 * t141 - t138;
t139 = t120 * t115;
t88 = t118 * t141 - t139;
t42 = Icges(5,5) * t88 + Icges(5,6) * t87 + Icges(5,3) * t145;
t46 = Icges(5,4) * t88 + Icges(5,2) * t87 + Icges(5,6) * t145;
t50 = Icges(5,1) * t88 + Icges(5,4) * t87 + Icges(5,5) * t145;
t11 = t42 * t145 + t87 * t46 + t88 * t50;
t143 = t116 * t120;
t89 = t117 * t118 - t119 * t139;
t142 = t117 * t115;
t90 = t119 * t138 + t142;
t43 = Icges(5,5) * t90 + Icges(5,6) * t89 + Icges(5,3) * t143;
t47 = Icges(5,4) * t90 + Icges(5,2) * t89 + Icges(5,6) * t143;
t51 = Icges(5,1) * t90 + Icges(5,4) * t89 + Icges(5,5) * t143;
t12 = t43 * t145 + t87 * t47 + t88 * t51;
t44 = Icges(4,5) * t88 + Icges(4,6) * t87 + Icges(4,3) * t145;
t48 = Icges(4,4) * t88 + Icges(4,2) * t87 + Icges(4,6) * t145;
t52 = Icges(4,1) * t88 + Icges(4,4) * t87 + Icges(4,5) * t145;
t13 = t44 * t145 + t87 * t48 + t88 * t52;
t45 = Icges(4,5) * t90 + Icges(4,6) * t89 + Icges(4,3) * t143;
t49 = Icges(4,4) * t90 + Icges(4,2) * t89 + Icges(4,6) * t143;
t53 = Icges(4,1) * t90 + Icges(4,4) * t89 + Icges(4,5) * t143;
t14 = t45 * t145 + t87 * t49 + t88 * t53;
t75 = -Icges(5,5) * t119 + (Icges(5,1) * t118 - Icges(5,4) * t115) * t116;
t26 = t67 * t145 + t87 * t71 + t88 * t75;
t76 = -Icges(4,5) * t119 + (Icges(4,1) * t118 - Icges(4,4) * t115) * t116;
t27 = t68 * t145 + t87 * t72 + t88 * t76;
t170 = (-t26 - t27) * t119 + ((t12 + t14) * t120 + (t11 + t13) * t117) * t116;
t15 = t42 * t143 + t89 * t46 + t90 * t50;
t16 = t43 * t143 + t89 * t47 + t90 * t51;
t17 = t44 * t143 + t89 * t48 + t90 * t52;
t18 = t45 * t143 + t89 * t49 + t90 * t53;
t28 = t67 * t143 + t89 * t71 + t90 * t75;
t29 = t68 * t143 + t89 * t72 + t90 * t76;
t169 = (-t28 - t29) * t119 + ((t16 + t18) * t120 + (t15 + t17) * t117) * t116;
t21 = -t119 * t42 + (-t115 * t46 + t118 * t50) * t116;
t23 = -t119 * t44 + (-t115 * t48 + t118 * t52) * t116;
t168 = -t21 - t23;
t22 = -t119 * t43 + (-t115 * t47 + t118 * t51) * t116;
t24 = -t119 * t45 + (-t115 * t49 + t118 * t53) * t116;
t167 = t22 + t24;
t166 = (t75 + t76) * t116 * t118;
t106 = t118 * pkin(3) + pkin(2);
t140 = t119 * t120;
t165 = t90 * rSges(5,1) + t89 * rSges(5,2) + pkin(3) * t142 + t106 * t140 + t176 * t143;
t164 = t88 * rSges(5,1) + t87 * rSges(5,2) - pkin(3) * t139;
t112 = t117 ^ 2;
t163 = t119 ^ 2;
t113 = t120 ^ 2;
t96 = t116 * rSges(3,1) + t119 * rSges(3,2);
t162 = m(3) * t96;
t161 = t117 / 0.2e1;
t160 = -t119 / 0.2e1;
t158 = pkin(2) * t119;
t157 = -pkin(2) + t106;
t156 = pkin(6) + t114;
t155 = t116 * t174 + t172 * t119 - t166;
t154 = (-t156 * t116 + t157 * t119) * t117 + rSges(5,3) * t145 + t164;
t137 = pkin(2) * t140 + pkin(6) * t143;
t153 = -t137 + t165;
t152 = (t156 - rSges(5,3)) * t119 + (rSges(5,1) * t118 - rSges(5,2) * t115 + t157) * t116;
t80 = -t119 * rSges(4,3) + (rSges(4,1) * t118 - rSges(4,2) * t115) * t116;
t99 = t116 * pkin(2) - t119 * pkin(6);
t151 = -t80 - t99;
t150 = t112 * (pkin(6) * t116 + t158) + t120 * t137;
t149 = t120 * rSges(3,3);
t147 = Icges(3,4) * t119;
t136 = t120 * pkin(1) + t117 * pkin(5);
t135 = t112 + t113;
t134 = -t99 - t152;
t59 = t90 * rSges(4,1) + t89 * rSges(4,2) + rSges(4,3) * t143;
t132 = -t88 * rSges(4,1) - t87 * rSges(4,2);
t130 = rSges(3,1) * t119 - rSges(3,2) * t116;
t126 = -Icges(3,2) * t116 + t147;
t125 = Icges(3,5) * t119 - Icges(3,6) * t116;
t124 = rSges(3,1) * t140 - rSges(3,2) * t143 + t117 * rSges(3,3);
t122 = t23 / 0.2e1 + t21 / 0.2e1 + t27 / 0.2e1 + t26 / 0.2e1;
t121 = t24 / 0.2e1 + t22 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1;
t109 = t120 * pkin(5);
t98 = t120 * rSges(2,1) - t117 * rSges(2,2);
t97 = -t117 * rSges(2,1) - t120 * rSges(2,2);
t93 = Icges(3,6) * t119 + t175;
t70 = Icges(3,3) * t117 + t125 * t120;
t69 = -Icges(3,3) * t120 + t125 * t117;
t63 = t124 + t136;
t62 = t149 + t109 + (-pkin(1) - t130) * t117;
t61 = t151 * t120;
t60 = t151 * t117;
t57 = rSges(4,3) * t145 - t132;
t41 = t120 * t124 + (t130 * t117 - t149) * t117;
t40 = t134 * t120;
t39 = t134 * t117;
t38 = t59 + t136 + t137;
t37 = t109 + (-t158 - pkin(1) + (-rSges(4,3) - pkin(6)) * t116) * t117 + t132;
t36 = -t119 * t59 - t80 * t143;
t35 = t119 * t57 + t80 * t145;
t32 = t136 + t165;
t31 = t109 + (-t106 * t119 - t176 * t116 - pkin(1)) * t117 - t164;
t30 = (-t117 * t59 + t120 * t57) * t116;
t25 = t117 * t57 + t120 * t59 + t150;
t20 = -t153 * t119 - t152 * t143;
t19 = t154 * t119 + t152 * t145;
t10 = (-t153 * t117 + t154 * t120) * t116;
t9 = t154 * t117 + t153 * t120 + t150;
t8 = t18 * t117 - t17 * t120;
t7 = t16 * t117 - t15 * t120;
t6 = t14 * t117 - t13 * t120;
t5 = -t11 * t120 + t12 * t117;
t1 = [Icges(2,3) + (Icges(3,4) * t116 + Icges(3,2) * t119 - t172) * t119 + (Icges(3,1) * t116 + t147 - t174) * t116 + m(5) * (t31 ^ 2 + t32 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(2) * (t97 ^ 2 + t98 ^ 2) + m(3) * (t62 ^ 2 + t63 ^ 2) + t166; m(5) * (t40 * t31 + t39 * t32) + m(4) * (t61 * t37 + t60 * t38) + (t126 * t117 * t160 - t62 * t162 - t122 + (t173 - Icges(3,6) * t160 + t93 / 0.2e1) * t120) * t120 + (t117 * t173 + t119 * (Icges(3,6) * t117 + t126 * t120) / 0.2e1 - t63 * t162 + t93 * t161 + t121) * t117; m(5) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(4) * (t25 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(3) * (t135 * t96 ^ 2 + t41 ^ 2) + (t112 * t70 + t7 + t8) * t117 + (-t113 * t69 - t5 - t6 + (-t117 * t69 + t120 * t70) * t117) * t120; t155 * t119 + m(5) * (t19 * t31 + t20 * t32) + m(4) * (t35 * t37 + t36 * t38) + (t122 * t117 + t121 * t120) * t116; m(5) * (t10 * t9 + t19 * t40 + t20 * t39) + m(4) * (t30 * t25 + t35 * t61 + t36 * t60) + ((t7 / 0.2e1 + t8 / 0.2e1) * t120 + (t6 / 0.2e1 + t5 / 0.2e1) * t117) * t116 + t169 * t161 + (t167 * t117 + t168 * t120) * t160 - t170 * t120 / 0.2e1; m(5) * (t10 ^ 2 + t19 ^ 2 + t20 ^ 2) + m(4) * (t30 ^ 2 + t35 ^ 2 + t36 ^ 2) - t155 * t163 + ((-t167 * t119 + t169) * t120 + (t168 * t119 + t170) * t117) * t116; m(5) * (t117 * t32 + t120 * t31) * t116; m(5) * (-t119 * t9 + (t117 * t39 + t120 * t40) * t116); m(5) * (-t119 * t10 + (t117 * t20 + t120 * t19) * t116); m(5) * (t135 * t116 ^ 2 + t163);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
