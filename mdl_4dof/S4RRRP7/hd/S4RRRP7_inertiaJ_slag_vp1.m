% Calculate joint inertia matrix for
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP7_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP7_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:09
% EndTime: 2019-12-31 17:20:12
% DurationCPUTime: 1.24s
% Computational Cost: add. (1674->241), mult. (4102->364), div. (0->0), fcn. (4484->6), ass. (0->123)
t112 = sin(qJ(2));
t167 = Icges(3,5) * t112;
t166 = t167 / 0.2e1;
t162 = rSges(5,3) + qJ(4);
t164 = rSges(5,1) + pkin(3);
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t116 = cos(qJ(1));
t132 = t116 * t114;
t113 = sin(qJ(1));
t115 = cos(qJ(2));
t135 = t113 * t115;
t88 = t111 * t135 + t132;
t133 = t116 * t111;
t89 = t114 * t135 - t133;
t165 = -t162 * t88 - t164 * t89;
t68 = -Icges(4,3) * t115 + (Icges(4,5) * t114 - Icges(4,6) * t111) * t112;
t71 = -Icges(5,2) * t115 + (Icges(5,4) * t114 + Icges(5,6) * t111) * t112;
t163 = -t68 - t71;
t138 = t112 * t113;
t42 = Icges(5,5) * t89 + Icges(5,6) * t138 + Icges(5,3) * t88;
t46 = Icges(5,4) * t89 + Icges(5,2) * t138 + Icges(5,6) * t88;
t50 = Icges(5,1) * t89 + Icges(5,4) * t138 + Icges(5,5) * t88;
t11 = t46 * t138 + t88 * t42 + t89 * t50;
t136 = t112 * t116;
t90 = -t113 * t114 + t115 * t133;
t91 = t113 * t111 + t115 * t132;
t43 = Icges(5,5) * t91 + Icges(5,6) * t136 + Icges(5,3) * t90;
t47 = Icges(5,4) * t91 + Icges(5,2) * t136 + Icges(5,6) * t90;
t51 = Icges(5,1) * t91 + Icges(5,4) * t136 + Icges(5,5) * t90;
t12 = t47 * t138 + t88 * t43 + t89 * t51;
t44 = Icges(4,5) * t89 - Icges(4,6) * t88 + Icges(4,3) * t138;
t48 = Icges(4,4) * t89 - Icges(4,2) * t88 + Icges(4,6) * t138;
t52 = Icges(4,1) * t89 - Icges(4,4) * t88 + Icges(4,5) * t138;
t13 = t44 * t138 - t88 * t48 + t89 * t52;
t45 = Icges(4,5) * t91 - Icges(4,6) * t90 + Icges(4,3) * t136;
t49 = Icges(4,4) * t91 - Icges(4,2) * t90 + Icges(4,6) * t136;
t53 = Icges(4,1) * t91 - Icges(4,4) * t90 + Icges(4,5) * t136;
t14 = t45 * t138 - t88 * t49 + t89 * t53;
t67 = -Icges(5,6) * t115 + (Icges(5,5) * t114 + Icges(5,3) * t111) * t112;
t75 = -Icges(5,4) * t115 + (Icges(5,1) * t114 + Icges(5,5) * t111) * t112;
t28 = t71 * t138 + t88 * t67 + t89 * t75;
t72 = -Icges(4,6) * t115 + (Icges(4,4) * t114 - Icges(4,2) * t111) * t112;
t76 = -Icges(4,5) * t115 + (Icges(4,1) * t114 - Icges(4,4) * t111) * t112;
t29 = t68 * t138 - t88 * t72 + t89 * t76;
t161 = (-t28 - t29) * t115 + ((t12 + t14) * t116 + (t11 + t13) * t113) * t112;
t15 = t46 * t136 + t90 * t42 + t91 * t50;
t16 = t47 * t136 + t90 * t43 + t91 * t51;
t17 = t44 * t136 - t90 * t48 + t91 * t52;
t18 = t45 * t136 - t90 * t49 + t91 * t53;
t30 = t71 * t136 + t90 * t67 + t91 * t75;
t31 = t68 * t136 - t90 * t72 + t91 * t76;
t160 = (-t30 - t31) * t115 + ((t16 + t18) * t116 + (t15 + t17) * t113) * t112;
t19 = -t115 * t46 + (t111 * t42 + t114 * t50) * t112;
t21 = -t115 * t44 + (-t111 * t48 + t114 * t52) * t112;
t159 = -t19 - t21;
t20 = -t115 * t47 + (t111 * t43 + t114 * t51) * t112;
t22 = -t115 * t45 + (-t111 * t49 + t114 * t53) * t112;
t158 = t20 + t22;
t139 = t111 * t112;
t157 = t67 * t139 + (t75 + t76) * t112 * t114;
t156 = t113 ^ 2;
t155 = t116 ^ 2;
t97 = t112 * rSges(3,1) + t115 * rSges(3,2);
t154 = m(3) * t97;
t153 = t113 / 0.2e1;
t152 = -t115 / 0.2e1;
t150 = pkin(2) * t115;
t149 = t163 * t115 - t72 * t139 + t157;
t148 = rSges(5,2) * t138 - t165;
t147 = rSges(5,2) * t136 + t162 * t90 + t164 * t91;
t145 = -t115 * rSges(5,2) + (t162 * t111 + t164 * t114) * t112;
t134 = t115 * t116;
t131 = pkin(2) * t134 + pkin(6) * t136;
t144 = t156 * (pkin(6) * t112 + t150) + t116 * t131;
t143 = t116 * rSges(3,3);
t100 = t112 * pkin(2) - t115 * pkin(6);
t80 = -t115 * rSges(4,3) + (rSges(4,1) * t114 - rSges(4,2) * t111) * t112;
t142 = -t100 - t80;
t140 = Icges(3,4) * t115;
t130 = t116 * pkin(1) + t113 * pkin(5);
t129 = -t100 - t145;
t57 = t91 * rSges(4,1) - t90 * rSges(4,2) + rSges(4,3) * t136;
t128 = -pkin(1) - t150;
t127 = t130 + t131;
t126 = -t89 * rSges(4,1) + t88 * rSges(4,2);
t125 = rSges(3,1) * t115 - rSges(3,2) * t112;
t121 = -Icges(3,2) * t112 + t140;
t120 = Icges(3,5) * t115 - Icges(3,6) * t112;
t119 = rSges(3,1) * t134 - rSges(3,2) * t136 + t113 * rSges(3,3);
t118 = t21 / 0.2e1 + t19 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1;
t117 = t22 / 0.2e1 + t20 / 0.2e1 + t31 / 0.2e1 + t30 / 0.2e1;
t109 = t116 * pkin(5);
t99 = t116 * rSges(2,1) - t113 * rSges(2,2);
t98 = -t113 * rSges(2,1) - t116 * rSges(2,2);
t94 = Icges(3,6) * t115 + t167;
t70 = Icges(3,3) * t113 + t120 * t116;
t69 = -Icges(3,3) * t116 + t120 * t113;
t63 = t119 + t130;
t62 = t143 + t109 + (-pkin(1) - t125) * t113;
t59 = t142 * t116;
t58 = t142 * t113;
t55 = rSges(4,3) * t138 - t126;
t41 = t116 * t119 + (t125 * t113 - t143) * t113;
t40 = t129 * t116;
t39 = t129 * t113;
t38 = t127 + t57;
t37 = t109 + ((-rSges(4,3) - pkin(6)) * t112 + t128) * t113 + t126;
t36 = -t115 * t57 - t80 * t136;
t35 = t115 * t55 + t80 * t138;
t32 = (-t113 * t57 + t116 * t55) * t112;
t27 = t127 + t147;
t26 = t109 + ((-rSges(5,2) - pkin(6)) * t112 + t128) * t113 + t165;
t25 = t113 * t55 + t116 * t57 + t144;
t24 = -t147 * t115 - t145 * t136;
t23 = t148 * t115 + t145 * t138;
t10 = (-t147 * t113 + t148 * t116) * t112;
t9 = t148 * t113 + t147 * t116 + t144;
t8 = t18 * t113 - t17 * t116;
t7 = t16 * t113 - t15 * t116;
t6 = t14 * t113 - t13 * t116;
t5 = -t11 * t116 + t12 * t113;
t1 = [Icges(2,3) + (Icges(3,1) * t112 - t111 * t72 + t140) * t112 + (Icges(3,4) * t112 + Icges(3,2) * t115 + t163) * t115 + m(5) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t62 ^ 2 + t63 ^ 2) + m(2) * (t98 ^ 2 + t99 ^ 2) + t157; m(5) * (t40 * t26 + t39 * t27) + m(4) * (t59 * t37 + t58 * t38) + (t121 * t113 * t152 - t62 * t154 - t118 + (t166 - Icges(3,6) * t152 + t94 / 0.2e1) * t116) * t116 + (-t63 * t154 + t113 * t166 + t115 * (Icges(3,6) * t113 + t121 * t116) / 0.2e1 + t94 * t153 + t117) * t113; m(5) * (t39 ^ 2 + t40 ^ 2 + t9 ^ 2) + m(4) * (t25 ^ 2 + t58 ^ 2 + t59 ^ 2) + m(3) * (t41 ^ 2 + (t155 + t156) * t97 ^ 2) + (t156 * t70 + t7 + t8) * t113 + (-t155 * t69 - t5 - t6 + (-t113 * t69 + t116 * t70) * t113) * t116; -t149 * t115 + m(5) * (t23 * t26 + t24 * t27) + m(4) * (t35 * t37 + t36 * t38) + (t118 * t113 + t117 * t116) * t112; m(5) * (t10 * t9 + t23 * t40 + t24 * t39) + m(4) * (t32 * t25 + t35 * t59 + t36 * t58) + ((t7 / 0.2e1 + t8 / 0.2e1) * t116 + (t6 / 0.2e1 + t5 / 0.2e1) * t113) * t112 + t160 * t153 + (t158 * t113 + t159 * t116) * t152 - t161 * t116 / 0.2e1; m(5) * (t10 ^ 2 + t23 ^ 2 + t24 ^ 2) + m(4) * (t32 ^ 2 + t35 ^ 2 + t36 ^ 2) + t149 * t115 ^ 2 + ((-t158 * t115 + t160) * t116 + (t159 * t115 + t161) * t113) * t112; m(5) * (t90 * t26 + t88 * t27); m(5) * (t9 * t139 + t88 * t39 + t90 * t40); m(5) * (t10 * t139 + t90 * t23 + t88 * t24); m(5) * (t112 ^ 2 * t111 ^ 2 + t88 ^ 2 + t90 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
