% Calculate joint inertia matrix for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_inertiaJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:39
% EndTime: 2019-12-05 18:08:45
% DurationCPUTime: 1.75s
% Computational Cost: add. (3334->275), mult. (8483->443), div. (0->0), fcn. (10100->8), ass. (0->130)
t122 = sin(qJ(3));
t160 = Icges(4,5) * t122;
t159 = t160 / 0.2e1;
t123 = sin(qJ(1));
t126 = cos(qJ(3));
t127 = cos(qJ(1));
t145 = t122 * t123;
t121 = sin(qJ(4));
t125 = cos(qJ(4));
t140 = t127 * t125;
t142 = t123 * t126;
t98 = t121 * t142 + t140;
t141 = t127 * t121;
t99 = t125 * t142 - t141;
t62 = Icges(5,5) * t99 - Icges(5,6) * t98 + Icges(5,3) * t145;
t65 = Icges(5,4) * t99 - Icges(5,2) * t98 + Icges(5,6) * t145;
t68 = Icges(5,1) * t99 - Icges(5,4) * t98 + Icges(5,5) * t145;
t26 = t62 * t145 - t98 * t65 + t99 * t68;
t100 = -t123 * t125 + t126 * t141;
t101 = t123 * t121 + t126 * t140;
t143 = t122 * t127;
t63 = Icges(5,5) * t101 - Icges(5,6) * t100 + Icges(5,3) * t143;
t66 = Icges(5,4) * t101 - Icges(5,2) * t100 + Icges(5,6) * t143;
t69 = Icges(5,1) * t101 - Icges(5,4) * t100 + Icges(5,5) * t143;
t27 = t63 * t145 - t98 * t66 + t99 * t69;
t120 = sin(qJ(5));
t124 = cos(qJ(5));
t75 = -t99 * t120 + t124 * t145;
t76 = t120 * t145 + t99 * t124;
t43 = Icges(6,5) * t76 + Icges(6,6) * t75 + Icges(6,3) * t98;
t45 = Icges(6,4) * t76 + Icges(6,2) * t75 + Icges(6,6) * t98;
t47 = Icges(6,1) * t76 + Icges(6,4) * t75 + Icges(6,5) * t98;
t14 = t98 * t43 + t75 * t45 + t76 * t47;
t77 = -t101 * t120 + t124 * t143;
t78 = t101 * t124 + t120 * t143;
t44 = Icges(6,5) * t78 + Icges(6,6) * t77 + Icges(6,3) * t100;
t46 = Icges(6,4) * t78 + Icges(6,2) * t77 + Icges(6,6) * t100;
t48 = Icges(6,1) * t78 + Icges(6,4) * t77 + Icges(6,5) * t100;
t15 = t98 * t44 + t75 * t46 + t76 * t48;
t146 = t121 * t122;
t144 = t122 * t125;
t96 = -t120 * t144 - t126 * t124;
t97 = -t126 * t120 + t124 * t144;
t61 = Icges(6,5) * t97 + Icges(6,6) * t96 + Icges(6,3) * t146;
t64 = Icges(6,4) * t97 + Icges(6,2) * t96 + Icges(6,6) * t146;
t67 = Icges(6,1) * t97 + Icges(6,4) * t96 + Icges(6,5) * t146;
t20 = t98 * t61 + t75 * t64 + t76 * t67;
t3 = -t20 * t126 + (t123 * t14 + t127 * t15) * t122;
t83 = -Icges(5,3) * t126 + (Icges(5,5) * t125 - Icges(5,6) * t121) * t122;
t86 = -Icges(5,6) * t126 + (Icges(5,4) * t125 - Icges(5,2) * t121) * t122;
t89 = -Icges(5,5) * t126 + (Icges(5,1) * t125 - Icges(5,4) * t121) * t122;
t37 = t83 * t145 - t98 * t86 + t99 * t89;
t158 = -t37 * t126 + (t123 * t26 + t127 * t27) * t122 + t3;
t28 = -t100 * t65 + t101 * t68 + t62 * t143;
t29 = -t100 * t66 + t101 * t69 + t63 * t143;
t38 = -t100 * t86 + t101 * t89 + t83 * t143;
t16 = t100 * t43 + t77 * t45 + t78 * t47;
t17 = t100 * t44 + t77 * t46 + t78 * t48;
t21 = t100 * t61 + t77 * t64 + t78 * t67;
t4 = -t21 * t126 + (t123 * t16 + t127 * t17) * t122;
t157 = -t38 * t126 + (t123 * t28 + t127 * t29) * t122 + t4;
t150 = rSges(4,1) * t126;
t156 = (-rSges(4,2) * t122 + t150) * t123 - t127 * rSges(4,3);
t118 = t123 ^ 2;
t119 = t127 ^ 2;
t155 = t98 / 0.2e1;
t154 = t100 / 0.2e1;
t153 = t123 / 0.2e1;
t152 = -t126 / 0.2e1;
t151 = -t127 / 0.2e1;
t147 = Icges(4,4) * t126;
t139 = t118 + t119;
t25 = t61 * t146 + t96 * t64 + t97 * t67;
t50 = t78 * rSges(6,1) + t77 * rSges(6,2) + t100 * rSges(6,3);
t72 = t101 * rSges(5,1) - t100 * rSges(5,2) + rSges(5,3) * t143;
t18 = t43 * t146 + t96 * t45 + t97 * t47;
t138 = t20 / 0.2e1 + t18 / 0.2e1;
t19 = t44 * t146 + t96 * t46 + t97 * t48;
t137 = t21 / 0.2e1 + t19 / 0.2e1;
t132 = -Icges(4,2) * t122 + t147;
t131 = Icges(4,5) * t126 - Icges(4,6) * t122;
t130 = -rSges(4,2) * t143 + t123 * rSges(4,3) + t127 * t150;
t31 = -t126 * t62 + (-t121 * t65 + t125 * t68) * t122;
t129 = t37 / 0.2e1 + t31 / 0.2e1 + t138;
t32 = -t126 * t63 + (-t121 * t66 + t125 * t69) * t122;
t128 = t32 / 0.2e1 + t38 / 0.2e1 + t137;
t49 = t76 * rSges(6,1) + t75 * rSges(6,2) + t98 * rSges(6,3);
t71 = t99 * rSges(5,1) - t98 * rSges(5,2) + rSges(5,3) * t145;
t116 = t127 * qJ(2);
t115 = t123 * qJ(2);
t110 = t127 * rSges(2,1) - t123 * rSges(2,2);
t109 = -t123 * rSges(2,1) - t127 * rSges(2,2);
t108 = t122 * rSges(4,1) + t126 * rSges(4,2);
t103 = -t123 * rSges(3,1) + t127 * rSges(3,3) + t116;
t102 = t127 * rSges(3,1) + t123 * rSges(3,3) + t115;
t92 = -t126 * rSges(5,3) + (rSges(5,1) * t125 - rSges(5,2) * t121) * t122;
t85 = Icges(4,3) * t123 + t131 * t127;
t84 = -Icges(4,3) * t127 + t131 * t123;
t81 = t116 - t156;
t80 = t115 + t130;
t79 = t89 * t144;
t70 = t97 * rSges(6,1) + t96 * rSges(6,2) + rSges(6,3) * t146;
t59 = t123 * t156 + t127 * t130;
t58 = t116 - t71;
t57 = t115 + t72;
t53 = -t126 * t72 - t92 * t143;
t52 = t126 * t71 + t92 * t145;
t51 = -t126 * t83 - t86 * t146 + t79;
t42 = t116 - t49;
t41 = t115 + t50;
t40 = t123 * t71 + t127 * t72;
t39 = (-t123 * t72 + t127 * t71) * t122;
t36 = -t126 * t50 - t70 * t143;
t35 = t126 * t49 + t70 * t145;
t34 = -t100 * t70 + t50 * t146;
t33 = -t49 * t146 + t98 * t70;
t30 = t123 * t49 + t127 * t50;
t24 = (-t123 * t50 + t127 * t49) * t122;
t23 = t25 * t146;
t22 = t100 * t49 - t98 * t50;
t13 = t29 * t123 - t28 * t127;
t12 = t27 * t123 - t26 * t127;
t9 = t19 * t123 - t18 * t127;
t8 = t17 * t123 - t16 * t127;
t7 = t15 * t123 - t14 * t127;
t6 = -t25 * t126 + (t18 * t123 + t19 * t127) * t122;
t5 = t19 * t100 + t18 * t98 + t23;
t2 = t17 * t100 + t21 * t146 + t16 * t98;
t1 = t15 * t100 + t14 * t98 + t20 * t146;
t10 = [Icges(3,2) + Icges(2,3) + t79 + (Icges(4,4) * t122 + Icges(4,2) * t126 - t83) * t126 + (Icges(4,1) * t122 - t121 * t86 + t147) * t122 + m(6) * (t41 ^ 2 + t42 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t80 ^ 2 + t81 ^ 2) + m(2) * (t109 ^ 2 + t110 ^ 2) + m(3) * (t102 ^ 2 + t103 ^ 2) + t25; m(6) * (t123 * t42 - t127 * t41) + m(5) * (t123 * t58 - t127 * t57) + m(4) * (t123 * t81 - t127 * t80) + m(3) * (-t127 * t102 + t123 * t103); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t139; (t127 * t159 + (-Icges(4,6) * t127 + t132 * t123) * t152 - t129) * t127 + (t123 * t159 + t126 * (Icges(4,6) * t123 + t132 * t127) / 0.2e1 + t128) * t123 + m(5) * (-t123 * t57 - t127 * t58) * t92 + m(6) * (-t123 * t41 - t127 * t42) * t70 + m(4) * (-t123 * t80 - t127 * t81) * t108 + (t118 / 0.2e1 + t119 / 0.2e1) * (Icges(4,6) * t126 + t160); 0; m(6) * (t139 * t70 ^ 2 + t30 ^ 2) + m(5) * (t139 * t92 ^ 2 + t40 ^ 2) + m(4) * (t139 * t108 ^ 2 + t59 ^ 2) + (t118 * t85 + t13 + t8) * t123 + (-t119 * t84 - t12 - t7 + (-t123 * t84 + t127 * t85) * t123) * t127; (-t25 - t51) * t126 + m(6) * (t35 * t42 + t36 * t41) + m(5) * (t52 * t58 + t53 * t57) + (t129 * t123 + t128 * t127) * t122; m(5) * (t52 * t123 - t53 * t127) + m(6) * (t35 * t123 - t36 * t127); m(6) * (t24 * t30 + (-t123 * t36 - t127 * t35) * t70) + m(5) * (t39 * t40 + (-t123 * t53 - t127 * t52) * t92) + ((t8 / 0.2e1 + t13 / 0.2e1) * t127 + (t7 / 0.2e1 + t12 / 0.2e1) * t123) * t122 + t157 * t153 + (t32 * t123 - t31 * t127 + t9) * t152 + t158 * t151; (t51 * t126 - t6) * t126 + m(6) * (t24 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t39 ^ 2 + t52 ^ 2 + t53 ^ 2) + ((-t126 * t32 + t157) * t127 + (-t126 * t31 + t158) * t123) * t122; t23 + m(6) * (t33 * t42 + t34 * t41) + t138 * t98 + t137 * t100; m(6) * (t33 * t123 - t34 * t127); t2 * t153 + t9 * t146 / 0.2e1 + m(6) * (t22 * t30 + (-t123 * t34 - t127 * t33) * t70) + t7 * t155 + t8 * t154 + t1 * t151; m(6) * (t22 * t24 + t33 * t35 + t34 * t36) + t5 * t152 + t3 * t155 + t4 * t154 + (t121 * t6 / 0.2e1 + t127 * t2 / 0.2e1 + t1 * t153) * t122; m(6) * (t22 ^ 2 + t33 ^ 2 + t34 ^ 2) + t100 * t2 + t98 * t1 + t5 * t146;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
