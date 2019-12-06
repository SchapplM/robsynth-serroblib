% Calculate joint inertia matrix for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_inertiaJ_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:02:52
% EndTime: 2019-12-05 17:02:56
% DurationCPUTime: 1.01s
% Computational Cost: add. (2681->215), mult. (3747->349), div. (0->0), fcn. (3943->8), ass. (0->122)
t112 = cos(qJ(2));
t105 = t112 ^ 2;
t109 = sin(qJ(2));
t106 = qJ(3) + qJ(4);
t100 = sin(t106);
t101 = cos(t106);
t139 = Icges(5,4) * t101;
t117 = -Icges(5,2) * t100 + t139;
t64 = Icges(5,6) * t109 + t117 * t112;
t140 = Icges(5,4) * t100;
t119 = Icges(5,1) * t101 - t140;
t66 = Icges(5,5) * t109 + t119 * t112;
t124 = -t100 * t64 + t101 * t66;
t63 = -Icges(5,6) * t112 + t117 * t109;
t65 = -Icges(5,5) * t112 + t119 * t109;
t125 = t100 * t63 - t101 * t65;
t115 = Icges(5,5) * t101 - Icges(5,6) * t100;
t61 = -Icges(5,3) * t112 + t115 * t109;
t62 = Icges(5,3) * t109 + t115 * t112;
t138 = t100 * t109;
t110 = cos(qJ(5));
t133 = t112 * t110;
t107 = sin(qJ(5));
t136 = t109 * t107;
t80 = -t101 * t136 - t133;
t134 = t112 * t107;
t135 = t109 * t110;
t81 = t101 * t135 - t134;
t39 = Icges(6,5) * t81 + Icges(6,6) * t80 + Icges(6,3) * t138;
t41 = Icges(6,4) * t81 + Icges(6,2) * t80 + Icges(6,6) * t138;
t43 = Icges(6,1) * t81 + Icges(6,4) * t80 + Icges(6,5) * t138;
t12 = t39 * t138 + t80 * t41 + t81 * t43;
t137 = t100 * t112;
t82 = -t101 * t134 + t135;
t83 = t101 * t133 + t136;
t40 = Icges(6,5) * t83 + Icges(6,6) * t82 + Icges(6,3) * t137;
t42 = Icges(6,4) * t83 + Icges(6,2) * t82 + Icges(6,6) * t137;
t44 = Icges(6,1) * t83 + Icges(6,4) * t82 + Icges(6,5) * t137;
t13 = t40 * t138 + t80 * t42 + t81 * t44;
t8 = t13 * t109 - t12 * t112;
t162 = -t105 * t61 - (t124 * t109 + (t125 - t62) * t112) * t109 - t8;
t108 = sin(qJ(3));
t111 = cos(qJ(3));
t161 = rSges(4,1) * t111 - rSges(4,2) * t108;
t104 = t109 ^ 2;
t93 = -t108 * rSges(4,1) - t111 * rSges(4,2);
t160 = m(4) * t93;
t88 = -t100 * rSges(5,1) - t101 * rSges(5,2);
t159 = m(5) * t88;
t58 = t101 * rSges(6,3) + (-rSges(6,1) * t110 + rSges(6,2) * t107) * t100;
t158 = m(6) * t58;
t157 = -t109 / 0.2e1;
t156 = t109 / 0.2e1;
t155 = -t112 / 0.2e1;
t154 = t112 / 0.2e1;
t14 = t39 * t137 + t82 * t41 + t83 * t43;
t15 = t40 * t137 + t82 * t42 + t83 * t44;
t9 = t15 * t109 - t14 * t112;
t153 = (t104 * t62 + t9 + (t125 * t112 + (t124 - t61) * t109) * t112) * t109;
t152 = pkin(2) * t108;
t151 = pkin(2) * t111;
t55 = Icges(6,3) * t101 + (-Icges(6,5) * t110 + Icges(6,6) * t107) * t100;
t56 = Icges(6,6) * t101 + (-Icges(6,4) * t110 + Icges(6,2) * t107) * t100;
t150 = t100 * t107 * t56 + t101 * t55;
t148 = rSges(5,1) * t101;
t57 = Icges(6,5) * t101 + (-Icges(6,1) * t110 + Icges(6,4) * t107) * t100;
t146 = t110 * t57;
t145 = t112 * rSges(5,3);
t19 = t101 * t39 + (t107 * t41 - t110 * t43) * t100;
t144 = t19 * t112;
t20 = t101 * t40 + (t107 * t42 - t110 * t44) * t100;
t143 = t20 * t109;
t142 = Icges(4,4) * t108;
t141 = Icges(4,4) * t111;
t132 = t104 + t105;
t46 = t83 * rSges(6,1) + t82 * rSges(6,2) + rSges(6,3) * t137;
t131 = t58 - t152;
t130 = t88 - t152;
t23 = t55 * t138 + t80 * t56 + t81 * t57;
t3 = t23 * t101 + (t109 * t12 + t112 * t13) * t100;
t24 = t55 * t137 + t82 * t56 + t83 * t57;
t4 = t24 * t101 + (t109 * t14 + t112 * t15) * t100;
t129 = t8 * t138 / 0.2e1 + t3 * t155 + t4 * t156 + t101 * (t143 - t144) / 0.2e1 + t9 * t137 / 0.2e1;
t128 = t132 * t151;
t127 = -t81 * rSges(6,1) - t80 * rSges(6,2);
t126 = -rSges(5,2) * t100 + t148;
t86 = -Icges(5,2) * t101 - t140;
t87 = -Icges(5,1) * t100 - t139;
t123 = t100 * t86 - t101 * t87;
t45 = rSges(6,3) * t138 - t127;
t27 = -t109 * t45 - t112 * t46;
t68 = -rSges(5,2) * t137 + t109 * rSges(5,3) + t112 * t148;
t38 = -t109 * (t126 * t109 - t145) - t112 * t68;
t120 = Icges(4,1) * t111 - t142;
t118 = -Icges(4,2) * t108 + t141;
t116 = Icges(4,5) * t111 - Icges(4,6) * t108;
t114 = t112 * t162 + t153;
t85 = -Icges(5,5) * t100 - Icges(5,6) * t101;
t113 = -t143 / 0.2e1 + t144 / 0.2e1 + (t123 * t109 + t112 * t85) * t155 + (-t109 * t85 + t123 * t112) * t156 + (-t100 * t66 - t101 * t64 + t24) * t157 + (-t100 * t65 - t101 * t63 + t23) * t154;
t99 = t112 * t151;
t95 = t112 * rSges(3,1) - t109 * rSges(3,2);
t94 = -t109 * rSges(3,1) - t112 * rSges(3,2);
t77 = t109 * rSges(4,3) + t112 * t161;
t75 = -t112 * rSges(4,3) + t109 * t161;
t70 = Icges(4,3) * t109 + t116 * t112;
t69 = -Icges(4,3) * t112 + t116 * t109;
t60 = t130 * t112;
t59 = t130 * t109;
t53 = t68 + t99;
t52 = t145 + (-t126 - t151) * t109;
t49 = t131 * t112;
t48 = t131 * t109;
t47 = -t109 * t75 - t112 * t77;
t37 = t99 + t46;
t36 = (-rSges(6,3) * t100 - t151) * t109 + t127;
t31 = -t128 + t38;
t30 = t101 * t45 - t58 * t138;
t29 = -t101 * t46 + t58 * t137;
t28 = t100 * t146 - t150;
t26 = (t109 * t46 - t112 * t45) * t100;
t25 = -t128 + t27;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t95 + m(4) * t77 + m(5) * t53 + m(6) * t37; -t101 * t86 - t108 * (-Icges(4,1) * t108 - t141) - t111 * (-Icges(4,2) * t111 - t142) + Icges(3,3) + (-t87 - t146) * t100 + m(6) * (t36 ^ 2 + t37 ^ 2) + m(5) * (t52 ^ 2 + t53 ^ 2) + m(4) * (t75 ^ 2 + t77 ^ 2) + m(3) * (t94 ^ 2 + t95 ^ 2) + t150; m(5) * t59 + m(6) * t48 + t109 * t160; (-t108 * (Icges(4,5) * t109 + t120 * t112) - t111 * (Icges(4,6) * t109 + t118 * t112)) * t157 + (-t108 * (-Icges(4,5) * t112 + t120 * t109) - t111 * (-Icges(4,6) * t112 + t118 * t109)) * t154 + m(6) * (t49 * t36 + t48 * t37) + m(5) * (t60 * t52 + t59 * t53) + (t109 * t77 - t112 * t75) * t160 + (-t104 / 0.2e1 - t105 / 0.2e1) * (-Icges(4,5) * t108 - Icges(4,6) * t111) + t113; m(6) * (t25 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t31 ^ 2 + t59 ^ 2 + t60 ^ 2) + t109 * t104 * t70 + m(4) * (t132 * t93 ^ 2 + t47 ^ 2) + t153 + (-t105 * t69 + (-t109 * t69 + t112 * t70) * t109 + t162) * t112; (t158 + t159) * t109; (t109 * t53 + t112 * t52) * t159 + (t109 * t37 + t112 * t36) * t158 + t113; m(6) * (t27 * t25 + (t109 * t48 + t112 * t49) * t58) + m(5) * (t38 * t31 + (t109 * t59 + t112 * t60) * t88) + t114; m(5) * (t132 * t88 ^ 2 + t38 ^ 2) + m(6) * (t132 * t58 ^ 2 + t27 ^ 2) + t114; m(6) * t29; t101 * t28 + m(6) * (t29 * t37 + t30 * t36) + ((-t20 / 0.2e1 - t24 / 0.2e1) * t112 + (-t19 / 0.2e1 - t23 / 0.2e1) * t109) * t100; m(6) * (t26 * t25 + t29 * t48 + t30 * t49) + t129; m(6) * (t26 * t27 + (t109 * t29 + t112 * t30) * t58) + t129; m(6) * (t26 ^ 2 + t29 ^ 2 + t30 ^ 2) - t101 ^ 2 * t28 + (t112 * t4 + t101 * (t109 * t19 + t112 * t20) + t109 * t3) * t100;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
