% Calculate joint inertia matrix for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:19
% DurationCPUTime: 0.81s
% Computational Cost: add. (2005->192), mult. (2004->281), div. (0->0), fcn. (1834->8), ass. (0->109)
t96 = sin(qJ(1));
t92 = t96 ^ 2;
t150 = t96 * pkin(5);
t98 = cos(qJ(1));
t93 = t98 ^ 2;
t135 = t92 + t93;
t94 = qJ(2) + qJ(3);
t85 = qJ(4) + t94;
t80 = sin(t85);
t81 = cos(t85);
t122 = rSges(5,1) * t81 - rSges(5,2) * t80;
t83 = sin(t94);
t84 = cos(t94);
t123 = rSges(4,1) * t84 - rSges(4,2) * t83;
t95 = sin(qJ(2));
t97 = cos(qJ(2));
t106 = Icges(3,5) * t97 - Icges(3,6) * t95;
t149 = t96 * (Icges(3,3) * t96 + t106 * t98);
t148 = t96 / 0.2e1;
t147 = -t98 / 0.2e1;
t99 = -pkin(6) - pkin(5);
t146 = pkin(2) * t95;
t82 = t97 * pkin(2) + pkin(1);
t145 = rSges(3,1) * t97;
t142 = rSges(3,2) * t95;
t139 = t98 * rSges(3,3);
t102 = t96 * rSges(5,3) + t122 * t98;
t16 = t96 * (-t98 * rSges(5,3) + t122 * t96) + t98 * t102;
t77 = t98 * t82;
t90 = t98 * pkin(5);
t138 = t96 * (t90 + (-pkin(1) + t82) * t96) + t98 * (-t98 * pkin(1) - t150 + t77);
t103 = t96 * rSges(4,3) + t123 * t98;
t17 = t96 * (-t98 * rSges(4,3) + t123 * t96) + t98 * t103;
t137 = t96 * rSges(3,3) + t98 * t145;
t134 = Icges(3,4) * t95;
t133 = Icges(3,4) * t97;
t132 = Icges(4,4) * t83;
t131 = Icges(4,4) * t84;
t130 = Icges(5,4) * t80;
t129 = Icges(5,4) * t81;
t107 = -Icges(5,2) * t80 + t129;
t40 = Icges(5,6) * t96 + t107 * t98;
t110 = Icges(5,1) * t81 - t130;
t42 = Icges(5,5) * t96 + t110 * t98;
t120 = -t40 * t80 + t42 * t81;
t39 = -Icges(5,6) * t98 + t107 * t96;
t41 = -Icges(5,5) * t98 + t110 * t96;
t121 = t39 * t80 - t41 * t81;
t104 = Icges(5,5) * t81 - Icges(5,6) * t80;
t37 = -Icges(5,3) * t98 + t104 * t96;
t38 = Icges(5,3) * t96 + t104 * t98;
t1 = t96 * (t92 * t38 + (t121 * t98 + (t120 - t37) * t96) * t98);
t2 = t93 * t37 + (t120 * t96 + (t121 - t38) * t98) * t96;
t128 = -t98 * t2 + t1;
t67 = t83 * rSges(4,1) + t84 * rSges(4,2);
t127 = -t67 - t146;
t61 = t80 * rSges(5,1) + t81 * rSges(5,2);
t126 = -pkin(3) * t83 - t61;
t59 = Icges(5,2) * t81 + t130;
t60 = Icges(5,1) * t80 + t129;
t115 = -t59 * t80 + t60 * t81;
t58 = Icges(5,5) * t80 + Icges(5,6) * t81;
t125 = (t115 * t98 + t81 * t40 + t80 * t42 + t96 * t58) * t148 + (t115 * t96 + t81 * t39 + t80 * t41 - t98 * t58) * t147;
t68 = pkin(3) * t84 + t82;
t62 = t98 * t68;
t6 = (t68 - t82) * t92 + t98 * (t62 - t77) + t16;
t124 = -t142 + t145;
t108 = -Icges(4,2) * t83 + t131;
t47 = -Icges(4,6) * t98 + t108 * t96;
t111 = Icges(4,1) * t84 - t132;
t49 = -Icges(4,5) * t98 + t111 * t96;
t119 = t47 * t83 - t49 * t84;
t48 = Icges(4,6) * t96 + t108 * t98;
t50 = Icges(4,5) * t96 + t111 * t98;
t118 = -t48 * t83 + t50 * t84;
t65 = Icges(4,2) * t84 + t132;
t66 = Icges(4,1) * t83 + t131;
t114 = -t65 * t83 + t66 * t84;
t105 = Icges(4,5) * t84 - Icges(4,6) * t83;
t45 = -Icges(4,3) * t98 + t105 * t96;
t46 = Icges(4,3) * t96 + t105 * t98;
t4 = t96 * (t92 * t46 + (t119 * t98 + (t118 - t45) * t96) * t98);
t5 = t93 * t45 + (t118 * t96 + (t119 - t46) * t98) * t96;
t113 = t1 + t4 + (-t5 - t2) * t98;
t112 = Icges(3,1) * t97 - t134;
t109 = -Icges(3,2) * t95 + t133;
t101 = t126 - t146;
t64 = Icges(4,5) * t83 + Icges(4,6) * t84;
t100 = t125 + (t114 * t98 + t84 * t48 + t83 * t50 + t96 * t64) * t148 + (t114 * t96 + t84 * t47 + t83 * t49 - t98 * t64) * t147;
t91 = -pkin(7) + t99;
t76 = t98 * rSges(2,1) - t96 * rSges(2,2);
t75 = -t96 * rSges(2,1) - t98 * rSges(2,2);
t74 = t95 * rSges(3,1) + t97 * rSges(3,2);
t44 = t127 * t98;
t43 = t127 * t96;
t32 = t150 + (pkin(1) - t142) * t98 + t137;
t31 = t139 + t90 + (-pkin(1) - t124) * t96;
t30 = t126 * t98;
t29 = t126 * t96;
t26 = t101 * t98;
t25 = t101 * t96;
t24 = -t96 * t99 + t103 + t77;
t23 = (rSges(4,3) - t99) * t98 + (-t123 - t82) * t96;
t20 = t98 * (-t98 * t142 + t137) + (t124 * t96 - t139) * t96;
t19 = -t96 * t91 + t102 + t62;
t18 = (rSges(5,3) - t91) * t98 + (-t122 - t68) * t96;
t7 = t17 + t138;
t3 = t6 + t138;
t8 = [t81 * t59 + t80 * t60 + t84 * t65 + t83 * t66 + t97 * (Icges(3,2) * t97 + t134) + t95 * (Icges(3,1) * t95 + t133) + Icges(2,3) + m(2) * (t75 ^ 2 + t76 ^ 2) + m(3) * (t31 ^ 2 + t32 ^ 2) + m(4) * (t23 ^ 2 + t24 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2); (t97 * (Icges(3,6) * t96 + t109 * t98) + t95 * (Icges(3,5) * t96 + t112 * t98)) * t148 + (t97 * (-Icges(3,6) * t98 + t109 * t96) + t95 * (-Icges(3,5) * t98 + t112 * t96)) * t147 + m(4) * (t44 * t23 + t43 * t24) + m(5) * (t26 * t18 + t25 * t19) + m(3) * (-t31 * t98 - t32 * t96) * t74 + (t92 / 0.2e1 + t93 / 0.2e1) * (Icges(3,5) * t95 + Icges(3,6) * t97) + t100; m(5) * (t25 ^ 2 + t26 ^ 2 + t3 ^ 2) + m(4) * (t43 ^ 2 + t44 ^ 2 + t7 ^ 2) + t4 + t92 * t149 + m(3) * (t135 * t74 ^ 2 + t20 ^ 2) + t128 + (t98 * t149 - t5 - t135 * (-Icges(3,3) * t98 + t106 * t96)) * t98; m(5) * (t30 * t18 + t29 * t19) + m(4) * (-t23 * t98 - t24 * t96) * t67 + t100; m(5) * (t29 * t25 + t30 * t26 + t6 * t3) + m(4) * (t17 * t7 + (-t43 * t96 - t44 * t98) * t67) + t113; m(4) * (t135 * t67 ^ 2 + t17 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2 + t6 ^ 2) + t113; m(5) * (-t18 * t98 - t19 * t96) * t61 + t125; m(5) * (t16 * t3 + (-t25 * t96 - t26 * t98) * t61) + t128; m(5) * (t16 * t6 + (-t29 * t96 - t30 * t98) * t61) + t128; m(5) * (t135 * t61 ^ 2 + t16 ^ 2) + t128;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t8(1), t8(2), t8(4), t8(7); t8(2), t8(3), t8(5), t8(8); t8(4), t8(5), t8(6), t8(9); t8(7), t8(8), t8(9), t8(10);];
Mq = res;
