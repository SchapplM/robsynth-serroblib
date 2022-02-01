% Calculate joint inertia matrix for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:02
% EndTime: 2022-01-20 10:48:03
% DurationCPUTime: 0.58s
% Computational Cost: add. (2310->158), mult. (1388->221), div. (0->0), fcn. (1200->10), ass. (0->91)
t89 = qJ(1) + qJ(2);
t82 = pkin(9) + t89;
t78 = sin(t82);
t79 = cos(t82);
t132 = t78 * t79;
t75 = t78 ^ 2;
t76 = t79 ^ 2;
t131 = t78 / 0.2e1;
t130 = -t79 / 0.2e1;
t90 = sin(qJ(4));
t92 = cos(qJ(4));
t64 = t90 * rSges(5,1) + t92 * rSges(5,2);
t129 = m(5) * t64;
t88 = qJ(4) + qJ(5);
t83 = sin(t88);
t85 = cos(t88);
t50 = t83 * rSges(6,1) + t85 * rSges(6,2);
t128 = m(6) * t50;
t84 = sin(t89);
t127 = pkin(2) * t84;
t91 = sin(qJ(1));
t126 = t91 * pkin(1);
t125 = rSges(5,1) * t92;
t124 = rSges(6,1) * t85;
t123 = rSges(5,2) * t90;
t122 = rSges(6,2) * t83;
t121 = t79 * rSges(6,3) + t78 * t122;
t97 = t78 * rSges(6,3) + (-t122 + t124) * t79;
t6 = t78 * (t78 * t124 - t121) + t79 * t97;
t120 = t79 * rSges(5,3) + t78 * t123;
t119 = -t79 * pkin(3) - t78 * pkin(7);
t118 = t75 + t76;
t117 = Icges(5,4) * t90;
t116 = Icges(5,4) * t92;
t115 = Icges(6,4) * t83;
t114 = Icges(6,4) * t85;
t101 = -Icges(6,2) * t83 + t114;
t103 = Icges(6,1) * t85 - t115;
t48 = Icges(6,2) * t85 + t115;
t49 = Icges(6,1) * t83 + t114;
t106 = -t48 * t83 + t49 * t85;
t47 = Icges(6,5) * t83 + Icges(6,6) * t85;
t113 = (t106 * t79 + t85 * (Icges(6,6) * t78 + t101 * t79) + t83 * (Icges(6,5) * t78 + t103 * t79) + t78 * t47) * t131 + (t106 * t78 + t85 * (-Icges(6,6) * t79 + t101 * t78) + t83 * (-Icges(6,5) * t79 + t103 * t78) - t79 * t47) * t130;
t99 = Icges(6,5) * t85 - Icges(6,6) * t83;
t24 = -Icges(6,3) * t79 + t99 * t78;
t25 = Icges(6,3) * t78 + t99 * t79;
t112 = -t79 * (-t25 * t132 + t76 * t24) + t78 * (-t24 * t132 + t75 * t25);
t111 = -pkin(4) * t90 - t50;
t86 = cos(t89);
t52 = t86 * rSges(3,1) - t84 * rSges(3,2);
t80 = pkin(2) * t86;
t41 = t79 * rSges(4,1) - t78 * rSges(4,2) + t80;
t51 = -t84 * rSges(3,1) - t86 * rSges(3,2);
t62 = Icges(5,2) * t92 + t117;
t63 = Icges(5,1) * t90 + t116;
t105 = -t62 * t90 + t63 * t92;
t104 = Icges(5,1) * t92 - t117;
t102 = -Icges(5,2) * t90 + t116;
t100 = Icges(5,5) * t92 - Icges(5,6) * t90;
t98 = t78 * rSges(5,3) + (-t123 + t125) * t79;
t96 = t85 * t48 + t83 * t49 + t92 * t62 + t90 * t63 + Icges(3,3) + Icges(4,3);
t61 = Icges(5,5) * t90 + Icges(5,6) * t92;
t95 = t113 + (t105 * t79 + t92 * (Icges(5,6) * t78 + t102 * t79) + t90 * (Icges(5,5) * t78 + t104 * t79) + t78 * t61) * t131 + (t105 * t78 + t92 * (-Icges(5,6) * t79 + t102 * t78) + t90 * (-Icges(5,5) * t79 + t104 * t78) - t79 * t61) * t130;
t40 = -t78 * rSges(4,1) - t79 * rSges(4,2) - t127;
t21 = t80 + t98 - t119;
t81 = t92 * pkin(4) + pkin(3);
t55 = t79 * t81;
t94 = -pkin(8) - pkin(7);
t17 = -t78 * t94 + t55 + t80 + t97;
t73 = t79 * pkin(7);
t20 = -t127 + t73 + (-pkin(3) - t125) * t78 + t120;
t16 = -t127 - t79 * t94 + (-t81 - t124) * t78 + t121;
t93 = cos(qJ(1));
t87 = t93 * pkin(1);
t66 = t93 * rSges(2,1) - t91 * rSges(2,2);
t65 = -t91 * rSges(2,1) - t93 * rSges(2,2);
t45 = t52 + t87;
t44 = t51 - t126;
t39 = t41 + t87;
t38 = t40 - t126;
t33 = Icges(5,3) * t78 + t100 * t79;
t32 = -Icges(5,3) * t79 + t100 * t78;
t31 = t111 * t79;
t30 = t111 * t78;
t19 = t21 + t87;
t18 = t20 - t126;
t13 = t17 + t87;
t12 = t16 - t126;
t11 = t78 * (t78 * t125 - t120) + t79 * t98;
t3 = t79 * (t55 + t119) + (t73 + (-pkin(3) + t81) * t78) * t78 + t6;
t1 = [Icges(2,3) + m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2) + m(3) * (t44 ^ 2 + t45 ^ 2) + m(2) * (t65 ^ 2 + t66 ^ 2) + t96; m(6) * (t16 * t12 + t17 * t13) + m(5) * (t20 * t18 + t21 * t19) + m(4) * (t40 * t38 + t41 * t39) + m(3) * (t51 * t44 + t52 * t45) + t96; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2) + m(3) * (t51 ^ 2 + t52 ^ 2) + t96; 0; 0; m(4) + m(5) + m(6); m(6) * (t31 * t12 + t30 * t13) + (-t18 * t79 - t19 * t78) * t129 + t95; m(6) * (t31 * t16 + t30 * t17) + (-t20 * t79 - t21 * t78) * t129 + t95; m(5) * t11 + m(6) * t3; m(5) * (t118 * t64 ^ 2 + t11 ^ 2) + t78 * (-t32 * t132 + t75 * t33) - t79 * (-t33 * t132 + t76 * t32) + m(6) * (t3 ^ 2 + t30 ^ 2 + t31 ^ 2) + t112; (-t12 * t79 - t13 * t78) * t128 + t113; (-t16 * t79 - t17 * t78) * t128 + t113; m(6) * t6; m(6) * (t6 * t3 + (-t30 * t78 - t31 * t79) * t50) + t112; m(6) * (t118 * t50 ^ 2 + t6 ^ 2) + t112;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
