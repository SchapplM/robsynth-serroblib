% Calculate joint inertia matrix for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:33
% EndTime: 2022-01-20 09:48:34
% DurationCPUTime: 0.59s
% Computational Cost: add. (2249->152), mult. (1344->215), div. (0->0), fcn. (1168->10), ass. (0->89)
t85 = qJ(1) + pkin(9);
t81 = qJ(3) + t85;
t76 = sin(t81);
t77 = cos(t81);
t130 = t76 * t77;
t73 = t76 ^ 2;
t74 = t77 ^ 2;
t129 = t76 / 0.2e1;
t128 = -t77 / 0.2e1;
t87 = sin(qJ(4));
t89 = cos(qJ(4));
t62 = t87 * rSges(5,1) + t89 * rSges(5,2);
t127 = m(5) * t62;
t86 = qJ(4) + qJ(5);
t82 = sin(t86);
t83 = cos(t86);
t50 = t82 * rSges(6,1) + t83 * rSges(6,2);
t126 = m(6) * t50;
t88 = sin(qJ(1));
t125 = t88 * pkin(1);
t124 = rSges(5,1) * t89;
t123 = rSges(6,1) * t83;
t122 = rSges(5,2) * t87;
t121 = rSges(6,2) * t82;
t120 = t77 * rSges(6,3) + t76 * t121;
t93 = t76 * rSges(6,3) + (-t121 + t123) * t77;
t6 = t76 * (t76 * t123 - t120) + t77 * t93;
t119 = t77 * rSges(5,3) + t76 * t122;
t118 = -t77 * pkin(3) - t76 * pkin(7);
t117 = t73 + t74;
t80 = cos(t85);
t90 = cos(qJ(1));
t84 = t90 * pkin(1);
t116 = pkin(2) * t80 + t84;
t115 = Icges(5,4) * t87;
t114 = Icges(5,4) * t89;
t113 = Icges(6,4) * t82;
t112 = Icges(6,4) * t83;
t48 = Icges(6,2) * t83 + t113;
t49 = Icges(6,1) * t82 + t112;
t102 = -t48 * t82 + t49 * t83;
t47 = Icges(6,5) * t82 + Icges(6,6) * t83;
t97 = -Icges(6,2) * t82 + t112;
t99 = Icges(6,1) * t83 - t113;
t111 = (t102 * t77 + t83 * (Icges(6,6) * t76 + t97 * t77) + t82 * (Icges(6,5) * t76 + t99 * t77) + t76 * t47) * t129 + (t102 * t76 + t83 * (-Icges(6,6) * t77 + t97 * t76) + t82 * (-Icges(6,5) * t77 + t99 * t76) - t77 * t47) * t128;
t95 = Icges(6,5) * t83 - Icges(6,6) * t82;
t24 = -Icges(6,3) * t77 + t95 * t76;
t25 = Icges(6,3) * t76 + t95 * t77;
t110 = -t77 * (-t25 * t130 + t74 * t24) + t76 * (-t24 * t130 + t73 * t25);
t109 = -pkin(4) * t87 - t50;
t45 = t77 * rSges(4,1) - t76 * rSges(4,2);
t79 = sin(t85);
t108 = -pkin(2) * t79 - t125;
t60 = Icges(5,2) * t89 + t115;
t61 = Icges(5,1) * t87 + t114;
t107 = t83 * t48 + t82 * t49 + t89 * t60 + t87 * t61 + Icges(4,3);
t44 = -t76 * rSges(4,1) - t77 * rSges(4,2);
t101 = -t60 * t87 + t61 * t89;
t100 = Icges(5,1) * t89 - t115;
t98 = -Icges(5,2) * t87 + t114;
t96 = Icges(5,5) * t89 - Icges(5,6) * t87;
t94 = t76 * rSges(5,3) + (-t122 + t124) * t77;
t59 = Icges(5,5) * t87 + Icges(5,6) * t89;
t92 = t111 + (t101 * t77 + t89 * (Icges(5,6) * t76 + t98 * t77) + t87 * (Icges(5,5) * t76 + t100 * t77) + t76 * t59) * t129 + (t101 * t76 + t89 * (-Icges(5,6) * t77 + t98 * t76) + t87 * (-Icges(5,5) * t77 + t100 * t76) - t77 * t59) * t128;
t21 = t94 - t118;
t71 = t77 * pkin(7);
t20 = t71 + (-pkin(3) - t124) * t76 + t119;
t78 = t89 * pkin(4) + pkin(3);
t53 = t77 * t78;
t91 = -pkin(8) - pkin(7);
t19 = -t76 * t91 + t53 + t93;
t18 = -t77 * t91 + (-t78 - t123) * t76 + t120;
t64 = t90 * rSges(2,1) - t88 * rSges(2,2);
t63 = -t88 * rSges(2,1) - t90 * rSges(2,2);
t43 = t80 * rSges(3,1) - t79 * rSges(3,2) + t84;
t42 = -t79 * rSges(3,1) - t80 * rSges(3,2) - t125;
t39 = t45 + t116;
t38 = t108 + t44;
t33 = Icges(5,3) * t76 + t96 * t77;
t32 = -Icges(5,3) * t77 + t96 * t76;
t31 = t109 * t77;
t30 = t109 * t76;
t17 = t21 + t116;
t16 = t108 + t20;
t13 = t19 + t116;
t12 = t108 + t18;
t11 = t76 * (t76 * t124 - t119) + t77 * t94;
t3 = t77 * (t53 + t118) + (t71 + (-pkin(3) + t78) * t76) * t76 + t6;
t1 = [Icges(2,3) + Icges(3,3) + m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(2) * (t63 ^ 2 + t64 ^ 2) + t107; 0; m(3) + m(4) + m(5) + m(6); m(6) * (t18 * t12 + t19 * t13) + m(5) * (t20 * t16 + t21 * t17) + m(4) * (t44 * t38 + t45 * t39) + t107; 0; m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + t107; m(6) * (t31 * t12 + t30 * t13) + (-t16 * t77 - t17 * t76) * t127 + t92; m(5) * t11 + m(6) * t3; m(6) * (t31 * t18 + t30 * t19) + (-t20 * t77 - t21 * t76) * t127 + t92; m(5) * (t117 * t62 ^ 2 + t11 ^ 2) + t76 * (-t32 * t130 + t73 * t33) - t77 * (-t33 * t130 + t74 * t32) + m(6) * (t3 ^ 2 + t30 ^ 2 + t31 ^ 2) + t110; (-t12 * t77 - t13 * t76) * t126 + t111; m(6) * t6; (-t18 * t77 - t19 * t76) * t126 + t111; m(6) * (t6 * t3 + (-t30 * t76 - t31 * t77) * t50) + t110; m(6) * (t117 * t50 ^ 2 + t6 ^ 2) + t110;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
