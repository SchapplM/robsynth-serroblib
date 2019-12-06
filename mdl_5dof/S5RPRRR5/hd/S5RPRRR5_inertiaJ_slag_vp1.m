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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:16:06
% EndTime: 2019-12-05 18:16:09
% DurationCPUTime: 0.55s
% Computational Cost: add. (2249->156), mult. (1344->219), div. (0->0), fcn. (1168->10), ass. (0->93)
t79 = qJ(1) + pkin(9);
t76 = qJ(3) + t79;
t71 = sin(t76);
t72 = cos(t76);
t126 = t71 * t72;
t69 = t71 ^ 2;
t70 = t72 ^ 2;
t125 = t71 / 0.2e1;
t124 = t72 / 0.2e1;
t80 = qJ(4) + qJ(5);
t77 = sin(t80);
t78 = cos(t80);
t87 = Icges(6,5) * t78 - Icges(6,6) * t77;
t23 = Icges(6,3) * t72 - t87 * t71;
t24 = Icges(6,3) * t71 + t87 * t72;
t123 = t71 * (t23 * t126 + t69 * t24) + t72 * (t24 * t126 + t70 * t23);
t81 = sin(qJ(4));
t83 = cos(qJ(4));
t61 = t81 * rSges(5,1) + t83 * rSges(5,2);
t122 = m(5) * t61;
t50 = t77 * rSges(6,1) + t78 * rSges(6,2);
t121 = m(6) * t50;
t82 = sin(qJ(1));
t120 = t82 * pkin(1);
t84 = cos(qJ(1));
t119 = t84 * pkin(1);
t73 = t83 * pkin(4) + pkin(3);
t118 = -pkin(3) + t73;
t117 = rSges(5,1) * t83;
t116 = rSges(6,1) * t78;
t115 = rSges(5,2) * t81;
t114 = rSges(6,2) * t77;
t113 = t72 * rSges(6,3) + t71 * t114;
t112 = t72 * rSges(5,3) + t71 * t115;
t111 = t69 + t70;
t110 = Icges(5,4) * t81;
t109 = Icges(5,4) * t83;
t108 = Icges(6,4) * t77;
t107 = Icges(6,4) * t78;
t47 = Icges(6,5) * t77 + Icges(6,6) * t78;
t89 = -Icges(6,2) * t77 + t107;
t91 = Icges(6,1) * t78 - t108;
t48 = Icges(6,2) * t78 + t108;
t49 = Icges(6,1) * t77 + t107;
t94 = t48 * t77 - t49 * t78;
t106 = (t78 * (Icges(6,6) * t71 + t89 * t72) + t77 * (Icges(6,5) * t71 + t91 * t72) + t71 * t47 - t94 * t72) * t125 + (t78 * (Icges(6,6) * t72 - t89 * t71) + t77 * (Icges(6,5) * t72 - t91 * t71) + t72 * t47 + t94 * t71) * t124;
t105 = -pkin(3) - t117;
t104 = pkin(4) * t81 + t50;
t45 = -t72 * rSges(4,1) + t71 * rSges(4,2);
t103 = -t73 - t116;
t102 = -t71 * rSges(6,3) + t72 * t114;
t74 = sin(t79);
t101 = -pkin(2) * t74 - t120;
t75 = cos(t79);
t100 = -pkin(2) * t75 - t119;
t59 = Icges(5,2) * t83 + t110;
t60 = Icges(5,1) * t81 + t109;
t99 = t78 * t48 + t77 * t49 + t83 * t59 + t81 * t60 + Icges(4,3);
t44 = -t71 * rSges(4,1) - t72 * rSges(4,2);
t93 = t59 * t81 - t60 * t83;
t92 = Icges(5,1) * t83 - t110;
t90 = -Icges(5,2) * t81 + t109;
t88 = Icges(5,5) * t83 - Icges(5,6) * t81;
t58 = Icges(5,5) * t81 + Icges(5,6) * t83;
t86 = t106 + (t83 * (Icges(5,6) * t71 + t90 * t72) + t81 * (Icges(5,5) * t71 + t92 * t72) + t71 * t58 - t93 * t72) * t125 + (t83 * (Icges(5,6) * t72 - t90 * t71) + t81 * (Icges(5,5) * t72 - t92 * t71) + t72 * t58 + t93 * t71) * t124;
t68 = t72 * pkin(7);
t20 = t105 * t71 + t112 + t68;
t85 = -pkin(8) - pkin(7);
t64 = t71 * t85;
t19 = t103 * t72 + t102 + t64;
t18 = t103 * t71 - t72 * t85 + t113;
t56 = t72 * t115;
t21 = t56 + t105 * t72 + (-rSges(5,3) - pkin(7)) * t71;
t63 = -t84 * rSges(2,1) + t82 * rSges(2,2);
t62 = -t82 * rSges(2,1) - t84 * rSges(2,2);
t43 = -t75 * rSges(3,1) + t74 * rSges(3,2) - t119;
t42 = -t74 * rSges(3,1) - t75 * rSges(3,2) - t120;
t39 = t100 + t45;
t38 = t101 + t44;
t33 = Icges(5,3) * t71 + t88 * t72;
t32 = Icges(5,3) * t72 - t88 * t71;
t31 = t104 * t72;
t30 = t104 * t71;
t29 = -t71 * t116 + t113;
t22 = t72 * (t72 * t116 - t102);
t17 = t100 + t21;
t16 = t101 + t20;
t13 = t100 + t19;
t12 = t101 + t18;
t11 = t72 * (t71 * rSges(5,3) + t72 * t117 - t56) - t71 * (-t71 * t117 + t112);
t6 = -t71 * t29 + t22;
t3 = t22 + (t118 * t72 - t64) * t72 + (-t29 + t68 + t118 * t71 + (-pkin(7) + t85) * t72) * t71;
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t62 ^ 2 + t63 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2) + t99; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t44 * t38 + t45 * t39) + m(5) * (t20 * t16 + t21 * t17) + m(6) * (t18 * t12 + t19 * t13) + t99; 0; m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + t99; m(6) * (-t31 * t12 + t30 * t13) + (-t16 * t72 + t17 * t71) * t122 + t86; m(5) * t11 + m(6) * t3; m(6) * (-t31 * t18 + t30 * t19) + (-t20 * t72 + t21 * t71) * t122 + t86; m(5) * (t111 * t61 ^ 2 + t11 ^ 2) + t72 * (t33 * t126 + t70 * t32) + t71 * (t32 * t126 + t69 * t33) + m(6) * (t3 ^ 2 + t30 ^ 2 + t31 ^ 2) + t123; (-t12 * t72 + t13 * t71) * t121 + t106; m(6) * t6; (-t18 * t72 + t19 * t71) * t121 + t106; m(6) * (t6 * t3 + (t30 * t71 + t31 * t72) * t50) + t123; m(6) * (t111 * t50 ^ 2 + t6 ^ 2) + t123;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
