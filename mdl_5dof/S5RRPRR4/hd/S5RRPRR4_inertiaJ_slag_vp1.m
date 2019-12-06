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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:31:50
% EndTime: 2019-12-05 18:31:51
% DurationCPUTime: 0.58s
% Computational Cost: add. (2310->163), mult. (1388->225), div. (0->0), fcn. (1200->10), ass. (0->95)
t83 = qJ(1) + qJ(2);
t77 = pkin(9) + t83;
t74 = sin(t77);
t75 = cos(t77);
t129 = t74 * t75;
t71 = t74 ^ 2;
t72 = t75 ^ 2;
t128 = t74 / 0.2e1;
t127 = t75 / 0.2e1;
t82 = qJ(4) + qJ(5);
t78 = sin(t82);
t80 = cos(t82);
t91 = Icges(6,5) * t80 - Icges(6,6) * t78;
t23 = Icges(6,3) * t75 - t91 * t74;
t24 = Icges(6,3) * t74 + t91 * t75;
t126 = t74 * (t23 * t129 + t71 * t24) + t75 * (t24 * t129 + t72 * t23);
t84 = sin(qJ(4));
t86 = cos(qJ(4));
t63 = t84 * rSges(5,1) + t86 * rSges(5,2);
t125 = m(5) * t63;
t50 = t78 * rSges(6,1) + t80 * rSges(6,2);
t124 = m(6) * t50;
t79 = sin(t83);
t123 = pkin(2) * t79;
t81 = cos(t83);
t122 = pkin(2) * t81;
t85 = sin(qJ(1));
t121 = t85 * pkin(1);
t87 = cos(qJ(1));
t120 = t87 * pkin(1);
t76 = t86 * pkin(4) + pkin(3);
t119 = -pkin(3) + t76;
t118 = rSges(5,1) * t86;
t117 = rSges(6,1) * t80;
t116 = rSges(5,2) * t84;
t115 = rSges(6,2) * t78;
t114 = t75 * rSges(6,3) + t74 * t115;
t113 = t75 * rSges(5,3) + t74 * t116;
t112 = t71 + t72;
t111 = Icges(5,4) * t84;
t110 = Icges(5,4) * t86;
t109 = Icges(6,4) * t78;
t108 = Icges(6,4) * t80;
t47 = Icges(6,5) * t78 + Icges(6,6) * t80;
t93 = -Icges(6,2) * t78 + t108;
t95 = Icges(6,1) * t80 - t109;
t48 = Icges(6,2) * t80 + t109;
t49 = Icges(6,1) * t78 + t108;
t98 = t48 * t78 - t49 * t80;
t107 = (t80 * (Icges(6,6) * t74 + t93 * t75) + t78 * (Icges(6,5) * t74 + t95 * t75) + t74 * t47 - t98 * t75) * t128 + (t80 * (Icges(6,6) * t75 - t93 * t74) + t78 * (Icges(6,5) * t75 - t95 * t74) + t75 * t47 + t98 * t74) * t127;
t106 = -pkin(3) - t118;
t105 = pkin(4) * t84 + t50;
t52 = -t81 * rSges(3,1) + t79 * rSges(3,2);
t104 = -t76 - t117;
t103 = -t74 * rSges(6,3) + t75 * t115;
t51 = -t79 * rSges(3,1) - t81 * rSges(3,2);
t61 = Icges(5,2) * t86 + t111;
t62 = Icges(5,1) * t84 + t110;
t97 = t61 * t84 - t62 * t86;
t96 = Icges(5,1) * t86 - t111;
t94 = -Icges(5,2) * t84 + t110;
t92 = Icges(5,5) * t86 - Icges(5,6) * t84;
t41 = -t75 * rSges(4,1) + t74 * rSges(4,2) - t122;
t90 = t80 * t48 + t78 * t49 + t86 * t61 + t84 * t62 + Icges(3,3) + Icges(4,3);
t60 = Icges(5,5) * t84 + Icges(5,6) * t86;
t89 = t107 + (t86 * (Icges(5,6) * t74 + t94 * t75) + t84 * (Icges(5,5) * t74 + t96 * t75) + t74 * t60 - t97 * t75) * t128 + (t86 * (Icges(5,6) * t75 - t94 * t74) + t84 * (Icges(5,5) * t75 - t96 * t74) + t75 * t60 + t97 * t74) * t127;
t40 = -t74 * rSges(4,1) - t75 * rSges(4,2) - t123;
t70 = t75 * pkin(7);
t20 = t106 * t74 + t113 - t123 + t70;
t88 = -pkin(8) - pkin(7);
t16 = t104 * t74 - t75 * t88 + t114 - t123;
t66 = t74 * t88;
t17 = t104 * t75 + t103 - t122 + t66;
t58 = t75 * t116;
t21 = -t122 + t58 + t106 * t75 + (-rSges(5,3) - pkin(7)) * t74;
t65 = -t87 * rSges(2,1) + t85 * rSges(2,2);
t64 = -t85 * rSges(2,1) - t87 * rSges(2,2);
t45 = t52 - t120;
t44 = t51 - t121;
t39 = t41 - t120;
t38 = t40 - t121;
t33 = Icges(5,3) * t74 + t92 * t75;
t32 = Icges(5,3) * t75 - t92 * t74;
t31 = t105 * t75;
t30 = t105 * t74;
t29 = -t74 * t117 + t114;
t22 = t75 * (t75 * t117 - t103);
t19 = t21 - t120;
t18 = t20 - t121;
t13 = t17 - t120;
t12 = t16 - t121;
t11 = t75 * (t74 * rSges(5,3) + t75 * t118 - t58) - t74 * (-t74 * t118 + t113);
t6 = -t74 * t29 + t22;
t3 = t22 + (t119 * t75 - t66) * t75 + (-t29 + t70 + t119 * t74 + (-pkin(7) + t88) * t75) * t74;
t1 = [Icges(2,3) + m(2) * (t64 ^ 2 + t65 ^ 2) + m(3) * (t44 ^ 2 + t45 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2) + t90; m(3) * (t51 * t44 + t52 * t45) + m(4) * (t40 * t38 + t41 * t39) + m(5) * (t20 * t18 + t21 * t19) + m(6) * (t16 * t12 + t17 * t13) + t90; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2) + m(3) * (t51 ^ 2 + t52 ^ 2) + t90; 0; 0; m(4) + m(5) + m(6); m(6) * (-t31 * t12 + t30 * t13) + (-t18 * t75 + t19 * t74) * t125 + t89; m(6) * (-t31 * t16 + t30 * t17) + (-t20 * t75 + t21 * t74) * t125 + t89; m(5) * t11 + m(6) * t3; m(5) * (t112 * t63 ^ 2 + t11 ^ 2) + t75 * (t33 * t129 + t72 * t32) + t74 * (t32 * t129 + t71 * t33) + m(6) * (t3 ^ 2 + t30 ^ 2 + t31 ^ 2) + t126; (-t12 * t75 + t13 * t74) * t124 + t107; (-t16 * t75 + t17 * t74) * t124 + t107; m(6) * t6; m(6) * (t6 * t3 + (t30 * t74 + t31 * t75) * t50) + t126; m(6) * (t112 * t50 ^ 2 + t6 ^ 2) + t126;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
