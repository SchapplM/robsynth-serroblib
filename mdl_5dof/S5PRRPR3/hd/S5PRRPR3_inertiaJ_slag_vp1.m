% Calculate joint inertia matrix for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR3_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:03
% EndTime: 2019-12-05 16:19:06
% DurationCPUTime: 0.83s
% Computational Cost: add. (1905->166), mult. (1310->239), div. (0->0), fcn. (1156->8), ass. (0->82)
t83 = qJ(3) + pkin(9);
t76 = sin(t83);
t78 = cos(t83);
t85 = sin(qJ(3));
t86 = cos(qJ(3));
t137 = Icges(4,5) * t86 + Icges(5,5) * t78 - Icges(4,6) * t85 - Icges(5,6) * t76;
t136 = Icges(4,3) + Icges(5,3);
t82 = pkin(8) + qJ(2);
t75 = sin(t82);
t135 = t75 * pkin(6);
t77 = cos(t82);
t134 = t75 * t77;
t79 = qJ(5) + t83;
t70 = sin(t79);
t71 = cos(t79);
t106 = rSges(6,1) * t71 - rSges(6,2) * t70;
t107 = rSges(5,1) * t78 - rSges(5,2) * t76;
t133 = t136 * t77 - t137 * t75;
t132 = t136 * t75 + t137 * t77;
t73 = t75 ^ 2;
t74 = t77 ^ 2;
t131 = t75 / 0.2e1;
t130 = -t77 / 0.2e1;
t129 = pkin(3) * t85;
t72 = t86 * pkin(3) + pkin(2);
t128 = rSges(4,1) * t86;
t125 = rSges(4,2) * t85;
t122 = t77 * rSges(4,3);
t84 = -qJ(4) - pkin(6);
t59 = t77 * t72;
t69 = t77 * pkin(6);
t121 = t75 * (t69 + (-pkin(2) + t72) * t75) + t77 * (-t77 * pkin(2) - t135 + t59);
t88 = t75 * rSges(6,3) + t106 * t77;
t9 = t75 * (-t77 * rSges(6,3) + t106 * t75) + t77 * t88;
t120 = t75 * rSges(4,3) + t77 * t128;
t119 = t73 + t74;
t118 = Icges(4,4) * t85;
t117 = Icges(4,4) * t86;
t116 = Icges(5,4) * t76;
t115 = Icges(5,4) * t78;
t114 = Icges(6,4) * t70;
t113 = Icges(6,4) * t71;
t45 = Icges(6,5) * t70 + Icges(6,6) * t71;
t93 = -Icges(6,2) * t70 + t113;
t96 = Icges(6,1) * t71 - t114;
t46 = Icges(6,2) * t71 + t114;
t47 = Icges(6,1) * t70 + t113;
t99 = -t46 * t70 + t47 * t71;
t112 = (t71 * (Icges(6,6) * t75 + t77 * t93) + t70 * (Icges(6,5) * t75 + t96 * t77) + t75 * t45 + t99 * t77) * t131 + (t71 * (-Icges(6,6) * t77 + t75 * t93) + t70 * (-Icges(6,5) * t77 + t96 * t75) - t77 * t45 + t99 * t75) * t130;
t90 = Icges(6,5) * t71 - Icges(6,6) * t70;
t23 = -Icges(6,3) * t77 + t75 * t90;
t24 = Icges(6,3) * t75 + t77 * t90;
t111 = -t77 * (-t24 * t134 + t74 * t23) + t75 * (-t23 * t134 + t73 * t24);
t110 = -t76 * rSges(5,1) - t78 * rSges(5,2) - t129;
t108 = -t125 + t128;
t98 = Icges(4,1) * t86 - t118;
t97 = Icges(5,1) * t78 - t116;
t95 = -Icges(4,2) * t85 + t117;
t94 = -Icges(5,2) * t76 + t115;
t89 = t75 * rSges(5,3) + t107 * t77;
t48 = t70 * rSges(6,1) + t71 * rSges(6,2);
t87 = -pkin(4) * t76 - t129 - t48;
t81 = -pkin(7) + t84;
t64 = t85 * rSges(4,1) + t86 * rSges(4,2);
t56 = pkin(4) * t78 + t72;
t54 = t77 * rSges(3,1) - t75 * rSges(3,2);
t52 = -t75 * rSges(3,1) - t77 * rSges(3,2);
t43 = t77 * t56;
t36 = t110 * t77;
t35 = t110 * t75;
t20 = t135 + (pkin(2) - t125) * t77 + t120;
t19 = t122 + t69 + (-pkin(2) - t108) * t75;
t16 = t87 * t77;
t15 = t87 * t75;
t14 = -t75 * t84 + t59 + t89;
t13 = (rSges(5,3) - t84) * t77 + (-t107 - t72) * t75;
t12 = -t75 * t81 + t43 + t88;
t11 = (rSges(6,3) - t81) * t77 + (-t106 - t56) * t75;
t10 = t77 * (-t77 * t125 + t120) + (t108 * t75 - t122) * t75;
t4 = t77 * t89 + (-t77 * rSges(5,3) + t107 * t75) * t75 + t121;
t3 = t77 * (t43 - t59) + (t56 - t72) * t73 + t9 + t121;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; t71 * t46 + t70 * t47 + t78 * (Icges(5,2) * t78 + t116) + t76 * (Icges(5,1) * t76 + t115) + t86 * (Icges(4,2) * t86 + t118) + t85 * (Icges(4,1) * t85 + t117) + Icges(3,3) + m(6) * (t11 ^ 2 + t12 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(4) * (t19 ^ 2 + t20 ^ 2) + m(3) * (t52 ^ 2 + t54 ^ 2); m(4) * t10 + m(5) * t4 + m(6) * t3; m(6) * (t16 * t11 + t15 * t12) + m(5) * (t36 * t13 + t35 * t14) + m(4) * (-t19 * t77 - t20 * t75) * t64 + t112 + (t78 * (Icges(5,6) * t75 + t94 * t77) + t76 * (Icges(5,5) * t75 + t97 * t77) + t86 * (Icges(4,6) * t75 + t95 * t77) + t85 * (Icges(4,5) * t75 + t98 * t77)) * t131 + (t78 * (-Icges(5,6) * t77 + t94 * t75) + t76 * (-Icges(5,5) * t77 + t97 * t75) + t86 * (-Icges(4,6) * t77 + t95 * t75) + t85 * (-Icges(4,5) * t77 + t98 * t75)) * t130 + (Icges(4,5) * t85 + Icges(5,5) * t76 + Icges(4,6) * t86 + Icges(5,6) * t78) * (t73 / 0.2e1 + t74 / 0.2e1); m(6) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2 + t4 ^ 2) + m(4) * (t119 * t64 ^ 2 + t10 ^ 2) + t111 + t133 * t77 * t74 + (t132 * t73 + (t132 * t77 + t133 * t75) * t77) * t75; 0; m(6) * (t75 * t11 - t77 * t12) + m(5) * (t75 * t13 - t77 * t14); m(6) * (-t77 * t15 + t75 * t16) + m(5) * (-t77 * t35 + t75 * t36); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t119; m(6) * t9; m(6) * (-t11 * t77 - t12 * t75) * t48 + t112; m(6) * (t9 * t3 + (-t15 * t75 - t16 * t77) * t48) + t111; 0; m(6) * (t119 * t48 ^ 2 + t9 ^ 2) + t111;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
