% Calculate joint inertia matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:31
% EndTime: 2019-12-05 15:30:35
% DurationCPUTime: 0.80s
% Computational Cost: add. (1637->196), mult. (2015->279), div. (0->0), fcn. (2101->6), ass. (0->95)
t77 = sin(pkin(8));
t78 = cos(pkin(8));
t80 = sin(qJ(4));
t81 = cos(qJ(4));
t48 = -Icges(6,3) * t78 + (Icges(6,5) * t81 - Icges(6,6) * t80) * t77;
t49 = -Icges(5,3) * t78 + (Icges(5,5) * t81 - Icges(5,6) * t80) * t77;
t123 = -t48 - t49;
t50 = -Icges(6,6) * t78 + (Icges(6,4) * t81 - Icges(6,2) * t80) * t77;
t51 = -Icges(5,6) * t78 + (Icges(5,4) * t81 - Icges(5,2) * t80) * t77;
t122 = (-t50 - t51) * t80;
t76 = pkin(7) + qJ(2);
t74 = cos(t76);
t111 = t74 * t77;
t73 = sin(t76);
t112 = t73 * t80;
t70 = t81 * pkin(4) + pkin(3);
t114 = t70 * t78;
t107 = t78 * t80;
t58 = -t74 * t107 + t73 * t81;
t106 = t78 * t81;
t59 = t74 * t106 + t112;
t121 = -t59 * rSges(6,1) - t58 * rSges(6,2) - rSges(6,3) * t111 - pkin(4) * t112 - t74 * t114;
t52 = -Icges(6,5) * t78 + (Icges(6,1) * t81 - Icges(6,4) * t80) * t77;
t53 = -Icges(5,5) * t78 + (Icges(5,1) * t81 - Icges(5,4) * t80) * t77;
t120 = (t52 + t53) * t77 * t81;
t110 = t74 * t80;
t56 = -t73 * t107 - t74 * t81;
t57 = t73 * t106 - t110;
t119 = t57 * rSges(6,1) + t56 * rSges(6,2) - pkin(4) * t110;
t118 = t78 ^ 2;
t117 = pkin(3) * t78;
t116 = -pkin(3) + t70;
t79 = -qJ(5) - pkin(6);
t115 = pkin(6) + t79;
t113 = t73 * t77;
t105 = t77 * t122 + t123 * t78 + t120;
t87 = Icges(6,3) * t77;
t21 = Icges(6,5) * t57 + Icges(6,6) * t56 + t73 * t87;
t88 = Icges(5,3) * t77;
t23 = Icges(5,5) * t57 + Icges(5,6) * t56 + t73 * t88;
t104 = t23 + t21;
t22 = Icges(6,5) * t59 + Icges(6,6) * t58 + t74 * t87;
t24 = Icges(5,5) * t59 + Icges(5,6) * t58 + t74 * t88;
t103 = t24 + t22;
t89 = Icges(6,6) * t77;
t26 = Icges(6,4) * t59 + Icges(6,2) * t58 + t74 * t89;
t90 = Icges(5,6) * t77;
t28 = Icges(5,4) * t59 + Icges(5,2) * t58 + t74 * t90;
t102 = t26 + t28;
t25 = Icges(6,4) * t57 + Icges(6,2) * t56 + t73 * t89;
t27 = Icges(5,4) * t57 + Icges(5,2) * t56 + t73 * t90;
t101 = t27 + t25;
t91 = Icges(6,5) * t77;
t30 = Icges(6,1) * t59 + Icges(6,4) * t58 + t74 * t91;
t92 = Icges(5,5) * t77;
t32 = Icges(5,1) * t59 + Icges(5,4) * t58 + t74 * t92;
t100 = t30 + t32;
t29 = Icges(6,1) * t57 + Icges(6,4) * t56 + t73 * t91;
t31 = Icges(5,1) * t57 + Icges(5,4) * t56 + t73 * t92;
t99 = t31 + t29;
t85 = t115 * t77;
t98 = (t116 * t78 - t85) * t73 + rSges(6,3) * t113 + t119;
t97 = -(-t85 - t117) * t74 + t121;
t96 = (t115 - rSges(6,3)) * t78 + (rSges(6,1) * t81 - rSges(6,2) * t80 + t116) * t77;
t94 = t74 * pkin(2) + t73 * qJ(3);
t93 = t73 ^ 2 + t74 ^ 2;
t38 = t59 * rSges(5,1) + t58 * rSges(5,2) + rSges(5,3) * t111;
t84 = rSges(4,1) * t78 - rSges(4,2) * t77;
t83 = -t57 * rSges(5,1) - t56 * rSges(5,2);
t68 = t74 * qJ(3);
t62 = t74 * rSges(3,1) - t73 * rSges(3,2);
t61 = -t73 * rSges(3,1) - t74 * rSges(3,2);
t55 = -t78 * rSges(5,3) + (rSges(5,1) * t81 - rSges(5,2) * t80) * t77;
t40 = t73 * rSges(4,3) + t84 * t74 + t94;
t39 = t74 * rSges(4,3) + t68 + (-pkin(2) - t84) * t73;
t36 = rSges(5,3) * t113 - t83;
t18 = (pkin(6) * t77 + t117) * t74 + t38 + t94;
t17 = t68 + (-t117 - pkin(2) + (-rSges(5,3) - pkin(6)) * t77) * t73 + t83;
t16 = -t55 * t111 - t78 * t38;
t15 = t55 * t113 + t78 * t36;
t14 = -t79 * t111 - t121 + t94;
t13 = t68 + (-t114 - pkin(2) + (-rSges(6,3) + t79) * t77) * t73 - t119;
t12 = t49 * t111 + t58 * t51 + t59 * t53;
t11 = t48 * t111 + t58 * t50 + t59 * t52;
t10 = t49 * t113 + t56 * t51 + t57 * t53;
t9 = t48 * t113 + t56 * t50 + t57 * t52;
t8 = (t36 * t74 - t38 * t73) * t77;
t7 = -t96 * t111 + t97 * t78;
t6 = t96 * t113 + t98 * t78;
t5 = -t78 * t24 + (-t28 * t80 + t32 * t81) * t77;
t4 = -t78 * t23 + (-t27 * t80 + t31 * t81) * t77;
t3 = -t78 * t22 + (-t26 * t80 + t30 * t81) * t77;
t2 = -t78 * t21 + (-t25 * t80 + t29 * t81) * t77;
t1 = (t97 * t73 + t98 * t74) * t77;
t19 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + (Icges(4,2) * t78 + t123) * t78 + (Icges(4,1) * t77 + 0.2e1 * Icges(4,4) * t78 + t122) * t77 + m(3) * (t61 ^ 2 + t62 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2) + t120; 0; m(4) * (t73 * t39 - t74 * t40) + m(5) * (t73 * t17 - t74 * t18) + m(6) * (t73 * t13 - t74 * t14); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1) * t93; m(5) * t8 + m(6) * t1; -t105 * t78 + m(5) * (t15 * t17 + t16 * t18) + m(6) * (t6 * t13 + t7 * t14) + ((t5 / 0.2e1 + t3 / 0.2e1 + t12 / 0.2e1 + t11 / 0.2e1) * t74 + (t4 / 0.2e1 + t2 / 0.2e1 + t10 / 0.2e1 + t9 / 0.2e1) * t73) * t77; m(5) * (t15 * t73 - t16 * t74) + m(6) * (t6 * t73 - t7 * t74); m(5) * (t15 ^ 2 + t16 ^ 2 + t8 ^ 2) + m(6) * (t1 ^ 2 + t6 ^ 2 + t7 ^ 2) + t105 * t118 + (((t100 * t59 + t102 * t58 + t103 * t111) * t111 + (-t11 - t12 - t3 - t5) * t78) * t74 + ((t101 * t56 + t104 * t113 + t99 * t57) * t113 + (-t10 - t4 - t9 - t2) * t78 + ((t103 * t73 + t104 * t74) * t77 + t99 * t59 + t101 * t58 + t100 * t57 + t102 * t56) * t111) * t73) * t77; -m(6) * t78; m(6) * (t13 * t74 + t14 * t73) * t77; 0; m(6) * (-t78 * t1 + (t6 * t74 + t7 * t73) * t77); m(6) * (t93 * t77 ^ 2 + t118);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
