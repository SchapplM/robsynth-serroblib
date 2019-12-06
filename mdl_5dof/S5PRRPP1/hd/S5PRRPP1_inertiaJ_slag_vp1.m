% Calculate joint inertia matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP1_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:12
% EndTime: 2019-12-05 16:06:16
% DurationCPUTime: 0.93s
% Computational Cost: add. (1288->117), mult. (1047->162), div. (0->0), fcn. (909->6), ass. (0->59)
t71 = qJ(3) + pkin(8);
t69 = cos(t71);
t133 = t69 ^ 2;
t70 = pkin(7) + qJ(2);
t66 = sin(t70);
t63 = t66 ^ 2;
t68 = cos(t70);
t65 = t68 ^ 2;
t105 = t63 + t65;
t131 = Icges(6,4) + Icges(5,5);
t130 = Icges(5,6) - Icges(6,6);
t67 = sin(t71);
t73 = sin(qJ(3));
t74 = cos(qJ(3));
t126 = Icges(4,5) * t74 - Icges(4,6) * t73 - t130 * t67 + t131 * t69;
t125 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t121 = t66 * pkin(6);
t118 = rSges(6,1) + pkin(4);
t120 = t118 * t69;
t119 = rSges(6,3) + qJ(5);
t117 = t125 * t68 - t126 * t66;
t116 = t125 * t66 + t126 * t68;
t113 = pkin(3) * t73;
t112 = rSges(4,1) * t74;
t111 = rSges(4,2) * t73;
t110 = t67 * t68;
t109 = t68 * rSges(4,3);
t108 = t68 * t69;
t62 = t74 * pkin(3) + pkin(2);
t51 = t68 * t62;
t61 = t68 * pkin(6);
t107 = t66 * (t61 + (-pkin(2) + t62) * t66) + t68 * (-t68 * pkin(2) - t121 + t51);
t106 = t66 * rSges(4,3) + t68 * t112;
t98 = qJ(5) * t67;
t97 = -t67 * rSges(5,1) - t69 * rSges(5,2) - t113;
t72 = -qJ(4) - pkin(6);
t96 = -t66 * t72 + t51;
t95 = -t118 * t67 + t119 * t69 - t113;
t93 = t66 * rSges(6,2) + rSges(6,3) * t110 + t108 * t118 + t68 * t98;
t92 = -t111 + t112;
t91 = rSges(5,1) * t69 - rSges(5,2) * t67;
t75 = rSges(5,1) * t108 - rSges(5,2) * t110 + t66 * rSges(5,3);
t56 = t73 * rSges(4,1) + t74 * rSges(4,2);
t44 = t68 * rSges(3,1) - t66 * rSges(3,2);
t40 = -t66 * rSges(3,1) - t68 * rSges(3,2);
t27 = t97 * t68;
t26 = t97 * t66;
t13 = t121 + (pkin(2) - t111) * t68 + t106;
t12 = t109 + t61 + (-pkin(2) - t92) * t66;
t9 = t95 * t68;
t8 = t95 * t66;
t7 = t75 + t96;
t6 = (rSges(5,3) - t72) * t68 + (-t62 - t91) * t66;
t5 = t68 * (-t68 * t111 + t106) + (t92 * t66 - t109) * t66;
t4 = t93 + t96;
t3 = (rSges(6,2) - t72) * t68 + (-t119 * t67 - t120 - t62) * t66;
t2 = t68 * t75 + (-t68 * rSges(5,3) + t91 * t66) * t66 + t107;
t1 = t93 * t68 + (-t68 * rSges(6,2) + (rSges(6,3) * t67 + t120 + t98) * t66) * t66 + t107;
t10 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(4,2) * t74 ^ 2 + Icges(3,3) + m(3) * (t40 ^ 2 + t44 ^ 2) + m(4) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t6 ^ 2 + t7 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2) + (Icges(5,2) + Icges(6,3)) * t133 + (Icges(4,1) * t73 + 0.2e1 * Icges(4,4) * t74) * t73 + (0.2e1 * (Icges(5,4) - Icges(6,5)) * t69 + (Icges(5,1) + Icges(6,1)) * t67) * t67; m(4) * t5 + m(5) * t2 + m(6) * t1; m(5) * (t26 * t7 + t27 * t6) + m(6) * (t9 * t3 + t8 * t4) + m(4) * (-t12 * t68 - t13 * t66) * t56 + t105 * (Icges(4,5) * t73 + Icges(4,6) * t74 + t130 * t69 + t131 * t67); m(6) * (t1 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t2 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(4) * (t105 * t56 ^ 2 + t5 ^ 2) + t116 * t66 * t63 + (t117 * t65 + (t116 * t68 + t117 * t66) * t66) * t68; 0; m(5) * (t66 * t6 - t68 * t7) + m(6) * (t66 * t3 - t68 * t4); m(6) * (t66 * t9 - t68 * t8) + m(5) * (-t68 * t26 + t66 * t27); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t105; -m(6) * t69; m(6) * (t3 * t68 + t4 * t66) * t67; m(6) * (-t69 * t1 + (t66 * t8 + t68 * t9) * t67); 0; m(6) * (t105 * t67 ^ 2 + t133);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
