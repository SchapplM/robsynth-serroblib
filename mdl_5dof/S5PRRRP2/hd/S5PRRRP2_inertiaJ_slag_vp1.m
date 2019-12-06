% Calculate joint inertia matrix for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP2_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:38
% EndTime: 2019-12-05 16:41:40
% DurationCPUTime: 0.61s
% Computational Cost: add. (1441->103), mult. (1029->143), div. (0->0), fcn. (881->6), ass. (0->61)
t77 = cos(qJ(4));
t126 = t77 ^ 2;
t125 = Icges(6,4) + Icges(5,5);
t124 = Icges(5,6) - Icges(6,6);
t76 = sin(qJ(4));
t122 = t125 * t76;
t123 = t124 * t77;
t117 = t122 + t123;
t74 = pkin(8) + qJ(2);
t73 = qJ(3) + t74;
t69 = sin(t73);
t66 = t69 ^ 2;
t70 = cos(t73);
t67 = t70 ^ 2;
t121 = -t124 * t76 + t125 * t77;
t119 = Icges(6,2) + Icges(5,3);
t113 = rSges(6,1) + pkin(4);
t115 = t113 * t77;
t114 = rSges(6,3) + qJ(5);
t112 = t119 * t70 - t121 * t69;
t111 = t119 * t69 + t121 * t70;
t57 = t76 * rSges(5,1) + t77 * rSges(5,2);
t108 = m(5) * t57;
t107 = m(6) * t76;
t71 = sin(t74);
t106 = pkin(2) * t71;
t105 = rSges(5,1) * t77;
t104 = t70 * t76;
t103 = t70 * t77;
t102 = t69 * t76 * rSges(5,2) + t70 * rSges(5,3);
t101 = -t113 * t76 + t114 * t77;
t100 = t70 * pkin(3) + t69 * pkin(7);
t99 = t66 + t67;
t94 = qJ(5) * t76;
t36 = t70 * rSges(4,1) - t69 * rSges(4,2);
t93 = t69 * rSges(6,2) + rSges(6,3) * t104 + t113 * t103 + t70 * t94;
t35 = -t69 * rSges(4,1) - t70 * rSges(4,2);
t80 = rSges(5,1) * t103 - rSges(5,2) * t104 + t69 * rSges(5,3);
t79 = Icges(4,3) + (Icges(5,2) + Icges(6,3)) * t126 + ((Icges(5,1) + Icges(6,1)) * t76 + (2 * Icges(5,4) - 2 * Icges(6,5)) * t77) * t76;
t10 = t93 + t100;
t78 = t117 * t67 + (t123 / 0.2e1 + t122 / 0.2e1 + t117 / 0.2e1) * t66;
t18 = t80 + t100;
t64 = t70 * pkin(7);
t17 = t64 + (-pkin(3) - t105) * t69 + t102;
t61 = t70 * rSges(6,2);
t9 = t61 + t64 + (-t114 * t76 - pkin(3) - t115) * t69;
t72 = cos(t74);
t68 = pkin(2) * t72;
t38 = t72 * rSges(3,1) - t71 * rSges(3,2);
t37 = -t71 * rSges(3,1) - t72 * rSges(3,2);
t34 = t36 + t68;
t33 = t35 - t106;
t20 = t101 * t70;
t19 = t101 * t69;
t16 = t18 + t68;
t15 = t17 - t106;
t8 = t68 + t10;
t7 = t9 - t106;
t6 = t69 * (t69 * t105 - t102) + t70 * t80;
t1 = t93 * t70 + (-t61 + (rSges(6,3) * t76 + t115 + t94) * t69) * t69;
t2 = [m(2) + m(3) + m(4) + m(5) + m(6); 0; Icges(3,3) + m(6) * (t7 ^ 2 + t8 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2) + m(3) * (t37 ^ 2 + t38 ^ 2) + t79; 0; m(6) * (t10 * t8 + t9 * t7) + m(5) * (t17 * t15 + t18 * t16) + m(4) * (t35 * t33 + t36 * t34) + t79; m(6) * (t10 ^ 2 + t9 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t35 ^ 2 + t36 ^ 2) + t79; m(5) * t6 + m(6) * t1; m(6) * (t19 * t8 + t20 * t7) + (-t15 * t70 - t16 * t69) * t108 + t78; m(6) * (t19 * t10 + t20 * t9) + (-t17 * t70 - t18 * t69) * t108 + t78; m(5) * (t99 * t57 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t19 ^ 2 + t20 ^ 2) + t111 * t69 * t66 + (t112 * t67 + (t111 * t70 + t112 * t69) * t69) * t70; -m(6) * t77; (t69 * t8 + t7 * t70) * t107; (t10 * t69 + t70 * t9) * t107; m(6) * (-t77 * t1 + (t19 * t69 + t20 * t70) * t76); m(6) * (t99 * t76 ^ 2 + t126);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
