% Calculate joint inertia matrix for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP3_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP3_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:03
% DurationCPUTime: 0.54s
% Computational Cost: add. (905->96), mult. (984->140), div. (0->0), fcn. (846->6), ass. (0->60)
t75 = cos(qJ(3));
t125 = t75 ^ 2;
t124 = Icges(5,4) + Icges(4,5);
t123 = Icges(4,6) - Icges(5,6);
t73 = sin(qJ(3));
t121 = t124 * t73;
t122 = t123 * t75;
t116 = t121 + t122;
t72 = qJ(1) + qJ(2);
t68 = sin(t72);
t66 = t68 ^ 2;
t69 = cos(t72);
t67 = t69 ^ 2;
t120 = -t123 * t73 + t124 * t75;
t118 = Icges(5,2) + Icges(4,3);
t112 = rSges(5,1) + pkin(3);
t114 = t112 * t75;
t113 = rSges(5,3) + qJ(4);
t111 = t118 * t69 - t120 * t68;
t110 = t118 * t68 + t120 * t69;
t49 = t73 * rSges(4,1) + t75 * rSges(4,2);
t107 = m(4) * t49;
t106 = m(5) * t73;
t74 = sin(qJ(1));
t105 = t74 * pkin(1);
t104 = rSges(4,1) * t75;
t103 = t69 * t73;
t102 = t69 * t75;
t101 = -t112 * t73 + t113 * t75;
t100 = t68 * t73 * rSges(4,2) + t69 * rSges(4,3);
t99 = t69 * pkin(2) + t68 * pkin(6);
t98 = t66 + t67;
t93 = qJ(4) * t73;
t36 = t69 * rSges(3,1) - t68 * rSges(3,2);
t92 = t68 * rSges(5,2) + rSges(5,3) * t103 + t112 * t102 + t69 * t93;
t35 = -t68 * rSges(3,1) - t69 * rSges(3,2);
t79 = rSges(4,1) * t102 - rSges(4,2) * t103 + t68 * rSges(4,3);
t78 = Icges(3,3) + (Icges(4,2) + Icges(5,3)) * t125 + ((Icges(4,1) + Icges(5,1)) * t73 + (2 * Icges(4,4) - 2 * Icges(5,5)) * t75) * t73;
t14 = t92 + t99;
t77 = t116 * t67 + (t122 / 0.2e1 + t121 / 0.2e1 + t116 / 0.2e1) * t66;
t18 = t79 + t99;
t64 = t69 * pkin(6);
t17 = t64 + (-pkin(2) - t104) * t68 + t100;
t61 = t69 * rSges(5,2);
t13 = t61 + t64 + (-t113 * t73 - pkin(2) - t114) * t68;
t76 = cos(qJ(1));
t70 = t76 * pkin(1);
t51 = t76 * rSges(2,1) - t74 * rSges(2,2);
t50 = -t74 * rSges(2,1) - t76 * rSges(2,2);
t34 = t36 + t70;
t33 = t35 - t105;
t20 = t101 * t69;
t19 = t101 * t68;
t16 = t18 + t70;
t15 = t17 - t105;
t8 = t70 + t14;
t7 = t13 - t105;
t6 = t68 * (t68 * t104 - t100) + t69 * t79;
t1 = t92 * t69 + (-t61 + (rSges(5,3) * t73 + t114 + t93) * t68) * t68;
t2 = [Icges(2,3) + m(2) * (t50 ^ 2 + t51 ^ 2) + m(3) * (t33 ^ 2 + t34 ^ 2) + m(4) * (t15 ^ 2 + t16 ^ 2) + m(5) * (t7 ^ 2 + t8 ^ 2) + t78; m(3) * (t35 * t33 + t36 * t34) + m(4) * (t17 * t15 + t18 * t16) + m(5) * (t13 * t7 + t14 * t8) + t78; m(3) * (t35 ^ 2 + t36 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + t78; m(5) * (t19 * t8 + t20 * t7) + (-t15 * t69 - t16 * t68) * t107 + t77; m(5) * (t20 * t13 + t19 * t14) + (-t17 * t69 - t18 * t68) * t107 + t77; m(4) * (t98 * t49 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t19 ^ 2 + t20 ^ 2) + t110 * t68 * t66 + (t111 * t67 + (t110 * t69 + t111 * t68) * t68) * t69; (t68 * t8 + t69 * t7) * t106; (t13 * t69 + t14 * t68) * t106; m(5) * (-t75 * t1 + (t19 * t68 + t20 * t69) * t73); m(5) * (t98 * t73 ^ 2 + t125);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
