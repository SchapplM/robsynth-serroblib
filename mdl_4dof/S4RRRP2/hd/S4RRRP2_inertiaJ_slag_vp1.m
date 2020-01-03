% Calculate joint inertia matrix for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP2_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:12:55
% DurationCPUTime: 0.54s
% Computational Cost: add. (852->96), mult. (894->136), div. (0->0), fcn. (760->6), ass. (0->58)
t122 = Icges(4,5) + Icges(5,5);
t121 = Icges(4,6) + Icges(5,6);
t73 = sin(qJ(3));
t119 = t122 * t73;
t75 = cos(qJ(3));
t120 = t121 * t75;
t114 = t119 + t120;
t71 = qJ(1) + qJ(2);
t68 = sin(t71);
t66 = t68 ^ 2;
t69 = cos(t71);
t67 = t69 ^ 2;
t118 = -t121 * t73 + t122 * t75;
t115 = Icges(4,3) + Icges(5,3);
t112 = t115 * t69 - t118 * t68;
t111 = t115 * t68 + t118 * t69;
t50 = t73 * rSges(4,1) + t75 * rSges(4,2);
t108 = m(4) * t50;
t74 = sin(qJ(1));
t107 = t74 * pkin(1);
t106 = rSges(4,1) * t75;
t105 = t68 * t73;
t104 = t69 * t73;
t103 = t69 * t75;
t102 = -rSges(5,2) * t105 - t69 * rSges(5,3);
t101 = rSges(4,2) * t105 + t69 * rSges(4,3);
t100 = t69 * pkin(2) + t68 * pkin(6);
t99 = t66 + t67;
t94 = -t75 * rSges(5,2) + (-rSges(5,1) - pkin(3)) * t73;
t65 = t75 * pkin(3) + pkin(2);
t93 = -rSges(5,1) * t75 - t65;
t36 = t69 * rSges(3,1) - t68 * rSges(3,2);
t92 = Icges(3,3) + (Icges(4,2) + Icges(5,2)) * t75 ^ 2 + ((Icges(4,1) + Icges(5,1)) * t73 + (2 * Icges(4,4) + 2 * Icges(5,4)) * t75) * t73;
t35 = -t68 * rSges(3,1) - t69 * rSges(3,2);
t79 = rSges(4,1) * t103 - rSges(4,2) * t104 + t68 * rSges(4,3);
t78 = rSges(5,1) * t103 - rSges(5,2) * t104 + t68 * rSges(5,3) + t69 * t65;
t77 = t114 * t67 + (t120 / 0.2e1 + t119 / 0.2e1 + t114 / 0.2e1) * t66;
t18 = t79 + t100;
t63 = t69 * pkin(6);
t17 = t63 + (-pkin(2) - t106) * t68 + t101;
t72 = -qJ(4) - pkin(6);
t14 = -t68 * t72 + t78;
t13 = t93 * t68 - t69 * t72 - t102;
t76 = cos(qJ(1));
t70 = t76 * pkin(1);
t52 = t76 * rSges(2,1) - t74 * rSges(2,2);
t51 = -t74 * rSges(2,1) - t76 * rSges(2,2);
t34 = t36 + t70;
t33 = t35 - t107;
t32 = t94 * t69;
t31 = t94 * t68;
t16 = t18 + t70;
t15 = t17 - t107;
t12 = t14 + t70;
t11 = t13 - t107;
t6 = t68 * (t68 * t106 - t101) + t69 * t79;
t1 = (t78 - t100) * t69 + (t63 + (-pkin(2) - t93) * t68 + t102) * t68;
t2 = [Icges(2,3) + m(2) * (t51 ^ 2 + t52 ^ 2) + m(3) * (t33 ^ 2 + t34 ^ 2) + m(4) * (t15 ^ 2 + t16 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2) + t92; m(3) * (t35 * t33 + t36 * t34) + m(4) * (t17 * t15 + t18 * t16) + m(5) * (t13 * t11 + t14 * t12) + t92; m(3) * (t35 ^ 2 + t36 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + t92; m(5) * (t32 * t11 + t31 * t12) + (-t15 * t69 - t16 * t68) * t108 + t77; m(5) * (t32 * t13 + t31 * t14) + (-t17 * t69 - t18 * t68) * t108 + t77; m(4) * (t99 * t50 ^ 2 + t6 ^ 2) + m(5) * (t1 ^ 2 + t31 ^ 2 + t32 ^ 2) + t111 * t68 * t66 + (t112 * t67 + (t111 * t69 + t112 * t68) * t68) * t69; m(5) * (t68 * t11 - t69 * t12); m(5) * (t68 * t13 - t69 * t14); m(5) * (-t69 * t31 + t68 * t32); m(5) * t99;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
