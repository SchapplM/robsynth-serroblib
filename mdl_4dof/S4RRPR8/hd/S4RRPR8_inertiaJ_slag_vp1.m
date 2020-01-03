% Calculate joint inertia matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_inertiaJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR8_inertiaJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR8_inertiaJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:48
% EndTime: 2019-12-31 17:07:51
% DurationCPUTime: 1.21s
% Computational Cost: add. (884->160), mult. (2160->249), div. (0->0), fcn. (2271->6), ass. (0->85)
t86 = cos(qJ(2));
t156 = t86 ^ 2;
t154 = Icges(4,4) + Icges(3,5);
t153 = Icges(3,6) - Icges(4,6);
t83 = sin(qJ(2));
t149 = -t153 * t83 + t154 * t86;
t148 = Icges(4,2) + Icges(3,3);
t146 = (Icges(3,6) / 0.2e1 - Icges(4,6) / 0.2e1) * t86 + (Icges(3,5) / 0.2e1 + Icges(4,4) / 0.2e1) * t83;
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t55 = -t86 * t82 + t83 * t85;
t145 = t55 / 0.2e1;
t98 = t83 * t82 + t86 * t85;
t142 = -t98 / 0.2e1;
t137 = Icges(5,5) * t145 + Icges(5,6) * t142;
t84 = sin(qJ(1));
t87 = cos(qJ(1));
t136 = t148 * t84 + t149 * t87;
t135 = t148 * t87 - t149 * t84;
t80 = t84 ^ 2;
t81 = t87 ^ 2;
t134 = m(4) / 0.2e1;
t133 = -rSges(5,3) - pkin(6);
t132 = t83 * t87;
t131 = t86 * t87;
t130 = t87 * rSges(4,2);
t129 = t87 * rSges(3,3);
t51 = t55 * t87;
t52 = t98 * t87;
t128 = t52 * rSges(5,1) + t51 * rSges(5,2);
t113 = qJ(3) * t83;
t125 = pkin(2) * t131 + t87 * t113;
t127 = t80 * (pkin(2) * t86 + t113) + t87 * t125;
t63 = t83 * pkin(2) - t86 * qJ(3);
t126 = -t83 * rSges(4,1) + t86 * rSges(4,3) - t63;
t124 = t87 * pkin(1) + t84 * pkin(5);
t123 = t80 + t81;
t50 = t98 * t84;
t122 = Icges(5,1) * t50;
t119 = Icges(5,4) * t50;
t116 = Icges(5,5) * t87;
t49 = t55 * t84;
t115 = Icges(5,2) * t49;
t114 = Icges(5,6) * t87;
t26 = Icges(5,4) * t55 - Icges(5,2) * t98;
t27 = Icges(5,1) * t55 - Icges(5,4) * t98;
t88 = t114 + t115 + t119;
t89 = Icges(5,4) * t49 + t116 + t122;
t112 = t88 * t142 + t89 * t145 + t87 * t137 + t49 * t26 / 0.2e1 + t50 * t27 / 0.2e1;
t18 = Icges(5,4) * t52 + Icges(5,2) * t51 - Icges(5,6) * t84;
t19 = Icges(5,1) * t52 + Icges(5,4) * t51 - Icges(5,5) * t84;
t111 = t98 * t18 / 0.2e1 - t55 * t19 / 0.2e1 + t84 * t137 - t51 * t26 / 0.2e1 - t52 * t27 / 0.2e1;
t110 = rSges(4,1) * t131 + t84 * rSges(4,2) + rSges(4,3) * t132;
t109 = t124 + t125;
t28 = t55 * rSges(5,1) - rSges(5,2) * t98;
t108 = -pkin(3) * t83 - t28 - t63;
t17 = Icges(5,5) * t52 + Icges(5,6) * t51 - Icges(5,3) * t84;
t106 = t87 * ((t87 * t17 + t49 * t18 + t50 * t19) * t84 - (Icges(5,3) * t81 + (0.2e1 * t116 + t122) * t50 + (0.2e1 * t114 + t115 + 0.2e1 * t119) * t49) * t87) - t84 * ((-t84 * t17 + t51 * t18 + t52 * t19) * t84 - (t52 * t89 + t51 * t88 - t84 * (Icges(5,5) * t50 + Icges(5,6) * t49 + Icges(5,3) * t87)) * t87);
t105 = rSges(3,1) * t86 - rSges(3,2) * t83;
t104 = -t50 * rSges(5,1) - t49 * rSges(5,2);
t14 = t108 * t84;
t15 = t108 * t87;
t103 = t14 * t84 + t15 * t87;
t91 = rSges(3,1) * t131 - rSges(3,2) * t132 + t84 * rSges(3,3);
t77 = t87 * pkin(5);
t12 = t77 + t133 * t87 + (-t113 - pkin(1) + (-pkin(2) - pkin(3)) * t86) * t84 + t104;
t72 = pkin(3) * t131;
t13 = t133 * t84 + t109 + t128 + t72;
t90 = m(5) * (t12 * t87 + t13 * t84);
t67 = t87 * rSges(2,1) - t84 * rSges(2,2);
t66 = -t84 * rSges(2,1) - t87 * rSges(2,2);
t65 = t83 * rSges(3,1) + t86 * rSges(3,2);
t32 = t126 * t87;
t31 = t126 * t84;
t30 = t91 + t124;
t29 = t129 + t77 + (-pkin(1) - t105) * t84;
t23 = t109 + t110;
t22 = t130 + t77 + (-pkin(1) + (-rSges(4,1) - pkin(2)) * t86 + (-rSges(4,3) - qJ(3)) * t83) * t84;
t21 = -t84 * rSges(5,3) + t128;
t20 = t87 * rSges(5,3) - t104;
t16 = t87 * t91 + (t105 * t84 - t129) * t84;
t11 = t87 * t110 + (-t130 + (rSges(4,1) * t86 + rSges(4,3) * t83) * t84) * t84 + t127;
t10 = -t84 * t20 - t87 * t21;
t5 = (t21 + t72) * t87 + (t84 * t86 * pkin(3) + t20) * t84 + t127;
t1 = [-t98 * t26 + t55 * t27 + Icges(2,3) + m(5) * (t12 ^ 2 + t13 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2) + m(3) * (t29 ^ 2 + t30 ^ 2) + m(2) * (t66 ^ 2 + t67 ^ 2) + (Icges(3,2) + Icges(4,3)) * t156 + (0.2e1 * (Icges(3,4) - Icges(4,5)) * t86 + (Icges(3,1) + Icges(4,1)) * t83) * t83; (t146 * t87 - t112) * t87 + (t146 * t84 - t111) * t84 + m(5) * (t15 * t12 + t14 * t13) + m(4) * (t32 * t22 + t31 * t23) + m(3) * (-t29 * t87 - t30 * t84) * t65 + (t153 * t86 + t154 * t83) * (t80 / 0.2e1 + t81 / 0.2e1); m(5) * (t14 ^ 2 + t15 ^ 2 + t5 ^ 2) + m(4) * (t11 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(3) * (t123 * t65 ^ 2 + t16 ^ 2) - t106 + t135 * t87 * t81 + (t136 * t80 + (t135 * t84 + t136 * t87) * t87) * t84; 0.2e1 * (t90 / 0.2e1 + (t22 * t87 + t23 * t84) * t134) * t83; m(5) * (t103 * t83 - t86 * t5) + m(4) * (-t86 * t11 + (t31 * t84 + t32 * t87) * t83); 0.2e1 * (t134 + m(5) / 0.2e1) * (t123 * t83 ^ 2 + t156); t111 * t84 + t112 * t87 + t28 * t90; m(5) * (t10 * t5 + t103 * t28) + t106; m(5) * (t123 * t83 * t28 - t10 * t86); m(5) * (t123 * t28 ^ 2 + t10 ^ 2) - t106;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
