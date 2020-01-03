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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:53:47
% EndTime: 2020-01-03 11:53:50
% DurationCPUTime: 0.63s
% Computational Cost: add. (2249->150), mult. (1344->213), div. (0->0), fcn. (1168->10), ass. (0->89)
t89 = qJ(1) + pkin(9);
t84 = qJ(3) + t89;
t79 = sin(t84);
t80 = cos(t84);
t134 = t79 * t80;
t90 = qJ(4) + qJ(5);
t85 = sin(t90);
t86 = cos(t90);
t136 = -rSges(6,1) * t86 + rSges(6,2) * t85;
t91 = sin(qJ(4));
t93 = cos(qJ(4));
t135 = -rSges(5,1) * t93 + rSges(5,2) * t91;
t29 = -t79 * rSges(6,3) + t136 * t80;
t81 = t93 * pkin(4) + pkin(3);
t133 = t80 * t81 - t29;
t75 = t79 ^ 2;
t76 = t80 ^ 2;
t132 = -t79 / 0.2e1;
t131 = -t80 / 0.2e1;
t63 = t91 * rSges(5,1) + t93 * rSges(5,2);
t130 = m(5) * t63;
t50 = t85 * rSges(6,1) + t86 * rSges(6,2);
t129 = m(6) * t50;
t95 = -pkin(8) - pkin(7);
t124 = t79 * t81 + t80 * t95;
t44 = t79 * rSges(4,1) + t80 * rSges(4,2);
t123 = -t80 * pkin(3) - t79 * pkin(7);
t122 = t75 + t76;
t82 = sin(t89);
t92 = sin(qJ(1));
t87 = t92 * pkin(1);
t121 = pkin(2) * t82 + t87;
t83 = cos(t89);
t94 = cos(qJ(1));
t88 = t94 * pkin(1);
t120 = pkin(2) * t83 + t88;
t119 = Icges(5,4) * t91;
t118 = Icges(5,4) * t93;
t117 = Icges(6,4) * t85;
t116 = Icges(6,4) * t86;
t101 = -Icges(6,2) * t85 + t116;
t103 = Icges(6,1) * t86 - t117;
t48 = Icges(6,2) * t86 + t117;
t49 = Icges(6,1) * t85 + t116;
t106 = t48 * t85 - t49 * t86;
t47 = Icges(6,5) * t85 + Icges(6,6) * t86;
t115 = (t106 * t80 + t86 * (-Icges(6,6) * t79 - t101 * t80) + t85 * (-Icges(6,5) * t79 - t103 * t80) - t79 * t47) * t132 + (-t106 * t79 + t86 * (-Icges(6,6) * t80 + t101 * t79) + t85 * (-Icges(6,5) * t80 + t103 * t79) - t80 * t47) * t131;
t114 = pkin(4) * t91 + t50;
t45 = t80 * rSges(4,1) - t79 * rSges(4,2);
t113 = t135 * t79;
t99 = Icges(6,5) * t86 - Icges(6,6) * t85;
t23 = -Icges(6,3) * t80 + t99 * t79;
t24 = -Icges(6,3) * t79 - t99 * t80;
t112 = -t80 * (t24 * t134 + t76 * t23) - t79 * (t23 * t134 + t75 * t24);
t61 = Icges(5,2) * t93 + t119;
t62 = Icges(5,1) * t91 + t118;
t111 = t86 * t48 + t85 * t49 + t93 * t61 + t91 * t62 + Icges(4,3);
t105 = t61 * t91 - t62 * t93;
t104 = Icges(5,1) * t93 - t119;
t102 = -Icges(5,2) * t91 + t118;
t100 = Icges(5,5) * t93 - Icges(5,6) * t91;
t98 = -t79 * rSges(5,3) + t135 * t80;
t60 = Icges(5,5) * t91 + Icges(5,6) * t93;
t97 = t115 + (t105 * t80 + t93 * (-Icges(5,6) * t79 - t102 * t80) + t91 * (-Icges(5,5) * t79 - t104 * t80) - t79 * t60) * t132 + (-t105 * t79 + t93 * (-Icges(5,6) * t80 + t102 * t79) + t91 * (-Icges(5,5) * t80 + t104 * t79) - t80 * t60) * t131;
t96 = -t80 * rSges(6,3) - t136 * t79;
t21 = -t98 - t123;
t18 = t96 + t124;
t19 = -t79 * t95 + t133;
t73 = t79 * pkin(3);
t20 = t73 + (-rSges(5,3) - pkin(7)) * t80 - t113;
t65 = t94 * rSges(2,1) - t92 * rSges(2,2);
t64 = t92 * rSges(2,1) + t94 * rSges(2,2);
t43 = t83 * rSges(3,1) - t82 * rSges(3,2) + t88;
t42 = t82 * rSges(3,1) + t83 * rSges(3,2) + t87;
t39 = t45 + t120;
t38 = t121 + t44;
t33 = -Icges(5,3) * t79 - t100 * t80;
t32 = -Icges(5,3) * t80 + t100 * t79;
t31 = t114 * t80;
t30 = t114 * t79;
t22 = t79 * t96;
t17 = t21 + t120;
t16 = t20 + t121;
t13 = t19 + t120;
t12 = t18 + t121;
t11 = -t80 * t98 + t79 * (-t80 * rSges(5,3) - t113);
t6 = -t80 * t29 + t22;
t3 = t79 * (-t73 + t124) + t22 + ((pkin(7) - t95) * t79 + t123 + t133) * t80;
t1 = [Icges(2,3) + Icges(3,3) + m(2) * (t64 ^ 2 + t65 ^ 2) + m(3) * (t42 ^ 2 + t43 ^ 2) + m(4) * (t38 ^ 2 + t39 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2) + t111; 0; m(3) + m(4) + m(5) + m(6); m(4) * (t44 * t38 + t45 * t39) + m(5) * (t20 * t16 + t21 * t17) + m(6) * (t18 * t12 + t19 * t13) + t111; 0; m(6) * (t18 ^ 2 + t19 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + t111; m(6) * (t31 * t12 - t30 * t13) + (t16 * t80 - t17 * t79) * t130 + t97; m(5) * t11 + m(6) * t3; m(6) * (t31 * t18 - t30 * t19) + (t20 * t80 - t21 * t79) * t130 + t97; m(5) * (t122 * t63 ^ 2 + t11 ^ 2) - t80 * (t33 * t134 + t76 * t32) - t79 * (t32 * t134 + t75 * t33) + m(6) * (t3 ^ 2 + t30 ^ 2 + t31 ^ 2) + t112; (t12 * t80 - t13 * t79) * t129 + t115; m(6) * t6; (t18 * t80 - t19 * t79) * t129 + t115; m(6) * (t6 * t3 + (t30 * t79 + t31 * t80) * t50) + t112; m(6) * (t122 * t50 ^ 2 + t6 ^ 2) + t112;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
