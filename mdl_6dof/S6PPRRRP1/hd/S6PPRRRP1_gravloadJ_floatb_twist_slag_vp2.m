% Calculate Gravitation load on the joints for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 14:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:50:44
% EndTime: 2018-11-23 14:50:45
% DurationCPUTime: 0.86s
% Computational Cost: add. (2734->106), mult. (2867->155), div. (0->0), fcn. (2814->22), ass. (0->84)
t132 = mrSges(6,1) + mrSges(7,1);
t124 = -mrSges(6,2) - mrSges(7,2);
t56 = cos(qJ(5));
t131 = m(6) * pkin(4) + m(7) * (pkin(5) * t56 + pkin(4)) + mrSges(5,1);
t130 = m(6) * pkin(10) - m(7) * (-qJ(6) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t99 = pkin(7) + qJ(3);
t84 = sin(t99) / 0.2e1;
t100 = pkin(7) - qJ(3);
t89 = sin(t100);
t123 = t84 - t89 / 0.2e1;
t98 = pkin(6) - pkin(12);
t80 = cos(t98) / 0.2e1;
t97 = pkin(6) + pkin(12);
t88 = cos(t97);
t47 = t80 - t88 / 0.2e1;
t58 = cos(qJ(3));
t79 = sin(t97) / 0.2e1;
t87 = sin(t98);
t66 = t79 + t87 / 0.2e1;
t129 = t123 * t66 + t47 * t58;
t102 = sin(pkin(11));
t105 = cos(pkin(11));
t46 = t79 - t87 / 0.2e1;
t51 = cos(pkin(12));
t41 = -t102 * t46 + t105 * t51;
t101 = sin(pkin(12));
t67 = t80 + t88 / 0.2e1;
t62 = t101 * t105 + t102 * t67;
t128 = -t123 * t62 + t41 * t58;
t40 = t102 * t51 + t105 * t46;
t61 = t101 * t102 - t105 * t67;
t127 = -t123 * t61 + t40 * t58;
t116 = m(7) * pkin(5);
t125 = mrSges(4,2) - mrSges(5,3);
t122 = -m(5) - m(6) - m(7);
t53 = sin(qJ(5));
t121 = t124 * t53 + t132 * t56 + t131;
t119 = -t116 - t132;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t118 = t130 * t54 + t131 * t57 + mrSges(4,1);
t117 = m(3) + m(4) - t122;
t104 = sin(pkin(6));
t90 = cos(t99);
t85 = t90 / 0.2e1;
t91 = cos(t100);
t70 = t85 - t91 / 0.2e1;
t65 = t70 * t104;
t23 = t105 * t65 + t127;
t115 = t23 * t53;
t26 = -t102 * t65 + t128;
t114 = t26 * t53;
t107 = cos(pkin(6));
t30 = -t107 * t70 + t129;
t113 = t30 * t53;
t109 = t53 * t57;
t108 = t56 * t57;
t106 = cos(pkin(7));
t103 = sin(pkin(7));
t86 = t91 / 0.2e1;
t82 = t105 * t104;
t81 = t104 * t102;
t78 = t86 - t90 / 0.2e1;
t74 = t78 * t104;
t71 = t86 + t85;
t68 = t84 + t89 / 0.2e1;
t64 = t68 * t104;
t63 = -t103 * t66 + t106 * t107;
t60 = t103 * t62 + t106 * t81;
t59 = t103 * t61 - t106 * t82;
t55 = sin(qJ(3));
t29 = t107 * t78 + t129;
t28 = -t107 * t68 + t47 * t55 - t66 * t71;
t25 = t102 * t74 + t128;
t24 = -t102 * t64 + t41 * t55 + t62 * t71;
t22 = -t105 * t74 + t127;
t21 = t105 * t64 + t40 * t55 + t61 * t71;
t18 = t29 * t57 + t54 * t63;
t17 = t29 * t54 - t57 * t63;
t14 = t25 * t57 + t54 * t60;
t13 = t25 * t54 - t57 * t60;
t12 = t22 * t57 + t54 * t59;
t11 = t22 * t54 - t57 * t59;
t1 = [(-m(2) - t117) * g(3) (-g(1) * t81 + g(2) * t82 - g(3) * t107) * t117 (-t113 * t116 + t122 * (-t28 * pkin(3) + t30 * pkin(9)) + t125 * t30 - t132 * (-t108 * t28 + t113) + t124 * (t109 * t28 + t30 * t56) + t118 * t28) * g(3) + (-t115 * t116 + t122 * (-t21 * pkin(3) + pkin(9) * t23) - t132 * (-t108 * t21 + t115) + t124 * (t109 * t21 + t23 * t56) + t125 * t23 + t118 * t21) * g(2) + (-t114 * t116 + t122 * (-t24 * pkin(3) + pkin(9) * t26) - t132 * (-t108 * t24 + t114) + t124 * (t109 * t24 + t26 * t56) + t125 * t26 + t118 * t24) * g(1) (t121 * t17 - t130 * t18) * g(3) + (t121 * t11 - t12 * t130) * g(2) + (t121 * t13 - t130 * t14) * g(1) (t124 * (-t18 * t56 - t28 * t53) + t119 * (-t18 * t53 + t28 * t56)) * g(3) + (t124 * (-t12 * t56 - t21 * t53) + t119 * (-t12 * t53 + t21 * t56)) * g(2) + (t124 * (-t14 * t56 - t24 * t53) + t119 * (-t14 * t53 + t24 * t56)) * g(1) (-g(1) * t13 - g(2) * t11 - g(3) * t17) * m(7)];
taug  = t1(:);
