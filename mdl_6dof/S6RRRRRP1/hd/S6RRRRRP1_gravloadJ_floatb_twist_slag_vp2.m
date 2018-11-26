% Calculate Gravitation load on the joints for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:26:34
% EndTime: 2018-11-23 18:26:35
% DurationCPUTime: 1.01s
% Computational Cost: add. (640->118), mult. (583->127), div. (0->0), fcn. (502->10), ass. (0->72)
t136 = mrSges(6,1) + mrSges(7,1);
t44 = cos(qJ(5));
t135 = t136 * t44;
t124 = mrSges(6,2) + mrSges(7,2);
t134 = -mrSges(6,3) - mrSges(7,3);
t39 = qJ(2) + qJ(3);
t36 = qJ(4) + t39;
t30 = sin(t36);
t133 = t135 * t30;
t132 = t30 * t124;
t31 = cos(t36);
t32 = pkin(5) * t44 + pkin(4);
t40 = -qJ(6) - pkin(10);
t62 = -t30 * t32 - t31 * t40;
t131 = -m(7) * t62 + t133;
t34 = sin(t39);
t35 = cos(t39);
t64 = mrSges(5,1) * t30 + mrSges(5,2) * t31;
t130 = mrSges(4,1) * t34 + mrSges(4,2) * t35 + t64;
t41 = sin(qJ(5));
t46 = cos(qJ(1));
t91 = t41 * t46;
t94 = t31 * t46;
t129 = -t91 * t132 + t134 * t94;
t43 = sin(qJ(1));
t92 = t41 * t43;
t95 = t31 * t43;
t128 = -t92 * t132 + t134 * t95;
t116 = g(1) * t46 + g(2) * t43;
t127 = t31 * mrSges(5,1) + (-mrSges(5,2) - t134) * t30;
t121 = t31 * pkin(4) + t30 * pkin(10);
t123 = -t30 * t40 + t31 * t32;
t126 = -m(6) * t121 - m(7) * t123;
t108 = m(7) * pkin(5);
t75 = t35 * mrSges(4,1) - t34 * mrSges(4,2);
t120 = t131 * t46 + t129;
t119 = t131 * t43 + t128;
t105 = pkin(4) * t30;
t106 = pkin(3) * t34;
t118 = -m(7) * (t62 - t106) - m(6) * (-t105 - t106) + t133;
t115 = m(5) + m(6) + m(7);
t114 = -t127 + (t124 * t41 - t135) * t31;
t112 = t108 + t136;
t47 = -pkin(8) - pkin(7);
t111 = -m(3) * pkin(7) + m(4) * t47 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t110 = -t75 + t114;
t45 = cos(qJ(2));
t37 = t45 * pkin(2);
t42 = sin(qJ(2));
t68 = t45 * mrSges(3,1) - t42 * mrSges(3,2);
t109 = mrSges(2,1) + m(4) * (t37 + pkin(1)) + t75 + m(3) * pkin(1) + t68 + t127;
t107 = pkin(2) * t42;
t29 = pkin(3) * t35;
t90 = t43 * t44;
t88 = t44 * t46;
t84 = t29 + t37;
t79 = t29 + t121;
t21 = pkin(10) * t95;
t72 = -t43 * t105 + t21;
t22 = pkin(10) * t94;
t71 = -t46 * t105 + t22;
t70 = t29 + t123;
t3 = -t31 * t91 + t90;
t1 = t31 * t92 + t88;
t38 = -pkin(9) + t47;
t14 = -t106 - t107;
t13 = pkin(1) + t84;
t7 = t46 * t14;
t6 = t43 * t14;
t4 = t31 * t88 + t92;
t2 = -t31 * t90 + t91;
t5 = [(-t92 * t108 - t115 * (t46 * t13 - t38 * t43) - t136 * t4 - t124 * t3 + t111 * t43 + (-t109 + t126) * t46) * g(2) + (-t136 * t2 - t124 * t1 + (-t108 * t41 + t115 * t38 + t111) * t46 + (m(5) * t13 - m(6) * (-t13 - t121) - m(7) * (-t13 - t123) + t109) * t43) * g(1) (-m(6) * (t6 + t72) - m(7) * t6 + t119) * g(2) + (-m(6) * (t7 + t71) - m(7) * t7 + t120) * g(1) + (-t68 - m(4) * t37 - m(5) * t84 - m(6) * (t37 + t79) - m(7) * (t37 + t70) + t110) * g(3) + t116 * (m(4) * t107 - m(5) * t14 + mrSges(3,1) * t42 + mrSges(3,2) * t45 + t130) (-m(6) * t21 + t118 * t43 + t128) * g(2) + (-m(6) * t22 + t118 * t46 + t129) * g(1) + (-m(5) * t29 - m(6) * t79 - m(7) * t70 + t110) * g(3) + t116 * (m(5) * t106 + t130) t116 * t64 + (-m(6) * t72 + t119) * g(2) + (-m(6) * t71 + t120) * g(1) + (t114 + t126) * g(3) (t112 * t41 + t124 * t44) * g(3) * t30 + (t112 * t1 - t124 * t2) * g(2) + (-t112 * t3 + t124 * t4) * g(1) (g(3) * t31 - t116 * t30) * m(7)];
taug  = t5(:);
