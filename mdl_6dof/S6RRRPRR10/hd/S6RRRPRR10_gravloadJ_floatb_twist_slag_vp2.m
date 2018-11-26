% Calculate Gravitation load on the joints for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2018-11-23 17:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:58:21
% EndTime: 2018-11-23 17:58:22
% DurationCPUTime: 1.35s
% Computational Cost: add. (506->153), mult. (1049->179), div. (0->0), fcn. (1126->10), ass. (0->80)
t123 = mrSges(4,1) + mrSges(5,1);
t121 = mrSges(4,2) - mrSges(5,3);
t131 = -mrSges(5,2) - mrSges(4,3);
t134 = -mrSges(6,3) - mrSges(7,3);
t113 = m(6) * pkin(9) - t134;
t133 = t113 + t131;
t44 = sin(qJ(3));
t48 = cos(qJ(3));
t50 = cos(qJ(1));
t46 = sin(qJ(1));
t49 = cos(qJ(2));
t90 = t46 * t49;
t21 = t44 * t90 + t48 * t50;
t87 = t50 * t44;
t89 = t48 * t49;
t22 = t46 * t89 - t87;
t42 = qJ(5) + qJ(6);
t35 = sin(t42);
t36 = cos(t42);
t119 = t21 * t36 - t22 * t35;
t67 = t21 * t35 + t22 * t36;
t105 = mrSges(7,1) * t119 - t67 * mrSges(7,2);
t43 = sin(qJ(5));
t100 = t21 * t43;
t47 = cos(qJ(5));
t65 = t22 * t47 + t100;
t132 = t65 * mrSges(6,2) - t105;
t97 = t43 * t44;
t61 = t47 * t48 + t97;
t96 = t43 * t48;
t62 = -t44 * t47 + t96;
t63 = t35 * t44 + t36 * t48;
t64 = t35 * t48 - t36 * t44;
t130 = t61 * mrSges(6,1) + t63 * mrSges(7,1) - t62 * mrSges(6,2) - t64 * mrSges(7,2) - t121 * t44 + t123 * t48;
t51 = -pkin(10) - pkin(9);
t129 = m(7) * t51;
t45 = sin(qJ(2));
t128 = t131 * t45;
t106 = m(7) * pkin(5);
t127 = -m(5) - m(6);
t126 = -m(6) - m(7);
t122 = mrSges(2,2) - mrSges(3,3);
t120 = t21 * t47 - t22 * t43;
t101 = (-mrSges(7,1) * t64 - mrSges(7,2) * t63) * t45;
t118 = t101 + (-mrSges(6,1) * t62 - mrSges(6,2) * t61) * t45;
t23 = -t46 * t48 + t49 * t87;
t88 = t49 * t50;
t24 = t44 * t46 + t48 * t88;
t5 = t23 * t36 - t24 * t35;
t6 = t23 * t35 + t24 * t36;
t104 = mrSges(7,1) * t5 - mrSges(7,2) * t6;
t7 = t23 * t47 - t24 * t43;
t8 = t23 * t43 + t24 * t47;
t117 = mrSges(6,1) * t7 - mrSges(6,2) * t8 + t104;
t116 = m(5) - t126;
t72 = t49 * mrSges(3,1) - mrSges(3,2) * t45;
t114 = -t129 * t45 - t72;
t112 = -m(7) * (pkin(5) * t43 + qJ(4)) + t121;
t111 = pkin(8) * (-m(4) - t116);
t102 = pkin(4) * t48;
t84 = qJ(4) * t44;
t58 = -pkin(3) * t48 - pkin(2) - t84;
t34 = pkin(5) * t47 + pkin(4);
t59 = pkin(5) * t97 + t34 * t48;
t110 = (-t129 + mrSges(3,2) + t133) * t49 + (-m(7) * (t58 - t59) - m(6) * (t58 - t102) - m(5) * t58 + m(4) * pkin(2) + mrSges(3,1) + t130) * t45;
t107 = m(6) * pkin(4) + m(7) * t34 + t123;
t103 = m(7) * t45;
t37 = t45 * pkin(8);
t39 = t49 * pkin(2);
t95 = t44 * t45;
t92 = t45 * t50;
t86 = t39 + t37;
t85 = pkin(1) * t50 + pkin(7) * t46;
t83 = -pkin(1) - t39;
t76 = pkin(2) * t88 + pkin(8) * t92 + t85;
t74 = pkin(3) * t24 + t76;
t40 = t50 * pkin(7);
t19 = t23 * pkin(3);
t17 = t21 * pkin(3);
t1 = [(-m(3) * t85 - m(4) * t76 - m(7) * t74 - t8 * mrSges(6,1) - t6 * mrSges(7,1) - t7 * mrSges(6,2) - t5 * mrSges(7,2) + t127 * (qJ(4) * t23 + t74) + (-mrSges(2,1) + t114) * t50 + t122 * t46 - t107 * t24 + t112 * t23 + t133 * t92) * g(2) + (t100 * t106 + t65 * mrSges(6,1) + t67 * mrSges(7,1) + t120 * mrSges(6,2) + t119 * mrSges(7,2) - t116 * (-pkin(3) * t22 - t21 * qJ(4) + t40) + t122 * t50 + (-m(3) - m(4)) * t40 + t107 * t22 - t121 * t21 + (m(3) * pkin(1) + mrSges(2,1) + t72 + t126 * t83 + (-m(4) - m(5)) * (t83 - t37) + (-m(6) * (-pkin(8) + pkin(9)) - m(7) * (-pkin(8) - t51) + t134) * t45 - t128) * t46) * g(1) (-m(4) * t86 - t116 * (pkin(3) * t89 + t49 * t84 + t86) + t113 * t45 + (-m(6) * t102 - m(7) * t59 - t130) * t49 + t114 + t128) * g(3) + (t110 * t46 + t111 * t90) * g(2) + (t110 * t50 + t111 * t88) * g(1) (-m(6) * (-pkin(3) - pkin(4)) * t95 - (pkin(5) * t96 + (-pkin(3) - t34) * t44) * t103 + t118 + ((m(5) * pkin(3) + t123) * t44 + (-qJ(4) * t116 + t121) * t48) * t45) * g(3) + (m(7) * t17 + mrSges(6,1) * t120 + t127 * (qJ(4) * t22 - t17) + t112 * t22 + t107 * t21 - t132) * g(2) + (m(7) * t19 + t127 * (qJ(4) * t24 - t19) + t112 * t24 + t107 * t23 + t117) * g(1), t116 * (-g(1) * t23 - g(2) * t21 - g(3) * t95) (pkin(5) * t103 * t62 - t118) * g(3) + ((-mrSges(6,1) - t106) * t120 + t132) * g(2) + (-t106 * t7 - t117) * g(1), -g(1) * t104 - g(2) * t105 - g(3) * t101];
taug  = t1(:);
