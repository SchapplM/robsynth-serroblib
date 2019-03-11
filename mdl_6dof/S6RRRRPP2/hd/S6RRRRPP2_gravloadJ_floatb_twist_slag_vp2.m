% Calculate Gravitation load on the joints for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:49:07
% EndTime: 2019-03-09 20:49:11
% DurationCPUTime: 1.15s
% Computational Cost: add. (560->134), mult. (735->150), div. (0->0), fcn. (689->8), ass. (0->64)
t121 = -mrSges(6,3) - mrSges(7,2);
t129 = mrSges(5,3) + mrSges(6,2);
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t128 = -t44 * mrSges(6,1) + t121 * t41;
t122 = m(6) + m(7);
t127 = -mrSges(5,1) - mrSges(7,1);
t40 = qJ(2) + qJ(3);
t37 = sin(t40);
t125 = t129 * t37;
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t118 = g(1) * t46 + g(2) * t43;
t67 = m(7) * (-pkin(4) - pkin(5)) - mrSges(7,1);
t83 = qJ(5) * t41;
t71 = -pkin(3) - t83;
t99 = t37 * t44;
t123 = mrSges(5,1) * t99 + (-m(7) * t71 - t44 * t67 - m(6) * (-pkin(4) * t44 + t71) - t128) * t37;
t38 = cos(t40);
t119 = t38 * mrSges(4,1) - t37 * mrSges(4,2);
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t66 = t45 * mrSges(3,1) - t42 * mrSges(3,2);
t117 = m(3) * pkin(1) + mrSges(2,1) + t119 + t66;
t116 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t115 = mrSges(5,2) + t121;
t105 = pkin(3) * t37;
t106 = pkin(2) * t42;
t114 = -m(7) * (-qJ(6) * t38 - t106) + t38 * mrSges(7,3) + m(6) * t106 - m(5) * (-t105 - t106) + t123;
t113 = -(-m(7) * qJ(6) - mrSges(7,3)) * t38 + t123;
t95 = t38 * t46;
t27 = pkin(9) * t95;
t94 = t41 * mrSges(5,2);
t81 = t37 * t94;
t112 = -t122 * t27 - t129 * t95 - t46 * t81;
t96 = t38 * t44;
t111 = t37 * mrSges(7,3) - t119 - t125 + t127 * t96 + (t94 + t128) * t38;
t110 = m(7) * pkin(5) + mrSges(6,1) - t127;
t109 = (-t81 + (-t129 + (-m(5) - t122) * pkin(9)) * t38) * t43;
t102 = g(3) * t37;
t32 = t37 * pkin(9);
t33 = t38 * pkin(3);
t39 = t45 * pkin(2);
t98 = t37 * t46;
t91 = t43 * t41;
t90 = t43 * t44;
t88 = t44 * t46;
t87 = t46 * t41;
t47 = -pkin(8) - pkin(7);
t86 = t46 * t47;
t84 = t33 + t32;
t82 = qJ(6) * t37;
t36 = t39 + pkin(1);
t75 = -t36 - t33;
t70 = t36 * t46 - t43 * t47;
t69 = pkin(4) * t96 + t38 * t83 + t84;
t63 = mrSges(4,1) * t37 + t38 * mrSges(4,2);
t61 = pkin(3) * t95 + pkin(9) * t98 + t70;
t57 = pkin(5) * t96 + t69 - t82;
t8 = t38 * t88 + t91;
t7 = t38 * t87 - t90;
t6 = t38 * t90 - t87;
t5 = t38 * t91 + t88;
t1 = [(-m(4) * t70 - m(5) * t61 - t122 * (pkin(4) * t8 + qJ(5) * t7 + t61) + (mrSges(7,3) - t129) * t98 - t110 * t8 + t115 * t7 + (m(7) * t82 - t117) * t46 + t116 * t43) * g(2) + (m(5) * t86 - t122 * (-pkin(4) * t6 - t5 * qJ(5) - t86) + t110 * t6 - t115 * t5 + (m(4) * t47 + t116) * t46 + (m(4) * t36 - m(7) * t75 - (m(7) * (-pkin(9) + qJ(6)) + mrSges(7,3)) * t37 + (-m(5) - m(6)) * (t75 - t32) + t117 + t125) * t43) * g(1) (t114 * t43 + t109) * g(2) + (-m(5) * t27 + t114 * t46 + t112) * g(1) + (-t66 - m(4) * t39 - m(5) * (t39 + t84) - m(6) * (t39 + t69) - m(7) * (t39 + t57) + t111) * g(3) + t118 * (m(4) * t106 + mrSges(3,1) * t42 + mrSges(3,2) * t45 + t63) t118 * t63 + ((m(5) * t105 + t113) * t43 + t109) * g(2) + (-m(5) * (-pkin(3) * t98 + t27) + t113 * t46 + t112) * g(1) + (-m(5) * t84 - m(6) * t69 - m(7) * t57 + t111) * g(3) -(-mrSges(5,1) * t41 - mrSges(5,2) * t44) * t102 + ((t121 * t44 + (m(6) * pkin(4) + mrSges(6,1) - t67) * t41) * t37 - t122 * qJ(5) * t99) * g(3) + (-t122 * (-pkin(4) * t5 + qJ(5) * t6) + t115 * t6 + t110 * t5) * g(2) + (-t122 * (-pkin(4) * t7 + qJ(5) * t8) + t115 * t8 + t110 * t7) * g(1), t122 * (-g(1) * t7 - g(2) * t5 - t102 * t41) (-g(3) * t38 + t118 * t37) * m(7)];
taug  = t1(:);
