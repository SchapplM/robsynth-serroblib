% Calculate Gravitation load on the joints for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:07
% EndTime: 2019-03-09 16:35:09
% DurationCPUTime: 0.88s
% Computational Cost: add. (579->111), mult. (554->125), div. (0->0), fcn. (492->10), ass. (0->64)
t128 = -mrSges(6,1) - mrSges(7,1);
t127 = mrSges(6,3) + mrSges(7,2);
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t126 = -t43 * mrSges(7,3) + t128 * t46;
t117 = m(6) + m(7);
t124 = m(5) + t117;
t42 = qJ(2) + qJ(3);
t37 = pkin(10) + t42;
t31 = sin(t37);
t123 = t127 * t31;
t32 = cos(t37);
t38 = sin(t42);
t39 = cos(t42);
t122 = mrSges(4,1) * t38 + mrSges(5,1) * t31 + mrSges(4,2) * t39 + mrSges(5,2) * t32;
t119 = pkin(5) * t46 + qJ(6) * t43;
t111 = (-m(7) * (-pkin(4) - t119) - t126) * t31;
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t120 = g(1) * t48 + g(2) * t45;
t118 = -t39 * mrSges(4,1) - t32 * mrSges(5,1) + t38 * mrSges(4,2) + t31 * mrSges(5,2);
t86 = t43 * mrSges(6,2);
t75 = t31 * t86;
t90 = t32 * t45;
t116 = -t127 * t90 - t45 * t75;
t89 = t32 * t48;
t115 = -t127 * t89 - t48 * t75;
t113 = -pkin(4) * t32 - pkin(9) * t31;
t101 = pkin(4) * t31;
t102 = pkin(3) * t38;
t112 = m(7) * t102 - m(6) * (-t101 - t102) + t111;
t47 = cos(qJ(2));
t40 = t47 * pkin(2);
t44 = sin(qJ(2));
t65 = t47 * mrSges(3,1) - t44 * mrSges(3,2);
t108 = -mrSges(2,1) - m(3) * pkin(1) - t65 - m(4) * (t40 + pkin(1)) + t118;
t49 = -pkin(8) - pkin(7);
t107 = -m(3) * pkin(7) + m(4) * t49 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t106 = t118 - t123 + (t86 + t126) * t32;
t105 = m(7) * pkin(5) - t128;
t104 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t103 = pkin(2) * t44;
t33 = pkin(3) * t39;
t97 = g(3) * t31;
t11 = -t102 - t103;
t23 = pkin(9) * t89;
t93 = t11 * t48 + t23;
t91 = t31 * t48;
t84 = t43 * t45;
t83 = t45 * t46;
t80 = t46 * t48;
t79 = t48 * t43;
t78 = t33 + t40;
t72 = t33 - t113;
t10 = pkin(1) + t78;
t41 = -qJ(4) + t49;
t69 = t10 * t48 - t41 * t45;
t66 = t119 * t32 + t72;
t20 = pkin(9) * t90;
t4 = t32 * t80 + t84;
t3 = t32 * t79 - t83;
t2 = t32 * t83 - t79;
t1 = t32 * t84 + t80;
t5 = [(-m(5) * t69 - t127 * t91 - t117 * (pkin(4) * t89 + pkin(9) * t91 + t69) - t105 * t4 - t104 * t3 + t108 * t48 + t107 * t45) * g(2) + (t105 * t2 + t104 * t1 + (m(5) * t10 - t117 * (-t10 + t113) - t108 + t123) * t45 + (t124 * t41 + t107) * t48) * g(1) (-t117 * (t11 * t45 + t20) + (m(6) * t101 + t111) * t45 + t116) * g(2) + (-m(6) * (-pkin(4) * t91 + t93) - m(7) * t93 + t111 * t48 + t115) * g(1) + (-t65 - m(4) * t40 - m(5) * t78 - m(6) * (t40 + t72) - m(7) * (t40 + t66) + t106) * g(3) + t120 * (m(4) * t103 - m(5) * t11 + mrSges(3,1) * t44 + mrSges(3,2) * t47 + t122) (t112 * t45 - t117 * t20 + t116) * g(2) + (t112 * t48 - t117 * t23 + t115) * g(1) + (-m(5) * t33 - m(6) * t72 - m(7) * t66 + t106) * g(3) + (m(5) * t102 + t122) * t120 (-g(1) * t45 + g(2) * t48) * t124 (-t104 * t46 + t105 * t43) * t97 + (t1 * t105 - t104 * t2) * g(2) + (-t104 * t4 + t105 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - t43 * t97) * m(7)];
taug  = t5(:);
