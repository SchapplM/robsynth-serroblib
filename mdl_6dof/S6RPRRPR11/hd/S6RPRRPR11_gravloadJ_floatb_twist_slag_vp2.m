% Calculate Gravitation load on the joints for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:11
% EndTime: 2019-03-09 05:39:14
% DurationCPUTime: 1.25s
% Computational Cost: add. (1175->115), mult. (3058->161), div. (0->0), fcn. (3874->16), ass. (0->65)
t49 = sin(pkin(13));
t50 = cos(pkin(13));
t109 = -m(7) * (pkin(5) * t50 + pkin(4)) - m(6) * pkin(4) - t50 * mrSges(6,1) + t49 * mrSges(6,2) - mrSges(5,1);
t48 = pkin(13) + qJ(6);
t45 = sin(t48);
t46 = cos(t48);
t140 = -t46 * mrSges(7,1) + t45 * mrSges(7,2) + t109;
t103 = cos(pkin(7));
t108 = cos(qJ(1));
t106 = sin(qJ(1));
t102 = cos(pkin(12));
t104 = cos(pkin(6));
t90 = t104 * t102;
t99 = sin(pkin(12));
t75 = t106 * t99 - t108 * t90;
t100 = sin(pkin(7));
t101 = sin(pkin(6));
t86 = t101 * t100;
t139 = t75 * t103 + t108 * t86;
t107 = cos(qJ(3));
t88 = t104 * t99;
t33 = t106 * t102 + t108 * t88;
t53 = sin(qJ(3));
t18 = -t107 * t33 + t139 * t53;
t52 = sin(qJ(4));
t54 = cos(qJ(4));
t87 = t103 * t101;
t61 = -t75 * t100 + t108 * t87;
t138 = t18 * t52 - t54 * t61;
t137 = t18 * t54 + t52 * t61;
t130 = -m(5) - m(6);
t111 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(11) - qJ(5)) - mrSges(7,3);
t69 = t106 * t90 + t108 * t99;
t55 = t69 * t100 + t106 * t87;
t129 = m(6) + m(7);
t118 = m(5) + t129;
t133 = pkin(3) * t118 - t111 * t52 - t140 * t54 + mrSges(4,1);
t110 = -t49 * mrSges(6,1) - t50 * mrSges(6,2) + mrSges(4,2) - mrSges(5,3) - m(7) * (pkin(5) * t49 + pkin(10));
t131 = -t45 * mrSges(7,1) - t46 * mrSges(7,2) + t110;
t15 = t139 * t107 + t33 * t53;
t121 = t69 * t103 - t106 * t86;
t120 = t104 * t100 + t102 * t87;
t112 = t130 * pkin(10) + t131;
t92 = t101 * t106;
t105 = t108 * pkin(1) + qJ(2) * t92;
t93 = t108 * t101;
t95 = -pkin(1) * t106 + qJ(2) * t93;
t85 = t101 * t99;
t68 = -t102 * t86 + t103 * t104;
t63 = -t33 * pkin(2) + t61 * pkin(9) + t95;
t60 = t18 * pkin(3) + t63;
t34 = t102 * t108 - t106 * t88;
t59 = t34 * pkin(2) + t55 * pkin(9) + t105;
t20 = t34 * t107 - t121 * t53;
t58 = t20 * pkin(3) + t59;
t26 = t107 * t85 + t120 * t53;
t25 = -t107 * t120 + t53 * t85;
t19 = t107 * t121 + t34 * t53;
t14 = t26 * t54 + t52 * t68;
t13 = t26 * t52 - t54 * t68;
t8 = t20 * t54 + t52 * t55;
t7 = t20 * t52 - t54 * t55;
t2 = t19 * t45 + t46 * t8;
t1 = t19 * t46 - t45 * t8;
t3 = [(-m(3) * t105 - m(4) * t59 - m(7) * t58 - mrSges(2,1) * t108 - t34 * mrSges(3,1) - t20 * mrSges(4,1) - t2 * mrSges(7,1) + t106 * mrSges(2,2) + mrSges(3,2) * t69 - t1 * mrSges(7,2) - mrSges(3,3) * t92 - mrSges(4,3) * t55 + t130 * (t19 * pkin(10) + t58) + t109 * t8 + t110 * t19 + t111 * t7) * g(2) + (-m(3) * t95 - m(4) * t63 - m(7) * t60 + t106 * mrSges(2,1) + t33 * mrSges(3,1) - t18 * mrSges(4,1) + mrSges(2,2) * t108 - mrSges(3,2) * t75 - mrSges(3,3) * t93 - mrSges(4,3) * t61 + t130 * (-pkin(10) * t15 + t60) + t140 * t137 - t131 * t15 + t111 * t138) * g(1) (-g(1) * t92 + g(2) * t93 - g(3) * t104) * (m(3) + m(4) + t118) (t112 * t26 + t133 * t25) * g(3) + (-t112 * t18 + t133 * t15) * g(2) + (t112 * t20 + t133 * t19) * g(1) (t111 * t14 - t13 * t140) * g(3) + (-t111 * t137 + t138 * t140) * g(2) + (t111 * t8 - t140 * t7) * g(1), t129 * (-g(1) * t7 + g(2) * t138 - g(3) * t13) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t137 * t45 + t15 * t46) * mrSges(7,1) + (t137 * t46 - t15 * t45) * mrSges(7,2)) - g(3) * ((-t14 * t45 + t25 * t46) * mrSges(7,1) + (-t14 * t46 - t25 * t45) * mrSges(7,2))];
taug  = t3(:);
