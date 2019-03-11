% Calculate Gravitation load on the joints for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:28:02
% EndTime: 2019-03-10 03:28:03
% DurationCPUTime: 0.68s
% Computational Cost: add. (710->121), mult. (517->120), div. (0->0), fcn. (428->12), ass. (0->67)
t38 = qJ(2) + qJ(3);
t34 = qJ(4) + t38;
t31 = qJ(5) + t34;
t24 = sin(t31);
t25 = cos(t31);
t39 = sin(qJ(6));
t93 = mrSges(7,2) * t39;
t118 = t24 * t93 + t25 * (m(7) * pkin(11) + mrSges(7,3));
t109 = t25 * pkin(5) + t24 * pkin(11);
t117 = m(7) * t109;
t116 = -t25 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t24;
t28 = sin(t34);
t29 = cos(t34);
t94 = mrSges(6,2) * t25;
t115 = mrSges(5,1) * t28 + mrSges(6,1) * t24 + mrSges(5,2) * t29 + t94;
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t114 = g(1) * t44 + g(2) * t41;
t32 = sin(t38);
t33 = cos(t38);
t113 = mrSges(4,1) * t32 + mrSges(4,2) * t33 + t115;
t110 = m(6) + m(7);
t71 = t29 * mrSges(5,1) - t28 * mrSges(5,2);
t72 = t33 * mrSges(4,1) - t32 * mrSges(4,2);
t42 = cos(qJ(6));
t84 = t42 * mrSges(7,1);
t108 = -(t84 - t93) * t25 + t116;
t105 = -t71 + t108;
t103 = -t72 + t105;
t45 = -pkin(8) - pkin(7);
t37 = -pkin(9) + t45;
t102 = -m(3) * pkin(7) + m(4) * t45 + m(5) * t37 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t43 = cos(qJ(2));
t36 = t43 * pkin(2);
t40 = sin(qJ(2));
t65 = t43 * mrSges(3,1) - t40 * mrSges(3,2);
t27 = pkin(3) * t33;
t81 = t27 + t36;
t101 = mrSges(2,1) + m(5) * (pkin(1) + t81) + t71 + m(4) * (t36 + pkin(1)) + t72 + m(3) * pkin(1) + t65 - t116;
t100 = pkin(2) * t40;
t99 = pkin(3) * t32;
t98 = pkin(4) * t28;
t23 = pkin(4) * t29;
t97 = pkin(5) * t24;
t87 = t39 * t41;
t86 = t39 * t44;
t85 = t41 * t42;
t83 = t42 * t44;
t79 = t24 * t84;
t78 = t23 + t109;
t77 = t23 + t81;
t69 = t27 + t78;
t68 = t118 * t41;
t67 = t118 * t44;
t8 = -t98 - t99;
t7 = t8 - t100;
t49 = m(7) * (t7 - t97) - t79;
t48 = m(7) * (t8 - t97) - t79;
t47 = m(7) * (-t97 - t98) - t79;
t46 = t94 + (m(7) * pkin(5) + mrSges(6,1) + t84) * t24;
t35 = -pkin(10) + t37;
t6 = pkin(1) + t77;
t5 = t25 * t83 + t87;
t4 = -t25 * t86 + t85;
t3 = -t25 * t85 + t86;
t2 = t25 * t87 + t83;
t1 = [(-t5 * mrSges(7,1) - t4 * mrSges(7,2) - t110 * (-t35 * t41 + t44 * t6) + t102 * t41 + (-t101 - t117) * t44) * g(2) + (-t3 * mrSges(7,1) - t2 * mrSges(7,2) + (t110 * t35 + t102) * t44 + (m(6) * t6 - m(7) * (-t6 - t109) + t101) * t41) * g(1), -g(1) * (t49 * t44 + t67) - g(2) * (t49 * t41 + t68) + (-t65 - m(4) * t36 - m(5) * t81 - m(6) * t77 - m(7) * (t36 + t69) + t103) * g(3) + t114 * (m(4) * t100 - m(5) * (-t99 - t100) - m(6) * t7 + mrSges(3,1) * t40 + mrSges(3,2) * t43 + t113) -g(1) * (t48 * t44 + t67) - g(2) * (t48 * t41 + t68) + (-m(5) * t27 - m(6) * (t23 + t27) - m(7) * t69 + t103) * g(3) + t114 * (m(5) * t99 - m(6) * t8 + t113) -g(1) * (t47 * t44 + t67) - g(2) * (t47 * t41 + t68) + (-m(6) * t23 - m(7) * t78 + t105) * g(3) + (m(6) * t98 + t115) * t114 (t108 - t117) * g(3) + (t46 * t41 - t68) * g(2) + (t46 * t44 - t67) * g(1), -g(1) * (mrSges(7,1) * t4 - mrSges(7,2) * t5) - g(2) * (-mrSges(7,1) * t2 + mrSges(7,2) * t3) - g(3) * (-mrSges(7,1) * t39 - mrSges(7,2) * t42) * t24];
taug  = t1(:);
