% Calculate Gravitation load on the joints for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:29
% EndTime: 2019-03-09 05:34:31
% DurationCPUTime: 0.89s
% Computational Cost: add. (288->122), mult. (653->144), div. (0->0), fcn. (662->8), ass. (0->64)
t92 = mrSges(5,1) + mrSges(6,1);
t91 = mrSges(5,2) - mrSges(6,3);
t90 = mrSges(5,3) + mrSges(6,2);
t29 = sin(qJ(4));
t30 = sin(qJ(3));
t35 = cos(qJ(1));
t68 = t35 * t30;
t31 = sin(qJ(1));
t33 = cos(qJ(4));
t73 = t31 * t33;
t12 = t29 * t68 + t73;
t67 = t35 * t33;
t75 = t31 * t29;
t13 = t30 * t67 - t75;
t28 = sin(qJ(6));
t32 = cos(qJ(6));
t46 = t12 * t28 + t13 * t32;
t95 = -t12 * t32 + t13 * t28;
t99 = t95 * mrSges(7,1) + t46 * mrSges(7,2);
t34 = cos(qJ(3));
t98 = t90 * t34;
t44 = t28 * t29 + t32 * t33;
t45 = t28 * t33 - t29 * t32;
t97 = t44 * mrSges(7,1) - t45 * mrSges(7,2) - t91 * t29 + t92 * t33;
t96 = m(6) + m(7);
t94 = t13 * pkin(4) + t12 * qJ(5);
t89 = -m(5) - t96;
t88 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t81 = pkin(4) * t33;
t43 = -qJ(5) * t29 - pkin(3) - t81;
t80 = pkin(5) * t33;
t87 = m(5) * pkin(3) - m(6) * t43 - m(7) * (t43 - t80) + t97;
t86 = m(7) * pkin(9) + mrSges(7,3);
t49 = t30 * mrSges(4,1) + t34 * mrSges(4,2);
t53 = m(7) * (-pkin(8) + pkin(9)) + mrSges(7,3);
t85 = -t53 * t34 + mrSges(2,2) - mrSges(3,3) - t49;
t84 = m(7) * pkin(5) + t92;
t82 = -pkin(1) - pkin(7);
t79 = g(1) * t31;
t78 = g(2) * t35;
t76 = t29 * t34;
t74 = t31 * t30;
t72 = t31 * t34;
t69 = t34 * t35;
t66 = pkin(3) * t72 + pkin(8) * t74;
t65 = t35 * pkin(1) + t31 * qJ(2);
t64 = qJ(5) * t34;
t62 = pkin(8) * t72;
t61 = t35 * pkin(7) + t65;
t24 = t35 * qJ(2);
t60 = -t31 * pkin(1) + t24;
t55 = pkin(3) * t74 + t61;
t10 = t29 * t74 - t67;
t11 = t29 * t35 + t30 * t73;
t1 = t10 * t32 - t11 * t28;
t2 = t10 * t28 + t11 * t32;
t52 = mrSges(7,1) * t1 - mrSges(7,2) * t2;
t51 = (-mrSges(7,1) * t45 - mrSges(7,2) * t44) * t34;
t50 = mrSges(4,1) * t34 - mrSges(4,2) * t30;
t40 = t11 * pkin(4) + qJ(5) * t10 + t55;
t20 = pkin(3) * t68;
t38 = -pkin(8) * t69 + t82 * t31 + t20 + t24;
t16 = t33 * t64;
t3 = [(-m(3) * t65 - m(4) * t61 - m(5) * (t55 - t62) - m(6) * (t40 - t62) - m(7) * t40 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t84 * t11 + t91 * t10 + t88 * t35 + (t85 + t98) * t31) * g(2) + (-m(3) * t60 - m(4) * t24 - m(5) * t38 - m(6) * (t38 + t94) - m(7) * (t20 + t60 + t94) - t46 * mrSges(7,1) + t95 * mrSges(7,2) + t90 * t69 - t84 * t13 + t91 * t12 + (-m(4) * t82 + m(7) * pkin(7) - t88) * t31 + t85 * t35) * g(1) (-t79 + t78) * (m(3) + m(4) - t89) -t50 * t79 + (-m(5) * t66 - t96 * (t64 * t75 + t72 * t81 + t66) + ((-m(7) * t80 - t97) * t34 + (t86 - t90) * t30) * t31) * g(1) + (t50 + t87 * t34 + (-t53 + (m(5) + m(6)) * pkin(8) + t90) * t30) * t78 + (t87 * t30 + t49 + (pkin(8) * t89 + t86) * t34 - t98) * g(3) (-m(6) * t16 - m(7) * (t16 + (-pkin(4) - pkin(5)) * t76) + t51 + (t91 * t33 + (m(6) * pkin(4) + t92) * t29) * t34) * g(3) + (-t96 * (t12 * pkin(4) - qJ(5) * t13) - t91 * t13 - t84 * t12 + t99) * g(2) + (t52 - t96 * (-t10 * pkin(4) + qJ(5) * t11) + t91 * t11 + t84 * t10) * g(1), t96 * (-g(1) * t10 + g(2) * t12 - g(3) * t76) -g(1) * t52 - g(2) * t99 - g(3) * t51];
taug  = t3(:);
