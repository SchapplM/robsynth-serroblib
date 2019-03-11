% Calculate Gravitation load on the joints for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:53:06
% EndTime: 2019-03-09 21:53:07
% DurationCPUTime: 0.64s
% Computational Cost: add. (630->113), mult. (475->116), div. (0->0), fcn. (391->12), ass. (0->62)
t38 = qJ(2) + qJ(3);
t35 = qJ(4) + t38;
t28 = pkin(11) + t35;
t23 = sin(t28);
t24 = cos(t28);
t39 = sin(qJ(6));
t91 = mrSges(7,2) * t39;
t114 = t23 * t91 + t24 * (m(7) * pkin(10) + mrSges(7,3));
t29 = sin(t35);
t30 = cos(t35);
t113 = mrSges(5,1) * t29 + mrSges(6,1) * t23 + mrSges(5,2) * t30 + mrSges(6,2) * t24;
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t112 = g(1) * t44 + g(2) * t41;
t32 = sin(t38);
t33 = cos(t38);
t111 = mrSges(4,1) * t32 + mrSges(4,2) * t33 + t113;
t110 = -t30 * mrSges(5,1) - t24 * mrSges(6,1) + t29 * mrSges(5,2) + (mrSges(6,2) - mrSges(7,3)) * t23;
t107 = m(6) + m(7);
t66 = t24 * pkin(5) + t23 * pkin(10);
t72 = t33 * mrSges(4,1) - t32 * mrSges(4,2);
t42 = cos(qJ(6));
t92 = mrSges(7,1) * t42;
t104 = -(-t91 + t92) * t24 + t110;
t102 = -t72 + t104;
t45 = -pkin(8) - pkin(7);
t37 = -pkin(9) + t45;
t100 = -m(3) * pkin(7) + m(4) * t45 + m(5) * t37 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t43 = cos(qJ(2));
t36 = t43 * pkin(2);
t40 = sin(qJ(2));
t64 = t43 * mrSges(3,1) - t40 * mrSges(3,2);
t27 = pkin(3) * t33;
t81 = t27 + t36;
t99 = mrSges(2,1) + m(5) * (pkin(1) + t81) + m(4) * (t36 + pkin(1)) + t72 + m(3) * pkin(1) + t64 - t110;
t98 = pkin(2) * t40;
t97 = pkin(3) * t32;
t96 = pkin(4) * t29;
t25 = pkin(4) * t30;
t95 = pkin(5) * t23;
t85 = t39 * t41;
t84 = t39 * t44;
t83 = t41 * t42;
t82 = t42 * t44;
t79 = t23 * t92;
t78 = t25 + t66;
t77 = t25 + t81;
t69 = t27 + t78;
t68 = t114 * t41;
t67 = t114 * t44;
t10 = -t96 - t97;
t7 = t10 - t98;
t48 = m(7) * (t7 - t95) - t79;
t47 = m(7) * (t10 - t95) - t79;
t46 = m(7) * (-t95 - t96) - t79;
t34 = -qJ(5) + t37;
t6 = pkin(1) + t77;
t5 = t24 * t82 + t85;
t4 = -t24 * t84 + t83;
t3 = -t24 * t83 + t84;
t2 = t24 * t85 + t82;
t1 = [(-t5 * mrSges(7,1) - t4 * mrSges(7,2) - t107 * (-t41 * t34 + t44 * t6) + t100 * t41 + (-m(7) * t66 - t99) * t44) * g(2) + (-t3 * mrSges(7,1) - t2 * mrSges(7,2) + (t107 * t34 + t100) * t44 + (m(6) * t6 - m(7) * (-t6 - t66) + t99) * t41) * g(1), -g(1) * (t48 * t44 + t67) - g(2) * (t48 * t41 + t68) + (-t64 - m(4) * t36 - m(5) * t81 - m(6) * t77 - m(7) * (t36 + t69) + t102) * g(3) + t112 * (m(4) * t98 - m(5) * (-t97 - t98) - m(6) * t7 + mrSges(3,1) * t40 + mrSges(3,2) * t43 + t111) -g(1) * (t47 * t44 + t67) - g(2) * (t47 * t41 + t68) + (-m(5) * t27 - m(6) * (t25 + t27) - m(7) * t69 + t102) * g(3) + t112 * (m(5) * t97 - m(6) * t10 + t111) -g(1) * (t46 * t44 + t67) - g(2) * (t46 * t41 + t68) + (-m(6) * t25 - m(7) * t78 + t104) * g(3) + (m(6) * t96 + t113) * t112, t107 * (-g(1) * t41 + g(2) * t44) -g(1) * (mrSges(7,1) * t4 - mrSges(7,2) * t5) - g(2) * (-mrSges(7,1) * t2 + mrSges(7,2) * t3) - g(3) * (-mrSges(7,1) * t39 - mrSges(7,2) * t42) * t23];
taug  = t1(:);
