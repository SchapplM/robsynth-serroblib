% Calculate Gravitation load on the joints for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 17:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:50:36
% EndTime: 2018-11-23 17:50:36
% DurationCPUTime: 0.57s
% Computational Cost: add. (596->111), mult. (453->113), div. (0->0), fcn. (373->12), ass. (0->62)
t38 = qJ(2) + qJ(3);
t32 = pkin(11) + t38;
t30 = qJ(5) + t32;
t24 = sin(t30);
t25 = cos(t30);
t39 = sin(qJ(6));
t92 = mrSges(7,2) * t39;
t115 = t24 * t92 + t25 * (m(7) * pkin(10) + mrSges(7,3));
t105 = t25 * pkin(5) + t24 * pkin(10);
t114 = m(7) * t105;
t27 = sin(t32);
t28 = cos(t32);
t33 = sin(t38);
t34 = cos(t38);
t113 = -t34 * mrSges(4,1) - t28 * mrSges(5,1) + t33 * mrSges(4,2) + t27 * mrSges(5,2);
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t112 = g(1) * t44 + g(2) * t41;
t111 = -t25 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t24;
t93 = mrSges(6,2) * t25;
t110 = mrSges(4,1) * t33 + mrSges(5,1) * t27 + mrSges(6,1) * t24 + mrSges(4,2) * t34 + mrSges(5,2) * t28 + t93;
t107 = m(6) + m(7);
t42 = cos(qJ(6));
t84 = t42 * mrSges(7,1);
t104 = -(t84 - t92) * t25 + t111;
t102 = t104 + t113;
t45 = -pkin(8) - pkin(7);
t37 = -qJ(4) + t45;
t100 = -m(3) * pkin(7) + m(4) * t45 + m(5) * t37 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t43 = cos(qJ(2));
t36 = t43 * pkin(2);
t40 = sin(qJ(2));
t63 = t43 * mrSges(3,1) - t40 * mrSges(3,2);
t29 = pkin(3) * t34;
t79 = t29 + t36;
t99 = mrSges(2,1) + m(5) * (pkin(1) + t79) + m(4) * (t36 + pkin(1)) + m(3) * pkin(1) + t63 - t111 - t113;
t98 = pkin(3) * t33;
t97 = pkin(5) * t24;
t94 = t40 * pkin(2);
t86 = t41 * t39;
t85 = t41 * t42;
t83 = t44 * t39;
t82 = t44 * t42;
t23 = pkin(4) * t28;
t80 = t23 + t29;
t76 = t24 * t84;
t75 = t23 + t79;
t67 = t80 + t105;
t66 = t115 * t41;
t65 = t115 * t44;
t8 = -pkin(4) * t27 - t98;
t7 = t8 - t94;
t48 = m(7) * (t7 - t97) - t76;
t47 = m(7) * (t8 - t97) - t76;
t46 = t93 + (m(7) * pkin(5) + mrSges(6,1) + t84) * t24;
t35 = -pkin(9) + t37;
t6 = pkin(1) + t75;
t5 = t25 * t82 + t86;
t4 = -t25 * t83 + t85;
t3 = -t25 * t85 + t83;
t2 = t25 * t86 + t82;
t1 = [(-t5 * mrSges(7,1) - t4 * mrSges(7,2) - t107 * (-t41 * t35 + t44 * t6) + t100 * t41 + (-t99 - t114) * t44) * g(2) + (-t3 * mrSges(7,1) - t2 * mrSges(7,2) + (t107 * t35 + t100) * t44 + (m(6) * t6 - m(7) * (-t6 - t105) + t99) * t41) * g(1), -g(1) * (t48 * t44 + t65) - g(2) * (t48 * t41 + t66) + (-t63 - m(4) * t36 - m(5) * t79 - m(6) * t75 - m(7) * (t36 + t67) + t102) * g(3) + t112 * (m(4) * t94 - m(5) * (-t94 - t98) - m(6) * t7 + mrSges(3,1) * t40 + mrSges(3,2) * t43 + t110) -g(1) * (t47 * t44 + t65) - g(2) * (t47 * t41 + t66) + (-m(5) * t29 - m(6) * t80 - m(7) * t67 + t102) * g(3) + t112 * (m(5) * t98 - m(6) * t8 + t110) (-g(1) * t41 + g(2) * t44) * (m(5) + t107) (t104 - t114) * g(3) + (t46 * t41 - t66) * g(2) + (t46 * t44 - t65) * g(1), -g(1) * (t4 * mrSges(7,1) - t5 * mrSges(7,2)) - g(2) * (-t2 * mrSges(7,1) + t3 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t39 - mrSges(7,2) * t42) * t24];
taug  = t1(:);
