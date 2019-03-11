% Calculate Gravitation load on the joints for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:33
% EndTime: 2019-03-09 07:00:34
% DurationCPUTime: 0.71s
% Computational Cost: add. (553->118), mult. (492->145), div. (0->0), fcn. (453->12), ass. (0->68)
t100 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t40 = qJ(4) + qJ(5);
t34 = cos(t40);
t44 = cos(qJ(4));
t36 = t44 * pkin(4);
t23 = pkin(5) * t34 + t36;
t21 = pkin(3) + t23;
t35 = qJ(6) + t40;
t28 = sin(t35);
t29 = cos(t35);
t30 = t36 + pkin(3);
t33 = sin(t40);
t41 = sin(qJ(4));
t99 = -m(5) * pkin(3) - m(6) * t30 - m(7) * t21 - mrSges(5,1) * t44 - mrSges(6,1) * t34 - mrSges(7,1) * t29 + mrSges(5,2) * t41 + mrSges(6,2) * t33 + mrSges(7,2) * t28;
t47 = -pkin(9) - pkin(8);
t39 = -pkin(10) + t47;
t98 = -m(5) * pkin(8) + m(6) * t47 + m(7) * t39 - t100;
t88 = m(6) * pkin(4);
t80 = t41 * pkin(4);
t84 = pkin(5) * t33;
t22 = t80 + t84;
t97 = m(7) * t22;
t96 = mrSges(5,1) + t88;
t56 = -mrSges(7,1) * t28 - mrSges(7,2) * t29;
t95 = mrSges(6,1) * t33 + mrSges(6,2) * t34 - t56;
t38 = qJ(1) + pkin(11);
t31 = sin(t38);
t32 = cos(t38);
t45 = cos(qJ(3));
t73 = t33 * t45;
t15 = t31 * t34 - t32 * t73;
t72 = t34 * t45;
t16 = t31 * t33 + t32 * t72;
t74 = t32 * t45;
t7 = -t28 * t74 + t29 * t31;
t8 = t28 * t31 + t29 * t74;
t85 = mrSges(7,1) * t7 - mrSges(7,2) * t8;
t94 = -mrSges(6,1) * t15 + mrSges(6,2) * t16 - t85;
t13 = t31 * t73 + t32 * t34;
t14 = -t31 * t72 + t32 * t33;
t76 = t31 * t45;
t5 = t28 * t76 + t29 * t32;
t6 = t28 * t32 - t29 * t76;
t86 = -mrSges(7,1) * t5 + mrSges(7,2) * t6;
t93 = mrSges(6,1) * t13 - mrSges(6,2) * t14 - t86;
t92 = -m(4) - m(5) - m(6) - m(7);
t90 = mrSges(3,2) - mrSges(4,3) - t97;
t42 = sin(qJ(3));
t59 = t45 * mrSges(4,1) - t42 * mrSges(4,2);
t89 = t100 * t42 + mrSges(3,1) + t59;
t87 = m(7) * pkin(5);
t81 = g(3) * t42;
t43 = sin(qJ(1));
t79 = t43 * pkin(1);
t46 = cos(qJ(1));
t37 = t46 * pkin(1);
t77 = t31 * t41;
t75 = t32 * t41;
t71 = t41 * t45;
t67 = t44 * t45;
t60 = pkin(3) * t45 + pkin(8) * t42;
t55 = t21 * t45 - t39 * t42;
t54 = t30 * t45 - t42 * t47;
t19 = t31 * t44 - t32 * t71;
t17 = t31 * t71 + t32 * t44;
t20 = t32 * t67 + t77;
t18 = -t31 * t67 + t75;
t1 = [(-t77 * t88 - m(3) * t37 - t46 * mrSges(2,1) - t20 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) + t43 * mrSges(2,2) - t19 * mrSges(5,2) - t15 * mrSges(6,2) - t7 * mrSges(7,2) + t92 * (pkin(2) * t32 + pkin(7) * t31 + t37) + t90 * t31 + (-m(5) * t60 - m(6) * t54 - m(7) * t55 - t89) * t32) * g(2) + (-t75 * t88 + m(3) * t79 + t43 * mrSges(2,1) - t18 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) + t46 * mrSges(2,2) - t17 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + t92 * (pkin(7) * t32 - t79) + t90 * t32 + (m(4) * pkin(2) - m(5) * (-pkin(2) - t60) - m(6) * (-pkin(2) - t54) - m(7) * (-pkin(2) - t55) + t89) * t31) * g(1) (-m(3) + t92) * g(3) (t42 * t98 + t45 * t99 - t59) * g(3) + (g(1) * t32 + g(2) * t31) * ((mrSges(4,2) + t98) * t45 + (mrSges(4,1) - t99) * t42) (m(6) * t80 + mrSges(5,1) * t41 + mrSges(5,2) * t44 + t95 + t97) * t81 + (-t18 * mrSges(5,2) - m(7) * (-t22 * t76 - t23 * t32) + t96 * t17 + t93) * g(2) + (t20 * mrSges(5,2) - m(7) * (-t22 * t74 + t23 * t31) - t96 * t19 + t94) * g(1) (m(7) * t84 + t95) * t81 + (t13 * t87 + t93) * g(2) + (-t15 * t87 + t94) * g(1), -g(1) * t85 - g(2) * t86 - t56 * t81];
taug  = t1(:);
