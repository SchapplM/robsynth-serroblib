% Calculate Gravitation load on the joints for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:31
% EndTime: 2019-03-09 07:12:32
% DurationCPUTime: 0.78s
% Computational Cost: add. (534->119), mult. (515->138), div. (0->0), fcn. (471->12), ass. (0->70)
t103 = mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t38 = qJ(4) + qJ(5);
t33 = cos(t38);
t44 = cos(qJ(4));
t35 = t44 * pkin(4);
t23 = pkin(5) * t33 + t35;
t21 = pkin(3) + t23;
t34 = qJ(6) + t38;
t27 = sin(t34);
t28 = cos(t34);
t29 = t35 + pkin(3);
t32 = sin(t38);
t42 = sin(qJ(4));
t102 = -m(5) * pkin(3) - m(6) * t29 - m(7) * t21 - mrSges(5,1) * t44 - mrSges(6,1) * t33 - mrSges(7,1) * t28 + mrSges(5,2) * t42 + mrSges(6,2) * t32 + mrSges(7,2) * t27;
t46 = -pkin(9) - pkin(8);
t37 = -pkin(10) + t46;
t101 = -m(5) * pkin(8) + m(6) * t46 + m(7) * t37 - t103;
t90 = m(6) * pkin(4);
t86 = pkin(4) * t42;
t100 = m(6) * t86;
t85 = pkin(5) * t32;
t22 = t85 + t86;
t99 = m(7) * t22;
t98 = mrSges(5,1) + t90;
t56 = -mrSges(7,1) * t27 - mrSges(7,2) * t28;
t97 = mrSges(6,1) * t32 + mrSges(6,2) * t33 - t56;
t36 = pkin(11) + qJ(3);
t31 = cos(t36);
t43 = sin(qJ(1));
t73 = t33 * t43;
t45 = cos(qJ(1));
t74 = t32 * t45;
t15 = -t31 * t74 + t73;
t72 = t33 * t45;
t75 = t32 * t43;
t16 = t31 * t72 + t75;
t76 = t31 * t45;
t7 = -t27 * t76 + t28 * t43;
t8 = t27 * t43 + t28 * t76;
t87 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t96 = -t15 * mrSges(6,1) + t16 * mrSges(6,2) - t87;
t13 = t31 * t75 + t72;
t14 = -t31 * t73 + t74;
t77 = t31 * t43;
t5 = t27 * t77 + t28 * t45;
t6 = t27 * t45 - t28 * t77;
t88 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t95 = t13 * mrSges(6,1) - t14 * mrSges(6,2) - t88;
t94 = m(4) + m(5) + m(6) + m(7);
t92 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - t99;
t30 = sin(t36);
t40 = cos(pkin(11));
t59 = t31 * mrSges(4,1) - t30 * mrSges(4,2);
t91 = m(3) * pkin(1) + t40 * mrSges(3,1) - sin(pkin(11)) * mrSges(3,2) + mrSges(2,1) + t59 + t103 * t30;
t89 = m(7) * pkin(5);
t82 = g(3) * t30;
t71 = t42 * t43;
t70 = t42 * t45;
t69 = t43 * t44;
t68 = t44 * t45;
t60 = pkin(3) * t31 + pkin(8) * t30;
t55 = t21 * t31 - t30 * t37;
t54 = t31 * t29 - t30 * t46;
t19 = -t31 * t70 + t69;
t17 = t31 * t71 + t68;
t41 = -pkin(7) - qJ(2);
t25 = pkin(2) * t40 + pkin(1);
t20 = t31 * t68 + t71;
t18 = -t31 * t69 + t70;
t1 = [(-t71 * t90 - t20 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) - t19 * mrSges(5,2) - t15 * mrSges(6,2) - t7 * mrSges(7,2) - t94 * (t45 * t25 - t43 * t41) + t92 * t43 + (-m(5) * t60 - m(6) * t54 - m(7) * t55 - t91) * t45) * g(2) + (-t18 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t17 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + (t94 * t41 - t100 + t92) * t45 + (m(4) * t25 - m(5) * (-t25 - t60) - m(6) * (-t25 - t54) - m(7) * (-t25 - t55) + t91) * t43) * g(1) (-g(1) * t43 + g(2) * t45) * (m(3) + t94) (t101 * t30 + t102 * t31 - t59) * g(3) + (g(1) * t45 + g(2) * t43) * ((mrSges(4,2) + t101) * t31 + (mrSges(4,1) - t102) * t30) (mrSges(5,1) * t42 + mrSges(5,2) * t44 + t100 + t97 + t99) * t82 + (-t18 * mrSges(5,2) - m(7) * (-t22 * t77 - t23 * t45) + t98 * t17 + t95) * g(2) + (t20 * mrSges(5,2) - m(7) * (-t22 * t76 + t23 * t43) - t98 * t19 + t96) * g(1) (m(7) * t85 + t97) * t82 + (t13 * t89 + t95) * g(2) + (-t15 * t89 + t96) * g(1), -g(1) * t87 - g(2) * t88 - t56 * t82];
taug  = t1(:);
