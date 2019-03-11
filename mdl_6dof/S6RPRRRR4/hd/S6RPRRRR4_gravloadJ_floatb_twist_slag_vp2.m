% Calculate Gravitation load on the joints for
% S6RPRRRR4
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:47
% EndTime: 2019-03-09 07:04:48
% DurationCPUTime: 0.55s
% Computational Cost: add. (552->103), mult. (404->106), div. (0->0), fcn. (332->12), ass. (0->58)
t35 = pkin(11) + qJ(3);
t31 = qJ(4) + t35;
t28 = qJ(5) + t31;
t22 = sin(t28);
t23 = cos(t28);
t39 = sin(qJ(6));
t77 = t39 * mrSges(7,2);
t102 = t22 * t77 + t23 * (m(7) * pkin(10) + mrSges(7,3));
t94 = t23 * pkin(5) + t22 * pkin(10);
t101 = m(7) * t94;
t100 = -t23 * mrSges(6,1) + (mrSges(6,2) - mrSges(7,3)) * t22;
t25 = sin(t31);
t26 = cos(t31);
t82 = mrSges(6,2) * t23;
t99 = mrSges(5,1) * t25 + mrSges(6,1) * t22 + mrSges(5,2) * t26 + t82;
t40 = sin(qJ(1));
t42 = cos(qJ(1));
t98 = g(1) * t42 + g(2) * t40;
t95 = m(6) + m(7);
t63 = t26 * mrSges(5,1) - t25 * mrSges(5,2);
t41 = cos(qJ(6));
t74 = t41 * mrSges(7,1);
t93 = -(t74 - t77) * t23 + t100;
t90 = -t63 + t93;
t38 = -pkin(7) - qJ(2);
t34 = -pkin(8) + t38;
t89 = -m(3) * qJ(2) + m(4) * t38 + m(5) * t34 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t30 = cos(t35);
t24 = pkin(3) * t30;
t37 = cos(pkin(11));
t27 = t37 * pkin(2) + pkin(1);
t10 = t24 + t27;
t29 = sin(t35);
t57 = t30 * mrSges(4,1) - t29 * mrSges(4,2);
t88 = m(4) * t27 + m(5) * t10 + mrSges(2,1) + m(3) * pkin(1) + t37 * mrSges(3,1) - sin(pkin(11)) * mrSges(3,2) + t57 + t63 - t100;
t87 = pkin(3) * t29;
t86 = pkin(4) * t25;
t21 = pkin(4) * t26;
t85 = pkin(5) * t22;
t76 = t40 * t39;
t75 = t40 * t41;
t73 = t42 * t39;
t72 = t42 * t41;
t69 = t22 * t74;
t68 = t21 + t94;
t60 = t102 * t40;
t59 = t102 * t42;
t7 = -t86 - t87;
t45 = m(7) * (t7 - t85) - t69;
t44 = m(7) * (-t85 - t86) - t69;
t43 = t82 + (m(7) * pkin(5) + mrSges(6,1) + t74) * t22;
t32 = -pkin(9) + t34;
t6 = t23 * t72 + t76;
t5 = -t23 * t73 + t75;
t4 = -t23 * t75 + t73;
t3 = t23 * t76 + t72;
t2 = t21 + t10;
t1 = [(-t6 * mrSges(7,1) - t5 * mrSges(7,2) - t95 * (t42 * t2 - t40 * t32) + t89 * t40 + (-t88 - t101) * t42) * g(2) + (-t4 * mrSges(7,1) - t3 * mrSges(7,2) + (t95 * t32 + t89) * t42 + (m(6) * t2 - m(7) * (-t2 - t94) + t88) * t40) * g(1) (-g(1) * t40 + g(2) * t42) * (m(3) + m(4) + m(5) + t95) -g(1) * (t45 * t42 + t59) - g(2) * (t45 * t40 + t60) + (-t57 - m(5) * t24 - m(6) * (t21 + t24) - m(7) * (t24 + t68) + t90) * g(3) + t98 * (m(5) * t87 - m(6) * t7 + mrSges(4,1) * t29 + mrSges(4,2) * t30 + t99) -g(1) * (t44 * t42 + t59) - g(2) * (t44 * t40 + t60) + (-m(6) * t21 - m(7) * t68 + t90) * g(3) + (m(6) * t86 + t99) * t98 (t93 - t101) * g(3) + (t43 * t40 - t60) * g(2) + (t43 * t42 - t59) * g(1), -g(1) * (t5 * mrSges(7,1) - t6 * mrSges(7,2)) - g(2) * (-t3 * mrSges(7,1) + t4 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t39 - mrSges(7,2) * t41) * t22];
taug  = t1(:);
