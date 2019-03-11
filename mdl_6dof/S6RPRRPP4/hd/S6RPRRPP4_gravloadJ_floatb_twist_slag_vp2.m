% Calculate Gravitation load on the joints for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:46
% EndTime: 2019-03-09 04:38:49
% DurationCPUTime: 0.80s
% Computational Cost: add. (459->101), mult. (500->114), div. (0->0), fcn. (457->10), ass. (0->55)
t90 = mrSges(6,1) + mrSges(7,1);
t89 = -mrSges(6,2) + mrSges(7,3);
t75 = m(6) + m(7);
t88 = mrSges(6,3) + mrSges(7,2);
t82 = pkin(4) * t75 + mrSges(5,1);
t25 = qJ(4) + pkin(10);
t21 = sin(t25);
t23 = cos(t25);
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t86 = m(5) * pkin(3) + t32 * mrSges(5,1) - t30 * mrSges(5,2) + t89 * t21 + t90 * t23;
t24 = pkin(9) + qJ(3);
t22 = cos(t24);
t31 = sin(qJ(1));
t59 = t31 * t32;
t33 = cos(qJ(1));
t63 = t30 * t33;
t7 = -t22 * t63 + t59;
t20 = sin(t24);
t84 = t88 * t20;
t83 = g(1) * t33 + g(2) * t31;
t81 = -m(4) - m(5);
t80 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t27 = cos(pkin(9));
t46 = t22 * mrSges(4,1) - t20 * mrSges(4,2);
t79 = mrSges(2,1) + m(3) * pkin(1) + t27 * mrSges(3,1) - sin(pkin(9)) * mrSges(3,2) + t46 + t20 * mrSges(5,3);
t77 = m(7) * pkin(5) + t90;
t76 = m(7) * qJ(6) + t89;
t71 = g(3) * t20;
t28 = -qJ(5) - pkin(8);
t67 = t20 * t28;
t66 = t20 * t33;
t19 = pkin(4) * t32 + pkin(3);
t10 = t22 * t19;
t65 = t23 * t33;
t29 = -pkin(7) - qJ(2);
t64 = t29 * t33;
t62 = t31 * t21;
t61 = t31 * t23;
t60 = t31 * t30;
t58 = t32 * t33;
t57 = t33 * t21;
t55 = m(5) * pkin(8) + mrSges(5,3);
t18 = pkin(2) * t27 + pkin(1);
t52 = t33 * t18 - t31 * t29;
t48 = pkin(3) * t22 + pkin(8) * t20;
t42 = pkin(5) * t23 + qJ(6) * t21;
t5 = t22 * t60 + t58;
t8 = t22 * t58 + t60;
t6 = -t22 * t59 + t63;
t4 = t22 * t65 + t62;
t3 = t22 * t57 - t61;
t2 = t22 * t61 - t57;
t1 = t22 * t62 + t65;
t9 = [(-t8 * mrSges(5,1) - t7 * mrSges(5,2) - t88 * t66 + t81 * t52 - t75 * (pkin(4) * t60 + t33 * t10 - t28 * t66 + t52) - t77 * t4 - t76 * t3 + t80 * t31 + (-m(5) * t48 - t79) * t33) * g(2) + (m(5) * t64 - t6 * mrSges(5,1) - t5 * mrSges(5,2) - t75 * (pkin(4) * t63 + t31 * t67 - t64) + t77 * t2 + t76 * t1 + (m(4) * t29 + t80) * t33 + (m(4) * t18 - m(5) * (-t18 - t48) - t75 * (-t18 - t10) + t79 + t84) * t31) * g(1) (-g(1) * t31 + g(2) * t33) * (m(3) + t75 - t81) (-t46 - t75 * (t10 - t67) - t84) * g(3) + ((-m(7) * t42 - t86) * g(3) + t83 * (t28 * t75 + mrSges(4,2) - t55 - t88)) * t22 + (-t55 * g(3) + t83 * (mrSges(4,1) + m(6) * t19 - m(7) * (-t19 - t42) + t86)) * t20 (mrSges(5,2) * t32 + t77 * t21 - t76 * t23 + t82 * t30) * t71 + (-t6 * mrSges(5,2) + t77 * t1 - t76 * t2 + t82 * t5) * g(2) + (t8 * mrSges(5,2) + t77 * t3 - t76 * t4 - t82 * t7) * g(1) (t22 * g(3) - t20 * t83) * t75 (-g(1) * t3 - g(2) * t1 - t21 * t71) * m(7)];
taug  = t9(:);
