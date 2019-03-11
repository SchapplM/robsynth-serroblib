% Calculate Gravitation load on the joints for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:19:54
% EndTime: 2019-03-09 05:19:55
% DurationCPUTime: 0.54s
% Computational Cost: add. (343->101), mult. (362->114), div. (0->0), fcn. (295->10), ass. (0->59)
t31 = qJ(3) + qJ(4);
t24 = pkin(10) + t31;
t22 = sin(t24);
t25 = sin(t31);
t90 = mrSges(5,2) * t25 + mrSges(6,2) * t22;
t23 = cos(t24);
t26 = cos(t31);
t95 = t25 * mrSges(5,1) + t22 * mrSges(6,1) + t26 * mrSges(5,2) + t23 * mrSges(6,2);
t94 = -m(4) - m(5);
t93 = m(6) + m(7);
t92 = t22 * pkin(5) - t23 * pkin(9);
t36 = cos(qJ(3));
t74 = t36 * pkin(3);
t79 = pkin(4) * t26;
t10 = t74 + t79;
t91 = m(5) * t74 + m(6) * t10;
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t89 = -g(1) * t34 + g(2) * t37;
t88 = -m(3) + t94;
t32 = sin(qJ(6));
t69 = mrSges(7,2) * t32;
t56 = t23 * t69;
t87 = (-t56 - t90) * t37;
t35 = cos(qJ(6));
t63 = t35 * mrSges(7,1);
t55 = t23 * t63;
t66 = t34 * t23;
t67 = t26 * t34;
t68 = t22 * t34;
t86 = -mrSges(5,1) * t67 - mrSges(6,1) * t66 - mrSges(7,3) * t68 - (t55 - t56) * t34;
t33 = sin(qJ(3));
t75 = t33 * pkin(3);
t85 = -m(5) * t75 - t33 * mrSges(4,1) - t36 * mrSges(4,2) - t95;
t84 = mrSges(5,1) * t26 + mrSges(6,1) * t23 + t22 * mrSges(7,3) + t55;
t83 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t19 = t23 * mrSges(7,3);
t82 = -m(7) * t92 + mrSges(2,2) - mrSges(3,3) + t19 + t85;
t38 = -pkin(8) - pkin(7);
t80 = pkin(4) * t25;
t65 = t34 * t32;
t64 = t34 * t35;
t62 = t37 * t32;
t61 = t37 * t35;
t60 = pkin(5) * t66 + pkin(9) * t68;
t59 = t37 * pkin(1) + t34 * qJ(2);
t57 = m(6) * t79;
t28 = t37 * qJ(2);
t54 = -t34 * pkin(1) + t28;
t52 = -pkin(5) * t23 - pkin(9) * t22;
t44 = t19 + (-t63 + t69) * t22;
t43 = -t80 - t92;
t30 = -qJ(5) + t38;
t9 = t75 + t80;
t4 = t22 * t61 - t65;
t3 = t22 * t62 + t64;
t2 = t22 * t64 + t62;
t1 = -t22 * t65 + t61;
t5 = [(-t2 * mrSges(7,1) - t1 * mrSges(7,2) - t93 * (-t37 * t30 + t34 * t9 + t59) + t88 * t59 + (-m(4) * pkin(7) + m(5) * t38 - t83) * t37 + t82 * t34) * g(2) + (-m(3) * t54 - t4 * mrSges(7,1) + t3 * mrSges(7,2) - t93 * (t34 * t30 + t37 * t9 + t54) + t94 * t28 + (-m(4) * (-pkin(1) - pkin(7)) - m(5) * (-pkin(1) + t38) + t83) * t34 + t82 * t37) * g(1), t89 * (-t88 + t93) t89 * (mrSges(4,1) * t36 - mrSges(4,2) * t33) + ((-m(7) * (-t10 + t52) + t84 + t91) * t37 + t87) * g(2) + (-m(7) * t60 + (-m(7) * t10 + t90 - t91) * t34 + t86) * g(1) + (m(6) * t9 - m(7) * (t43 - t75) - t44 - t85) * g(3) (m(6) * t80 - m(7) * t43 - t44 + t95) * g(3) + ((t57 - m(7) * (t52 - t79) + t84) * t37 + t87) * g(2) + (-m(7) * (pkin(4) * t67 + t60) + (-t57 + t90) * t34 + t86) * g(1), t93 * (-g(1) * t37 - g(2) * t34) -g(1) * (t1 * mrSges(7,1) - t2 * mrSges(7,2)) - g(2) * (t3 * mrSges(7,1) + t4 * mrSges(7,2)) - g(3) * (-mrSges(7,1) * t32 - mrSges(7,2) * t35) * t23];
taug  = t5(:);
