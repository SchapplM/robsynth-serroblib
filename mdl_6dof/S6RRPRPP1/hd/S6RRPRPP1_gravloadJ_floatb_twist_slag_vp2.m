% Calculate Gravitation load on the joints for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:09
% EndTime: 2019-03-09 09:45:12
% DurationCPUTime: 0.90s
% Computational Cost: add. (475->107), mult. (532->120), div. (0->0), fcn. (485->10), ass. (0->61)
t98 = mrSges(6,1) + mrSges(7,1);
t97 = -mrSges(6,2) + mrSges(7,3);
t96 = m(4) + m(5);
t80 = m(6) + m(7);
t26 = qJ(2) + pkin(9);
t21 = sin(t26);
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t89 = g(1) * t34 + g(2) * t31;
t95 = t21 * t89;
t88 = pkin(4) * t80 + mrSges(5,1);
t94 = mrSges(4,2) - mrSges(5,3);
t93 = -mrSges(6,3) - mrSges(7,2);
t25 = qJ(4) + pkin(10);
t20 = sin(t25);
t22 = cos(t25);
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t92 = m(5) * pkin(3) + t32 * mrSges(5,1) - t29 * mrSges(5,2) + t20 * t97 + t22 * t98;
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t49 = t33 * mrSges(3,1) - t30 * mrSges(3,2);
t90 = m(3) * pkin(1) + mrSges(2,1) + t49;
t23 = cos(t26);
t62 = t31 * t32;
t66 = t29 * t34;
t7 = -t23 * t66 + t62;
t86 = -t23 * mrSges(4,1) + t21 * t94;
t85 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t84 = -t21 * t93 - t86;
t83 = m(7) * pkin(5) + t98;
t82 = m(7) * qJ(6) + t97;
t79 = pkin(2) * t30;
t77 = pkin(8) * t21;
t74 = g(3) * t21;
t24 = t33 * pkin(2);
t27 = -qJ(5) - pkin(8);
t70 = t21 * t27;
t69 = t21 * t34;
t18 = pkin(4) * t32 + pkin(3);
t10 = t23 * t18;
t68 = t23 * t34;
t28 = -qJ(3) - pkin(7);
t67 = t28 * t34;
t65 = t31 * t20;
t64 = t31 * t22;
t63 = t31 * t29;
t61 = t32 * t34;
t60 = t34 * t20;
t19 = t24 + pkin(1);
t56 = t19 * t34 - t28 * t31;
t52 = pkin(3) * t23 + t77;
t44 = pkin(5) * t22 + qJ(6) * t20;
t5 = t23 * t63 + t61;
t8 = t23 * t61 + t63;
t6 = -t23 * t62 + t66;
t4 = t22 * t68 + t65;
t3 = t23 * t60 - t64;
t2 = t23 * t64 - t60;
t1 = t22 * t34 + t23 * t65;
t9 = [(-t8 * mrSges(5,1) - t7 * mrSges(5,2) + t93 * t69 - t96 * t56 - t80 * (pkin(4) * t63 + t18 * t68 - t27 * t69 + t56) - t83 * t4 - t82 * t3 + t85 * t31 + (-m(5) * t52 + t86 - t90) * t34) * g(2) + (m(5) * t67 - t6 * mrSges(5,1) - t5 * mrSges(5,2) - t80 * (pkin(4) * t66 + t31 * t70 - t67) + t83 * t2 + t82 * t1 + (m(4) * t28 + t85) * t34 + (m(4) * t19 - m(5) * (-t19 - t52) - t80 * (-t19 - t10) + t84 + t90) * t31) * g(1) (-t49 - m(4) * t24 - m(5) * (t24 + t77) - t80 * (t10 + t24 - t70) - t84) * g(3) + t89 * (mrSges(3,1) * t30 + mrSges(3,2) * t33 + t96 * t79 - t80 * (-t23 * t27 - t79)) + ((-m(7) * t44 - t92) * g(3) + t89 * (-m(5) * pkin(8) + t93 + t94)) * t23 + (mrSges(4,1) + m(6) * t18 - m(7) * (-t18 - t44) + t92) * t95 (-g(1) * t31 + g(2) * t34) * (t96 + t80) (mrSges(5,2) * t32 + t83 * t20 - t82 * t22 + t29 * t88) * t74 + (-t6 * mrSges(5,2) + t83 * t1 - t82 * t2 + t5 * t88) * g(2) + (t8 * mrSges(5,2) + t83 * t3 - t82 * t4 - t7 * t88) * g(1) (t23 * g(3) - t95) * t80 (-g(1) * t3 - g(2) * t1 - t20 * t74) * m(7)];
taug  = t9(:);
