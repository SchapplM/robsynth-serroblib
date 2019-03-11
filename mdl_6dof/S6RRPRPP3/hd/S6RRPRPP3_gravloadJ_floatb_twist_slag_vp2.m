% Calculate Gravitation load on the joints for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:20
% EndTime: 2019-03-09 09:54:22
% DurationCPUTime: 0.92s
% Computational Cost: add. (451->111), mult. (641->120), div. (0->0), fcn. (610->8), ass. (0->56)
t95 = mrSges(5,1) - mrSges(6,2);
t94 = -m(4) * qJ(3) - mrSges(7,1) - mrSges(4,3);
t83 = -mrSges(6,3) - mrSges(7,2);
t79 = mrSges(5,2) + t83;
t28 = -pkin(8) - qJ(3);
t93 = -m(7) * (pkin(5) - t28) + t94;
t29 = sin(qJ(2));
t30 = sin(qJ(1));
t32 = cos(qJ(1));
t88 = g(1) * t32 + g(2) * t30;
t92 = t29 * t88;
t91 = mrSges(5,3) + mrSges(6,1);
t25 = pkin(9) + qJ(4);
t20 = sin(t25);
t21 = cos(t25);
t90 = -t79 * t20 + t95 * t21;
t31 = cos(qJ(2));
t26 = sin(pkin(9));
t27 = cos(pkin(9));
t39 = m(4) * pkin(2) + t27 * mrSges(4,1) - t26 * mrSges(4,2);
t84 = t39 * t31;
t89 = t93 * t29 - t84;
t87 = m(3) + m(4);
t86 = m(6) + m(7);
t62 = t32 * t31;
t7 = t20 * t62 - t30 * t21;
t64 = t30 * t20;
t8 = t21 * t62 + t64;
t85 = t8 * pkin(4) + t7 * qJ(5);
t59 = qJ(5) * t20;
t63 = t31 * t21;
t82 = pkin(4) * t63 + t31 * t59;
t81 = -m(5) - t86;
t80 = -mrSges(4,1) * t26 - mrSges(4,2) * t27 + mrSges(2,2) - mrSges(3,3);
t47 = t31 * mrSges(3,1) - t29 * mrSges(3,2);
t78 = t91 * t29 + t47;
t77 = m(7) * qJ(6) + mrSges(7,3);
t75 = t77 + t95;
t74 = pkin(3) * t26;
t71 = g(3) * t29;
t69 = t21 * t29;
t68 = t28 * t29;
t65 = t29 * t32;
t19 = pkin(3) * t27 + pkin(2);
t12 = t31 * t19;
t60 = t32 * pkin(1) + t30 * pkin(7);
t23 = t32 * pkin(7);
t58 = t30 * t68 + t32 * t74 + t23;
t52 = t12 - t68;
t51 = t19 * t62 + t30 * t74 + t60;
t50 = -t19 - t59;
t48 = m(7) * (-pkin(4) - qJ(6)) - mrSges(7,3);
t40 = -t28 * t65 + t51;
t6 = -t32 * t20 + t30 * t63;
t5 = t21 * t32 + t31 * t64;
t1 = [(-m(5) * t40 - m(6) * (t40 + t85) - m(7) * (t51 + t85) - t91 * t65 - t87 * t60 - t75 * t8 + t79 * t7 + t80 * t30 + (-mrSges(2,1) - t47 + t89) * t32) * g(2) + (-m(5) * t58 - t86 * (-t6 * pkin(4) - qJ(5) * t5 + t58) - t87 * t23 + t75 * t6 - t79 * t5 + t80 * t32 + (t84 + mrSges(2,1) + t87 * pkin(1) + t81 * (-pkin(1) - t12) + (m(7) * pkin(5) - t94) * t29 + t78) * t30) * g(1) (-m(5) * t52 - m(6) * (t52 + t82) - m(7) * (t12 + t82) - t78 + t89) * g(3) + ((-t77 * t21 - t90) * g(3) + t88 * (mrSges(3,2) - t91 + (m(5) + m(6)) * t28 + t93)) * t31 + (mrSges(3,1) + t39 + m(5) * t19 - m(6) * (-pkin(4) * t21 + t50) - m(7) * t50 - t48 * t21 + t90) * t92 (t31 * g(3) - t92) * (m(4) - t81) -(-mrSges(5,1) * t20 - mrSges(5,2) * t21) * t71 + ((t83 * t21 + (m(6) * pkin(4) - mrSges(6,2) - t48) * t20) * t29 - t86 * qJ(5) * t69) * g(3) + (-t86 * (-t5 * pkin(4) + qJ(5) * t6) + t79 * t6 + t75 * t5) * g(2) + (-t86 * (-t7 * pkin(4) + qJ(5) * t8) + t79 * t8 + t75 * t7) * g(1), t86 * (-g(1) * t7 - g(2) * t5 - t20 * t71) (-g(1) * t8 - g(2) * t6 - g(3) * t69) * m(7)];
taug  = t1(:);
