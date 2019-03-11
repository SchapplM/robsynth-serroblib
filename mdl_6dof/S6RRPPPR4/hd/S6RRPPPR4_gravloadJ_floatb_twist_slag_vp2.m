% Calculate Gravitation load on the joints for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:23
% EndTime: 2019-03-09 08:17:25
% DurationCPUTime: 0.86s
% Computational Cost: add. (265->98), mult. (603->122), div. (0->0), fcn. (578->8), ass. (0->53)
t87 = m(6) + m(7);
t62 = -m(5) - t87;
t85 = mrSges(5,2) - mrSges(6,3);
t78 = m(7) * pkin(8) - mrSges(6,2) - mrSges(5,3) + mrSges(7,3);
t79 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1);
t27 = sin(pkin(9));
t28 = cos(pkin(9));
t29 = sin(qJ(6));
t32 = cos(qJ(6));
t46 = t27 * t29 + t28 * t32;
t47 = t27 * t32 - t28 * t29;
t90 = -t47 * mrSges(7,1) + t46 * mrSges(7,2) + (t87 * qJ(5) - t85) * t28 + t79 * t27;
t89 = t62 * (-pkin(2) - qJ(4)) - t78;
t88 = -m(4) - m(5);
t86 = m(4) - t62;
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t84 = (mrSges(3,1) - mrSges(4,2)) * t33 + (-mrSges(3,2) + mrSges(4,3)) * t30;
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t83 = g(1) * t34 + g(2) * t31;
t82 = -mrSges(2,1) - t84;
t81 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t80 = (-mrSges(4,3) + t90) * t33 + (m(4) * pkin(2) - mrSges(4,2) + t89) * t30;
t77 = pkin(4) * t27;
t74 = g(3) * t33;
t23 = t33 * pkin(2);
t72 = t27 * t31;
t71 = t30 * t34;
t70 = t31 * t28;
t69 = t33 * t34;
t19 = t30 * qJ(3);
t66 = t23 + t19;
t24 = t34 * pkin(7);
t65 = t34 * pkin(3) + t24;
t64 = t34 * pkin(1) + t31 * pkin(7);
t63 = qJ(3) * t33;
t20 = t33 * qJ(4);
t61 = t33 * t77;
t60 = t20 + t66;
t56 = pkin(2) * t69 + t34 * t19 + t64;
t7 = t27 * t34 + t30 * t70;
t8 = t28 * t34 - t30 * t72;
t55 = t29 * t8 - t32 * t7;
t54 = t29 * t7 + t32 * t8;
t43 = t31 * pkin(3) + t34 * t20 + t56;
t15 = t34 * t63;
t14 = t31 * t63;
t6 = t27 * t71 + t70;
t5 = -t28 * t71 + t72;
t2 = t29 * t5 + t32 * t6;
t1 = -t29 * t6 + t32 * t5;
t3 = [(-m(3) * t64 - m(4) * t56 - m(5) * t43 - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t87 * (t6 * pkin(4) + t5 * qJ(5) + t43) + t79 * t6 + t85 * t5 + t78 * t69 + t82 * t34 + t81 * t31) * g(2) + (-m(5) * t65 - t54 * mrSges(7,1) + t55 * mrSges(7,2) - t87 * (t8 * pkin(4) + qJ(5) * t7 + t65) + t79 * t8 + t85 * t7 + (-m(3) - m(4)) * t24 + t81 * t34 + (m(3) * pkin(1) + m(4) * t23 - t86 * (-pkin(1) - t19) + t89 * t33 - t82) * t31) * g(1), t83 * (mrSges(3,1) * t30 + mrSges(3,2) * t33) + (-t87 * (t31 * t61 + t14) + t88 * t14 + t80 * t31) * g(2) + (-t87 * (t34 * t61 + t15) + t88 * t15 + t80 * t34) * g(1) + (-m(4) * t66 - m(5) * t60 - t87 * (t30 * t77 + t60) + t78 * t33 + t90 * t30 - t84) * g(3) (-t83 * t30 + t74) * t86 (t30 * g(3) + t83 * t33) * t62, t87 * (-g(1) * t5 + g(2) * t7 - t28 * t74) -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (t55 * mrSges(7,1) + t54 * mrSges(7,2)) - (t46 * mrSges(7,1) + t47 * mrSges(7,2)) * t74];
taug  = t3(:);
