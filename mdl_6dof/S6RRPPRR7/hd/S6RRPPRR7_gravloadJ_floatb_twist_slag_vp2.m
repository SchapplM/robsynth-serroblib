% Calculate Gravitation load on the joints for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:11
% EndTime: 2019-03-09 09:17:13
% DurationCPUTime: 0.81s
% Computational Cost: add. (434->114), mult. (1020->150), div. (0->0), fcn. (1149->10), ass. (0->60)
t83 = m(6) + m(7);
t91 = t83 * pkin(9);
t39 = sin(qJ(6));
t43 = cos(qJ(6));
t88 = m(7) * pkin(5) + t43 * mrSges(7,1) - t39 * mrSges(7,2) + mrSges(6,1);
t75 = m(5) + t83;
t55 = t39 * mrSges(7,1) + t43 * mrSges(7,2);
t62 = mrSges(5,2) - mrSges(3,1) - mrSges(4,1) - mrSges(6,3);
t90 = -t55 + t62;
t40 = sin(qJ(5));
t44 = cos(qJ(5));
t74 = -mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t89 = t40 * mrSges(6,2) - t88 * t44 + t74;
t68 = m(7) * pkin(10) + mrSges(7,3);
t60 = mrSges(6,2) - t68;
t87 = t62 - t91;
t69 = m(4) + t75;
t86 = t69 * qJ(3) - t74;
t85 = -m(6) * qJ(3) - m(7) * (pkin(10) * t40 + qJ(3)) - t40 * mrSges(7,3) + t89;
t84 = -t90 + t91;
t38 = sin(pkin(6));
t41 = sin(qJ(2));
t82 = t38 * t41;
t42 = sin(qJ(1));
t81 = t38 * t42;
t45 = cos(qJ(2));
t80 = t38 * t45;
t46 = cos(qJ(1));
t79 = t38 * t46;
t78 = pkin(2) * t80 + qJ(3) * t82;
t77 = t46 * pkin(1) + pkin(8) * t81;
t76 = cos(pkin(6));
t66 = t42 * t76;
t26 = -t41 * t66 + t45 * t46;
t71 = t26 * pkin(2) + t77;
t70 = pkin(3) * t80 + t78;
t67 = -t42 * pkin(1) + pkin(8) * t79;
t65 = t46 * t76;
t23 = t41 * t42 - t45 * t65;
t11 = t23 * pkin(2);
t24 = t41 * t65 + t42 * t45;
t64 = qJ(3) * t24 - t11;
t25 = t46 * t41 + t45 * t66;
t17 = t25 * pkin(2);
t63 = qJ(3) * t26 - t17;
t61 = t26 * pkin(3) + t71;
t59 = -t24 * pkin(2) + t67;
t58 = t25 * pkin(4) + t61;
t54 = -t24 * pkin(3) + t59;
t3 = -t23 * t40 + t44 * t79;
t4 = t23 * t44 + t40 * t79;
t48 = mrSges(2,2) + (t75 * qJ(4) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3)) * t38;
t22 = t76 * t40 + t44 * t80;
t16 = t25 * pkin(3);
t10 = t23 * pkin(3);
t8 = t25 * t44 - t40 * t81;
t7 = t25 * t40 + t44 * t81;
t2 = t26 * t39 + t43 * t8;
t1 = t26 * t43 - t39 * t8;
t5 = [(-t46 * mrSges(2,1) - m(3) * t77 - m(4) * t71 - m(5) * t61 - m(6) * t58 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t58) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t60 * t7 - t86 * t25 + t48 * t42 + t87 * t26) * g(2) + (t42 * mrSges(2,1) - m(3) * t67 - m(4) * t59 - m(5) * t54 + t60 * t3 + t88 * t4 + t86 * t23 + t48 * t46 + (t55 - t87) * t24 + t83 * (t23 * pkin(4) - t54)) * g(1) (-m(4) * t78 - m(5) * t70 - t83 * (pkin(4) * t82 + pkin(9) * t80 + t70) + (t90 * t45 + (-t68 * t40 + t89) * t41) * t38) * g(3) + (-m(4) * t64 - m(5) * (-t10 + t64) - t83 * (t24 * pkin(4) - t10 - t11) + t85 * t24 + t84 * t23) * g(2) + (-m(4) * t63 - m(5) * (-t16 + t63) - t83 * (t26 * pkin(4) - t16 - t17) + t85 * t26 + t84 * t25) * g(1) (-g(1) * t25 - g(2) * t23 + g(3) * t80) * t69 ((t42 * g(1) - t46 * g(2)) * t38 + g(3) * t76) * t75 (-t60 * t22 - t88 * (t40 * t80 - t76 * t44)) * g(3) + (-t88 * t3 + t60 * t4) * g(2) + (t60 * t8 + t88 * t7) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * ((t24 * t43 - t39 * t4) * mrSges(7,1) + (-t24 * t39 - t4 * t43) * mrSges(7,2)) - g(3) * ((t22 * t39 + t43 * t82) * mrSges(7,1) + (t22 * t43 - t39 * t82) * mrSges(7,2))];
taug  = t5(:);
