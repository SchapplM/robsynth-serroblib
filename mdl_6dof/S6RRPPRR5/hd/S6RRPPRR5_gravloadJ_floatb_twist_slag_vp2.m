% Calculate Gravitation load on the joints for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:07:19
% EndTime: 2019-03-09 09:07:21
% DurationCPUTime: 1.13s
% Computational Cost: add. (434->117), mult. (1020->154), div. (0->0), fcn. (1149->10), ass. (0->61)
t90 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t96 = m(6) + m(7);
t100 = t96 * (qJ(3) - pkin(9)) + t90;
t39 = sin(qJ(6));
t43 = cos(qJ(6));
t94 = m(7) * pkin(5) + t43 * mrSges(7,1) - t39 * mrSges(7,2) + mrSges(6,1);
t88 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t45 = cos(qJ(2));
t46 = cos(qJ(1));
t72 = cos(pkin(6));
t62 = t46 * t72;
t23 = t41 * t42 - t45 * t62;
t24 = t41 * t62 + t42 * t45;
t40 = sin(qJ(5));
t44 = cos(qJ(5));
t38 = sin(pkin(6));
t79 = t38 * t46;
t4 = t24 * t44 + t40 * t79;
t98 = t23 * t43 + t39 * t4;
t97 = t23 * t39 - t4 * t43;
t71 = m(5) + t96;
t93 = mrSges(3,1) + mrSges(4,1) + mrSges(5,1);
t91 = -mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t87 = -t88 * t40 + t94 * t44 + t93;
t55 = -t39 * mrSges(7,1) - t43 * mrSges(7,2);
t86 = -t55 - t100;
t83 = t38 * t41;
t82 = t38 * t42;
t81 = t38 * t44;
t80 = t38 * t45;
t77 = pkin(2) * t80 + qJ(3) * t83;
t76 = t46 * pkin(1) + pkin(8) * t82;
t63 = t42 * t72;
t25 = t46 * t41 + t45 * t63;
t75 = qJ(3) * t25;
t74 = qJ(4) * t38;
t73 = t23 * qJ(3);
t26 = -t41 * t63 + t45 * t46;
t68 = t26 * pkin(2) + t76;
t67 = pkin(3) * t80 + t77;
t64 = -t42 * pkin(1) + pkin(8) * t79;
t11 = t23 * pkin(2);
t61 = qJ(3) * t24 - t11;
t17 = t25 * pkin(2);
t60 = qJ(3) * t26 - t17;
t58 = -t24 * pkin(2) + t64;
t54 = -t24 * t40 + t44 * t79;
t52 = t26 * pkin(3) - t42 * t74 + t68;
t51 = t26 * pkin(4) + t52;
t50 = -t24 * pkin(3) - t46 * t74 + t58;
t48 = -t24 * pkin(4) + t50;
t22 = -t72 * t40 + t41 * t81;
t16 = t25 * pkin(3);
t10 = t23 * pkin(3);
t8 = t26 * t44 - t40 * t82;
t7 = t26 * t40 + t42 * t81;
t2 = -t25 * t39 + t43 * t8;
t1 = -t25 * t43 - t39 * t8;
t3 = [(-t46 * mrSges(2,1) + t42 * mrSges(2,2) - m(3) * t76 - m(4) * (t68 + t75) - m(5) * (t52 + t75) - m(6) * t51 - t8 * mrSges(6,1) - m(7) * (pkin(5) * t8 + t51) - t2 * mrSges(7,1) - t1 * mrSges(7,2) + t88 * t7 + t91 * t82 - t93 * t26 - t100 * t25) * g(2) + (t42 * mrSges(2,1) + t46 * mrSges(2,2) - m(3) * t64 - m(4) * (t58 - t73) - m(5) * (t50 - t73) - m(6) * t48 + t4 * mrSges(6,1) - m(7) * (-pkin(5) * t4 + t48) - t97 * mrSges(7,1) - t98 * mrSges(7,2) + t88 * t54 + t91 * t79 + t93 * t24 + t100 * t23) * g(1) (-m(4) * t77 - m(5) * t67 - t96 * (pkin(4) * t80 + t67) + (-t87 * t45 + (t96 * pkin(9) - t55 - t90) * t41) * t38) * g(3) + (-m(4) * t61 - m(5) * (-t10 + t61) - t96 * (-t23 * pkin(4) - t10 - t11) + t86 * t24 + t87 * t23) * g(2) + (-m(4) * t60 - m(5) * (-t16 + t60) - t96 * (-t25 * pkin(4) - t16 - t17) + t86 * t26 + t87 * t25) * g(1) (-g(1) * t25 - g(2) * t23 + g(3) * t80) * (m(4) + t71) ((g(1) * t42 - g(2) * t46) * t38 + g(3) * t72) * t71 (t88 * t22 - t94 * (-t40 * t83 - t72 * t44)) * g(3) + (t88 * t4 - t94 * t54) * g(2) + (t94 * t7 + t88 * t8) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (-t98 * mrSges(7,1) + t97 * mrSges(7,2)) - g(3) * ((-t22 * t39 + t43 * t80) * mrSges(7,1) + (-t22 * t43 - t39 * t80) * mrSges(7,2))];
taug  = t3(:);
