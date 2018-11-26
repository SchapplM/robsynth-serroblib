% Calculate Gravitation load on the joints for
% S6RRPPRR9
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:55
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:54:43
% EndTime: 2018-11-23 16:54:44
% DurationCPUTime: 1.03s
% Computational Cost: add. (1023->117), mult. (1206->146), div. (0->0), fcn. (1134->14), ass. (0->66)
t106 = m(6) + m(7);
t103 = m(5) + t106;
t112 = t106 * (qJ(3) - pkin(9)) - mrSges(6,3) + mrSges(5,2) + mrSges(4,3) - mrSges(3,2);
t42 = sin(qJ(6));
t46 = cos(qJ(6));
t104 = m(7) * pkin(5) + t46 * mrSges(7,1) - t42 * mrSges(7,2) + mrSges(6,1);
t111 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t102 = mrSges(3,1) - mrSges(4,2) + mrSges(5,3);
t43 = sin(qJ(5));
t47 = cos(qJ(5));
t110 = -qJ(4) * t103 - t104 * t43 + t111 * t47 - t102;
t107 = -m(4) - m(5);
t101 = -mrSges(4,1) - mrSges(5,1) - mrSges(3,3);
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t82 = pkin(6) - qJ(2);
t65 = cos(t82) / 0.2e1;
t81 = pkin(6) + qJ(2);
t70 = cos(t81);
t50 = t65 + t70 / 0.2e1;
t92 = cos(qJ(1));
t18 = t44 * t45 - t92 * t50;
t68 = sin(t81);
t63 = t68 / 0.2e1;
t69 = sin(t82);
t61 = t63 - t69 / 0.2e1;
t48 = cos(qJ(2));
t88 = t45 * t48;
t19 = t92 * t61 + t88;
t40 = sin(pkin(6));
t77 = t40 * t92;
t56 = -t19 * t43 + t47 * t77;
t99 = -t18 * t46 + t42 * t56;
t98 = t18 * t42 + t46 * t56;
t95 = t42 * mrSges(7,1) + t46 * mrSges(7,2) - t112;
t89 = t40 * t45;
t85 = t92 * pkin(1) + pkin(8) * t89;
t84 = qJ(3) * t18;
t21 = t92 * t44 + t45 * t50;
t83 = qJ(3) * t21;
t76 = t92 * t48;
t22 = -t45 * t61 + t76;
t80 = t22 * pkin(2) + t85;
t75 = -pkin(1) * t45 + pkin(8) * t77;
t67 = -t19 * pkin(2) + t75;
t64 = t69 / 0.2e1;
t62 = t64 - t68 / 0.2e1;
t60 = pkin(3) * t89 + qJ(4) * t22 + t80;
t57 = t19 * t47 + t43 * t77;
t55 = pkin(4) * t89 + t60;
t53 = pkin(3) * t77 - qJ(4) * t19 + t67;
t52 = pkin(4) * t77 + t53;
t41 = cos(pkin(6));
t29 = t65 - t70 / 0.2e1;
t28 = t63 + t64;
t27 = t28 * pkin(2);
t23 = t45 * t62 + t76;
t20 = -t92 * t62 + t88;
t17 = t29 * t43 + t41 * t47;
t14 = t21 * pkin(2);
t12 = t18 * pkin(2);
t4 = t22 * t43 + t47 * t89;
t3 = -t22 * t47 + t43 * t89;
t2 = -t21 * t42 + t4 * t46;
t1 = -t21 * t46 - t4 * t42;
t5 = [(-t92 * mrSges(2,1) + t45 * mrSges(2,2) - m(3) * t85 - m(4) * (t80 + t83) - m(5) * (t60 + t83) - m(6) * t55 - t4 * mrSges(6,1) - m(7) * (pkin(5) * t4 + t55) - t2 * mrSges(7,1) - t1 * mrSges(7,2) - t111 * t3 + t101 * t89 - t102 * t22 - t112 * t21) * g(2) + (t45 * mrSges(2,1) + t92 * mrSges(2,2) - m(3) * t75 - m(4) * (t67 - t84) - m(5) * (t53 - t84) - m(6) * t52 - t56 * mrSges(6,1) - m(7) * (pkin(5) * t56 + t52) - t98 * mrSges(7,1) + t99 * mrSges(7,2) - t111 * t57 + t101 * t77 + t102 * t19 + t112 * t18) * g(1) (t107 * (qJ(3) * t29 + t27) - t106 * t27 + t95 * t29 + t110 * t28) * g(3) + (t107 * (qJ(3) * t20 - t12) + t106 * t12 + t95 * t20 - t110 * t18) * g(2) + (t107 * (qJ(3) * t23 - t14) + t106 * t14 + t95 * t23 - t110 * t21) * g(1) (-g(1) * t21 - g(2) * t18 + g(3) * t28) * (m(4) + t103) t103 * (-g(1) * t22 - g(2) * t19 - g(3) * t29) (-t111 * t17 - t104 * (t29 * t47 - t41 * t43)) * g(3) + (-t104 * t57 + t111 * t56) * g(2) + (t104 * t3 - t111 * t4) * g(1), -g(1) * (mrSges(7,1) * t1 - mrSges(7,2) * t2) - g(2) * (t99 * mrSges(7,1) + t98 * mrSges(7,2)) - g(3) * ((-t17 * t42 + t28 * t46) * mrSges(7,1) + (-t17 * t46 - t28 * t42) * mrSges(7,2))];
taug  = t5(:);
