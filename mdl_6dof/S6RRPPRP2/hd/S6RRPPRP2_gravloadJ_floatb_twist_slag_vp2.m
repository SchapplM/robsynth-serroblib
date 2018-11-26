% Calculate Gravitation load on the joints for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2018-11-23 16:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6RRPPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:45:36
% EndTime: 2018-11-23 16:45:37
% DurationCPUTime: 0.77s
% Computational Cost: add. (343->93), mult. (444->93), div. (0->0), fcn. (382->8), ass. (0->48)
t71 = -mrSges(6,1) - mrSges(7,1);
t70 = m(7) * pkin(5) - t71;
t69 = mrSges(6,2) + mrSges(7,2);
t73 = m(5) + m(7);
t46 = m(6) + t73;
t82 = m(4) + t46;
t81 = mrSges(4,1) - mrSges(5,2);
t80 = -mrSges(4,2) + mrSges(5,3);
t20 = sin(qJ(5));
t23 = cos(qJ(5));
t79 = -t70 * t20 - t69 * t23;
t18 = -qJ(6) - pkin(8);
t78 = -m(6) * (-pkin(3) - pkin(8)) + mrSges(6,3) - m(7) * (-pkin(3) + t18) + mrSges(7,3);
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t67 = g(1) * t25 + g(2) * t22;
t17 = qJ(2) + pkin(9);
t14 = sin(t17);
t15 = cos(t17);
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t75 = -t24 * mrSges(3,1) + t21 * mrSges(3,2) - t80 * t14 - t81 * t15;
t10 = t14 * qJ(4);
t52 = t15 * t25;
t72 = pkin(3) * t52 + t25 * t10;
t64 = -m(3) * pkin(1) - mrSges(2,1) + t75;
t19 = -qJ(3) - pkin(7);
t62 = -m(6) * (pkin(4) - t19) - m(7) * (pkin(5) * t23 + pkin(4)) - mrSges(5,1) + mrSges(2,2) - mrSges(4,3) - m(3) * pkin(7) - mrSges(3,3);
t59 = pkin(5) * t20;
t56 = g(3) * t15;
t11 = t15 * pkin(3);
t16 = t24 * pkin(2);
t54 = t15 * mrSges(7,3);
t53 = t15 * t18;
t51 = t20 * t25;
t50 = t22 * t20;
t49 = t22 * t23;
t48 = t23 * t25;
t45 = t11 + t10 + t16;
t13 = t16 + pkin(1);
t9 = t25 * t13;
t41 = -t22 * t19 + t9;
t40 = -t13 - t10;
t1 = t14 * t48 - t50;
t3 = t14 * t49 + t51;
t4 = -t14 * t50 + t48;
t2 = t14 * t51 + t49;
t5 = [(-m(4) * t41 - m(6) * (pkin(8) * t52 + t72 + t9) - mrSges(6,3) * t52 - t73 * (t41 + t72) + t71 * t2 - t69 * t1 + t62 * t22 + (-m(7) * (t14 * t59 - t53) - t54 + t64) * t25) * g(2) + (t71 * t4 + t69 * t3 + ((m(4) + t73) * t19 + t62) * t25 + (m(4) * t13 - m(5) * (t40 - t11) - m(6) * t40 - m(7) * (-t13 + (-qJ(4) - t59) * t14) + t78 * t15 - t64) * t22) * g(1) (-m(4) * t16 - m(5) * t45 - m(6) * (pkin(8) * t15 + t45) - t15 * mrSges(6,3) - m(7) * (t45 - t53) - t54 + t79 * t14 + t75) * g(3) + (mrSges(3,2) * t24 + (m(5) * pkin(3) + t78 + t81) * t14 + (-qJ(4) * t46 + t79 - t80) * t15 + (t82 * pkin(2) + mrSges(3,1)) * t21) * t67 (-g(1) * t22 + g(2) * t25) * t82 (-t67 * t14 + t56) * t46 (-t69 * t20 + t70 * t23) * t56 + (-t3 * t70 - t69 * t4) * g(2) + (-t1 * t70 + t69 * t2) * g(1) (-g(3) * t14 - t67 * t15) * m(7)];
taug  = t5(:);
