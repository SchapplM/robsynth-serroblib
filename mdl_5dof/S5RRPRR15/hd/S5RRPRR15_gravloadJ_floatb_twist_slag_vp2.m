% Calculate Gravitation load on the joints for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR15_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:41:14
% DurationCPUTime: 0.63s
% Computational Cost: add. (211->94), mult. (381->105), div. (0->0), fcn. (341->8), ass. (0->51)
t84 = mrSges(3,1) - mrSges(4,2);
t83 = -mrSges(3,2) + mrSges(4,3);
t26 = qJ(4) + qJ(5);
t19 = sin(t26);
t20 = cos(t26);
t27 = sin(qJ(4));
t30 = cos(qJ(4));
t66 = pkin(4) * t27;
t82 = -m(6) * t66 - t27 * mrSges(5,1) - t19 * mrSges(6,1) - t30 * mrSges(5,2) - t20 * mrSges(6,2);
t33 = -pkin(8) - pkin(7);
t81 = -m(5) * (-pkin(2) - pkin(7)) + mrSges(5,3) - m(6) * (-pkin(2) + t33) + mrSges(6,3);
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t75 = g(1) * t32 + g(2) * t29;
t77 = -m(6) * pkin(4) - mrSges(5,1);
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t76 = t83 * t28 + t84 * t31;
t74 = m(4) + m(5) + m(6);
t72 = -t31 * mrSges(6,3) - t76;
t62 = t30 * pkin(4);
t70 = -m(5) * pkin(3) - m(6) * (pkin(3) + t62) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t53 = t32 * t20;
t60 = t29 * t19;
t5 = t28 * t53 - t60;
t54 = t32 * t19;
t59 = t29 * t20;
t6 = t28 * t54 + t59;
t68 = t5 * mrSges(6,1) - t6 * mrSges(6,2);
t7 = t28 * t59 + t54;
t8 = -t28 * t60 + t53;
t67 = t7 * mrSges(6,1) + t8 * mrSges(6,2);
t63 = g(3) * t31;
t23 = t31 * pkin(2);
t61 = mrSges(6,1) * t20;
t58 = t29 * t27;
t57 = t29 * t30;
t55 = t31 * t33;
t52 = t32 * t27;
t51 = t32 * t30;
t50 = t32 * t31;
t21 = t28 * qJ(3);
t49 = t23 + t21;
t48 = t32 * pkin(1) + t29 * pkin(6);
t45 = -pkin(1) - t21;
t9 = t28 * t51 - t58;
t11 = t28 * t57 + t52;
t13 = t31 * t19 * mrSges(6,2);
t12 = -t28 * t58 + t51;
t10 = t28 * t52 + t57;
t1 = [(-m(3) * t48 - t10 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (-m(5) * pkin(7) - mrSges(5,3)) * t50 - t74 * (pkin(2) * t50 + t32 * t21 + t48) + t70 * t29 + (-mrSges(2,1) - m(6) * (t28 * t66 - t55) + t72) * t32) * g(2) + (-t12 * mrSges(5,1) - t8 * mrSges(6,1) + t11 * mrSges(5,2) + t7 * mrSges(6,2) + (mrSges(2,1) + m(3) * pkin(1) - m(4) * (t45 - t23) - m(5) * t45 - m(6) * (-pkin(1) + (-qJ(3) - t66) * t28) + t81 * t31 + t76) * t29 + ((-m(3) - t74) * pkin(6) + t70) * t32) * g(1), (-m(4) * t49 - m(5) * (t31 * pkin(7) + t49) - t31 * mrSges(5,3) - m(6) * (t49 - t55) + t82 * t28 + t72) * g(3) + ((m(4) * pkin(2) + t81 + t84) * t28 + (-qJ(3) * t74 + t82 - t83) * t31) * t75, (-t75 * t28 + t63) * t74, -(-mrSges(5,1) * t30 + mrSges(5,2) * t27) * t63 - g(3) * (t13 + (-m(6) * t62 - t61) * t31) + (-t12 * mrSges(5,2) + t77 * t11 - t67) * g(2) + (t10 * mrSges(5,2) + t77 * t9 - t68) * g(1), -g(1) * t68 - g(2) * t67 - g(3) * (-t31 * t61 + t13)];
taug = t1(:);
