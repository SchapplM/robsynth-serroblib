% Calculate Gravitation load on the joints for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:43
% EndTime: 2019-12-31 19:09:45
% DurationCPUTime: 0.44s
% Computational Cost: add. (273->81), mult. (316->93), div. (0->0), fcn. (283->10), ass. (0->48)
t70 = mrSges(5,3) + mrSges(6,3);
t27 = cos(qJ(4));
t15 = t27 * pkin(4) + pkin(3);
t21 = qJ(4) + qJ(5);
t18 = sin(t21);
t19 = cos(t21);
t25 = sin(qJ(4));
t69 = -m(5) * pkin(3) - m(6) * t15 - t27 * mrSges(5,1) - t19 * mrSges(6,1) + t25 * mrSges(5,2) + t18 * mrSges(6,2);
t29 = -pkin(8) - pkin(7);
t68 = -m(5) * pkin(7) + m(6) * t29 - t70;
t61 = m(6) * pkin(4);
t67 = t25 * t61;
t66 = mrSges(5,1) + t61;
t65 = m(4) + m(5) + m(6);
t64 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t20 = pkin(9) + qJ(3);
t16 = sin(t20);
t23 = cos(pkin(9));
t17 = cos(t20);
t38 = t17 * mrSges(4,1) - t16 * mrSges(4,2);
t62 = m(3) * pkin(1) + t23 * mrSges(3,1) - sin(pkin(9)) * mrSges(3,2) + mrSges(2,1) + t38 + t70 * t16;
t28 = cos(qJ(1));
t47 = t28 * t19;
t26 = sin(qJ(1));
t52 = t26 * t18;
t5 = t17 * t52 + t47;
t48 = t28 * t18;
t51 = t26 * t19;
t6 = -t17 * t51 + t48;
t60 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = -t17 * t48 + t51;
t8 = t17 * t47 + t52;
t59 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t55 = g(3) * t16;
t50 = t26 * t25;
t49 = t26 * t27;
t46 = t28 * t25;
t45 = t28 * t27;
t39 = t17 * pkin(3) + t16 * pkin(7);
t36 = -mrSges(6,1) * t18 - mrSges(6,2) * t19;
t35 = t17 * t15 - t16 * t29;
t11 = -t17 * t46 + t49;
t9 = t17 * t50 + t45;
t24 = -pkin(6) - qJ(2);
t14 = t23 * pkin(2) + pkin(1);
t12 = t17 * t45 + t50;
t10 = -t17 * t49 + t46;
t1 = [(-t50 * t61 - t12 * mrSges(5,1) - t8 * mrSges(6,1) - t11 * mrSges(5,2) - t7 * mrSges(6,2) - t65 * (t28 * t14 - t26 * t24) + t64 * t26 + (-m(5) * t39 - m(6) * t35 - t62) * t28) * g(2) + (-t10 * mrSges(5,1) - t6 * mrSges(6,1) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + (t65 * t24 + t64 - t67) * t28 + (m(4) * t14 - m(5) * (-t14 - t39) - m(6) * (-t14 - t35) + t62) * t26) * g(1), (-g(1) * t26 + g(2) * t28) * (m(3) + t65), (t68 * t16 + t69 * t17 - t38) * g(3) + (g(1) * t28 + g(2) * t26) * ((mrSges(4,2) + t68) * t17 + (mrSges(4,1) - t69) * t16), (mrSges(5,1) * t25 + mrSges(5,2) * t27 - t36 + t67) * t55 + (-t10 * mrSges(5,2) + t66 * t9 - t60) * g(2) + (t12 * mrSges(5,2) - t66 * t11 - t59) * g(1), -g(1) * t59 - g(2) * t60 - t36 * t55];
taug = t1(:);
