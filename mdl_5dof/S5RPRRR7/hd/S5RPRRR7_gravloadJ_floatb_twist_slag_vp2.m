% Calculate Gravitation load on the joints for
% S5RPRRR7
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:12
% EndTime: 2019-12-31 19:03:13
% DurationCPUTime: 0.44s
% Computational Cost: add. (282->80), mult. (296->99), div. (0->0), fcn. (267->10), ass. (0->46)
t66 = mrSges(5,3) + mrSges(6,3);
t27 = cos(qJ(4));
t16 = t27 * pkin(4) + pkin(3);
t23 = qJ(4) + qJ(5);
t19 = sin(t23);
t20 = cos(t23);
t24 = sin(qJ(4));
t65 = -m(5) * pkin(3) - m(6) * t16 - t27 * mrSges(5,1) - t20 * mrSges(6,1) + t24 * mrSges(5,2) + t19 * mrSges(6,2);
t30 = -pkin(8) - pkin(7);
t64 = -m(5) * pkin(7) + m(6) * t30 - t66;
t58 = m(6) * pkin(4);
t63 = mrSges(5,1) + t58;
t62 = mrSges(3,2) - mrSges(4,3);
t61 = -m(4) - m(5) - m(6);
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t38 = t28 * mrSges(4,1) - t25 * mrSges(4,2);
t59 = t66 * t25 + mrSges(3,1) + t38;
t22 = qJ(1) + pkin(9);
t17 = sin(t22);
t18 = cos(t22);
t46 = t28 * t19;
t5 = t17 * t46 + t18 * t20;
t45 = t28 * t20;
t6 = -t17 * t45 + t18 * t19;
t57 = -t5 * mrSges(6,1) + t6 * mrSges(6,2);
t7 = t17 * t20 - t18 * t46;
t8 = t17 * t19 + t18 * t45;
t56 = t7 * mrSges(6,1) - t8 * mrSges(6,2);
t53 = g(3) * t25;
t26 = sin(qJ(1));
t52 = t26 * pkin(1);
t29 = cos(qJ(1));
t21 = t29 * pkin(1);
t51 = t17 * t24;
t50 = t18 * t24;
t49 = t24 * t28;
t44 = t28 * t27;
t39 = t28 * pkin(3) + t25 * pkin(7);
t36 = -mrSges(6,1) * t19 - mrSges(6,2) * t20;
t35 = t28 * t16 - t25 * t30;
t11 = t17 * t27 - t18 * t49;
t9 = t17 * t49 + t18 * t27;
t12 = t18 * t44 + t51;
t10 = -t17 * t44 + t50;
t1 = [(-t51 * t58 - m(3) * t21 - t29 * mrSges(2,1) - t12 * mrSges(5,1) - t8 * mrSges(6,1) + t26 * mrSges(2,2) - t11 * mrSges(5,2) - t7 * mrSges(6,2) + t61 * (t18 * pkin(2) + t17 * pkin(6) + t21) + t62 * t17 + (-m(5) * t39 - m(6) * t35 - t59) * t18) * g(2) + (-t50 * t58 + m(3) * t52 + t26 * mrSges(2,1) - t10 * mrSges(5,1) - t6 * mrSges(6,1) + t29 * mrSges(2,2) - t9 * mrSges(5,2) - t5 * mrSges(6,2) + t61 * (t18 * pkin(6) - t52) + t62 * t18 + (m(4) * pkin(2) - m(5) * (-pkin(2) - t39) - m(6) * (-pkin(2) - t35) + t59) * t17) * g(1), (-m(3) + t61) * g(3), (t64 * t25 + t65 * t28 - t38) * g(3) + (g(1) * t18 + g(2) * t17) * ((mrSges(4,2) + t64) * t28 + (mrSges(4,1) - t65) * t25), (mrSges(5,2) * t27 + t63 * t24 - t36) * t53 + (-t10 * mrSges(5,2) + t63 * t9 - t57) * g(2) + (t12 * mrSges(5,2) - t63 * t11 - t56) * g(1), -g(1) * t56 - g(2) * t57 - t36 * t53];
taug = t1(:);
