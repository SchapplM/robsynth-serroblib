% Calculate Gravitation load on the joints for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:48
% EndTime: 2019-12-05 17:11:49
% DurationCPUTime: 0.41s
% Computational Cost: add. (260->74), mult. (326->97), div. (0->0), fcn. (297->10), ass. (0->39)
t21 = qJ(3) + qJ(4);
t17 = cos(t21);
t26 = cos(qJ(3));
t19 = t26 * pkin(3);
t11 = pkin(4) * t17 + t19;
t18 = qJ(5) + t21;
t13 = sin(t18);
t14 = cos(t18);
t16 = sin(t21);
t24 = sin(qJ(3));
t69 = mrSges(3,1) + m(5) * (t19 + pkin(2)) + t17 * mrSges(5,1) - t16 * mrSges(5,2) + m(6) * (pkin(2) + t11) + t14 * mrSges(6,1) - t13 * mrSges(6,2) + m(4) * pkin(2) + t26 * mrSges(4,1) - t24 * mrSges(4,2);
t28 = -pkin(7) - pkin(6);
t68 = mrSges(3,2) + m(6) * (-pkin(8) + t28) - mrSges(6,3) + m(5) * t28 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t67 = -m(5) * pkin(3) - mrSges(4,1);
t39 = -mrSges(6,1) * t13 - mrSges(6,2) * t14;
t66 = mrSges(5,1) * t16 + mrSges(5,2) * t17 - t39;
t22 = sin(pkin(9));
t23 = cos(pkin(9));
t27 = cos(qJ(2));
t50 = t23 * t27;
t38 = -t16 * t50 + t22 * t17;
t45 = t27 * t17;
t46 = t27 * t14;
t47 = t27 * t13;
t59 = (t22 * t14 - t23 * t47) * mrSges(6,1) + (-t22 * t13 - t23 * t46) * mrSges(6,2);
t65 = -t38 * mrSges(5,1) - (-t22 * t16 - t23 * t45) * mrSges(5,2) - t59;
t51 = t22 * t27;
t37 = -t16 * t51 - t23 * t17;
t60 = (-t23 * t14 - t22 * t47) * mrSges(6,1) + (t23 * t13 - t22 * t46) * mrSges(6,2);
t64 = -t37 * mrSges(5,1) - (t23 * t16 - t22 * t45) * mrSges(5,2) - t60;
t61 = m(6) * pkin(4);
t56 = pkin(4) * t16;
t25 = sin(qJ(2));
t53 = g(3) * t25;
t52 = t24 * pkin(3);
t49 = t24 * t27;
t48 = t26 * t27;
t10 = -t52 - t56;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), (t68 * t25 - t69 * t27) * g(3) + (g(1) * t23 + g(2) * t22) * (t69 * t25 + t68 * t27), (m(5) * t52 - m(6) * t10 + mrSges(4,1) * t24 + mrSges(4,2) * t26 + t66) * t53 + (-(-t22 * t48 + t23 * t24) * mrSges(4,2) - m(6) * (t10 * t51 - t23 * t11) + t67 * (-t22 * t49 - t23 * t26) + t64) * g(2) + (-(-t22 * t24 - t23 * t48) * mrSges(4,2) - m(6) * (t10 * t50 + t22 * t11) + t67 * (t22 * t26 - t23 * t49) + t65) * g(1), (m(6) * t56 + t66) * t53 + (-t37 * t61 + t64) * g(2) + (-t38 * t61 + t65) * g(1), -g(1) * t59 - g(2) * t60 - t39 * t53];
taug = t1(:);
