% Calculate Gravitation load on the joints for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:20
% DurationCPUTime: 0.29s
% Computational Cost: add. (215->59), mult. (215->73), div. (0->0), fcn. (171->10), ass. (0->36)
t21 = qJ(2) + qJ(3);
t17 = pkin(9) + t21;
t14 = sin(t17);
t15 = cos(t17);
t24 = sin(qJ(5));
t46 = t24 * mrSges(6,2);
t66 = t14 * t46 + t15 * (m(6) * pkin(7) + mrSges(6,3));
t18 = sin(t21);
t19 = cos(t21);
t65 = mrSges(4,1) * t18 + mrSges(5,1) * t14 + mrSges(4,2) * t19 + mrSges(5,2) * t15;
t22 = sin(pkin(8));
t23 = cos(pkin(8));
t64 = g(1) * t23 + g(2) * t22;
t61 = m(5) + m(6);
t26 = cos(qJ(5));
t45 = t26 * mrSges(6,1);
t59 = -t19 * mrSges(4,1) + t18 * mrSges(4,2) + (-t45 + t46 - mrSges(5,1)) * t15 + (-mrSges(6,3) + mrSges(5,2)) * t14;
t25 = sin(qJ(2));
t57 = pkin(2) * t25;
t56 = pkin(3) * t18;
t16 = pkin(3) * t19;
t55 = pkin(4) * t14;
t27 = cos(qJ(2));
t20 = t27 * pkin(2);
t50 = t22 * t24;
t49 = t22 * t26;
t48 = t23 * t24;
t47 = t23 * t26;
t43 = t14 * t45;
t42 = t15 * pkin(4) + t14 * pkin(7) + t16;
t39 = t66 * t22;
t38 = t66 * t23;
t3 = -t56 - t57;
t29 = m(6) * (t3 - t55) - t43;
t28 = m(6) * (-t55 - t56) - t43;
t1 = [(-m(2) - m(3) - m(4) - t61) * g(3), -g(1) * (t29 * t23 + t38) - g(2) * (t29 * t22 + t39) + (-mrSges(3,1) * t27 + t25 * mrSges(3,2) - m(4) * t20 - m(5) * (t16 + t20) - m(6) * (t20 + t42) + t59) * g(3) + t64 * (m(4) * t57 - m(5) * t3 + mrSges(3,1) * t25 + mrSges(3,2) * t27 + t65), -g(1) * (t28 * t23 + t38) - g(2) * (t28 * t22 + t39) + (-m(5) * t16 - m(6) * t42 + t59) * g(3) + (m(5) * t56 + t65) * t64, t61 * (-g(1) * t22 + g(2) * t23), -g(1) * ((-t15 * t48 + t49) * mrSges(6,1) + (-t15 * t47 - t50) * mrSges(6,2)) - g(2) * ((-t15 * t50 - t47) * mrSges(6,1) + (-t15 * t49 + t48) * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t24 - mrSges(6,2) * t26) * t14];
taug = t1(:);
