% Calculate Gravitation load on the joints for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:04
% EndTime: 2019-12-05 18:27:05
% DurationCPUTime: 0.30s
% Computational Cost: add. (266->59), mult. (234->56), div. (0->0), fcn. (177->10), ass. (0->30)
t21 = qJ(2) + pkin(9);
t15 = sin(t21);
t16 = cos(t21);
t23 = sin(qJ(2));
t17 = qJ(4) + t21;
t11 = sin(t17);
t12 = cos(t17);
t13 = qJ(5) + t17;
t8 = sin(t13);
t9 = cos(t13);
t40 = t9 * mrSges(6,1) - t8 * mrSges(6,2);
t34 = -t12 * mrSges(5,1) + t11 * mrSges(5,2) - t40;
t56 = t16 * mrSges(4,1) - t23 * mrSges(3,2) - t15 * mrSges(4,2) - t34;
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t53 = g(1) * t26 + g(2) * t24;
t25 = cos(qJ(2));
t19 = t25 * pkin(2);
t45 = pkin(3) * t16 + t19;
t7 = pkin(4) * t12;
t43 = t7 + t45;
t52 = mrSges(2,1) + m(4) * (t19 + pkin(1)) + m(5) * (pkin(1) + t45) + m(3) * pkin(1) + t25 * mrSges(3,1) + m(6) * (pkin(1) + t43) + t56;
t22 = -qJ(3) - pkin(6);
t20 = -pkin(7) + t22;
t51 = mrSges(2,2) + m(6) * (-pkin(8) + t20) - mrSges(6,3) + m(5) * t20 - mrSges(5,3) + m(4) * t22 - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3);
t41 = m(4) * pkin(2) + mrSges(3,1);
t36 = mrSges(6,1) * t8 + mrSges(6,2) * t9;
t4 = -t23 * pkin(2) - pkin(3) * t15;
t30 = mrSges(5,2) * t12 + t36;
t1 = [(t51 * t24 - t52 * t26) * g(2) + (t52 * t24 + t51 * t26) * g(1), (-m(5) * t45 - m(6) * t43 - t41 * t25 - t56) * g(3) + t53 * (-m(5) * t4 - m(6) * (-pkin(4) * t11 + t4) + mrSges(4,1) * t15 + mrSges(5,1) * t11 + mrSges(3,2) * t25 + mrSges(4,2) * t16 + t41 * t23 + t30), (-g(1) * t24 + g(2) * t26) * (m(4) + m(5) + m(6)), (-m(6) * t7 + t34) * g(3) + t53 * ((m(6) * pkin(4) + mrSges(5,1)) * t11 + t30), -g(3) * t40 + t53 * t36];
taug = t1(:);
