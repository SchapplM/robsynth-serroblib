% Calculate Gravitation load on the joints for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:35
% DurationCPUTime: 0.25s
% Computational Cost: add. (137->48), mult. (163->51), div. (0->0), fcn. (119->8), ass. (0->27)
t13 = pkin(8) + qJ(4);
t8 = qJ(5) + t13;
t4 = sin(t8);
t5 = cos(t8);
t22 = t4 * mrSges(6,1) + t5 * mrSges(6,2);
t7 = cos(t13);
t39 = -t7 * mrSges(5,2) - t22;
t25 = -m(4) - m(5) - m(6);
t38 = m(3) - t25;
t23 = m(6) * pkin(4) + mrSges(5,1);
t6 = sin(t13);
t36 = mrSges(5,2) * t6 - t23 * t7;
t35 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t14 = sin(pkin(8));
t30 = pkin(3) * t14;
t34 = -m(5) * t30 - t6 * mrSges(5,1) + mrSges(2,2) - mrSges(3,3) - m(6) * (pkin(4) * t6 + t30) - t14 * mrSges(4,1) - cos(pkin(8)) * mrSges(4,2) + t39;
t33 = mrSges(6,1) * t5;
t31 = mrSges(6,2) * t4;
t17 = sin(qJ(1));
t29 = g(1) * t17;
t18 = cos(qJ(1));
t28 = g(2) * t18;
t16 = -pkin(6) - qJ(3);
t12 = -pkin(7) + t16;
t3 = t18 * t31;
t2 = t17 * t33;
t1 = [(-t38 * (t18 * pkin(1) + t17 * qJ(2)) + (-m(4) * qJ(3) + m(5) * t16 + m(6) * t12 - t35) * t18 + t34 * t17) * g(2) + ((m(3) * pkin(1) - m(4) * (-pkin(1) - qJ(3)) - m(5) * (-pkin(1) + t16) - m(6) * (-pkin(1) + t12) + t35) * t17 + (-qJ(2) * t38 + t34) * t18) * g(1), (t28 - t29) * t38, (g(1) * t18 + g(2) * t17) * t25, -g(1) * t2 - g(2) * t3 + (t23 * t6 - t39) * g(3) + (t33 - t36) * t28 + (t31 + t36) * t29, -g(1) * (-t17 * t31 + t2) - g(2) * (-t18 * t33 + t3) + g(3) * t22];
taug = t1(:);
