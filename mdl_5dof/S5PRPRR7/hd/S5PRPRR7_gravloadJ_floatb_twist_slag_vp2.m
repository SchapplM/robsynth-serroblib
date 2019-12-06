% Calculate Gravitation load on the joints for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:59:44
% EndTime: 2019-12-05 15:59:46
% DurationCPUTime: 0.39s
% Computational Cost: add. (130->50), mult. (244->68), div. (0->0), fcn. (209->8), ass. (0->26)
t41 = m(6) * pkin(4);
t45 = -mrSges(5,1) - t41;
t52 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t12 = qJ(4) + qJ(5);
t8 = sin(t12);
t9 = cos(t12);
t51 = -t8 * mrSges(6,1) - t17 * mrSges(5,2) - t9 * mrSges(6,2) + t45 * t15 + mrSges(3,2) - mrSges(4,3);
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t44 = g(1) * t14 + g(2) * t13;
t28 = m(4) + m(5) + m(6);
t16 = sin(qJ(2));
t33 = t14 * t16;
t40 = (-t13 * t8 + t9 * t33) * mrSges(6,1) + (-t13 * t9 - t8 * t33) * mrSges(6,2);
t34 = t13 * t16;
t39 = (t14 * t8 + t9 * t34) * mrSges(6,1) + (t14 * t9 - t8 * t34) * mrSges(6,2);
t38 = mrSges(6,1) * t9;
t18 = cos(qJ(2));
t35 = g(3) * t18;
t32 = t15 * t16;
t31 = t16 * t17;
t19 = -pkin(7) - pkin(6);
t5 = t18 * t8 * mrSges(6,2);
t1 = [(-m(2) - m(3) - t28) * g(3), (-t28 * (t18 * pkin(2) + t16 * qJ(3)) + (-m(5) * pkin(6) + m(6) * t19 - t52) * t18 + t51 * t16) * g(3) + ((-m(6) * (-pkin(2) + t19) - m(5) * (-pkin(2) - pkin(6)) + m(4) * pkin(2) + t52) * t16 + (-t28 * qJ(3) + t51) * t18) * t44, (-t44 * t16 + t35) * t28, -(-mrSges(5,1) * t17 + mrSges(5,2) * t15) * t35 - g(3) * (t5 + (-t17 * t41 - t38) * t18) + (-(-t13 * t32 + t14 * t17) * mrSges(5,2) - t39 + t45 * (t13 * t31 + t14 * t15)) * g(2) + (-(-t13 * t17 - t14 * t32) * mrSges(5,2) - t40 + t45 * (-t13 * t15 + t14 * t31)) * g(1), -g(1) * t40 - g(2) * t39 - g(3) * (-t18 * t38 + t5)];
taug = t1(:);
