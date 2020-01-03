% Calculate Gravitation load on the joints for
% S5RPPRR8
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:05
% DurationCPUTime: 0.28s
% Computational Cost: add. (182->56), mult. (210->63), div. (0->0), fcn. (208->8), ass. (0->30)
t17 = sin(qJ(5));
t19 = cos(qJ(5));
t23 = -mrSges(6,1) * t19 + mrSges(6,2) * t17;
t41 = m(6) * pkin(4) + mrSges(5,1) - t23;
t40 = mrSges(2,1) + mrSges(3,1);
t39 = mrSges(2,2) - mrSges(3,3);
t18 = sin(qJ(1));
t20 = cos(qJ(1));
t32 = pkin(8) + qJ(4);
t26 = sin(t32);
t27 = cos(t32);
t1 = -t18 * t26 - t20 * t27;
t2 = -t18 * t27 + t20 * t26;
t38 = t2 * mrSges(5,2) + t1 * t41;
t37 = t1 * mrSges(5,2) - t2 * t41;
t15 = sin(pkin(8));
t36 = t15 * t20;
t35 = t18 * t15;
t34 = t20 * pkin(1) + t18 * qJ(2);
t33 = m(4) + m(5) + m(6);
t30 = -m(6) * pkin(7) - mrSges(6,3);
t13 = t20 * qJ(2);
t29 = -pkin(1) * t18 + t13;
t16 = cos(pkin(8));
t11 = pkin(3) * t16 + pkin(2);
t28 = pkin(3) * t35 + t20 * t11 + t34;
t10 = pkin(3) * t36;
t4 = t16 * t20 + t35;
t3 = -t18 * t16 + t36;
t5 = [(-t4 * mrSges(4,1) + t3 * mrSges(4,2) - m(5) * t28 - m(6) * (pkin(7) * t2 + t28) - t2 * mrSges(6,3) + (-m(3) - m(4)) * t34 + (-m(4) * pkin(2) - t40) * t20 + t39 * t18 + t38) * g(2) + (-m(3) * t29 - m(4) * t13 - t3 * mrSges(4,1) - t4 * mrSges(4,2) - m(5) * (t10 + t13) - m(6) * (pkin(7) * t1 + t10 + t29) - t1 * mrSges(6,3) + t39 * t20 + (-m(4) * (-pkin(1) - pkin(2)) - m(5) * (-pkin(1) - t11) + m(6) * t11 + t40) * t18 + t37) * g(1), (-g(1) * t18 + g(2) * t20) * (m(3) + t33), t33 * g(3), (-t30 * t2 - t38) * g(2) + (-t30 * t1 - t37) * g(1), -g(3) * t23 + (-g(1) * t1 - g(2) * t2) * (mrSges(6,1) * t17 + mrSges(6,2) * t19)];
taug = t5(:);
