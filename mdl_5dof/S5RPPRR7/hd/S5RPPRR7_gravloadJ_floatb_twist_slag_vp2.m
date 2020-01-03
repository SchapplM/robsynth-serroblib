% Calculate Gravitation load on the joints for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:32
% DurationCPUTime: 0.29s
% Computational Cost: add. (164->53), mult. (181->67), div. (0->0), fcn. (149->8), ass. (0->27)
t12 = qJ(1) + pkin(8);
t10 = cos(t12);
t9 = sin(t12);
t37 = -g(1) * t9 + g(2) * t10;
t38 = -m(5) - m(6);
t27 = m(4) - t38;
t36 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t14 = sin(qJ(4));
t17 = cos(qJ(4));
t21 = t14 * mrSges(5,1) + t17 * mrSges(5,2);
t35 = mrSges(3,2) - mrSges(4,3) - m(6) * (pkin(4) * t14 - pkin(7) * t17) + t17 * mrSges(6,3) - t21;
t15 = sin(qJ(1));
t31 = pkin(1) * t15;
t18 = cos(qJ(1));
t11 = t18 * pkin(1);
t13 = sin(qJ(5));
t29 = t14 * t13;
t16 = cos(qJ(5));
t28 = t14 * t16;
t26 = t10 * pkin(2) + t9 * qJ(3) + t11;
t25 = m(6) * pkin(7) + mrSges(6,3);
t20 = m(6) * pkin(4) + t16 * mrSges(6,1) - t13 * mrSges(6,2);
t4 = t10 * t28 - t13 * t9;
t3 = t10 * t29 + t16 * t9;
t2 = t10 * t13 + t9 * t28;
t1 = t10 * t16 - t9 * t29;
t5 = [(-m(3) * t11 - m(4) * t26 - mrSges(2,1) * t18 - t2 * mrSges(6,1) + t15 * mrSges(2,2) - t1 * mrSges(6,2) + t38 * (t10 * pkin(6) + t26) + t36 * t10 + t35 * t9) * g(2) + (m(3) * t31 + t15 * mrSges(2,1) - t4 * mrSges(6,1) + mrSges(2,2) * t18 + t3 * mrSges(6,2) - t27 * (t10 * qJ(3) - t31) + (m(4) * pkin(2) + t38 * (-pkin(2) - pkin(6)) - t36) * t9 + t35 * t10) * g(1), (-m(3) - t27) * g(3), t37 * t27, (t20 * t14 - t25 * t17 + t21) * g(3) + ((mrSges(5,1) + t20) * t17 + (-mrSges(5,2) + t25) * t14) * t37, -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * (mrSges(6,1) * t3 + mrSges(6,2) * t4) - g(3) * (-mrSges(6,1) * t13 - mrSges(6,2) * t16) * t17];
taug = t5(:);
