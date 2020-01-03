% Calculate Gravitation load on the joints for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:51
% EndTime: 2019-12-31 17:43:52
% DurationCPUTime: 0.33s
% Computational Cost: add. (153->52), mult. (171->64), div. (0->0), fcn. (144->8), ass. (0->27)
t33 = m(5) + m(6);
t29 = m(4) + t33;
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t35 = -mrSges(3,1) + (-mrSges(4,1) - mrSges(5,1)) * t15 + (mrSges(4,2) - mrSges(5,3)) * t14;
t34 = pkin(6) * m(6) + mrSges(3,2) - mrSges(5,2) - mrSges(4,3) + mrSges(6,3);
t17 = sin(qJ(1));
t32 = pkin(1) * t17;
t19 = cos(qJ(1));
t12 = t19 * pkin(1);
t13 = qJ(1) + pkin(7);
t11 = cos(t13);
t31 = t11 * t15;
t30 = qJ(4) * t14;
t10 = sin(t13);
t28 = t11 * pkin(2) + t10 * qJ(3) + t12;
t26 = -pkin(2) - t30;
t25 = pkin(3) * t31 + t11 * t30 + t28;
t16 = sin(qJ(5));
t18 = cos(qJ(5));
t21 = t14 * t18 - t15 * t16;
t20 = t14 * t16 + t15 * t18;
t4 = t20 * t11;
t3 = t21 * t11;
t2 = t20 * t10;
t1 = t21 * t10;
t5 = [(-mrSges(2,1) * t19 + t17 * mrSges(2,2) - m(3) * t12 - m(4) * t28 - m(5) * t25 - m(6) * (pkin(4) * t31 + t25) - t4 * mrSges(6,1) - t3 * mrSges(6,2) + t35 * t11 + t34 * t10) * g(2) + (m(3) * t32 + t17 * mrSges(2,1) + t2 * mrSges(6,1) + mrSges(2,2) * t19 + t1 * mrSges(6,2) - t29 * (t11 * qJ(3) - t32) + (m(4) * pkin(2) - m(5) * (-pkin(3) * t15 + t26) - m(6) * ((-pkin(3) - pkin(4)) * t15 + t26) - t35) * t10 + t34 * t11) * g(1), (-m(3) - t29) * g(3), (-g(1) * t10 + g(2) * t11) * t29, (t15 * g(3) + t14 * (-g(1) * t11 - g(2) * t10)) * t33, -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(3) * (-t20 * mrSges(6,1) - t21 * mrSges(6,2))];
taug = t5(:);
