% Calculate Gravitation load on the joints for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:29
% DurationCPUTime: 0.13s
% Computational Cost: add. (138->36), mult. (101->38), div. (0->0), fcn. (72->8), ass. (0->20)
t31 = mrSges(4,2) - mrSges(5,3);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t30 = t20 * mrSges(5,1) - t18 * mrSges(5,2);
t29 = -mrSges(4,1) - t30;
t28 = m(4) + m(5);
t17 = qJ(1) + qJ(2);
t15 = cos(t17);
t12 = pkin(2) * t15;
t13 = pkin(7) + t17;
t10 = sin(t13);
t11 = cos(t13);
t25 = t11 * pkin(3) + t10 * pkin(6) + t12;
t14 = sin(t17);
t23 = -t15 * mrSges(3,1) + t14 * mrSges(3,2) + t31 * t10 + t29 * t11;
t22 = mrSges(3,2) * t15 + (t28 * pkin(2) + mrSges(3,1)) * t14 + (m(5) * pkin(3) - t29) * t10 + (-m(5) * pkin(6) + t31) * t11;
t21 = cos(qJ(1));
t19 = sin(qJ(1));
t16 = t21 * pkin(1);
t1 = [(t19 * mrSges(2,2) - m(4) * (t12 + t16) - m(5) * (t16 + t25) + (-m(3) * pkin(1) - mrSges(2,1)) * t21 + t23) * g(2) + (mrSges(2,2) * t21 + (mrSges(2,1) + (m(3) + t28) * pkin(1)) * t19 + t22) * g(1), (-m(4) * t12 - m(5) * t25 + t23) * g(2) + t22 * g(1), -t28 * g(3), -g(3) * t30 + (g(1) * t11 + g(2) * t10) * (mrSges(5,1) * t18 + mrSges(5,2) * t20)];
taug = t1(:);
