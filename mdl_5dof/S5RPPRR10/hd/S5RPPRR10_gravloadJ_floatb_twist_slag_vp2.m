% Calculate Gravitation load on the joints for
% S5RPPRR10
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:01
% DurationCPUTime: 0.47s
% Computational Cost: add. (164->70), mult. (292->88), div. (0->0), fcn. (271->8), ass. (0->37)
t41 = m(4) + m(5) + m(6);
t52 = m(3) + t41;
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t51 = -mrSges(2,1) + (-mrSges(3,1) - mrSges(4,1)) * t25 + (mrSges(3,2) - mrSges(4,3)) * t24;
t50 = m(5) * pkin(6) - m(6) * (-pkin(7) - pkin(6)) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t49 = m(6) * pkin(4);
t27 = sin(qJ(1));
t23 = qJ(4) + qJ(5);
t18 = sin(t23);
t19 = cos(t23);
t35 = t25 * t18 - t24 * t19;
t5 = t35 * t27;
t34 = t24 * t18 + t25 * t19;
t6 = t34 * t27;
t48 = -t5 * mrSges(6,1) - t6 * mrSges(6,2);
t29 = cos(qJ(1));
t7 = t35 * t29;
t8 = t34 * t29;
t47 = -t7 * mrSges(6,1) - t8 * mrSges(6,2);
t46 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t26 = sin(qJ(4));
t45 = t24 * t26;
t44 = t25 * t29;
t43 = t29 * pkin(1) + t27 * qJ(2);
t42 = qJ(3) * t24;
t39 = -pkin(1) - t42;
t28 = cos(qJ(4));
t33 = t24 * t28 - t25 * t26;
t32 = t25 * t28 + t45;
t31 = t33 * t49;
t17 = t28 * pkin(4) + pkin(3);
t14 = t32 * t29;
t13 = t33 * t29;
t12 = t32 * t27;
t11 = t33 * t27;
t1 = [(-m(5) * pkin(3) * t44 - m(3) * t43 - t14 * mrSges(5,1) - t8 * mrSges(6,1) - t13 * mrSges(5,2) + t7 * mrSges(6,2) - t41 * (pkin(2) * t44 + t29 * t42 + t43) + (-m(6) * (pkin(4) * t45 + t17 * t25) + t51) * t29 + t50 * t27) * g(2) + (t12 * mrSges(5,1) + t6 * mrSges(6,1) + t11 * mrSges(5,2) - t5 * mrSges(6,2) + (m(3) * pkin(1) - m(4) * (-pkin(2) * t25 + t39) - m(5) * ((-pkin(2) - pkin(3)) * t25 + t39) - m(6) * (-pkin(1) + (-pkin(2) - t17) * t25 + (-pkin(4) * t26 - qJ(3)) * t24) - t51) * t27 + (-t52 * qJ(2) + t50) * t29) * g(1), (-g(1) * t27 + g(2) * t29) * t52, (g(3) * t25 + (-g(1) * t29 - g(2) * t27) * t24) * t41, (t33 * mrSges(5,2) - t46 + (mrSges(5,1) + t49) * t32) * g(3) + (-mrSges(5,1) * t11 + mrSges(5,2) * t12 - t27 * t31 - t48) * g(2) + (-mrSges(5,1) * t13 + mrSges(5,2) * t14 - t29 * t31 - t47) * g(1), -g(1) * t47 - g(2) * t48 - g(3) * t46];
taug = t1(:);
