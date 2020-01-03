% Calculate Gravitation load on the joints for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:42
% EndTime: 2019-12-31 16:53:43
% DurationCPUTime: 0.24s
% Computational Cost: add. (124->49), mult. (167->60), div. (0->0), fcn. (143->8), ass. (0->26)
t37 = m(4) + m(5);
t35 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t11 = cos(pkin(7));
t9 = pkin(7) + qJ(3);
t7 = sin(t9);
t8 = cos(t9);
t21 = mrSges(4,1) * t8 - t7 * mrSges(4,2);
t34 = mrSges(2,1) + m(3) * pkin(1) + t11 * mrSges(3,1) - sin(pkin(7)) * mrSges(3,2) + t21 + t7 * mrSges(5,3);
t13 = sin(qJ(4));
t16 = cos(qJ(1));
t30 = t13 * t16;
t14 = sin(qJ(1));
t29 = t14 * t13;
t15 = cos(qJ(4));
t28 = t14 * t15;
t27 = t15 * t16;
t25 = m(5) * pkin(6) + mrSges(5,3);
t22 = pkin(3) * t8 + pkin(6) * t7;
t18 = m(5) * pkin(3) + t15 * mrSges(5,1) - t13 * mrSges(5,2);
t12 = -pkin(5) - qJ(2);
t6 = pkin(2) * t11 + pkin(1);
t4 = t8 * t27 + t29;
t3 = -t8 * t30 + t28;
t2 = -t8 * t28 + t30;
t1 = t8 * t29 + t27;
t5 = [(-t4 * mrSges(5,1) - t3 * mrSges(5,2) - t37 * (-t14 * t12 + t16 * t6) + t35 * t14 + (-m(5) * t22 - t34) * t16) * g(2) + (-t2 * mrSges(5,1) - t1 * mrSges(5,2) + (t37 * t12 + t35) * t16 + (m(4) * t6 - m(5) * (-t22 - t6) + t34) * t14) * g(1), (-g(1) * t14 + g(2) * t16) * (m(3) + t37), (-t18 * t8 - t25 * t7 - t21) * g(3) + ((mrSges(4,2) - t25) * t8 + (mrSges(4,1) + t18) * t7) * (g(1) * t16 + g(2) * t14), -g(1) * (mrSges(5,1) * t3 - mrSges(5,2) * t4) - g(2) * (-mrSges(5,1) * t1 + mrSges(5,2) * t2) - g(3) * (-mrSges(5,1) * t13 - mrSges(5,2) * t15) * t7];
taug = t5(:);
