% Calculate Gravitation load on the joints for
% S4RRPR9
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR9_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:21
% EndTime: 2019-12-31 17:09:22
% DurationCPUTime: 0.31s
% Computational Cost: add. (134->59), mult. (231->71), div. (0->0), fcn. (205->8), ass. (0->31)
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t44 = g(1) * t18 + g(2) * t16;
t43 = m(4) + m(5);
t17 = cos(qJ(2));
t12 = sin(pkin(7));
t13 = cos(pkin(7));
t22 = m(4) * pkin(2) + t13 * mrSges(4,1) - t12 * mrSges(4,2);
t42 = t22 * t17;
t41 = -m(3) - t43;
t15 = sin(qJ(2));
t26 = t17 * mrSges(3,1) - t15 * mrSges(3,2);
t40 = t15 * mrSges(5,3) + mrSges(2,1) + t26;
t38 = -t13 * mrSges(4,2) + mrSges(2,2) - mrSges(3,3) + (-m(5) * pkin(3) - mrSges(4,1)) * t12;
t32 = t16 * t17;
t31 = t17 * t18;
t14 = -pkin(6) - qJ(3);
t30 = -m(5) * t14 + mrSges(5,3);
t29 = m(4) * qJ(3) + mrSges(4,3);
t5 = pkin(3) * t13 + pkin(2);
t27 = -t14 * t15 + t17 * t5;
t11 = pkin(7) + qJ(4);
t6 = sin(t11);
t7 = cos(t11);
t23 = m(5) * t5 + t7 * mrSges(5,1) - t6 * mrSges(5,2);
t20 = t29 * t15 + t42;
t4 = t16 * t6 + t7 * t31;
t3 = t16 * t7 - t6 * t31;
t2 = t18 * t6 - t7 * t32;
t1 = t18 * t7 + t6 * t32;
t8 = [(-t4 * mrSges(5,1) - t3 * mrSges(5,2) + t41 * (t18 * pkin(1) + t16 * pkin(5)) + t38 * t16 + (-m(5) * t27 - t20 - t40) * t18) * g(2) + (-t2 * mrSges(5,1) - t1 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * (-qJ(3) * t15 - pkin(1)) + t15 * mrSges(4,3) + t42 - m(5) * (-pkin(1) - t27) + t40) * t16 + (t41 * pkin(5) + t38) * t18) * g(1), (-t20 - t26) * g(3) + (-t23 * g(3) + t44 * (mrSges(3,2) - t29 - t30)) * t17 + (-t30 * g(3) + t44 * (mrSges(3,1) + t22 + t23)) * t15, (t17 * g(3) - t15 * t44) * t43, -g(1) * (mrSges(5,1) * t3 - mrSges(5,2) * t4) - g(2) * (-mrSges(5,1) * t1 + t2 * mrSges(5,2)) - g(3) * (-mrSges(5,1) * t6 - mrSges(5,2) * t7) * t15];
taug = t8(:);
