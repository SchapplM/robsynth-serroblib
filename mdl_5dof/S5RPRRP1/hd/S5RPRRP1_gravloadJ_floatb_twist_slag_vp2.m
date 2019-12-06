% Calculate Gravitation load on the joints for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:14
% EndTime: 2019-12-05 17:59:16
% DurationCPUTime: 0.30s
% Computational Cost: add. (147->51), mult. (200->56), div. (0->0), fcn. (148->6), ass. (0->28)
t16 = cos(qJ(3));
t28 = mrSges(5,2) + mrSges(6,2);
t13 = qJ(3) + qJ(4);
t8 = cos(t13);
t24 = t28 * t8;
t29 = mrSges(5,1) + mrSges(6,1);
t14 = sin(qJ(3));
t36 = pkin(3) * t14;
t7 = sin(t13);
t48 = t24 + t16 * mrSges(4,2) + m(6) * (pkin(4) * t7 + t36) + t29 * t7;
t25 = t28 * t7;
t46 = t29 * t8;
t44 = m(3) + m(4) + m(5) + m(6);
t26 = m(5) * pkin(3) + mrSges(4,1);
t39 = pkin(4) * t8;
t43 = -t26 * t16 - m(6) * (pkin(3) * t16 + t39) + mrSges(4,2) * t14;
t42 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3) - mrSges(3,2);
t41 = -m(5) * t36 - t14 * mrSges(4,1) + mrSges(2,2) - mrSges(3,3) - t48;
t18 = -pkin(7) - pkin(6);
t15 = sin(qJ(1));
t38 = t15 * t46;
t17 = cos(qJ(1));
t37 = t17 * t25;
t35 = g(1) * t15;
t34 = g(2) * t17;
t23 = m(6) * pkin(4) + t29;
t12 = -qJ(5) + t18;
t1 = [(-t44 * (t17 * pkin(1) + t15 * qJ(2)) + (-m(4) * pkin(6) + m(5) * t18 + m(6) * t12 - t42) * t17 + t41 * t15) * g(2) + ((m(3) * pkin(1) - m(4) * (-pkin(1) - pkin(6)) - m(5) * (-pkin(1) + t18) - m(6) * (-pkin(1) + t12) + t42) * t15 + (-t44 * qJ(2) + t41) * t17) * g(1), (t34 - t35) * t44, -t37 * g(2) - t38 * g(1) + (t26 * t14 + t48) * g(3) + (-t43 + t46) * t34 + (t25 + t43) * t35, (t23 * t7 + t24) * g(3) + (t23 * t8 * t17 - t37) * g(2) + ((-m(6) * t39 + t25) * t15 - t38) * g(1), (-g(1) * t17 - g(2) * t15) * m(6)];
taug = t1(:);
