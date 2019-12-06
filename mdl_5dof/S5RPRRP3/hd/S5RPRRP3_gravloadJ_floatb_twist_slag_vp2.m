% Calculate Gravitation load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:10
% EndTime: 2019-12-05 18:03:11
% DurationCPUTime: 0.30s
% Computational Cost: add. (203->51), mult. (182->52), div. (0->0), fcn. (134->8), ass. (0->30)
t19 = sin(qJ(3));
t18 = qJ(3) + qJ(4);
t13 = sin(t18);
t14 = cos(t18);
t32 = mrSges(5,2) + mrSges(6,2);
t33 = mrSges(5,1) + mrSges(6,1);
t27 = t32 * t13 - t33 * t14;
t50 = -mrSges(4,2) * t19 - t27;
t29 = t32 * t14;
t46 = m(3) + m(4) + m(5) + m(6);
t47 = t46 * pkin(1) + mrSges(2,1);
t21 = cos(qJ(3));
t30 = m(5) * pkin(3) + mrSges(4,1);
t45 = -t19 * t30 + m(6) * (-pkin(3) * t19 - pkin(4) * t13) - mrSges(4,2) * t21;
t15 = t21 * pkin(3);
t9 = pkin(4) * t14;
t37 = t9 + t15;
t44 = mrSges(3,1) + m(4) * pkin(2) + mrSges(4,1) * t21 + m(5) * (t15 + pkin(2)) + m(6) * (pkin(2) + t37) + t50;
t23 = -pkin(7) - pkin(6);
t42 = m(4) * pkin(6) - m(5) * t23 - m(6) * (-qJ(5) + t23) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t41 = m(6) * pkin(4);
t17 = qJ(1) + pkin(8);
t12 = cos(t17);
t40 = g(3) * t12;
t11 = sin(t17);
t35 = t11 * t13;
t31 = -t11 * t29 - t33 * t35;
t22 = cos(qJ(1));
t20 = sin(qJ(1));
t1 = [(mrSges(2,2) * t22 + t44 * t11 - t42 * t12 + t47 * t20) * g(3) + (-mrSges(2,2) * t20 + t42 * t11 + t44 * t12 + t47 * t22) * g(2), -t46 * g(1), (t45 * t11 + t31) * g(2) + (-m(6) * t37 - t30 * t21 - t50) * g(1) + (t13 * t33 + t29 - t45) * t40, (-t35 * t41 + t31) * g(2) + (-m(6) * t9 + t27) * g(1) + (t29 + (t33 + t41) * t13) * t40, (-g(2) * t12 - g(3) * t11) * m(6)];
taug = t1(:);
