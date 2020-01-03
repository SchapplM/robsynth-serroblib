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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:46:54
% EndTime: 2020-01-03 11:46:56
% DurationCPUTime: 0.30s
% Computational Cost: add. (203->50), mult. (182->52), div. (0->0), fcn. (134->8), ass. (0->30)
t21 = sin(qJ(3));
t20 = qJ(3) + qJ(4);
t13 = sin(t20);
t14 = cos(t20);
t34 = mrSges(5,2) + mrSges(6,2);
t35 = mrSges(5,1) + mrSges(6,1);
t29 = t34 * t13 - t35 * t14;
t50 = mrSges(4,2) * t21 + t29;
t31 = t34 * t14;
t46 = -m(3) - m(4) - m(5) - m(6);
t47 = pkin(1) * t46 - mrSges(2,1);
t23 = cos(qJ(3));
t32 = m(5) * pkin(3) + mrSges(4,1);
t45 = -t32 * t21 + m(6) * (-pkin(3) * t21 - pkin(4) * t13) - mrSges(4,2) * t23;
t16 = t23 * pkin(3);
t9 = pkin(4) * t14;
t39 = t9 + t16;
t44 = -mrSges(3,1) - m(4) * pkin(2) - mrSges(4,1) * t23 - m(5) * (t16 + pkin(2)) - m(6) * (pkin(2) + t39) + t50;
t25 = -pkin(7) - pkin(6);
t42 = m(4) * pkin(6) - m(5) * t25 - m(6) * (-qJ(5) + t25) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t41 = m(6) * pkin(4);
t19 = qJ(1) + pkin(8);
t11 = sin(t19);
t40 = g(2) * t11;
t12 = cos(t19);
t37 = t12 * t13;
t33 = -t12 * t31 - t35 * t37;
t24 = cos(qJ(1));
t22 = sin(qJ(1));
t1 = [(-mrSges(2,2) * t24 + t44 * t11 + t42 * t12 + t47 * t22) * g(3) + (mrSges(2,2) * t22 - t42 * t11 + t44 * t12 + t47 * t24) * g(2), t46 * g(1), (t45 * t12 + t33) * g(3) + (-m(6) * t39 - t32 * t23 + t50) * g(1) + (t35 * t13 + t31 - t45) * t40, (-t37 * t41 + t33) * g(3) + (-m(6) * t9 + t29) * g(1) + (t31 + (t35 + t41) * t13) * t40, (g(2) * t12 + g(3) * t11) * m(6)];
taug = t1(:);
