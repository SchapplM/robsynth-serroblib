% Calculate Gravitation load on the joints for
% S5RPRRP10
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:55
% EndTime: 2019-12-31 18:50:57
% DurationCPUTime: 0.47s
% Computational Cost: add. (229->64), mult. (301->74), div. (0->0), fcn. (264->8), ass. (0->35)
t49 = -mrSges(5,1) - mrSges(6,1);
t48 = mrSges(5,2) + mrSges(6,2);
t53 = mrSges(5,3) + mrSges(6,3);
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t7 = t17 * pkin(4) + pkin(3);
t52 = -m(5) * pkin(3) - m(6) * t7 + t48 * t15 + t49 * t17;
t13 = -qJ(5) - pkin(7);
t51 = -m(5) * pkin(7) + m(6) * t13 - t53;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t50 = g(1) * t18 + g(2) * t16;
t42 = m(6) * pkin(4);
t47 = m(4) + m(5) + m(6);
t46 = t42 - t49;
t45 = -m(3) * qJ(2) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t12 = cos(pkin(8));
t10 = pkin(8) + qJ(3);
t8 = sin(t10);
t9 = cos(t10);
t25 = t9 * mrSges(4,1) - t8 * mrSges(4,2);
t43 = m(3) * pkin(1) + t12 * mrSges(3,1) - sin(pkin(8)) * mrSges(3,2) + mrSges(2,1) + t25 + t53 * t8;
t36 = t16 * t15;
t35 = t16 * t17;
t34 = t18 * t15;
t33 = t18 * t17;
t27 = t9 * pkin(3) + t8 * pkin(7);
t26 = -t8 * t13 + t9 * t7;
t3 = -t9 * t34 + t35;
t1 = t9 * t36 + t33;
t14 = -pkin(6) - qJ(2);
t6 = t12 * pkin(2) + pkin(1);
t4 = t9 * t33 + t36;
t2 = -t9 * t35 + t34;
t5 = [(-t36 * t42 + t49 * t4 - t48 * t3 - t47 * (-t16 * t14 + t18 * t6) + t45 * t16 + (-m(5) * t27 - m(6) * t26 - t43) * t18) * g(2) + (t49 * t2 - t48 * t1 + (t47 * t14 - t42 * t15 + t45) * t18 + (m(4) * t6 - m(5) * (-t27 - t6) - m(6) * (-t26 - t6) + t43) * t16) * g(1), (-g(1) * t16 + g(2) * t18) * (m(3) + t47), -g(3) * t25 + (t52 * g(3) + t50 * (mrSges(4,2) + t51)) * t9 + (t51 * g(3) + t50 * (mrSges(4,1) - t52)) * t8, (t46 * t15 + t48 * t17) * g(3) * t8 + (t46 * t1 - t48 * t2) * g(2) + (-t46 * t3 + t48 * t4) * g(1), (g(3) * t9 - t50 * t8) * m(6)];
taug = t5(:);
