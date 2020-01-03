% Calculate Gravitation load on the joints for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR11_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:12
% DurationCPUTime: 0.38s
% Computational Cost: add. (215->74), mult. (265->84), div. (0->0), fcn. (231->8), ass. (0->40)
t49 = m(5) + m(6);
t18 = pkin(8) + qJ(3);
t16 = sin(t18);
t13 = t16 * qJ(4);
t25 = cos(qJ(1));
t17 = cos(t18);
t44 = t17 * t25;
t54 = pkin(3) * t44 + t25 * t13;
t53 = (-mrSges(4,1) - mrSges(5,1)) * t17 + (mrSges(4,2) - mrSges(5,3)) * t16;
t23 = sin(qJ(1));
t52 = g(1) * t25 + g(2) * t23;
t20 = cos(pkin(8));
t51 = -mrSges(2,1) - m(3) * pkin(1) - t20 * mrSges(3,1) + sin(pkin(8)) * mrSges(3,2) + t53;
t21 = -pkin(6) - qJ(2);
t50 = mrSges(2,2) - m(3) * qJ(2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,2) - m(6) * (-pkin(7) - t21) + mrSges(6,3);
t48 = -pkin(3) - pkin(4);
t14 = t17 * pkin(3);
t24 = cos(qJ(5));
t45 = t16 * t24;
t43 = t14 + t13;
t42 = qJ(4) * t17;
t41 = t16 * t48;
t15 = t20 * pkin(2) + pkin(1);
t11 = t25 * t15;
t38 = -t23 * t21 + t11;
t37 = -t15 - t13;
t22 = sin(qJ(5));
t6 = -t17 * t22 + t45;
t1 = t6 * t23;
t28 = t16 * t22 + t17 * t24;
t2 = t28 * t23;
t35 = t1 * mrSges(6,1) - t2 * mrSges(6,2);
t3 = t22 * t44 - t25 * t45;
t4 = t28 * t25;
t34 = -t3 * mrSges(6,1) - t4 * mrSges(6,2);
t33 = -mrSges(6,1) * t28 - t6 * mrSges(6,2);
t26 = t17 * mrSges(5,3) + (-m(5) * pkin(3) - mrSges(5,1)) * t16;
t10 = t25 * t42;
t8 = t23 * t42;
t5 = [(-m(4) * t38 - m(5) * (t38 + t54) - m(6) * (pkin(4) * t44 + t11 + t54) - t4 * mrSges(6,1) + t3 * mrSges(6,2) + t51 * t25 + t50 * t23) * g(2) + (t2 * mrSges(6,1) + t1 * mrSges(6,2) + ((m(4) + m(5)) * t21 + t50) * t25 + (m(4) * t15 - m(5) * (t37 - t14) - m(6) * (t48 * t17 + t37) - t51) * t23) * g(1), (-g(1) * t23 + g(2) * t25) * (m(3) + m(4) + t49), t52 * (mrSges(4,1) * t16 + mrSges(4,2) * t17) + (-m(5) * t8 - t23 * t26 - m(6) * (t23 * t41 + t8) + t35) * g(2) + (-m(5) * t10 - t25 * t26 - m(6) * (t25 * t41 + t10) + t34) * g(1) + (-m(5) * t43 - m(6) * (t17 * pkin(4) + t43) + t33 + t53) * g(3), (g(3) * t17 - t16 * t52) * t49, -g(1) * t34 - g(2) * t35 - g(3) * t33];
taug = t5(:);
