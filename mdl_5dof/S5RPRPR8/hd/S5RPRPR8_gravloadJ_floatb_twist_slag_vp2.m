% Calculate Gravitation load on the joints for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:18
% EndTime: 2019-12-31 18:21:19
% DurationCPUTime: 0.36s
% Computational Cost: add. (237->69), mult. (251->81), div. (0->0), fcn. (217->10), ass. (0->36)
t44 = m(5) + m(6);
t48 = -m(4) - t44;
t15 = qJ(1) + pkin(8);
t10 = sin(t15);
t12 = cos(t15);
t50 = g(1) * t12 + g(2) * t10;
t21 = cos(qJ(3));
t16 = sin(pkin(9));
t17 = cos(pkin(9));
t26 = m(5) * pkin(3) + t17 * mrSges(5,1) - t16 * mrSges(5,2);
t49 = t26 * t21;
t19 = sin(qJ(3));
t30 = t21 * mrSges(4,1) - t19 * mrSges(4,2);
t47 = t19 * mrSges(6,3) + mrSges(3,1) + t30;
t45 = -t17 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3) + (-m(6) * pkin(4) - mrSges(5,1)) * t16;
t20 = sin(qJ(1));
t40 = t20 * pkin(1);
t14 = pkin(9) + qJ(5);
t9 = sin(t14);
t39 = t21 * t9;
t22 = cos(qJ(1));
t13 = t22 * pkin(1);
t11 = cos(t14);
t37 = t21 * t11;
t18 = -pkin(7) - qJ(4);
t34 = -m(6) * t18 + mrSges(6,3);
t33 = m(5) * qJ(4) + mrSges(5,3);
t8 = t17 * pkin(4) + pkin(3);
t31 = -t19 * t18 + t21 * t8;
t27 = m(6) * t8 + t11 * mrSges(6,1) - t9 * mrSges(6,2);
t24 = t19 * t33 + t49;
t4 = t10 * t9 + t12 * t37;
t3 = t10 * t11 - t12 * t39;
t2 = -t10 * t37 + t12 * t9;
t1 = t10 * t39 + t12 * t11;
t5 = [(-m(3) * t13 - t22 * mrSges(2,1) - t4 * mrSges(6,1) + t20 * mrSges(2,2) - t3 * mrSges(6,2) + t48 * (t12 * pkin(2) + t10 * pkin(6) + t13) + t45 * t10 + (-m(6) * t31 - t24 - t47) * t12) * g(2) + (m(3) * t40 + t20 * mrSges(2,1) - t2 * mrSges(6,1) + t22 * mrSges(2,2) - t1 * mrSges(6,2) + t48 * (t12 * pkin(6) - t40) + t45 * t12 + (m(4) * pkin(2) - m(5) * (-t19 * qJ(4) - pkin(2)) + t19 * mrSges(5,3) + t49 - m(6) * (-pkin(2) - t31) + t47) * t10) * g(1), (-m(3) + t48) * g(3), (-t24 - t30) * g(3) + (-t27 * g(3) + t50 * (mrSges(4,2) - t33 - t34)) * t21 + (-t34 * g(3) + t50 * (mrSges(4,1) + t26 + t27)) * t19, (t21 * g(3) - t19 * t50) * t44, -g(1) * (mrSges(6,1) * t3 - mrSges(6,2) * t4) - g(2) * (-mrSges(6,1) * t1 + t2 * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t9 - mrSges(6,2) * t11) * t19];
taug = t5(:);
