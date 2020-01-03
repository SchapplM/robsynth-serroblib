% Calculate Gravitation load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:08
% EndTime: 2020-01-03 11:25:09
% DurationCPUTime: 0.30s
% Computational Cost: add. (181->50), mult. (203->63), div. (0->0), fcn. (179->8), ass. (0->27)
t34 = m(6) * pkin(4);
t39 = -mrSges(5,1) - mrSges(6,1);
t38 = mrSges(3,2) - mrSges(4,3);
t37 = mrSges(5,2) + mrSges(6,2);
t28 = m(4) + m(5) + m(6);
t36 = t34 - t39;
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t21 = cos(qJ(4));
t35 = -mrSges(3,1) + (-mrSges(4,1) - m(5) * pkin(3) - m(6) * (pkin(4) * t21 + pkin(3))) * t17 + (mrSges(4,2) - m(5) * pkin(6) - mrSges(5,3) - m(6) * (qJ(5) + pkin(6)) - mrSges(6,3)) * t16;
t20 = sin(qJ(1));
t13 = t20 * pkin(1);
t22 = cos(qJ(1));
t14 = t22 * pkin(1);
t15 = qJ(1) + pkin(7);
t11 = sin(t15);
t19 = sin(qJ(4));
t32 = t11 * t19;
t12 = cos(t15);
t31 = t12 * t19;
t30 = t17 * t19;
t29 = t17 * t21;
t3 = -t11 * t21 + t12 * t30;
t1 = -t11 * t30 - t12 * t21;
t4 = t12 * t29 + t32;
t2 = t11 * t29 - t31;
t5 = [(t31 * t34 - m(3) * t13 - t20 * mrSges(2,1) - mrSges(2,2) * t22 - t28 * (t11 * pkin(2) - qJ(3) * t12 + t13) + t39 * t2 - t38 * t12 - t37 * t1 + t35 * t11) * g(3) + (-t32 * t34 - m(3) * t14 - mrSges(2,1) * t22 + t20 * mrSges(2,2) + t39 * t4 + t37 * t3 - t28 * (t12 * pkin(2) + t11 * qJ(3) + t14) + t38 * t11 + t35 * t12) * g(2), (-m(3) - t28) * g(1), (g(2) * t12 + g(3) * t11) * t28, (t19 * t36 + t21 * t37) * g(1) * t16 + (-t3 * t36 - t37 * t4) * g(3) + (-t1 * t36 + t2 * t37) * g(2), (g(1) * t17 + (-g(2) * t11 + g(3) * t12) * t16) * m(6)];
taug = t5(:);
