% Calculate Gravitation load on the joints for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:38
% EndTime: 2019-12-31 18:08:39
% DurationCPUTime: 0.30s
% Computational Cost: add. (191->50), mult. (173->51), div. (0->0), fcn. (129->8), ass. (0->25)
t35 = m(5) + m(6);
t40 = -mrSges(5,1) - mrSges(6,1);
t39 = mrSges(5,2) - mrSges(6,3);
t10 = qJ(1) + pkin(7);
t4 = sin(t10);
t6 = cos(t10);
t38 = g(1) * t6 + g(2) * t4;
t12 = sin(qJ(3));
t14 = cos(qJ(3));
t9 = qJ(3) + pkin(8);
t3 = sin(t9);
t5 = cos(t9);
t37 = -t14 * mrSges(4,1) + t12 * mrSges(4,2) + t39 * t3 + t40 * t5;
t36 = m(3) + m(4);
t32 = m(4) * pkin(2) + mrSges(3,1) - t37;
t31 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t13 = sin(qJ(1));
t28 = pkin(1) * t13;
t7 = t14 * pkin(3);
t15 = cos(qJ(1));
t8 = t15 * pkin(1);
t21 = pkin(4) * t5 + qJ(5) * t3;
t11 = -qJ(4) - pkin(6);
t2 = t7 + pkin(2);
t1 = [(-mrSges(2,1) * t15 + t13 * mrSges(2,2) - t36 * t8 - t35 * (-t11 * t4 + t6 * t2 + t8) + (-m(6) * t21 - t32) * t6 + t31 * t4) * g(2) + (t13 * mrSges(2,1) + mrSges(2,2) * t15 + t36 * t28 - t35 * (-t11 * t6 - t28) + t31 * t6 + (m(5) * t2 - m(6) * (-t2 - t21) + t32) * t4) * g(1), (-t35 - t36) * g(3), (-m(5) * t7 - m(6) * (t21 + t7) + t37) * g(3) + t38 * (mrSges(4,2) * t14 + (-m(6) * qJ(5) + t39) * t5 + (m(6) * pkin(4) - t40) * t3 + (t35 * pkin(3) + mrSges(4,1)) * t12), t35 * (-g(1) * t4 + g(2) * t6), (g(3) * t5 - t38 * t3) * m(6)];
taug = t1(:);
