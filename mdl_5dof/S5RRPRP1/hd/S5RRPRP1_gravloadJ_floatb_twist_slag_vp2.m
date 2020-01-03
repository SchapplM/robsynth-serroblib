% Calculate Gravitation load on the joints for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:58:54
% EndTime: 2020-01-03 11:58:55
% DurationCPUTime: 0.23s
% Computational Cost: add. (243->56), mult. (173->56), div. (0->0), fcn. (127->8), ass. (0->31)
t33 = cos(qJ(4));
t31 = sin(qJ(4));
t44 = mrSges(5,2) + mrSges(6,2);
t39 = t44 * t31;
t51 = -mrSges(6,1) - mrSges(5,1);
t52 = t33 * t51 - mrSges(4,1) + t39;
t49 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t29 = qJ(1) + qJ(2);
t25 = sin(t29);
t21 = pkin(2) * t25;
t26 = cos(t29);
t22 = pkin(2) * t26;
t32 = sin(qJ(1));
t27 = t32 * pkin(1);
t43 = t21 + t27;
t24 = pkin(8) + t29;
t19 = sin(t24);
t20 = cos(t24);
t23 = pkin(4) * t33 + pkin(3);
t30 = -qJ(5) - pkin(7);
t42 = t19 * t23 + t20 * t30 + t21;
t41 = t20 * pkin(3) + t19 * pkin(7) + t22;
t40 = -m(3) * pkin(1) - mrSges(2,1);
t38 = m(6) * pkin(4) - t51;
t37 = -t19 * t30 + t20 * t23 + t22;
t36 = -t26 * mrSges(3,1) + t25 * mrSges(3,2) + t49 * t19 + t52 * t20;
t35 = -t25 * mrSges(3,1) - t26 * mrSges(3,2) + (m(5) * pkin(7) - t49) * t20 + t52 * t19;
t34 = cos(qJ(1));
t28 = t34 * pkin(1);
t14 = t19 * pkin(3);
t1 = [(-mrSges(2,2) * t34 - m(4) * t43 - m(5) * (t14 + t43) - m(6) * (t27 + t42) + t40 * t32 + t35) * g(3) + (t32 * mrSges(2,2) - m(4) * (t22 + t28) - m(5) * (t28 + t41) - m(6) * (t28 + t37) + t40 * t34 + t36) * g(2), (-m(4) * t21 - m(5) * (t14 + t21) - m(6) * t42 + t35) * g(3) + (-m(4) * t22 - m(5) * t41 - m(6) * t37 + t36) * g(2), (-m(4) - m(5) - m(6)) * g(1), (-t38 * t33 + t39) * g(1) + (g(2) * t19 - g(3) * t20) * (t38 * t31 + t44 * t33), (g(2) * t20 + g(3) * t19) * m(6)];
taug = t1(:);
