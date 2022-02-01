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
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 10:19:33
% EndTime: 2022-01-20 10:19:34
% DurationCPUTime: 0.24s
% Computational Cost: add. (243->50), mult. (173->50), div. (0->0), fcn. (127->8), ass. (0->26)
t27 = cos(qJ(4));
t25 = sin(qJ(4));
t37 = mrSges(5,2) + mrSges(6,2);
t34 = t37 * t25;
t38 = mrSges(5,1) + mrSges(6,1);
t43 = -t38 * t27 - mrSges(4,1) + t34;
t41 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t23 = qJ(1) + qJ(2);
t21 = cos(t23);
t17 = pkin(2) * t21;
t36 = m(4) + m(5) + m(6);
t19 = pkin(8) + t23;
t15 = sin(t19);
t16 = cos(t19);
t35 = t16 * pkin(3) + t15 * pkin(7) + t17;
t33 = m(6) * pkin(4) + t38;
t18 = t27 * pkin(4) + pkin(3);
t24 = -qJ(5) - pkin(7);
t32 = -t15 * t24 + t16 * t18 + t17;
t20 = sin(t23);
t30 = -t21 * mrSges(3,1) + t20 * mrSges(3,2) + t41 * t15 + t43 * t16;
t29 = t21 * mrSges(3,2) + (t36 * pkin(2) + mrSges(3,1)) * t20 + (-m(5) * pkin(7) + m(6) * t24 + t41) * t16 + (m(5) * pkin(3) + m(6) * t18 - t43) * t15;
t28 = cos(qJ(1));
t26 = sin(qJ(1));
t22 = t28 * pkin(1);
t1 = [(t26 * mrSges(2,2) - m(4) * (t17 + t22) - m(5) * (t22 + t35) - m(6) * (t22 + t32) + (-m(3) * pkin(1) - mrSges(2,1)) * t28 + t30) * g(2) + (t28 * mrSges(2,2) + (mrSges(2,1) + (m(3) + t36) * pkin(1)) * t26 + t29) * g(1), (-m(4) * t17 - m(5) * t35 - m(6) * t32 + t30) * g(2) + t29 * g(1), -t36 * g(3), (-t33 * t27 + t34) * g(3) + (g(1) * t16 + g(2) * t15) * (t33 * t25 + t37 * t27), (-g(1) * t15 + g(2) * t16) * m(6)];
taug = t1(:);
