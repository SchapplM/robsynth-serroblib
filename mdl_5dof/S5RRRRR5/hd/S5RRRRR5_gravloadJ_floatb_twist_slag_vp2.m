% Calculate Gravitation load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:40
% EndTime: 2022-01-20 12:01:41
% DurationCPUTime: 0.28s
% Computational Cost: add. (348->55), mult. (223->56), div. (0->0), fcn. (169->10), ass. (0->34)
t26 = qJ(4) + qJ(5);
t20 = sin(t26);
t22 = cos(t26);
t56 = t22 * mrSges(6,1) - t20 * mrSges(6,2);
t28 = sin(qJ(4));
t55 = t28 * mrSges(5,2) - t56;
t30 = cos(qJ(4));
t54 = -t30 * mrSges(5,1) - mrSges(4,1) + t55;
t52 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t27 = qJ(1) + qJ(2);
t24 = qJ(3) + t27;
t17 = sin(t24);
t18 = cos(t24);
t51 = g(1) * t18 + g(2) * t17;
t23 = cos(t27);
t16 = pkin(2) * t23;
t45 = t18 * pkin(3) + t17 * pkin(8);
t44 = m(4) + m(5) + m(6);
t43 = t16 + t45;
t42 = m(6) * pkin(4) + mrSges(5,1);
t19 = t30 * pkin(4) + pkin(3);
t32 = -pkin(9) - pkin(8);
t41 = -t17 * t32 + t18 * t19;
t40 = t16 + t41;
t39 = mrSges(6,1) * t20 + mrSges(6,2) * t22;
t36 = t52 * t17 + t54 * t18;
t21 = sin(t27);
t35 = -t23 * mrSges(3,1) + t21 * mrSges(3,2) + t36;
t34 = (m(5) * pkin(3) + m(6) * t19 - t54) * t17 + (-m(5) * pkin(8) + m(6) * t32 + t52) * t18;
t33 = t23 * mrSges(3,2) + (t44 * pkin(2) + mrSges(3,1)) * t21 + t34;
t31 = cos(qJ(1));
t29 = sin(qJ(1));
t25 = t31 * pkin(1);
t1 = [(t29 * mrSges(2,2) - m(4) * (t16 + t25) - m(5) * (t25 + t43) - m(6) * (t25 + t40) + (-m(3) * pkin(1) - mrSges(2,1)) * t31 + t35) * g(2) + (t31 * mrSges(2,2) + (mrSges(2,1) + (m(3) + t44) * pkin(1)) * t29 + t33) * g(1), (-m(4) * t16 - m(5) * t43 - m(6) * t40 + t35) * g(2) + t33 * g(1), (-m(5) * t45 - m(6) * t41 + t36) * g(2) + t34 * g(1), (-t42 * t30 + t55) * g(3) + t51 * (mrSges(5,2) * t30 + t42 * t28 + t39), -g(3) * t56 + t51 * t39];
taug = t1(:);
