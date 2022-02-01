% Calculate Gravitation load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:03
% EndTime: 2022-01-20 10:48:04
% DurationCPUTime: 0.24s
% Computational Cost: add. (271->53), mult. (182->53), div. (0->0), fcn. (135->10), ass. (0->30)
t26 = qJ(4) + qJ(5);
t21 = sin(t26);
t23 = cos(t26);
t52 = t23 * mrSges(6,1) - t21 * mrSges(6,2);
t28 = sin(qJ(4));
t51 = t28 * mrSges(5,2) - t52;
t30 = cos(qJ(4));
t50 = -t30 * mrSges(5,1) - mrSges(4,1) + t51;
t48 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t27 = qJ(1) + qJ(2);
t20 = pkin(9) + t27;
t16 = sin(t20);
t17 = cos(t20);
t47 = g(1) * t17 + g(2) * t16;
t24 = cos(t27);
t18 = pkin(2) * t24;
t41 = m(4) + m(5) + m(6);
t40 = t17 * pkin(3) + t16 * pkin(7) + t18;
t39 = m(6) * pkin(4) + mrSges(5,1);
t19 = t30 * pkin(4) + pkin(3);
t32 = -pkin(8) - pkin(7);
t38 = -t16 * t32 + t17 * t19 + t18;
t37 = mrSges(6,1) * t21 + mrSges(6,2) * t23;
t22 = sin(t27);
t34 = -t24 * mrSges(3,1) + t22 * mrSges(3,2) + t48 * t16 + t50 * t17;
t33 = t24 * mrSges(3,2) + (t41 * pkin(2) + mrSges(3,1)) * t22 + (m(5) * pkin(3) + m(6) * t19 - t50) * t16 + (-m(5) * pkin(7) + m(6) * t32 + t48) * t17;
t31 = cos(qJ(1));
t29 = sin(qJ(1));
t25 = t31 * pkin(1);
t1 = [(t29 * mrSges(2,2) - m(4) * (t18 + t25) - m(5) * (t25 + t40) - m(6) * (t25 + t38) + (-m(3) * pkin(1) - mrSges(2,1)) * t31 + t34) * g(2) + (t31 * mrSges(2,2) + (mrSges(2,1) + (m(3) + t41) * pkin(1)) * t29 + t33) * g(1), (-m(4) * t18 - m(5) * t40 - m(6) * t38 + t34) * g(2) + t33 * g(1), -t41 * g(3), (-t39 * t30 + t51) * g(3) + t47 * (mrSges(5,2) * t30 + t39 * t28 + t37), -g(3) * t52 + t47 * t37];
taug = t1(:);
