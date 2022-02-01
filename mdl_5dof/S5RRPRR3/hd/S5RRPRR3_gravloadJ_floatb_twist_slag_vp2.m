% Calculate Gravitation load on the joints for
% S5RRPRR3
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:33:58
% EndTime: 2022-01-20 10:33:59
% DurationCPUTime: 0.19s
% Computational Cost: add. (284->49), mult. (156->49), div. (0->0), fcn. (112->10), ass. (0->30)
t42 = mrSges(5,2) - mrSges(6,3);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t41 = mrSges(6,1) * t25 - mrSges(6,2) * t23;
t40 = -mrSges(5,1) - t41;
t39 = m(5) + m(6);
t22 = qJ(1) + qJ(2);
t18 = pkin(9) + t22;
t17 = qJ(4) + t18;
t11 = sin(t17);
t12 = cos(t17);
t38 = t12 * pkin(4) + t11 * pkin(8);
t20 = cos(t22);
t16 = pkin(2) * t20;
t15 = cos(t18);
t10 = pkin(3) * t15;
t35 = t10 + t16;
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t34 = t16 + t21;
t33 = m(4) + t39;
t32 = t35 + t38;
t30 = t42 * t11 + t40 * t12;
t29 = (m(6) * pkin(4) - t40) * t11 + (-m(6) * pkin(8) + t42) * t12;
t14 = sin(t18);
t19 = sin(t22);
t28 = -t20 * mrSges(3,1) - t15 * mrSges(4,1) + t19 * mrSges(3,2) + t14 * mrSges(4,2) + t30;
t27 = mrSges(3,2) * t20 + t15 * mrSges(4,2) + (t39 * pkin(3) + mrSges(4,1)) * t14 + (t33 * pkin(2) + mrSges(3,1)) * t19 + t29;
t24 = sin(qJ(1));
t1 = [(t24 * mrSges(2,2) - m(4) * t34 - m(5) * (t10 + t34) - m(6) * (t21 + t32) + (-m(3) * pkin(1) - mrSges(2,1)) * t26 + t28) * g(2) + (mrSges(2,2) * t26 + (mrSges(2,1) + (m(3) + t33) * pkin(1)) * t24 + t27) * g(1), (-m(4) * t16 - m(5) * t35 - m(6) * t32 + t28) * g(2) + t27 * g(1), -t33 * g(3), (-m(6) * t38 + t30) * g(2) + t29 * g(1), -g(3) * t41 + (g(1) * t12 + g(2) * t11) * (mrSges(6,1) * t23 + mrSges(6,2) * t25)];
taug = t1(:);
