% Calculate Gravitation load on the joints for
% S5RRPRR5
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:12
% EndTime: 2022-01-20 11:02:13
% DurationCPUTime: 0.26s
% Computational Cost: add. (275->53), mult. (207->56), div. (0->0), fcn. (157->10), ass. (0->31)
t33 = pkin(9) + qJ(4);
t27 = qJ(5) + t33;
t21 = sin(t27);
t22 = cos(t27);
t62 = t22 * mrSges(6,1) - t21 * mrSges(6,2);
t25 = sin(t33);
t61 = t25 * mrSges(5,2) - t62;
t26 = cos(t33);
t36 = cos(pkin(9));
t60 = -t36 * mrSges(4,1) - mrSges(3,1) + sin(pkin(9)) * mrSges(4,2) - t26 * mrSges(5,1) + t61;
t59 = mrSges(3,2) - mrSges(6,3) - mrSges(5,3) - mrSges(4,3);
t34 = qJ(1) + qJ(2);
t28 = sin(t34);
t29 = cos(t34);
t57 = g(1) * t29 + g(2) * t28;
t23 = t36 * pkin(3) + pkin(2);
t37 = -pkin(7) - qJ(3);
t49 = t29 * pkin(2) + t28 * qJ(3);
t48 = m(4) + m(5) + m(6);
t47 = m(6) * pkin(4) + mrSges(5,1);
t2 = pkin(4) * t26 + t23;
t32 = pkin(8) - t37;
t46 = t29 * t2 + t28 * t32;
t45 = t29 * t23 - t28 * t37;
t44 = mrSges(6,1) * t21 + mrSges(6,2) * t22;
t41 = t59 * t28 + t60 * t29;
t40 = (m(4) * pkin(2) + m(5) * t23 + m(6) * t2 - t60) * t28 + (-m(4) * qJ(3) + m(5) * t37 - m(6) * t32 + t59) * t29;
t39 = cos(qJ(1));
t38 = sin(qJ(1));
t31 = t39 * pkin(1);
t1 = [(t38 * mrSges(2,2) - m(4) * (t31 + t49) - m(5) * (t31 + t45) - m(6) * (t31 + t46) + (-m(3) * pkin(1) - mrSges(2,1)) * t39 + t41) * g(2) + (t39 * mrSges(2,2) + (mrSges(2,1) + (m(3) + t48) * pkin(1)) * t38 + t40) * g(1), (-m(4) * t49 - m(5) * t45 - m(6) * t46 + t41) * g(2) + t40 * g(1), (-g(1) * t28 + g(2) * t29) * t48, (-t26 * t47 + t61) * g(3) + t57 * (mrSges(5,2) * t26 + t25 * t47 + t44), -g(3) * t62 + t57 * t44];
taug = t1(:);
