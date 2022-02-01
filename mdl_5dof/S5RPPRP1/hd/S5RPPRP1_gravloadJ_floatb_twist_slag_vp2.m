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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:12:13
% EndTime: 2022-01-23 09:12:14
% DurationCPUTime: 0.33s
% Computational Cost: add. (181->54), mult. (203->67), div. (0->0), fcn. (179->8), ass. (0->29)
t34 = m(6) * pkin(4);
t39 = -mrSges(5,1) - mrSges(6,1);
t38 = mrSges(3,2) - mrSges(4,3);
t37 = mrSges(5,2) + mrSges(6,2);
t25 = m(4) + m(5) + m(6);
t36 = t34 - t39;
t13 = sin(pkin(8));
t14 = cos(pkin(8));
t35 = t14 * mrSges(4,1) + mrSges(3,1) + (-mrSges(4,2) + mrSges(5,3) + mrSges(6,3)) * t13;
t17 = sin(qJ(1));
t33 = pkin(1) * t17;
t16 = sin(qJ(4));
t12 = qJ(1) + pkin(7);
t9 = sin(t12);
t31 = t16 * t9;
t19 = cos(qJ(1));
t11 = t19 * pkin(1);
t10 = cos(t12);
t30 = t10 * t16;
t27 = t14 * t16;
t18 = cos(qJ(4));
t26 = t14 * t18;
t22 = pkin(3) * t14 + pkin(6) * t13;
t21 = -t13 * (-qJ(5) - pkin(6)) + t14 * (pkin(4) * t18 + pkin(3));
t1 = t10 * t18 + t27 * t9;
t3 = -t10 * t27 + t18 * t9;
t4 = t10 * t26 + t31;
t2 = -t26 * t9 + t30;
t5 = [(-t31 * t34 - m(3) * t11 - mrSges(2,1) * t19 + t17 * mrSges(2,2) + t38 * t9 + t39 * t4 - t37 * t3 - t25 * (t10 * pkin(2) + t9 * qJ(3) + t11) + (-m(5) * t22 - m(6) * t21 - t35) * t10) * g(2) + (-t30 * t34 + m(3) * t33 + t17 * mrSges(2,1) + mrSges(2,2) * t19 - t25 * (t10 * qJ(3) - t33) + t39 * t2 + t38 * t10 - t37 * t1 + (m(4) * pkin(2) - m(5) * (-pkin(2) - t22) - m(6) * (-pkin(2) - t21) + t35) * t9) * g(1), (-m(3) - t25) * g(3), (-g(1) * t9 + g(2) * t10) * t25, (t16 * t36 + t18 * t37) * g(3) * t13 + (t1 * t36 - t2 * t37) * g(2) + (-t3 * t36 + t37 * t4) * g(1), (g(3) * t14 + (-g(1) * t10 - g(2) * t9) * t13) * m(6)];
taug = t5(:);
