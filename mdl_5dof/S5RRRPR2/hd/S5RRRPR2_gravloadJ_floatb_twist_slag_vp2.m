% Calculate Gravitation load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:33
% EndTime: 2022-01-20 11:30:34
% DurationCPUTime: 0.21s
% Computational Cost: add. (307->51), mult. (167->50), div. (0->0), fcn. (120->10), ass. (0->29)
t41 = mrSges(5,2) - mrSges(6,3);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t40 = mrSges(6,1) * t25 - mrSges(6,2) * t23;
t39 = -mrSges(5,1) - t40;
t38 = m(5) + m(6);
t22 = qJ(1) + qJ(2);
t19 = cos(t22);
t14 = pkin(2) * t19;
t20 = qJ(3) + t22;
t17 = cos(t20);
t12 = pkin(3) * t17;
t26 = cos(qJ(1));
t21 = t26 * pkin(1);
t35 = t14 + t21;
t34 = m(4) + t38;
t15 = pkin(9) + t20;
t10 = sin(t15);
t11 = cos(t15);
t33 = t11 * pkin(4) + t10 * pkin(8) + t12;
t32 = t14 + t33;
t16 = sin(t20);
t30 = -t17 * mrSges(4,1) + t16 * mrSges(4,2) + t41 * t10 + t39 * t11;
t18 = sin(t22);
t29 = -t19 * mrSges(3,1) + t18 * mrSges(3,2) + t30;
t28 = t17 * mrSges(4,2) + (t38 * pkin(3) + mrSges(4,1)) * t16 + (m(6) * pkin(4) - t39) * t10 + (-m(6) * pkin(8) + t41) * t11;
t27 = mrSges(3,2) * t19 + (t34 * pkin(2) + mrSges(3,1)) * t18 + t28;
t24 = sin(qJ(1));
t1 = [(t24 * mrSges(2,2) - m(4) * t35 - m(5) * (t12 + t35) - m(6) * (t21 + t32) + (-m(3) * pkin(1) - mrSges(2,1)) * t26 + t29) * g(2) + (mrSges(2,2) * t26 + (mrSges(2,1) + (m(3) + t34) * pkin(1)) * t24 + t27) * g(1), (-m(4) * t14 - m(5) * (t12 + t14) - m(6) * t32 + t29) * g(2) + t27 * g(1), (-m(5) * t12 - m(6) * t33 + t30) * g(2) + t28 * g(1), -t38 * g(3), -g(3) * t40 + (g(1) * t11 + g(2) * t10) * (mrSges(6,1) * t23 + mrSges(6,2) * t25)];
taug = t1(:);
