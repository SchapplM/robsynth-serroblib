% Calculate Gravitation load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:28
% EndTime: 2022-01-23 09:34:29
% DurationCPUTime: 0.19s
% Computational Cost: add. (265->47), mult. (143->48), div. (0->0), fcn. (102->10), ass. (0->28)
t41 = mrSges(5,2) - mrSges(6,3);
t22 = sin(qJ(5));
t24 = cos(qJ(5));
t40 = mrSges(6,1) * t24 - mrSges(6,2) * t22;
t39 = -mrSges(5,1) - t40;
t38 = m(5) + m(6);
t21 = qJ(1) + pkin(9);
t19 = qJ(3) + t21;
t16 = qJ(4) + t19;
t11 = sin(t16);
t12 = cos(t16);
t37 = t12 * pkin(4) + t11 * pkin(8);
t15 = cos(t19);
t10 = pkin(3) * t15;
t18 = cos(t21);
t25 = cos(qJ(1));
t34 = t25 * pkin(1) + pkin(2) * t18;
t33 = m(4) + t38;
t32 = t10 + t34;
t31 = m(3) + t33;
t29 = t41 * t11 + t39 * t12;
t14 = sin(t19);
t28 = -t15 * mrSges(4,1) + t14 * mrSges(4,2) + t29;
t27 = (m(6) * pkin(4) - t39) * t11 + (-m(6) * pkin(8) + t41) * t12;
t26 = t15 * mrSges(4,2) + (t38 * pkin(3) + mrSges(4,1)) * t14 + t27;
t23 = sin(qJ(1));
t17 = sin(t21);
t1 = [(t23 * mrSges(2,2) - t18 * mrSges(3,1) + t17 * mrSges(3,2) - m(4) * t34 - m(5) * t32 - m(6) * (t32 + t37) + (-m(3) * pkin(1) - mrSges(2,1)) * t25 + t28) * g(2) + (mrSges(2,2) * t25 + mrSges(3,2) * t18 + (t33 * pkin(2) + mrSges(3,1)) * t17 + (t31 * pkin(1) + mrSges(2,1)) * t23 + t26) * g(1), -t31 * g(3), (-m(5) * t10 - m(6) * (t10 + t37) + t28) * g(2) + t26 * g(1), (-m(6) * t37 + t29) * g(2) + t27 * g(1), -g(3) * t40 + (g(1) * t12 + g(2) * t11) * (mrSges(6,1) * t22 + mrSges(6,2) * t24)];
taug = t1(:);
