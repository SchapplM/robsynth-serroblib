% Calculate Gravitation load on the joints for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [4x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:40
% EndTime: 2019-12-31 16:34:41
% DurationCPUTime: 0.25s
% Computational Cost: add. (103->42), mult. (178->58), div. (0->0), fcn. (159->8), ass. (0->21)
t11 = sin(qJ(3));
t13 = cos(qJ(3));
t8 = qJ(3) + qJ(4);
t6 = sin(t8);
t7 = cos(t8);
t39 = m(4) * pkin(2) + t13 * mrSges(4,1) - t11 * mrSges(4,2) + mrSges(3,1) + m(5) * (pkin(3) * t13 + pkin(2)) + t7 * mrSges(5,1) - t6 * mrSges(5,2);
t38 = mrSges(3,2) + m(5) * (-pkin(6) - pkin(5)) - mrSges(5,3) - m(4) * pkin(5) - mrSges(4,3);
t37 = m(5) * pkin(3) + mrSges(4,1);
t10 = cos(pkin(7));
t14 = cos(qJ(2));
t9 = sin(pkin(7));
t29 = t14 * t9;
t33 = (-t10 * t7 - t6 * t29) * mrSges(5,1) + (t10 * t6 - t7 * t29) * mrSges(5,2);
t28 = t10 * t14;
t32 = (-t6 * t28 + t7 * t9) * mrSges(5,1) + (-t7 * t28 - t6 * t9) * mrSges(5,2);
t12 = sin(qJ(2));
t30 = g(3) * t12;
t27 = t11 * t14;
t26 = t13 * t14;
t23 = -mrSges(5,1) * t6 - mrSges(5,2) * t7;
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), (t38 * t12 - t39 * t14) * g(3) + (g(1) * t10 + g(2) * t9) * (t39 * t12 + t38 * t14), (mrSges(4,2) * t13 + t37 * t11 - t23) * t30 + (-(t10 * t11 - t9 * t26) * mrSges(4,2) - t33 - t37 * (-t10 * t13 - t9 * t27)) * g(2) + (-(-t10 * t26 - t11 * t9) * mrSges(4,2) - t32 - t37 * (-t10 * t27 + t13 * t9)) * g(1), -g(1) * t32 - g(2) * t33 - t23 * t30];
taug = t1(:);
