% Calculate Gravitation load on the joints for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:03
% EndTime: 2019-12-05 16:19:03
% DurationCPUTime: 0.21s
% Computational Cost: add. (188->44), mult. (150->42), div. (0->0), fcn. (109->8), ass. (0->23)
t16 = qJ(3) + pkin(9);
t11 = cos(t16);
t18 = sin(qJ(3));
t12 = qJ(5) + t16;
t5 = sin(t12);
t6 = cos(t12);
t29 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t9 = sin(t16);
t41 = t11 * mrSges(5,1) - t18 * mrSges(4,2) - t9 * mrSges(5,2) + t29;
t40 = m(5) + m(6);
t15 = pkin(8) + qJ(2);
t10 = cos(t15);
t8 = sin(t15);
t39 = g(1) * t10 + g(2) * t8;
t19 = cos(qJ(3));
t13 = t19 * pkin(3);
t33 = pkin(4) * t11 + t13;
t38 = m(4) * pkin(2) + t19 * mrSges(4,1) + mrSges(3,1) + m(5) * (t13 + pkin(2)) + m(6) * (pkin(2) + t33) + t41;
t17 = -qJ(4) - pkin(6);
t37 = mrSges(3,2) + m(6) * (-pkin(7) + t17) - mrSges(6,3) + m(5) * t17 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t30 = m(5) * pkin(3) + mrSges(4,1);
t25 = mrSges(6,1) * t5 + mrSges(6,2) * t6;
t1 = [(-m(2) - m(3) - m(4) - t40) * g(3), (-t38 * t10 + t37 * t8) * g(2) + (t37 * t10 + t38 * t8) * g(1), (-m(6) * t33 - t30 * t19 - t41) * g(3) + t39 * (-m(6) * (-pkin(3) * t18 - pkin(4) * t9) + mrSges(5,1) * t9 + mrSges(4,2) * t19 + mrSges(5,2) * t11 + t30 * t18 + t25), t40 * (-g(1) * t8 + g(2) * t10), -g(3) * t29 + t39 * t25];
taug = t1(:);
