% Calculate Gravitation load on the joints for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:18
% EndTime: 2018-11-14 13:53:18
% DurationCPUTime: 0.13s
% Computational Cost: add. (124->34), mult. (78->34), div. (0->0), fcn. (50->8), ass. (0->20)
t23 = m(4) + m(5);
t15 = qJ(1) + qJ(2);
t13 = cos(t15);
t9 = pkin(2) * t13;
t17 = cos(qJ(1));
t22 = t17 * pkin(1) + t9;
t11 = pkin(7) + t15;
t10 = qJ(4) + t11;
t4 = sin(t10);
t5 = cos(t10);
t21 = t5 * mrSges(5,1) - t4 * mrSges(5,2);
t20 = -t4 * mrSges(5,1) - t5 * mrSges(5,2);
t12 = sin(t15);
t7 = sin(t11);
t8 = cos(t11);
t19 = -t13 * mrSges(3,1) - t8 * mrSges(4,1) + t12 * mrSges(3,2) + t7 * mrSges(4,2) - t21;
t18 = t13 * mrSges(3,2) + t8 * mrSges(4,2) + (m(5) * pkin(3) + mrSges(4,1)) * t7 + (t23 * pkin(2) + mrSges(3,1)) * t12 - t20;
t16 = sin(qJ(1));
t3 = pkin(3) * t8;
t1 = [(t16 * mrSges(2,2) - m(4) * t22 - m(5) * (t3 + t22) + (-m(3) * pkin(1) - mrSges(2,1)) * t17 + t19) * g(2) + (mrSges(2,2) * t17 + (mrSges(2,1) + (m(3) + t23) * pkin(1)) * t16 + t18) * g(1) (-m(4) * t9 - m(5) * (t3 + t9) + t19) * g(2) + t18 * g(1), -t23 * g(3), -g(1) * t20 - g(2) * t21];
taug  = t1(:);
