% Calculate Gravitation load on the joints for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S4RPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:34
% EndTime: 2018-11-14 13:50:34
% DurationCPUTime: 0.12s
% Computational Cost: add. (109->32), mult. (67->33), div. (0->0), fcn. (42->8), ass. (0->20)
t13 = qJ(1) + pkin(7);
t11 = qJ(3) + t13;
t7 = cos(t11);
t23 = pkin(3) * t7;
t22 = m(4) + m(5);
t10 = cos(t13);
t15 = cos(qJ(1));
t21 = t15 * pkin(1) + pkin(2) * t10;
t20 = m(3) + t22;
t8 = qJ(4) + t11;
t3 = sin(t8);
t4 = cos(t8);
t19 = t4 * mrSges(5,1) - t3 * mrSges(5,2);
t18 = -t3 * mrSges(5,1) - t4 * mrSges(5,2);
t6 = sin(t11);
t17 = -t7 * mrSges(4,1) + t6 * mrSges(4,2) - t19;
t16 = t7 * mrSges(4,2) + (m(5) * pkin(3) + mrSges(4,1)) * t6 - t18;
t14 = sin(qJ(1));
t9 = sin(t13);
t1 = [(t14 * mrSges(2,2) - t10 * mrSges(3,1) + t9 * mrSges(3,2) - m(4) * t21 - m(5) * (t21 + t23) + (-m(3) * pkin(1) - mrSges(2,1)) * t15 + t17) * g(2) + (mrSges(2,2) * t15 + mrSges(3,2) * t10 + (pkin(2) * t22 + mrSges(3,1)) * t9 + (pkin(1) * t20 + mrSges(2,1)) * t14 + t16) * g(1), -t20 * g(3) (-m(5) * t23 + t17) * g(2) + t16 * g(1), -g(1) * t18 - g(2) * t19];
taug  = t1(:);
