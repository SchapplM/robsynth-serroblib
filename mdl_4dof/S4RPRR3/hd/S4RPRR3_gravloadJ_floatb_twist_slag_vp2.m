% Calculate Gravitation load on the joints for
% S4RPRR3
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (112->32), mult. (108->33), div. (0->0), fcn. (79->8), ass. (0->20)
t10 = sin(qJ(3));
t9 = qJ(3) + qJ(4);
t5 = sin(t9);
t6 = cos(t9);
t19 = t6 * mrSges(5,1) - t5 * mrSges(5,2);
t31 = -mrSges(4,2) * t10 + t19;
t28 = m(3) + m(4) + m(5);
t30 = pkin(1) * t28 + mrSges(2,1);
t8 = qJ(1) + pkin(7);
t3 = sin(t8);
t4 = cos(t8);
t29 = g(1) * t4 + g(2) * t3;
t12 = cos(qJ(3));
t27 = m(4) * pkin(2) + mrSges(4,1) * t12 + mrSges(3,1) + m(5) * (t12 * pkin(3) + pkin(2)) + t31;
t26 = -m(4) * pkin(5) + m(5) * (-pkin(6) - pkin(5)) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t20 = m(5) * pkin(3) + mrSges(4,1);
t18 = mrSges(5,1) * t5 + mrSges(5,2) * t6;
t13 = cos(qJ(1));
t11 = sin(qJ(1));
t1 = [(mrSges(2,2) * t11 - t30 * t13 + t26 * t3 - t27 * t4) * g(2) + (mrSges(2,2) * t13 + t30 * t11 + t26 * t4 + t27 * t3) * g(1), -t28 * g(3), (-t20 * t12 - t31) * g(3) + t29 * (mrSges(4,2) * t12 + t20 * t10 + t18), -g(3) * t19 + t29 * t18];
taug = t1(:);
