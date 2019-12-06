% Calculate Gravitation load on the joints for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:59
% EndTime: 2019-12-05 17:38:00
% DurationCPUTime: 0.25s
% Computational Cost: add. (100->37), mult. (155->42), div. (0->0), fcn. (111->6), ass. (0->21)
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t9 = qJ(4) + qJ(5);
t3 = sin(t9);
t4 = cos(t9);
t17 = t3 * mrSges(6,1) + t4 * mrSges(6,2);
t18 = m(6) * pkin(4) + mrSges(5,1);
t37 = t12 * mrSges(5,2) + t18 * t10 + t17;
t21 = -m(4) - m(5) - m(6);
t33 = m(3) - t21;
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t32 = g(1) * t13 + g(2) * t11;
t31 = mrSges(2,1) - mrSges(3,2) + mrSges(4,3) + t37;
t30 = m(5) * pkin(6) - m(6) * (-pkin(7) - pkin(6)) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3) + mrSges(6,3);
t29 = t13 * pkin(1) + t11 * qJ(2);
t28 = mrSges(6,1) * t4;
t27 = mrSges(6,2) * t3;
t2 = t13 * t28;
t1 = t11 * t28;
t5 = [(-m(3) * t29 + t21 * (t13 * qJ(3) + t29) - t31 * t13 + t30 * t11) * g(2) + ((m(3) * pkin(1) + t21 * (-pkin(1) - qJ(3)) + t31) * t11 + (-t33 * qJ(2) + t30) * t13) * g(1), (-g(1) * t11 + g(2) * t13) * t33, t32 * t21, -g(1) * t2 - g(2) * t1 + t37 * g(3) + t32 * (mrSges(5,2) * t10 - t18 * t12 + t27), -g(1) * (-t13 * t27 + t2) - g(2) * (-t11 * t27 + t1) + g(3) * t17];
taug = t5(:);
