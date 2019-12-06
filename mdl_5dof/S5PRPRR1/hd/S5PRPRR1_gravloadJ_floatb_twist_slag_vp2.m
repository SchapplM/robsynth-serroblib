% Calculate Gravitation load on the joints for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:38
% EndTime: 2019-12-05 15:42:39
% DurationCPUTime: 0.20s
% Computational Cost: add. (174->37), mult. (130->36), div. (0->0), fcn. (93->8), ass. (0->21)
t13 = pkin(9) + qJ(4);
t10 = qJ(5) + t13;
t3 = sin(t10);
t4 = cos(t10);
t26 = t4 * mrSges(6,1) - t3 * mrSges(6,2);
t6 = sin(t13);
t36 = -t6 * mrSges(5,2) + t26;
t14 = pkin(8) + qJ(2);
t7 = sin(t14);
t9 = cos(t14);
t35 = g(1) * t9 + g(2) * t7;
t16 = cos(pkin(9));
t5 = t16 * pkin(3) + pkin(2);
t8 = cos(t13);
t34 = mrSges(3,1) + m(4) * pkin(2) + t16 * mrSges(4,1) - sin(pkin(9)) * mrSges(4,2) + m(6) * (pkin(4) * t8 + t5) + m(5) * t5 + t8 * mrSges(5,1) + t36;
t17 = -pkin(6) - qJ(3);
t33 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) + m(6) * (-pkin(7) + t17) - mrSges(6,3) + m(5) * t17 - mrSges(5,3);
t28 = m(4) + m(5) + m(6);
t27 = m(6) * pkin(4) + mrSges(5,1);
t22 = mrSges(6,1) * t3 + mrSges(6,2) * t4;
t1 = [(-m(2) - m(3) - t28) * g(3), (t33 * t7 - t34 * t9) * g(2) + (t33 * t9 + t34 * t7) * g(1), (-g(1) * t7 + g(2) * t9) * t28, (-t27 * t8 - t36) * g(3) + t35 * (mrSges(5,2) * t8 + t27 * t6 + t22), -g(3) * t26 + t35 * t22];
taug = t1(:);
