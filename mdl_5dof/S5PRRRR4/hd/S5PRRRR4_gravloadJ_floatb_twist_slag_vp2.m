% Calculate Gravitation load on the joints for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:40
% EndTime: 2019-12-05 17:07:41
% DurationCPUTime: 0.19s
% Computational Cost: add. (241->43), mult. (154->45), div. (0->0), fcn. (113->8), ass. (0->27)
t25 = qJ(4) + qJ(5);
t22 = sin(t25);
t23 = cos(t25);
t48 = t23 * mrSges(6,1) - t22 * mrSges(6,2);
t26 = sin(qJ(4));
t47 = t26 * mrSges(5,2) - t48;
t27 = cos(qJ(4));
t46 = -t27 * mrSges(5,1) - mrSges(4,1) + t47;
t44 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t24 = pkin(9) + qJ(2);
t21 = qJ(3) + t24;
t16 = sin(t21);
t17 = cos(t21);
t43 = g(1) * t17 + g(2) * t16;
t37 = t17 * pkin(3) + t16 * pkin(7);
t36 = m(4) + m(5) + m(6);
t35 = m(6) * pkin(4) + mrSges(5,1);
t18 = pkin(4) * t27 + pkin(3);
t28 = -pkin(8) - pkin(7);
t34 = -t16 * t28 + t17 * t18;
t33 = mrSges(6,1) * t22 + mrSges(6,2) * t23;
t30 = t44 * t16 + t46 * t17;
t29 = (m(5) * pkin(3) + m(6) * t18 - t46) * t16 + (-m(5) * pkin(7) + m(6) * t28 + t44) * t17;
t20 = cos(t24);
t19 = sin(t24);
t14 = pkin(2) * t20;
t1 = [(-m(2) - m(3) - t36) * g(3), (mrSges(3,2) * t19 - m(5) * (t14 + t37) - m(6) * (t14 + t34) + (-m(4) * pkin(2) - mrSges(3,1)) * t20 + t30) * g(2) + (mrSges(3,2) * t20 + (t36 * pkin(2) + mrSges(3,1)) * t19 + t29) * g(1), (-m(5) * t37 - m(6) * t34 + t30) * g(2) + t29 * g(1), (-t35 * t27 + t47) * g(3) + t43 * (mrSges(5,2) * t27 + t35 * t26 + t33), -g(3) * t48 + t43 * t33];
taug = t1(:);
