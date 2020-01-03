% Calculate Gravitation load on the joints for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:05
% EndTime: 2019-12-31 16:55:06
% DurationCPUTime: 0.19s
% Computational Cost: add. (79->34), mult. (121->39), div. (0->0), fcn. (89->6), ass. (0->22)
t11 = cos(qJ(3));
t8 = qJ(3) + qJ(4);
t3 = sin(t8);
t4 = cos(t8);
t16 = mrSges(5,1) * t3 + t4 * mrSges(5,2);
t17 = m(5) * pkin(3) + mrSges(4,1);
t9 = sin(qJ(3));
t33 = mrSges(4,2) * t11 + t17 * t9 + t16;
t30 = m(3) + m(4) + m(5);
t29 = mrSges(4,2) * t9 - t17 * t11;
t28 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t27 = mrSges(2,2) - mrSges(3,3) - t33;
t24 = mrSges(5,1) * t4;
t22 = mrSges(5,2) * t3;
t10 = sin(qJ(1));
t21 = g(1) * t10;
t12 = cos(qJ(1));
t20 = g(2) * t12;
t13 = -pkin(6) - pkin(5);
t2 = t12 * t22;
t1 = t10 * t24;
t5 = [(-t30 * (t12 * pkin(1) + t10 * qJ(2)) + (-m(4) * pkin(5) + m(5) * t13 - t28) * t12 + t27 * t10) * g(2) + ((m(3) * pkin(1) - m(4) * (-pkin(1) - pkin(5)) - m(5) * (-pkin(1) + t13) + t28) * t10 + (-t30 * qJ(2) + t27) * t12) * g(1), (t20 - t21) * t30, -g(1) * t1 - g(2) * t2 + t33 * g(3) + (t24 - t29) * t20 + (t22 + t29) * t21, -g(1) * (-t10 * t22 + t1) - g(2) * (-t12 * t24 + t2) + g(3) * t16];
taug = t5(:);
