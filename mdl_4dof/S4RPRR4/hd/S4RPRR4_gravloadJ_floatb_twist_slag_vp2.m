% Calculate Gravitation load on the joints for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:14
% EndTime: 2019-12-31 16:50:15
% DurationCPUTime: 0.21s
% Computational Cost: add. (126->49), mult. (150->64), div. (0->0), fcn. (129->8), ass. (0->25)
t35 = -m(4) - m(5);
t34 = mrSges(3,2) - mrSges(4,3);
t13 = sin(qJ(3));
t16 = cos(qJ(3));
t21 = mrSges(4,1) * t16 - mrSges(4,2) * t13;
t32 = t13 * mrSges(5,3) + mrSges(3,1) + t21;
t14 = sin(qJ(1));
t29 = pkin(1) * t14;
t17 = cos(qJ(1));
t10 = t17 * pkin(1);
t12 = sin(qJ(4));
t27 = t16 * t12;
t15 = cos(qJ(4));
t26 = t16 * t15;
t24 = m(5) * pkin(6) + mrSges(5,3);
t22 = pkin(3) * t16 + pkin(6) * t13;
t19 = m(5) * pkin(3) + mrSges(5,1) * t15 - mrSges(5,2) * t12;
t11 = qJ(1) + pkin(7);
t9 = cos(t11);
t8 = sin(t11);
t4 = t12 * t8 + t9 * t26;
t3 = t15 * t8 - t9 * t27;
t2 = t12 * t9 - t8 * t26;
t1 = t15 * t9 + t8 * t27;
t5 = [(-m(3) * t10 - mrSges(2,1) * t17 - t4 * mrSges(5,1) + t14 * mrSges(2,2) - t3 * mrSges(5,2) + t34 * t8 + t35 * (t9 * pkin(2) + t8 * pkin(5) + t10) + (-m(5) * t22 - t32) * t9) * g(2) + (m(3) * t29 + t14 * mrSges(2,1) - t2 * mrSges(5,1) + mrSges(2,2) * t17 - t1 * mrSges(5,2) + t34 * t9 + t35 * (t9 * pkin(5) - t29) + (m(4) * pkin(2) - m(5) * (-pkin(2) - t22) + t32) * t8) * g(1), (-m(3) + t35) * g(3), (-t24 * t13 - t19 * t16 - t21) * g(3) + ((mrSges(4,2) - t24) * t16 + (mrSges(4,1) + t19) * t13) * (g(1) * t9 + g(2) * t8), -g(1) * (mrSges(5,1) * t3 - mrSges(5,2) * t4) - g(2) * (-mrSges(5,1) * t1 + mrSges(5,2) * t2) - g(3) * (-mrSges(5,1) * t12 - mrSges(5,2) * t15) * t13];
taug = t5(:);
