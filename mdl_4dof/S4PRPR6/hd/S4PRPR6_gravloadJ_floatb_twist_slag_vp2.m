% Calculate Gravitation load on the joints for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:20
% DurationCPUTime: 0.18s
% Computational Cost: add. (74->30), mult. (133->41), div. (0->0), fcn. (109->8), ass. (0->15)
t4 = pkin(7) + qJ(4);
t2 = sin(t4);
t3 = cos(t4);
t7 = cos(pkin(7));
t28 = mrSges(3,1) + m(5) * (pkin(3) * t7 + pkin(2)) + t3 * mrSges(5,1) - t2 * mrSges(5,2) + m(4) * pkin(2) + t7 * mrSges(4,1) - sin(pkin(7)) * mrSges(4,2);
t27 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) + m(5) * (-pkin(5) - qJ(3)) - mrSges(5,3);
t6 = sin(pkin(6));
t8 = cos(pkin(6));
t26 = g(1) * t8 + g(2) * t6;
t22 = m(4) + m(5);
t11 = cos(qJ(2));
t21 = t11 * t6;
t20 = t11 * t8;
t10 = sin(qJ(2));
t1 = [(-m(2) - m(3) - t22) * g(3), (-t28 * g(3) + t26 * t27) * t11 + (t27 * g(3) + t26 * t28) * t10, (g(3) * t11 - t10 * t26) * t22, -g(1) * ((-t2 * t20 + t6 * t3) * mrSges(5,1) + (-t6 * t2 - t3 * t20) * mrSges(5,2)) - g(2) * ((-t2 * t21 - t8 * t3) * mrSges(5,1) + (t8 * t2 - t3 * t21) * mrSges(5,2)) - g(3) * (-mrSges(5,1) * t2 - mrSges(5,2) * t3) * t10];
taug = t1(:);
