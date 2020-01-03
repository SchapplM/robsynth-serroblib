% Calculate Gravitation load on the joints for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (89->29), mult. (82->30), div. (0->0), fcn. (58->8), ass. (0->16)
t23 = m(4) + m(5);
t22 = m(3) + t23;
t24 = pkin(1) * t22 + mrSges(2,1);
t10 = cos(pkin(7));
t7 = pkin(7) + qJ(4);
t2 = sin(t7);
t4 = cos(t7);
t17 = t4 * mrSges(5,1) - t2 * mrSges(5,2);
t21 = mrSges(3,1) + m(4) * pkin(2) + t10 * mrSges(4,1) - sin(pkin(7)) * mrSges(4,2) + m(5) * (pkin(3) * t10 + pkin(2)) + t17;
t20 = -m(4) * qJ(3) + m(5) * (-pkin(5) - qJ(3)) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t13 = cos(qJ(1));
t12 = sin(qJ(1));
t8 = qJ(1) + pkin(6);
t5 = cos(t8);
t3 = sin(t8);
t1 = [(t12 * mrSges(2,2) - t24 * t13 + t20 * t3 - t21 * t5) * g(2) + (mrSges(2,2) * t13 + t24 * t12 + t20 * t5 + t21 * t3) * g(1), -t22 * g(3), t23 * (-g(1) * t3 + g(2) * t5), -g(3) * t17 + (g(1) * t5 + g(2) * t3) * (mrSges(5,1) * t2 + mrSges(5,2) * t4)];
taug = t1(:);
