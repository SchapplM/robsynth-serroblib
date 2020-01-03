% Calculate Gravitation load on the joints for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (70->42), mult. (151->54), div. (0->0), fcn. (132->6), ass. (0->22)
t27 = m(4) + m(5);
t30 = m(3) + t27;
t10 = sin(pkin(6));
t11 = cos(pkin(6));
t29 = -mrSges(2,1) + (-mrSges(3,1) - mrSges(4,1)) * t11 + (mrSges(3,2) - mrSges(4,3)) * t10;
t28 = m(5) * pkin(5) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3) + mrSges(5,3);
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t26 = t15 * pkin(1) + t13 * qJ(2);
t25 = t15 * t11;
t24 = qJ(3) * t10;
t22 = pkin(2) * t25 + t15 * t24 + t26;
t21 = -pkin(1) - t24;
t12 = sin(qJ(4));
t14 = cos(qJ(4));
t17 = t10 * t14 - t11 * t12;
t16 = t10 * t12 + t11 * t14;
t4 = t16 * t15;
t3 = t17 * t15;
t2 = t16 * t13;
t1 = t17 * t13;
t5 = [(-m(3) * t26 - m(4) * t22 - m(5) * (pkin(3) * t25 + t22) - t4 * mrSges(5,1) - t3 * mrSges(5,2) + t29 * t15 + t28 * t13) * g(2) + (t2 * mrSges(5,1) + t1 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * (-pkin(2) * t11 + t21) - m(5) * ((-pkin(2) - pkin(3)) * t11 + t21) - t29) * t13 + (-t30 * qJ(2) + t28) * t15) * g(1), (-g(1) * t13 + g(2) * t15) * t30, (g(3) * t11 + t10 * (-g(1) * t15 - g(2) * t13)) * t27, -g(1) * (mrSges(5,1) * t3 - mrSges(5,2) * t4) - g(2) * (mrSges(5,1) * t1 - mrSges(5,2) * t2) - g(3) * (-t16 * mrSges(5,1) - t17 * mrSges(5,2))];
taug = t5(:);
