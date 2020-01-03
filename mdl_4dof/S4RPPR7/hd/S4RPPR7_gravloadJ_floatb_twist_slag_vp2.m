% Calculate Gravitation load on the joints for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR7_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:36
% DurationCPUTime: 0.17s
% Computational Cost: add. (62->28), mult. (95->30), div. (0->0), fcn. (68->6), ass. (0->13)
t25 = m(4) + m(5);
t10 = sin(qJ(1));
t11 = cos(qJ(1));
t24 = -g(1) * t10 + g(2) * t11;
t23 = m(3) + t25;
t22 = mrSges(2,1) + mrSges(4,3) + mrSges(5,3) - mrSges(3,2);
t6 = pkin(6) + qJ(4);
t1 = sin(t6);
t2 = cos(t6);
t14 = -t1 * mrSges(5,1) - t2 * mrSges(5,2);
t21 = -cos(pkin(6)) * mrSges(4,2) + mrSges(2,2) - mrSges(3,3) + t14 + (-m(5) * pkin(3) - mrSges(4,1)) * sin(pkin(6));
t9 = -pkin(5) - qJ(3);
t3 = [(-t23 * (t11 * pkin(1) + t10 * qJ(2)) + (-m(4) * qJ(3) + m(5) * t9 - t22) * t11 + t21 * t10) * g(2) + ((m(3) * pkin(1) - m(4) * (-pkin(1) - qJ(3)) - m(5) * (-pkin(1) + t9) + t22) * t10 + (-t23 * qJ(2) + t21) * t11) * g(1), t24 * t23, t25 * (-g(1) * t11 - g(2) * t10), -g(3) * t14 + t24 * (mrSges(5,1) * t2 - mrSges(5,2) * t1)];
taug = t3(:);
