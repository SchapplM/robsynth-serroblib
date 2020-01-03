% Calculate Gravitation load on the joints for
% S4RPPR4
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (75->27), mult. (76->32), div. (0->0), fcn. (52->6), ass. (0->15)
t23 = m(4) + m(5);
t7 = qJ(1) + pkin(6);
t4 = sin(t7);
t5 = cos(t7);
t22 = -g(1) * t4 + g(2) * t5;
t21 = mrSges(3,1) + mrSges(5,3) - mrSges(4,2);
t10 = cos(qJ(4));
t8 = sin(qJ(4));
t13 = t8 * mrSges(5,1) + t10 * mrSges(5,2);
t20 = mrSges(3,2) - mrSges(4,3) - t13;
t9 = sin(qJ(1));
t19 = pkin(1) * t9;
t11 = cos(qJ(1));
t6 = t11 * pkin(1);
t1 = [(-m(3) * t6 - mrSges(2,1) * t11 + t9 * mrSges(2,2) - t23 * (t5 * pkin(2) + t4 * qJ(3) + t6) + (-m(5) * pkin(5) - t21) * t5 + t20 * t4) * g(2) + (m(3) * t19 + t9 * mrSges(2,1) + mrSges(2,2) * t11 - t23 * (t5 * qJ(3) - t19) + t20 * t5 + (m(4) * pkin(2) - m(5) * (-pkin(2) - pkin(5)) + t21) * t4) * g(1), (-m(3) - t23) * g(3), t23 * t22, g(3) * t13 + t22 * (mrSges(5,1) * t10 - mrSges(5,2) * t8)];
taug = t1(:);
