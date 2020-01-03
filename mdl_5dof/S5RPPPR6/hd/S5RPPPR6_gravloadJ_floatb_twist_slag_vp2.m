% Calculate Gravitation load on the joints for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:37
% EndTime: 2019-12-31 17:47:38
% DurationCPUTime: 0.42s
% Computational Cost: add. (120->65), mult. (261->86), div. (0->0), fcn. (245->8), ass. (0->30)
t50 = m(5) + m(6);
t32 = m(4) + t50;
t22 = sin(qJ(1));
t24 = cos(qJ(1));
t49 = -g(1) * t24 - g(2) * t22;
t18 = sin(pkin(7));
t20 = cos(pkin(7));
t48 = -mrSges(2,1) + (-mrSges(3,1) + mrSges(4,2)) * t20 + (mrSges(3,2) - mrSges(4,3)) * t18;
t47 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t46 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t17 = sin(pkin(8));
t42 = t17 * t22;
t21 = sin(qJ(5));
t41 = t20 * t21;
t23 = cos(qJ(5));
t40 = t20 * t23;
t19 = cos(pkin(8));
t39 = t22 * t19;
t38 = t24 * t18;
t37 = t24 * t20;
t34 = t24 * pkin(1) + t22 * qJ(2);
t33 = qJ(3) * t18;
t29 = pkin(2) * t37 + t24 * t33 + t34;
t25 = t22 * pkin(3) + qJ(4) * t37 + t29;
t13 = t24 * qJ(2);
t6 = -t18 * t42 + t19 * t24;
t4 = t17 * t38 + t39;
t2 = t21 * t37 + t4 * t23;
t1 = -t4 * t21 + t23 * t37;
t3 = [(-m(3) * t34 - m(4) * t29 - m(5) * t25 - t4 * mrSges(5,1) - mrSges(5,3) * t37 - m(6) * (pkin(4) * t4 + t25) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t46 * (-t19 * t38 + t42) + t48 * t24 + t47 * t22) * g(2) + ((-m(6) * pkin(4) - t23 * mrSges(6,1) + t21 * mrSges(6,2) - mrSges(5,1)) * t6 - t50 * (t24 * pkin(3) + t13) + t46 * (t17 * t24 + t18 * t39) + (-m(3) - m(4)) * t13 + t47 * t24 + (m(3) * pkin(1) - t32 * (-pkin(1) - t33) + (m(4) * pkin(2) + t21 * mrSges(6,1) + t23 * mrSges(6,2) + mrSges(5,3) - t50 * (-pkin(2) - qJ(4))) * t20 - t48) * t22) * g(1), (-g(1) * t22 + g(2) * t24) * (m(3) + t32), (g(3) * t20 + t18 * t49) * t32, (-g(3) * t18 + t20 * t49) * t50, -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((t21 * t6 + t22 * t40) * mrSges(6,1) + (-t22 * t41 + t23 * t6) * mrSges(6,2)) - g(3) * ((t17 * t41 + t18 * t23) * mrSges(6,1) + (t17 * t40 - t18 * t21) * mrSges(6,2))];
taug = t3(:);
