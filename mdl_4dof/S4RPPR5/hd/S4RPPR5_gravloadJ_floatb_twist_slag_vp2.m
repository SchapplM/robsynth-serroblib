% Calculate Gravitation load on the joints for
% S4RPPR5
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RPPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (66->33), mult. (119->39), div. (0->0), fcn. (112->6), ass. (0->17)
t23 = m(4) + m(5);
t27 = -mrSges(2,1) - mrSges(3,1);
t10 = sin(qJ(4));
t12 = cos(qJ(4));
t16 = -t12 * mrSges(5,1) + t10 * mrSges(5,2);
t26 = m(5) * pkin(3) + mrSges(4,1) - t16;
t25 = mrSges(2,2) - mrSges(3,3);
t24 = -m(5) * pkin(5) + mrSges(4,2) - mrSges(5,3);
t11 = sin(qJ(1));
t13 = cos(qJ(1));
t22 = t13 * pkin(1) + t11 * qJ(2);
t21 = cos(pkin(6));
t20 = sin(pkin(6));
t7 = t13 * qJ(2);
t2 = -t11 * t21 + t13 * t20;
t1 = -t11 * t20 - t13 * t21;
t3 = [(-m(3) * t22 - t23 * (t13 * pkin(2) + t22) + t24 * t2 + t27 * t13 + t25 * t11 + t26 * t1) * g(2) + (-m(4) * t7 - t26 * t2 + (-m(3) - m(5)) * (-pkin(1) * t11 + t7) + t25 * t13 + (-m(4) * (-pkin(1) - pkin(2)) + m(5) * pkin(2) - t27) * t11 + t24 * t1) * g(1), (-g(1) * t11 + g(2) * t13) * (m(3) + t23), t23 * g(3), -g(3) * t16 + (-g(1) * t1 - g(2) * t2) * (mrSges(5,1) * t10 + mrSges(5,2) * t12)];
taug = t3(:);
