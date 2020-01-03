% Calculate Gravitation load on the joints for
% S4PRPR5
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:22:59
% DurationCPUTime: 0.20s
% Computational Cost: add. (74->28), mult. (111->39), div. (0->0), fcn. (89->8), ass. (0->18)
t24 = m(4) + m(5);
t27 = t24 * pkin(2) + mrSges(3,1);
t6 = sin(qJ(4));
t8 = cos(qJ(4));
t26 = m(5) * pkin(3) + t8 * mrSges(5,1) - t6 * mrSges(5,2) + mrSges(4,1);
t25 = -m(5) * pkin(5) + mrSges(4,2) - mrSges(5,3);
t4 = sin(pkin(6));
t18 = t4 * t6;
t17 = t4 * t8;
t5 = cos(pkin(6));
t16 = t5 * t6;
t15 = t5 * t8;
t9 = cos(qJ(2));
t7 = sin(qJ(2));
t3 = qJ(2) + pkin(7);
t2 = cos(t3);
t1 = sin(t3);
t10 = [(-m(2) - m(3) - t24) * g(3), (t7 * mrSges(3,2) + t25 * t1 - t26 * t2 - t27 * t9) * g(3) + (g(1) * t5 + g(2) * t4) * (mrSges(3,2) * t9 + t26 * t1 + t25 * t2 + t27 * t7), t24 * (-g(1) * t4 + g(2) * t5), -g(1) * ((-t2 * t16 + t17) * mrSges(5,1) + (-t2 * t15 - t18) * mrSges(5,2)) - g(2) * ((-t2 * t18 - t15) * mrSges(5,1) + (-t2 * t17 + t16) * mrSges(5,2)) - g(3) * (-mrSges(5,1) * t6 - mrSges(5,2) * t8) * t1];
taug = t10(:);
