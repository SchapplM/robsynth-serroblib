% Calculate Gravitation load on the joints for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:28
% DurationCPUTime: 0.16s
% Computational Cost: add. (120->35), mult. (111->39), div. (0->0), fcn. (80->6), ass. (0->21)
t43 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t41 = mrSges(5,1) * t18 + mrSges(5,2) * t20;
t42 = -t41 + mrSges(3,2) - mrSges(4,3);
t17 = qJ(1) + qJ(2);
t14 = sin(t17);
t15 = cos(t17);
t40 = -g(1) * t14 + g(2) * t15;
t39 = t42 * t15 + (-m(5) * (-pkin(2) - pkin(6)) - t43) * t14;
t38 = t42 * t14 + t43 * t15;
t19 = sin(qJ(1));
t37 = pkin(1) * t19;
t21 = cos(qJ(1));
t16 = t21 * pkin(1);
t34 = t15 * pkin(2) + t14 * qJ(3);
t29 = t16 + t34;
t6 = t15 * qJ(3);
t28 = -pkin(2) * t14 + t6;
t12 = t15 * pkin(6);
t1 = [(-mrSges(2,1) * t21 + t19 * mrSges(2,2) - m(3) * t16 - m(4) * t29 - m(5) * (t12 + t29) + t38) * g(2) + (t19 * mrSges(2,1) + mrSges(2,2) * t21 + m(3) * t37 - m(4) * (t28 - t37) - m(5) * (t6 - t37) + t39) * g(1), (-m(4) * t34 - m(5) * (t12 + t34) + t38) * g(2) + (-m(4) * t28 - m(5) * t6 + t39) * g(1), (m(4) + m(5)) * t40, g(3) * t41 + t40 * (mrSges(5,1) * t20 - mrSges(5,2) * t18)];
taug = t1(:);
