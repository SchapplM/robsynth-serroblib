% Calculate Gravitation load on the joints for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:11
% DurationCPUTime: 0.23s
% Computational Cost: add. (126->39), mult. (162->40), div. (0->0), fcn. (122->6), ass. (0->21)
t12 = sin(qJ(2));
t28 = mrSges(4,2) + mrSges(5,2);
t29 = mrSges(4,1) + mrSges(5,1);
t11 = qJ(2) + qJ(3);
t7 = sin(t11);
t8 = cos(t11);
t22 = t28 * t7 - t29 * t8;
t38 = t12 * mrSges(3,2) + t22;
t13 = sin(qJ(1));
t15 = cos(qJ(1));
t35 = g(1) * t15 + g(2) * t13;
t14 = cos(qJ(2));
t5 = pkin(3) * t8;
t9 = t14 * pkin(2);
t6 = t9 + pkin(1);
t34 = m(3) * pkin(1) + t14 * mrSges(3,1) + mrSges(2,1) + m(5) * (t5 + t6) + m(4) * t6 - t38;
t16 = -pkin(6) - pkin(5);
t33 = mrSges(2,2) + m(5) * (-qJ(4) + t16) - mrSges(5,3) + m(4) * t16 - mrSges(4,3) - m(3) * pkin(5) - mrSges(3,3);
t26 = m(4) * pkin(2) + mrSges(3,1);
t25 = t28 * t8;
t1 = [(t33 * t13 - t34 * t15) * g(2) + (t34 * t13 + t33 * t15) * g(1), (-m(5) * (t5 + t9) - t26 * t14 + t38) * g(3) + t35 * (-m(5) * (-pkin(2) * t12 - pkin(3) * t7) + mrSges(3,2) * t14 + t26 * t12 + t29 * t7 + t25), (-m(5) * t5 + t22) * g(3) + t35 * (t25 + (m(5) * pkin(3) + t29) * t7), (-g(1) * t13 + g(2) * t15) * m(5)];
taug = t1(:);
