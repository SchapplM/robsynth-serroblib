% Calculate Gravitation load on the joints for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:06
% DurationCPUTime: 0.25s
% Computational Cost: add. (175->49), mult. (151->48), div. (0->0), fcn. (111->8), ass. (0->24)
t32 = m(3) + m(4);
t31 = -m(5) - m(6);
t9 = qJ(1) + pkin(7);
t4 = sin(t9);
t6 = cos(t9);
t30 = g(1) * t6 + g(2) * t4;
t11 = cos(pkin(8));
t8 = pkin(8) + qJ(4);
t3 = sin(t8);
t5 = cos(t8);
t22 = t5 * mrSges(5,1) - t3 * mrSges(5,2);
t29 = -mrSges(3,1) - m(4) * pkin(2) - t11 * mrSges(4,1) + sin(pkin(8)) * mrSges(4,2) - t22;
t28 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t13 = sin(qJ(1));
t25 = pkin(1) * t13;
t14 = cos(qJ(1));
t7 = t14 * pkin(1);
t24 = m(4) - t31;
t20 = t5 * mrSges(6,1) + t3 * mrSges(6,3);
t18 = pkin(4) * t5 + qJ(5) * t3;
t16 = m(6) * t18 + t20;
t12 = -pkin(6) - qJ(3);
t2 = pkin(3) * t11 + pkin(2);
t1 = [(-mrSges(2,1) * t14 + t13 * mrSges(2,2) - t32 * t7 + t31 * (-t12 * t4 + t6 * t2 + t7) + (-t16 + t29) * t6 + t28 * t4) * g(2) + (t13 * mrSges(2,1) + mrSges(2,2) * t14 + t32 * t25 + t31 * (-t12 * t6 - t25) + t28 * t6 + (m(5) * t2 - m(6) * (-t18 - t2) + t20 - t29) * t4) * g(1), (-m(3) - t24) * g(3), (-g(1) * t4 + g(2) * t6) * t24, (-t16 - t22) * g(3) + ((-m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3)) * t5 + (m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1)) * t3) * t30, (g(3) * t5 - t30 * t3) * m(6)];
taug = t1(:);
