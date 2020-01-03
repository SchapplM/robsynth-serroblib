% Calculate Gravitation load on the joints for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:23
% DurationCPUTime: 0.18s
% Computational Cost: add. (126->40), mult. (145->41), div. (0->0), fcn. (109->8), ass. (0->21)
t15 = sin(qJ(2));
t13 = qJ(2) + pkin(7);
t10 = qJ(4) + t13;
t5 = sin(t10);
t6 = cos(t10);
t28 = t6 * mrSges(5,1) - t5 * mrSges(5,2);
t8 = sin(t13);
t9 = cos(t13);
t39 = t9 * mrSges(4,1) - t15 * mrSges(3,2) - t8 * mrSges(4,2) + t28;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t38 = g(1) * t18 + g(2) * t16;
t17 = cos(qJ(2));
t11 = t17 * pkin(2);
t32 = pkin(3) * t9 + t11;
t37 = m(3) * pkin(1) + t17 * mrSges(3,1) + mrSges(2,1) + m(5) * (pkin(1) + t32) + m(4) * (t11 + pkin(1)) + t39;
t14 = -qJ(3) - pkin(5);
t36 = mrSges(2,2) + m(5) * (-pkin(6) + t14) - mrSges(5,3) + m(4) * t14 - mrSges(4,3) - m(3) * pkin(5) - mrSges(3,3);
t29 = m(4) * pkin(2) + mrSges(3,1);
t24 = mrSges(5,1) * t5 + mrSges(5,2) * t6;
t1 = [(t36 * t16 - t37 * t18) * g(2) + (t37 * t16 + t36 * t18) * g(1), (-m(5) * t32 - t29 * t17 - t39) * g(3) + t38 * (-m(5) * (-pkin(2) * t15 - pkin(3) * t8) + mrSges(4,1) * t8 + mrSges(3,2) * t17 + mrSges(4,2) * t9 + t29 * t15 + t24), (m(4) + m(5)) * (-g(1) * t16 + g(2) * t18), -g(3) * t28 + t38 * t24];
taug = t1(:);
