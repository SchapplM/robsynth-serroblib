% Calculate Gravitation load on the joints for
% S4RRPR4
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRPR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:27
% DurationCPUTime: 0.18s
% Computational Cost: add. (142->37), mult. (123->41), div. (0->0), fcn. (92->8), ass. (0->21)
t20 = pkin(7) + qJ(4);
t15 = sin(t20);
t16 = cos(t20);
t40 = t16 * mrSges(5,1) - t15 * mrSges(5,2);
t39 = mrSges(3,2) - mrSges(5,3) - mrSges(4,3);
t23 = cos(pkin(7));
t38 = -t23 * mrSges(4,1) - mrSges(3,1) + sin(pkin(7)) * mrSges(4,2) - t40;
t37 = m(4) + m(5);
t21 = qJ(1) + qJ(2);
t17 = sin(t21);
t18 = cos(t21);
t36 = t18 * pkin(2) + t17 * qJ(3);
t13 = pkin(3) * t23 + pkin(2);
t24 = -pkin(6) - qJ(3);
t31 = t18 * t13 - t17 * t24;
t28 = t39 * t17 + t38 * t18;
t27 = (m(4) * pkin(2) + m(5) * t13 - t38) * t17 + (-m(4) * qJ(3) + m(5) * t24 + t39) * t18;
t26 = cos(qJ(1));
t25 = sin(qJ(1));
t19 = t26 * pkin(1);
t1 = [(t25 * mrSges(2,2) - m(4) * (t19 + t36) - m(5) * (t19 + t31) + (-m(3) * pkin(1) - mrSges(2,1)) * t26 + t28) * g(2) + (mrSges(2,2) * t26 + (mrSges(2,1) + (m(3) + t37) * pkin(1)) * t25 + t27) * g(1), (-m(4) * t36 - m(5) * t31 + t28) * g(2) + t27 * g(1), t37 * (-g(1) * t17 + g(2) * t18), -g(3) * t40 + (g(1) * t18 + g(2) * t17) * (mrSges(5,1) * t15 + mrSges(5,2) * t16)];
taug = t1(:);
