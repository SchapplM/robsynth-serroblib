% Calculate Gravitation load on the joints for
% S4RRRP5
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4RRRP5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:16:39
% DurationCPUTime: 0.31s
% Computational Cost: add. (144->50), mult. (184->55), div. (0->0), fcn. (143->6), ass. (0->27)
t53 = -mrSges(4,1) - mrSges(5,1);
t16 = sin(qJ(2));
t18 = cos(qJ(2));
t15 = qJ(2) + qJ(3);
t12 = sin(t15);
t13 = cos(t15);
t46 = t53 * t13 + (mrSges(4,2) - mrSges(5,3)) * t12;
t52 = -t18 * mrSges(3,1) + t16 * mrSges(3,2) + t46;
t47 = t13 * pkin(3) + t12 * qJ(4);
t51 = m(5) * t47;
t50 = (-m(5) * qJ(4) - mrSges(5,3)) * t13;
t49 = m(4) + m(5);
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t45 = g(1) * t19 + g(2) * t17;
t44 = m(3) * pkin(1) + mrSges(2,1) - t52;
t43 = -m(3) * pkin(5) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t42 = pkin(2) * t16;
t14 = t18 * pkin(2);
t38 = mrSges(4,2) * t13;
t34 = t50 * t17;
t33 = t50 * t19;
t22 = t38 + (m(5) * pkin(3) - t53) * t12;
t21 = m(5) * (-pkin(3) * t12 - t42) - t12 * mrSges(5,1);
t20 = -pkin(6) - pkin(5);
t11 = t14 + pkin(1);
t1 = [(-t49 * (t19 * t11 - t17 * t20) + (-t44 - t51) * t19 + t43 * t17) * g(2) + ((t49 * t20 + t43) * t19 + (m(4) * t11 - m(5) * (-t11 - t47) + t44) * t17) * g(1), -g(1) * (t21 * t19 - t33) - g(2) * (t21 * t17 - t34) + (-m(4) * t14 - m(5) * (t14 + t47) + t52) * g(3) + (m(4) * t42 + mrSges(3,1) * t16 + mrSges(4,1) * t12 + mrSges(3,2) * t18 + t38) * t45, (t46 - t51) * g(3) + (t22 * t17 + t34) * g(2) + (t22 * t19 + t33) * g(1), (g(3) * t13 - t45 * t12) * m(5)];
taug = t1(:);
