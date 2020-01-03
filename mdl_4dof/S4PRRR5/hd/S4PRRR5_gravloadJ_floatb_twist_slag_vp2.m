% Calculate Gravitation load on the joints for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S4PRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:36
% DurationCPUTime: 0.22s
% Computational Cost: add. (112->44), mult. (153->60), div. (0->0), fcn. (126->8), ass. (0->27)
t19 = cos(qJ(4));
t39 = mrSges(5,1) * t19;
t49 = -mrSges(4,1) - t39;
t14 = qJ(2) + qJ(3);
t12 = sin(t14);
t13 = cos(t14);
t17 = sin(qJ(4));
t37 = mrSges(5,2) * t17;
t48 = -t12 * t37 + t13 * (-m(5) * pkin(6) - mrSges(5,3));
t45 = (t37 + t49) * t13 + (mrSges(4,2) - mrSges(5,3)) * t12;
t18 = sin(qJ(2));
t43 = pkin(2) * t18;
t20 = cos(qJ(2));
t42 = pkin(2) * t20;
t38 = mrSges(4,2) * t13;
t15 = sin(pkin(7));
t34 = t15 * t17;
t33 = t15 * t19;
t16 = cos(pkin(7));
t32 = t16 * t17;
t31 = t16 * t19;
t30 = t13 * pkin(3) + t12 * pkin(6);
t27 = t48 * t15;
t26 = t48 * t16;
t22 = m(5) * (-pkin(3) * t12 - t43) - t12 * t39;
t21 = t38 + (m(5) * pkin(3) - t49) * t12;
t1 = [(-m(2) - m(3) - m(4) - m(5)) * g(3), -g(1) * (t22 * t16 - t26) - g(2) * (t22 * t15 - t27) + (-mrSges(3,1) * t20 + t18 * mrSges(3,2) - m(4) * t42 - m(5) * (t30 + t42) + t45) * g(3) + (m(4) * t43 + mrSges(3,1) * t18 + mrSges(4,1) * t12 + mrSges(3,2) * t20 + t38) * (g(1) * t16 + g(2) * t15), (-m(5) * t30 + t45) * g(3) + (t21 * t15 + t27) * g(2) + (t21 * t16 + t26) * g(1), -g(1) * ((-t13 * t32 + t33) * mrSges(5,1) + (-t13 * t31 - t34) * mrSges(5,2)) - g(2) * ((-t13 * t34 - t31) * mrSges(5,1) + (-t13 * t33 + t32) * mrSges(5,2)) - g(3) * (-mrSges(5,1) * t17 - mrSges(5,2) * t19) * t12];
taug = t1(:);
