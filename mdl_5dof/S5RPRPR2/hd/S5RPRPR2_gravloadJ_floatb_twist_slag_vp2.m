% Calculate Gravitation load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:33:37
% EndTime: 2020-01-03 11:33:38
% DurationCPUTime: 0.22s
% Computational Cost: add. (223->50), mult. (143->54), div. (0->0), fcn. (104->10), ass. (0->28)
t28 = pkin(9) + qJ(5);
t21 = sin(t28);
t23 = cos(t28);
t52 = mrSges(6,1) * t23 - mrSges(6,2) * t21;
t51 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t31 = cos(pkin(9));
t50 = sin(pkin(9)) * mrSges(5,2) - mrSges(4,1) - mrSges(5,1) * t31 - t52;
t49 = m(5) + m(6);
t29 = qJ(1) + pkin(8);
t25 = qJ(3) + t29;
t18 = sin(t25);
t19 = cos(t25);
t20 = pkin(4) * t31 + pkin(3);
t32 = -pkin(7) - qJ(4);
t48 = t18 * t20 + t19 * t32;
t47 = t19 * pkin(3) + t18 * qJ(4);
t22 = sin(t29);
t33 = sin(qJ(1));
t43 = t33 * pkin(1) + pkin(2) * t22;
t24 = cos(t29);
t34 = cos(qJ(1));
t42 = t34 * pkin(1) + pkin(2) * t24;
t41 = -m(3) * pkin(1) - mrSges(2,1);
t40 = -t18 * t32 + t19 * t20;
t36 = t51 * t18 + t50 * t19;
t35 = (m(5) * qJ(4) - t51) * t19 + t50 * t18;
t14 = t18 * pkin(3);
t1 = [(-mrSges(2,2) * t34 - mrSges(3,1) * t22 - mrSges(3,2) * t24 - m(4) * t43 - m(5) * (t14 + t43) - m(6) * (t43 + t48) + t41 * t33 + t35) * g(3) + (t33 * mrSges(2,2) - t24 * mrSges(3,1) + t22 * mrSges(3,2) - m(4) * t42 - m(5) * (t42 + t47) - m(6) * (t40 + t42) + t41 * t34 + t36) * g(2), (-m(3) - m(4) - t49) * g(1), (-m(5) * t14 - m(6) * t48 + t35) * g(3) + (-m(5) * t47 - m(6) * t40 + t36) * g(2), t49 * (g(2) * t19 + g(3) * t18), -g(1) * t52 + (g(2) * t18 - g(3) * t19) * (mrSges(6,1) * t21 + mrSges(6,2) * t23)];
taug = t1(:);
