% Calculate Gravitation load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:07
% EndTime: 2020-01-03 11:45:08
% DurationCPUTime: 0.27s
% Computational Cost: add. (224->51), mult. (160->54), div. (0->0), fcn. (117->8), ass. (0->28)
t30 = cos(qJ(4));
t28 = sin(qJ(4));
t41 = mrSges(5,2) + mrSges(6,2);
t35 = t41 * t28;
t49 = -mrSges(6,1) - mrSges(5,1);
t50 = t30 * t49 - mrSges(4,1) + t35;
t47 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t26 = qJ(1) + pkin(8);
t23 = qJ(3) + t26;
t18 = sin(t23);
t19 = cos(t23);
t20 = pkin(4) * t30 + pkin(3);
t27 = -qJ(5) - pkin(7);
t44 = t18 * t20 + t19 * t27;
t40 = t19 * pkin(3) + t18 * pkin(7);
t21 = sin(t26);
t29 = sin(qJ(1));
t39 = t29 * pkin(1) + pkin(2) * t21;
t22 = cos(t26);
t31 = cos(qJ(1));
t38 = t31 * pkin(1) + pkin(2) * t22;
t37 = -m(3) * pkin(1) - mrSges(2,1);
t36 = -t18 * t27 + t19 * t20;
t34 = m(6) * pkin(4) - t49;
t33 = t47 * t18 + t19 * t50;
t32 = (m(5) * pkin(7) - t47) * t19 + t50 * t18;
t14 = t18 * pkin(3);
t1 = [(-mrSges(2,2) * t31 - mrSges(3,1) * t21 - mrSges(3,2) * t22 - m(4) * t39 - m(5) * (t14 + t39) - m(6) * (t39 + t44) + t37 * t29 + t32) * g(3) + (t29 * mrSges(2,2) - t22 * mrSges(3,1) + t21 * mrSges(3,2) - m(4) * t38 - m(5) * (t38 + t40) - m(6) * (t36 + t38) + t37 * t31 + t33) * g(2), (-m(3) - m(4) - m(5) - m(6)) * g(1), (-m(5) * t14 - m(6) * t44 + t32) * g(3) + (-m(5) * t40 - m(6) * t36 + t33) * g(2), (-t30 * t34 + t35) * g(1) + (g(2) * t18 - g(3) * t19) * (t28 * t34 + t30 * t41), (g(2) * t19 + g(3) * t18) * m(6)];
taug = t1(:);
