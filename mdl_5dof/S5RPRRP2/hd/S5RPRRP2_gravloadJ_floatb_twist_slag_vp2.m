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
% m [6x1]
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:27:52
% EndTime: 2022-01-23 09:27:52
% DurationCPUTime: 0.20s
% Computational Cost: add. (224->47), mult. (160->49), div. (0->0), fcn. (117->8), ass. (0->26)
t26 = cos(qJ(4));
t24 = sin(qJ(4));
t38 = mrSges(5,2) + mrSges(6,2);
t32 = t38 * t24;
t39 = mrSges(5,1) + mrSges(6,1);
t44 = -t39 * t26 - mrSges(4,1) + t32;
t42 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t22 = qJ(1) + pkin(8);
t20 = qJ(3) + t22;
t15 = sin(t20);
t16 = cos(t20);
t37 = t16 * pkin(3) + t15 * pkin(7);
t19 = cos(t22);
t27 = cos(qJ(1));
t36 = t27 * pkin(1) + pkin(2) * t19;
t35 = m(4) + m(5) + m(6);
t34 = m(3) + t35;
t17 = pkin(4) * t26 + pkin(3);
t23 = -qJ(5) - pkin(7);
t33 = -t15 * t23 + t16 * t17;
t31 = m(6) * pkin(4) + t39;
t29 = t42 * t15 + t44 * t16;
t28 = (-m(5) * pkin(7) + m(6) * t23 + t42) * t16 + (m(5) * pkin(3) + m(6) * t17 - t44) * t15;
t25 = sin(qJ(1));
t18 = sin(t22);
t1 = [(t25 * mrSges(2,2) - t19 * mrSges(3,1) + t18 * mrSges(3,2) - m(4) * t36 - m(5) * (t36 + t37) - m(6) * (t33 + t36) + (-m(3) * pkin(1) - mrSges(2,1)) * t27 + t29) * g(2) + (mrSges(2,2) * t27 + mrSges(3,2) * t19 + (t35 * pkin(2) + mrSges(3,1)) * t18 + (t34 * pkin(1) + mrSges(2,1)) * t25 + t28) * g(1), -t34 * g(3), (-m(5) * t37 - m(6) * t33 + t29) * g(2) + t28 * g(1), (-t31 * t26 + t32) * g(3) + (g(1) * t16 + g(2) * t15) * (t31 * t24 + t38 * t26), (-g(1) * t15 + g(2) * t16) * m(6)];
taug = t1(:);
