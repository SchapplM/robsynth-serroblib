% Calculate Gravitation load on the joints for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:34
% EndTime: 2019-12-05 16:15:35
% DurationCPUTime: 0.15s
% Computational Cost: add. (212->40), mult. (128->42), div. (0->0), fcn. (92->8), ass. (0->23)
t23 = pkin(9) + qJ(5);
t18 = sin(t23);
t20 = cos(t23);
t42 = t20 * mrSges(6,1) - t18 * mrSges(6,2);
t41 = mrSges(4,2) - mrSges(6,3) - mrSges(5,3);
t26 = cos(pkin(9));
t40 = -t26 * mrSges(5,1) - mrSges(4,1) + sin(pkin(9)) * mrSges(5,2) - t42;
t39 = m(5) + m(6);
t24 = pkin(8) + qJ(2);
t22 = qJ(3) + t24;
t15 = sin(t22);
t16 = cos(t22);
t38 = t16 * pkin(3) + t15 * qJ(4);
t33 = m(4) + t39;
t17 = pkin(4) * t26 + pkin(3);
t27 = -pkin(7) - qJ(4);
t32 = -t15 * t27 + t16 * t17;
t29 = t41 * t15 + t40 * t16;
t28 = (m(5) * pkin(3) + m(6) * t17 - t40) * t15 + (-m(5) * qJ(4) + m(6) * t27 + t41) * t16;
t21 = cos(t24);
t19 = sin(t24);
t14 = pkin(2) * t21;
t1 = [(-m(2) - m(3) - t33) * g(3), (mrSges(3,2) * t19 - m(5) * (t14 + t38) - m(6) * (t14 + t32) + (-m(4) * pkin(2) - mrSges(3,1)) * t21 + t29) * g(2) + (mrSges(3,2) * t21 + (t33 * pkin(2) + mrSges(3,1)) * t19 + t28) * g(1), (-m(5) * t38 - m(6) * t32 + t29) * g(2) + t28 * g(1), t39 * (-g(1) * t15 + g(2) * t16), -g(3) * t42 + (g(1) * t16 + g(2) * t15) * (mrSges(6,1) * t18 + mrSges(6,2) * t20)];
taug = t1(:);
