% Calculate potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RPPPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:08
% EndTime: 2022-01-20 09:12:08
% DurationCPUTime: 0.37s
% Computational Cost: add. (175->56), mult. (153->44), div. (0->0), fcn. (125->10), ass. (0->22)
t14 = sin(pkin(9));
t16 = cos(pkin(9));
t12 = pkin(9) + qJ(5);
t5 = sin(t12);
t7 = cos(t12);
t41 = -m(5) * pkin(3) - t16 * mrSges(5,1) + t14 * mrSges(5,2) - mrSges(4,1) - m(6) * (t16 * pkin(4) + pkin(3)) - t7 * mrSges(6,1) + t5 * mrSges(6,2);
t40 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3) + m(6) * (-pkin(6) - qJ(4)) - mrSges(6,3);
t39 = -m(1) - m(2);
t38 = m(4) + m(5) + m(6);
t15 = sin(pkin(8));
t17 = cos(pkin(8));
t35 = t40 * t15 + t41 * t17 - mrSges(3,1);
t34 = t5 * mrSges(6,1) + t16 * mrSges(5,2) + t7 * mrSges(6,2) - mrSges(3,2) + mrSges(4,3) + (m(6) * pkin(4) + mrSges(5,1)) * t14;
t32 = pkin(5) + r_base(3);
t19 = sin(qJ(1));
t31 = t19 * pkin(1) + r_base(2);
t20 = cos(qJ(1));
t30 = t20 * pkin(1) + r_base(1);
t13 = qJ(1) + pkin(7);
t8 = cos(t13);
t6 = sin(t13);
t1 = (-m(1) * r_base(3) - m(2) * t32 - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + (-m(3) - t38) * (qJ(2) + t32) - t40 * t17 + t41 * t15) * g(3) + (-m(3) * t31 - t19 * mrSges(2,1) - t20 * mrSges(2,2) - mrSges(1,2) + t39 * r_base(2) - t38 * (t6 * pkin(2) + t31) + (t38 * qJ(3) + t34) * t8 + t35 * t6) * g(2) + (-m(3) * t30 - t20 * mrSges(2,1) + t19 * mrSges(2,2) - mrSges(1,1) + t39 * r_base(1) - t38 * (t8 * pkin(2) + t6 * qJ(3) + t30) + t35 * t8 - t34 * t6) * g(1);
U = t1;
