% Calculate potential energy for
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PPRRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:17
% EndTime: 2019-12-05 15:16:18
% DurationCPUTime: 0.72s
% Computational Cost: add. (168->65), mult. (266->61), div. (0->0), fcn. (266->10), ass. (0->29)
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t58 = -(m(6) * pkin(4) + mrSges(5,1)) * t23 - t25 * mrSges(5,2) + mrSges(3,2);
t19 = sin(pkin(9));
t21 = cos(pkin(9));
t51 = m(4) + m(5) + m(6);
t54 = -mrSges(2,1) + (-t51 * pkin(2) - mrSges(3,1)) * t21 + t58 * t19;
t53 = -m(1) - m(2);
t52 = mrSges(2,2) - mrSges(3,3);
t18 = qJ(4) + qJ(5);
t12 = sin(t18);
t13 = cos(t18);
t47 = t12 * mrSges(6,1) + t13 * mrSges(6,2) + mrSges(4,3);
t46 = -m(5) * pkin(6) + m(6) * (-pkin(7) - pkin(6)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t45 = -m(5) * pkin(3) - m(6) * (pkin(4) * t25 + pkin(3)) - t25 * mrSges(5,1) - t13 * mrSges(6,1) + t23 * mrSges(5,2) + t12 * mrSges(6,2) - mrSges(4,1);
t20 = sin(pkin(8));
t43 = t19 * t20;
t22 = cos(pkin(8));
t42 = t19 * t22;
t24 = sin(qJ(3));
t39 = t20 * t24;
t26 = cos(qJ(3));
t38 = t20 * t26;
t37 = t22 * t24;
t36 = t22 * t26;
t17 = qJ(1) + r_base(3);
t33 = t22 * pkin(1) + t20 * qJ(2) + r_base(1);
t29 = t20 * pkin(1) - t22 * qJ(2) + r_base(2);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) - m(3)) * t17 - t51 * (t19 * pkin(2) + t17) + (t51 * pkin(5) + t47 - t58) * t21 + (t46 * t24 + t45 * t26 - mrSges(3,1)) * t19) * g(3) + (-m(3) * t29 - mrSges(1,2) + t53 * r_base(2) - t51 * (pkin(5) * t43 + t29) - t52 * t22 - t47 * t43 + t45 * (t21 * t38 - t37) + t46 * (t21 * t39 + t36) + t54 * t20) * g(2) + (-m(3) * t33 - mrSges(1,1) + t53 * r_base(1) - t51 * (pkin(5) * t42 + t33) - t47 * t42 + t45 * (t21 * t36 + t39) + t52 * t20 + t46 * (t21 * t37 - t38) + t54 * t22) * g(1);
U = t1;
