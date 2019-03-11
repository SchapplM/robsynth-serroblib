% Calculate potential energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:38
% EndTime: 2019-03-09 02:02:38
% DurationCPUTime: 0.47s
% Computational Cost: add. (212->72), mult. (188->60), div. (0->0), fcn. (158->8), ass. (0->29)
t52 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t51 = -m(1) - m(2);
t50 = -m(6) - m(7);
t49 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t48 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t47 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t22 = sin(qJ(4));
t25 = cos(qJ(4));
t46 = -t22 * mrSges(5,1) - t52 * t25 + mrSges(3,2) - mrSges(4,3);
t45 = pkin(4) * t22;
t44 = pkin(8) * t25;
t21 = sin(qJ(5));
t43 = t21 * t22;
t24 = cos(qJ(5));
t42 = t22 * t24;
t39 = pkin(6) + r_base(3);
t23 = sin(qJ(1));
t38 = t23 * pkin(1) + r_base(2);
t26 = cos(qJ(1));
t37 = t26 * pkin(1) + r_base(1);
t15 = qJ(2) + t39;
t20 = qJ(1) + pkin(9);
t13 = sin(t20);
t35 = t13 * pkin(2) + t38;
t34 = pkin(3) + t15;
t33 = t13 * pkin(7) + t35;
t14 = cos(t20);
t32 = t14 * pkin(2) + t13 * qJ(3) + t37;
t1 = (-m(1) * r_base(3) - m(2) * t39 - m(5) * t34 - mrSges(4,1) - mrSges(1,3) - mrSges(2,3) - mrSges(3,3) + t50 * (t25 * pkin(4) + t22 * pkin(8) + t34) + (-m(3) - m(4)) * t15 + (-t47 * t21 + t48 * t24 - mrSges(5,1)) * t25 + t52 * t22) * g(3) + (-m(3) * t38 - m(4) * t35 - m(5) * t33 - t23 * mrSges(2,1) - t26 * mrSges(2,2) - mrSges(1,2) + t51 * r_base(2) + t50 * (t14 * t44 + t33) + t48 * (t13 * t21 - t14 * t42) + t47 * (t13 * t24 + t14 * t43) + t49 * t13 + (t50 * (-qJ(3) - t45) + (m(4) + m(5)) * qJ(3) - t46) * t14) * g(2) + (-m(3) * t37 - m(4) * t32 - t26 * mrSges(2,1) + t23 * mrSges(2,2) + t51 * r_base(1) - mrSges(1,1) + (-m(5) + t50) * (t14 * pkin(7) + t32) + (t48 * t21 + t47 * t24 + t49) * t14 + (t50 * (-t44 + t45) + t48 * t42 - t47 * t43 + t46) * t13) * g(1);
U  = t1;
