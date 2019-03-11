% Calculate potential energy for
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RPPRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:04
% EndTime: 2019-03-09 02:05:05
% DurationCPUTime: 0.48s
% Computational Cost: add. (202->76), mult. (312->72), div. (0->0), fcn. (340->8), ass. (0->33)
t57 = -m(1) - m(2);
t56 = m(6) + m(7);
t55 = -mrSges(2,1) - mrSges(3,1);
t27 = sin(qJ(4));
t29 = cos(qJ(4));
t54 = t29 * mrSges(5,1) - mrSges(5,2) * t27 + mrSges(4,1);
t53 = mrSges(2,2) - mrSges(3,3);
t52 = mrSges(4,2) - mrSges(5,3);
t51 = mrSges(6,3) + mrSges(7,2);
t50 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t49 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t48 = pkin(4) * t29;
t47 = sin(qJ(1));
t30 = cos(qJ(1));
t41 = sin(pkin(9));
t42 = cos(pkin(9));
t13 = -t30 * t42 - t47 * t41;
t46 = t13 * t27;
t14 = t30 * t41 - t47 * t42;
t45 = t14 * t27;
t26 = sin(qJ(5));
t44 = t26 * t29;
t28 = cos(qJ(5));
t43 = t28 * t29;
t25 = pkin(6) + r_base(3);
t18 = -qJ(3) + t25;
t40 = t30 * pkin(1) + t47 * qJ(2) + r_base(1);
t38 = t30 * pkin(2) + t40;
t36 = t47 * pkin(1) - qJ(2) * t30 + r_base(2);
t35 = t47 * pkin(2) + t36;
t34 = -t13 * pkin(3) + pkin(7) * t14 + t38;
t32 = -t14 * pkin(3) - t13 * pkin(7) + t35;
t1 = (-m(1) * r_base(3) - mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + mrSges(4,3) - t56 * (t29 * pkin(8) + t18) + (-m(2) - m(3)) * t25 + (-m(4) - m(5)) * t18 + (mrSges(5,2) - t51) * t29 + (pkin(4) * t56 + t26 * t49 + t28 * t50 + mrSges(5,1)) * t27) * g(3) + (-m(3) * t36 - m(4) * t35 - m(5) * t32 - mrSges(1,2) + t57 * r_base(2) + t55 * t47 + t51 * t45 - t56 * (-pkin(8) * t45 - t14 * t48 + t32) - t53 * t30 - t50 * (-t13 * t26 - t14 * t43) + t54 * t14 - t52 * t13 - t49 * (t13 * t28 - t14 * t44)) * g(2) + (-m(3) * t40 - m(4) * t38 - m(5) * t34 - mrSges(1,1) + t57 * r_base(1) + t53 * t47 + t51 * t46 - t56 * (-pkin(8) * t46 - t13 * t48 + t34) - t50 * (-t13 * t43 + t14 * t26) + t55 * t30 - t49 * (-t13 * t44 - t14 * t28) + t52 * t14 + t54 * t13) * g(1);
U  = t1;
