% Calculate potential energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S4RRRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:18
% EndTime: 2019-12-31 17:29:19
% DurationCPUTime: 0.39s
% Computational Cost: add. (145->71), mult. (282->83), div. (0->0), fcn. (311->10), ass. (0->34)
t50 = -m(1) - m(2);
t49 = -m(4) - m(5);
t48 = -m(5) * pkin(8) + mrSges(4,2) - mrSges(5,3);
t23 = sin(qJ(4));
t27 = cos(qJ(4));
t47 = -t23 * mrSges(5,1) - t27 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t46 = -m(5) * pkin(3) - t27 * mrSges(5,1) + t23 * mrSges(5,2) - mrSges(4,1);
t21 = sin(pkin(4));
t25 = sin(qJ(2));
t45 = t21 * t25;
t29 = cos(qJ(2));
t44 = t21 * t29;
t26 = sin(qJ(1));
t43 = t26 * t21;
t42 = t26 * t25;
t41 = t26 * t29;
t30 = cos(qJ(1));
t40 = t30 * t21;
t39 = t30 * t25;
t38 = t30 * t29;
t37 = pkin(5) + r_base(3);
t22 = cos(pkin(4));
t36 = t22 * pkin(6) + t37;
t35 = t30 * pkin(1) + pkin(6) * t43 + r_base(1);
t34 = t26 * pkin(1) - pkin(6) * t40 + r_base(2);
t32 = pkin(2) * t45 - pkin(7) * t44 + t36;
t28 = cos(qJ(3));
t24 = sin(qJ(3));
t12 = -t22 * t42 + t38;
t11 = t22 * t41 + t39;
t10 = t22 * t39 + t41;
t9 = -t22 * t38 + t42;
t8 = t22 * t24 + t28 * t45;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t37 - mrSges(2,3) - m(3) * t36 - t22 * mrSges(3,3) - (t25 * mrSges(3,1) + t29 * mrSges(3,2)) * t21 - m(4) * t32 - t8 * mrSges(4,1) + mrSges(4,3) * t44 - m(5) * (t8 * pkin(3) + t32) - (-t23 * t44 + t8 * t27) * mrSges(5,1) - (-t8 * t23 - t27 * t44) * mrSges(5,2) + t48 * (-t22 * t28 + t24 * t45)) * g(3) + (-m(3) * t34 - t26 * mrSges(2,1) - t10 * mrSges(3,1) - t30 * mrSges(2,2) + mrSges(3,3) * t40 - mrSges(1,2) + t50 * r_base(2) + t49 * (t10 * pkin(2) + t9 * pkin(7) + t34) + t46 * (t10 * t28 - t24 * t40) + t47 * t9 + t48 * (t10 * t24 + t28 * t40)) * g(2) + (-m(3) * t35 - t30 * mrSges(2,1) - t12 * mrSges(3,1) + t26 * mrSges(2,2) - mrSges(3,3) * t43 - mrSges(1,1) + t50 * r_base(1) + t49 * (t12 * pkin(2) + t11 * pkin(7) + t35) + t48 * (t12 * t24 - t28 * t43) + t46 * (t12 * t28 + t24 * t43) + t47 * t11) * g(1);
U = t1;
