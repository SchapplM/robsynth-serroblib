% Calculate potential energy for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:45
% EndTime: 2019-03-09 08:13:46
% DurationCPUTime: 0.49s
% Computational Cost: add. (173->80), mult. (240->62), div. (0->0), fcn. (206->8), ass. (0->33)
t58 = -m(6) * qJ(5) - mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(8) - qJ(5)) - mrSges(7,3);
t18 = sin(pkin(9));
t19 = cos(pkin(9));
t16 = pkin(9) + qJ(6);
t8 = sin(t16);
t9 = cos(t16);
t57 = t19 * mrSges(6,1) + t9 * mrSges(7,1) - t18 * mrSges(6,2) - t8 * mrSges(7,2) + mrSges(5,1) - mrSges(3,2) + mrSges(4,3);
t56 = -m(1) - m(2);
t55 = -m(6) - m(7);
t22 = sin(qJ(1));
t21 = sin(qJ(2));
t44 = qJ(3) * t21;
t23 = cos(qJ(2));
t46 = t22 * t23;
t54 = pkin(2) * t46 + t22 * t44;
t24 = cos(qJ(1));
t53 = pkin(3) * t46 + t24 * qJ(4);
t52 = m(5) - t55;
t7 = pkin(5) * t19 + pkin(4);
t49 = -mrSges(2,1) + t58 * t23 + (-m(6) * pkin(4) - m(7) * t7 - t57) * t21;
t48 = -t18 * mrSges(6,1) - t8 * mrSges(7,1) - t19 * mrSges(6,2) - t9 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3);
t47 = pkin(5) * t18;
t45 = t23 * t24;
t17 = pkin(6) + r_base(3);
t43 = t22 * pkin(1) + r_base(2);
t42 = t21 * pkin(2) + t17;
t41 = t24 * pkin(1) + t22 * pkin(7) + r_base(1);
t32 = -t24 * pkin(7) + t43;
t31 = pkin(2) * t45 + t24 * t44 + t41;
t28 = -qJ(3) * t23 + t42;
t27 = t32 + t54;
t11 = t21 * pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t28 - m(5) * (t11 + t28) + t55 * (t11 + t42) + (-m(2) - m(3)) * t17 + (-m(6) * (-pkin(4) - qJ(3)) - m(7) * (-qJ(3) - t7) + t57) * t23 + t58 * t21) * g(3) + (-mrSges(1,2) - m(3) * t32 - m(4) * t27 - m(5) * (t27 + t53) + t56 * r_base(2) + t55 * (t43 + t53 + t54) + (m(6) * pkin(7) - m(7) * (-pkin(7) + t47) + t48) * t24 + t49 * t22) * g(2) + (-m(3) * t41 - m(4) * t31 - mrSges(1,1) + t56 * r_base(1) - t52 * (pkin(3) * t45 + t31) + t49 * t24 + (m(7) * t47 + t52 * qJ(4) - t48) * t22) * g(1);
U  = t1;
