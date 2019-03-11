% Calculate potential energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:28
% EndTime: 2019-03-09 08:33:29
% DurationCPUTime: 0.50s
% Computational Cost: add. (163->75), mult. (240->66), div. (0->0), fcn. (206->6), ass. (0->33)
t56 = mrSges(5,1) + mrSges(4,3) - mrSges(3,2);
t55 = -m(6) * pkin(8) + m(7) * (-qJ(6) - pkin(8)) - mrSges(3,1) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t54 = m(7) * pkin(5);
t53 = -m(1) - m(2);
t52 = -m(6) - m(7);
t51 = mrSges(6,1) + mrSges(7,1);
t50 = -mrSges(6,2) - mrSges(7,2);
t49 = -m(5) + t52;
t48 = mrSges(2,2) - mrSges(3,3) - mrSges(4,2) + mrSges(5,3);
t23 = cos(qJ(5));
t11 = pkin(5) * t23 + pkin(4);
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t47 = -mrSges(2,1) + (-m(6) * pkin(4) - m(7) * t11 - t56) * t21 + t55 * t24;
t20 = sin(qJ(5));
t25 = cos(qJ(1));
t46 = t20 * t25;
t22 = sin(qJ(1));
t45 = t22 * t20;
t44 = t22 * t23;
t43 = t22 * t24;
t42 = t23 * t25;
t41 = t24 * t25;
t40 = qJ(3) * t21;
t18 = pkin(6) + r_base(3);
t39 = t21 * pkin(2) + t18;
t38 = t25 * pkin(1) + t22 * pkin(7) + r_base(1);
t33 = t22 * pkin(1) - pkin(7) * t25 + r_base(2);
t32 = pkin(2) * t41 + t25 * t40 + t38;
t31 = -qJ(3) * t24 + t39;
t30 = pkin(2) * t43 + t22 * t40 + t33;
t13 = t21 * pkin(3);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t31 - m(5) * (t13 + t31) + t52 * (t13 + t39) + (-m(2) - m(3)) * t18 + (-m(6) * (-pkin(4) - qJ(3)) - m(7) * (-qJ(3) - t11) + t51 * t23 + t50 * t20 + t56) * t24 + t55 * t21) * g(3) + (-t46 * t54 - m(3) * t33 - m(4) * t30 - mrSges(1,2) + t53 * r_base(2) + t49 * (pkin(3) * t43 + t25 * qJ(4) + t30) - t51 * (t21 * t44 + t46) + t50 * (-t21 * t45 + t42) - t48 * t25 + t47 * t22) * g(2) + (t45 * t54 - m(3) * t38 - m(4) * t32 - mrSges(1,1) + t53 * r_base(1) - t51 * (t21 * t42 - t45) + t50 * (-t21 * t46 - t44) + t49 * (pkin(3) * t41 - t22 * qJ(4) + t32) + t48 * t22 + t47 * t25) * g(1);
U  = t1;
