% Calculate potential energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:02
% EndTime: 2019-03-10 01:18:02
% DurationCPUTime: 0.58s
% Computational Cost: add. (256->77), mult. (254->65), div. (0->0), fcn. (228->10), ass. (0->32)
t25 = qJ(3) + qJ(4);
t17 = qJ(5) + t25;
t12 = cos(t17);
t29 = cos(qJ(3));
t13 = t29 * pkin(3) + pkin(2);
t14 = sin(t25);
t15 = cos(t25);
t26 = sin(qJ(3));
t7 = pkin(4) * t15 + t13;
t58 = -m(6) * t7 - m(7) * (pkin(5) * t12 + t7) - mrSges(3,1) - m(5) * t13 - t15 * mrSges(5,1) + t14 * mrSges(5,2) - m(4) * pkin(2) - t29 * mrSges(4,1) + t26 * mrSges(4,2);
t32 = -pkin(9) - pkin(8);
t24 = -pkin(10) + t32;
t57 = m(6) * t24 + m(7) * (-qJ(6) + t24) + mrSges(3,2) - mrSges(6,3) - mrSges(7,3) + m(5) * t32 - mrSges(5,3) - m(4) * pkin(8) - mrSges(4,3);
t56 = -m(1) - m(2);
t55 = -m(4) - m(5);
t54 = -mrSges(6,1) - mrSges(7,1);
t53 = mrSges(6,2) + mrSges(7,2);
t52 = -m(3) - m(6) - m(7);
t49 = t52 + t55;
t27 = sin(qJ(2));
t30 = cos(qJ(2));
t48 = t57 * t27 + t58 * t30 - mrSges(2,1);
t11 = sin(t17);
t18 = t26 * pkin(3);
t8 = pkin(4) * t14 + t18;
t47 = m(6) * t8 + m(7) * (pkin(5) * t11 + t8) - mrSges(2,2) + mrSges(3,3) + t14 * mrSges(5,1) + t15 * mrSges(5,2) + t26 * mrSges(4,1) + t29 * mrSges(4,2);
t28 = sin(qJ(1));
t46 = t28 * t30;
t31 = cos(qJ(1));
t45 = t30 * t31;
t44 = t28 * pkin(1) + r_base(2);
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + (-m(2) + t49) * (pkin(6) + r_base(3)) - t57 * t30 + (t53 * t11 + t54 * t12 + t58) * t27) * g(3) + (-mrSges(1,2) + t56 * r_base(2) + t55 * t44 + t54 * (-t11 * t31 + t12 * t46) - t53 * (-t11 * t46 - t12 * t31) + t52 * (-t31 * pkin(7) + t44) + (m(4) * pkin(7) - m(5) * (-pkin(7) - t18) + t47) * t31 + t48 * t28) * g(2) + (-mrSges(1,1) + t56 * r_base(1) + t54 * (t11 * t28 + t12 * t45) - t53 * (-t11 * t45 + t12 * t28) + t49 * (t31 * pkin(1) + t28 * pkin(7) + r_base(1)) + (-m(5) * t18 - t47) * t28 + t48 * t31) * g(1);
U  = t1;
