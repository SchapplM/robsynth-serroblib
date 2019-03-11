% Calculate potential energy for
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:22
% EndTime: 2019-03-09 20:58:23
% DurationCPUTime: 0.58s
% Computational Cost: add. (268->79), mult. (269->71), div. (0->0), fcn. (247->10), ass. (0->36)
t29 = cos(qJ(3));
t15 = t29 * pkin(3) + pkin(2);
t25 = qJ(3) + qJ(4);
t17 = sin(t25);
t18 = cos(t25);
t26 = sin(qJ(3));
t61 = -m(4) * pkin(2) - m(5) * t15 - t29 * mrSges(4,1) - t18 * mrSges(5,1) + t26 * mrSges(4,2) + t17 * mrSges(5,2) - mrSges(3,1);
t32 = -pkin(9) - pkin(8);
t60 = -m(4) * pkin(8) + m(5) * t32 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t59 = -m(1) - m(2);
t58 = -m(6) - m(7);
t57 = -mrSges(6,3) - mrSges(7,2);
t56 = m(3) + m(4) + m(5);
t27 = sin(qJ(2));
t30 = cos(qJ(2));
t53 = t60 * t27 + t30 * t61 - mrSges(2,1);
t52 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t51 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t49 = pkin(3) * t26;
t50 = m(5) * t49 + t26 * mrSges(4,1) + t17 * mrSges(5,1) + t29 * mrSges(4,2) + t18 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3);
t28 = sin(qJ(1));
t48 = t27 * t28;
t31 = cos(qJ(1));
t47 = t27 * t31;
t46 = t28 * t30;
t45 = t30 * t31;
t24 = pkin(6) + r_base(3);
t44 = t28 * pkin(1) + r_base(2);
t43 = t31 * pkin(1) + t28 * pkin(7) + r_base(1);
t23 = -qJ(5) + t32;
t16 = pkin(10) + t25;
t14 = cos(t16);
t13 = sin(t16);
t10 = pkin(4) * t17 + t49;
t9 = pkin(4) * t18 + t15;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t58 * (t30 * t23 + t27 * t9 + t24) + (-m(2) - t56) * t24 + (-t57 - t60) * t30 + (t51 * t13 + t52 * t14 + t61) * t27) * g(3) + (-mrSges(1,2) + t59 * r_base(2) + t57 * t48 + t58 * (-t23 * t48 + t9 * t46 + (-pkin(7) - t10) * t31 + t44) + t52 * (-t13 * t31 + t14 * t46) + t51 * (t13 * t46 + t14 * t31) - t56 * t44 + (t56 * pkin(7) + t50) * t31 + t53 * t28) * g(2) + (-mrSges(1,1) + t59 * r_base(1) + t57 * t47 + t58 * (t28 * t10 - t23 * t47 + t9 * t45 + t43) + t52 * (t13 * t28 + t14 * t45) + t51 * (t13 * t45 - t28 * t14) - t56 * t43 + t53 * t31 - t50 * t28) * g(1);
U  = t1;
