% Calculate potential energy for
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:38
% EndTime: 2019-03-09 16:47:38
% DurationCPUTime: 0.59s
% Computational Cost: add. (268->79), mult. (269->71), div. (0->0), fcn. (247->10), ass. (0->36)
t30 = cos(qJ(3));
t15 = t30 * pkin(3) + pkin(2);
t24 = qJ(3) + pkin(10);
t16 = sin(t24);
t17 = cos(t24);
t27 = sin(qJ(3));
t61 = -m(4) * pkin(2) - m(5) * t15 - t30 * mrSges(4,1) - t17 * mrSges(5,1) + t27 * mrSges(4,2) + t16 * mrSges(5,2) - mrSges(3,1);
t26 = -qJ(4) - pkin(8);
t60 = -m(4) * pkin(8) + m(5) * t26 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t59 = -m(1) - m(2);
t58 = -m(6) - m(7);
t57 = -mrSges(6,3) - mrSges(7,2);
t56 = m(3) + m(4) + m(5);
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t53 = t60 * t28 + t31 * t61 - mrSges(2,1);
t52 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t51 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t49 = pkin(3) * t27;
t50 = m(5) * t49 + t27 * mrSges(4,1) + t16 * mrSges(5,1) + t30 * mrSges(4,2) + t17 * mrSges(5,2) - mrSges(2,2) + mrSges(3,3);
t29 = sin(qJ(1));
t48 = t28 * t29;
t32 = cos(qJ(1));
t47 = t28 * t32;
t46 = t29 * t31;
t45 = t31 * t32;
t25 = pkin(6) + r_base(3);
t44 = t29 * pkin(1) + r_base(2);
t43 = t32 * pkin(1) + t29 * pkin(7) + r_base(1);
t23 = -pkin(9) + t26;
t18 = qJ(5) + t24;
t14 = cos(t18);
t13 = sin(t18);
t10 = pkin(4) * t16 + t49;
t9 = pkin(4) * t17 + t15;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) + t58 * (t31 * t23 + t28 * t9 + t25) + (-m(2) - t56) * t25 + (-t57 - t60) * t31 + (t51 * t13 + t52 * t14 + t61) * t28) * g(3) + (-mrSges(1,2) + t59 * r_base(2) + t57 * t48 + t58 * (-t23 * t48 + t9 * t46 + (-pkin(7) - t10) * t32 + t44) + t52 * (-t13 * t32 + t14 * t46) + t51 * (t13 * t46 + t14 * t32) - t56 * t44 + (t56 * pkin(7) + t50) * t32 + t53 * t29) * g(2) + (-mrSges(1,1) + t59 * r_base(1) + t57 * t47 + t58 * (t29 * t10 - t23 * t47 + t9 * t45 + t43) + t52 * (t29 * t13 + t14 * t45) + t51 * (t13 * t45 - t29 * t14) - t56 * t43 + t53 * t32 - t50 * t29) * g(1);
U  = t1;
