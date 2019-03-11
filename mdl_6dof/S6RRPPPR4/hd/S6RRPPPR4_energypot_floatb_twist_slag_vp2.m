% Calculate potential energy for
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:10
% EndTime: 2019-03-09 08:17:10
% DurationCPUTime: 0.67s
% Computational Cost: add. (186->81), mult. (307->76), div. (0->0), fcn. (297->8), ass. (0->39)
t60 = -mrSges(3,1) + mrSges(4,2);
t59 = mrSges(3,2) - mrSges(4,3);
t58 = -m(1) - m(2);
t57 = -m(6) - m(7);
t26 = sin(qJ(1));
t25 = sin(qJ(2));
t45 = qJ(3) * t25;
t28 = cos(qJ(2));
t47 = t26 * t28;
t56 = pkin(2) * t47 + t26 * t45;
t55 = t59 * t25 + t60 * t28 - mrSges(2,1);
t54 = mrSges(2,2) - mrSges(3,3) - mrSges(4,1);
t53 = m(7) * pkin(8) - mrSges(6,2) - mrSges(5,3) + mrSges(7,3);
t24 = sin(qJ(6));
t27 = cos(qJ(6));
t52 = t24 * mrSges(7,1) + t27 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t51 = m(7) * pkin(5) + t27 * mrSges(7,1) - t24 * mrSges(7,2) + mrSges(5,1) + mrSges(6,1);
t22 = sin(pkin(9));
t50 = t22 * t26;
t29 = cos(qJ(1));
t49 = t25 * t29;
t23 = cos(pkin(9));
t48 = t26 * t23;
t46 = t28 * t29;
t44 = qJ(4) * t28;
t21 = pkin(6) + r_base(3);
t43 = t26 * pkin(1) + r_base(2);
t42 = t25 * pkin(2) + t21;
t40 = t29 * pkin(1) + t26 * pkin(7) + r_base(1);
t39 = t25 * qJ(4) + t42;
t36 = -pkin(7) * t29 + t43;
t34 = pkin(2) * t46 + t29 * t45 + t40;
t33 = t26 * pkin(3) + t29 * t44 + t34;
t32 = t26 * t44 + (-pkin(3) - pkin(7)) * t29 + t43 + t56;
t6 = -t23 * t29 + t25 * t50;
t5 = t22 * t29 + t25 * t48;
t4 = t22 * t49 + t48;
t3 = -t23 * t49 + t50;
t1 = (-m(1) * r_base(3) - m(4) * t42 - m(5) * t39 - mrSges(1,3) - mrSges(2,3) + t57 * (t28 * t23 * qJ(5) + t39) + (-m(2) - m(3)) * t21 + (t57 * (-pkin(4) * t22 - qJ(3)) - t52 * t23 + t51 * t22 + (m(4) + m(5)) * qJ(3) - t59) * t28 + (t53 + t60) * t25) * g(3) + (-mrSges(1,2) - m(3) * t36 - m(4) * (t36 + t56) - m(5) * t32 + t58 * r_base(2) + t57 * (t6 * pkin(4) - t5 * qJ(5) + t32) - t51 * t6 + t52 * t5 + t53 * t47 - t54 * t29 + t55 * t26) * g(2) + (-m(3) * t40 - m(4) * t34 - m(5) * t33 - mrSges(1,1) + t58 * r_base(1) + t57 * (t4 * pkin(4) + t3 * qJ(5) + t33) - t51 * t4 - t52 * t3 + t53 * t46 + t55 * t29 + t54 * t26) * g(1);
U  = t1;
