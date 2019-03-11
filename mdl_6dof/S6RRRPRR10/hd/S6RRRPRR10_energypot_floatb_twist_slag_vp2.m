% Calculate potential energy for
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:31
% EndTime: 2019-03-09 19:15:31
% DurationCPUTime: 0.70s
% Computational Cost: add. (227->88), mult. (376->82), div. (0->0), fcn. (388->10), ass. (0->40)
t60 = -m(1) - m(2);
t59 = -m(5) - m(6);
t24 = sin(qJ(3));
t25 = sin(qJ(2));
t28 = cos(qJ(3));
t58 = (pkin(3) * t28 + qJ(4) * t24) * t25;
t57 = mrSges(2,2) - mrSges(3,3);
t56 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3) - mrSges(7,3);
t29 = cos(qJ(2));
t31 = -pkin(10) - pkin(9);
t55 = -t29 * mrSges(3,1) - mrSges(2,1) + (-m(7) * t31 + mrSges(3,2)) * t25;
t54 = m(6) * pkin(9) - t56;
t22 = qJ(5) + qJ(6);
t15 = sin(t22);
t16 = cos(t22);
t27 = cos(qJ(5));
t53 = -t15 * mrSges(7,1) - mrSges(6,2) * t27 - t16 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t23 = sin(qJ(5));
t52 = -m(7) * (pkin(5) * t23 + qJ(4)) - mrSges(6,1) * t23 + t53;
t51 = -m(6) * pkin(4) - m(7) * (pkin(5) * t27 + pkin(4)) - t27 * mrSges(6,1) - t16 * mrSges(7,1) + t23 * mrSges(6,2) + t15 * mrSges(7,2) - mrSges(4,1) - mrSges(5,1);
t26 = sin(qJ(1));
t49 = t25 * t26;
t30 = cos(qJ(1));
t48 = t25 * t30;
t46 = t26 * t29;
t45 = t29 * t30;
t21 = pkin(6) + r_base(3);
t44 = t25 * pkin(2) + t21;
t42 = t30 * pkin(1) + t26 * pkin(7) + r_base(1);
t39 = t26 * pkin(1) - t30 * pkin(7) + r_base(2);
t38 = pkin(2) * t45 + pkin(8) * t48 + t42;
t37 = -t29 * pkin(8) + t44;
t6 = t24 * t26 + t28 * t45;
t36 = t6 * pkin(3) + t38;
t35 = pkin(2) * t46 + pkin(8) * t49 + t39;
t4 = -t24 * t30 + t28 * t46;
t34 = t4 * pkin(3) + t35;
t5 = t24 * t45 - t26 * t28;
t3 = t24 * t46 + t28 * t30;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t37 - m(5) * (t37 + t58) + (-m(6) - m(7)) * (t44 + t58) + (-m(2) - m(3)) * t21 + (-mrSges(3,2) - m(6) * (-pkin(8) + pkin(9)) - m(7) * (-pkin(8) - t31) + t56) * t29 + (t51 * t28 - mrSges(3,1) + ((-m(7) * pkin(5) - mrSges(6,1)) * t23 + t53) * t24) * t25) * g(3) + (-m(3) * t39 - m(4) * t35 - m(7) * t34 - mrSges(1,2) + t60 * r_base(2) + t59 * (t3 * qJ(4) + t34) - t57 * t30 + t51 * t4 + t52 * t3 + t55 * t26 + t54 * t49) * g(2) + (-m(3) * t42 - m(4) * t38 - m(7) * t36 - mrSges(1,1) + t60 * r_base(1) + t59 * (t5 * qJ(4) + t36) + t51 * t6 + t52 * t5 + t55 * t30 + t57 * t26 + t54 * t48) * g(1);
U  = t1;
