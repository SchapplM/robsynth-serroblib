% Calculate potential energy for
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:26
% EndTime: 2019-03-09 15:56:27
% DurationCPUTime: 0.71s
% Computational Cost: add. (227->88), mult. (376->82), div. (0->0), fcn. (388->10), ass. (0->40)
t60 = -m(1) - m(2);
t59 = -m(5) - m(6);
t26 = sin(qJ(3));
t27 = sin(qJ(2));
t29 = cos(qJ(3));
t58 = (pkin(3) * t29 + qJ(4) * t26) * t27;
t57 = mrSges(2,2) - mrSges(3,3);
t56 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3) - mrSges(7,3);
t30 = cos(qJ(2));
t55 = -t30 * mrSges(3,1) - mrSges(2,1) + (m(6) * qJ(5) + mrSges(3,2)) * t27;
t25 = -pkin(9) - qJ(5);
t54 = -m(7) * t25 - t56;
t21 = pkin(10) + qJ(6);
t15 = sin(t21);
t16 = cos(t21);
t24 = cos(pkin(10));
t53 = -t15 * mrSges(7,1) - t24 * mrSges(6,2) - t16 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t23 = sin(pkin(10));
t52 = -m(7) * (pkin(5) * t23 + qJ(4)) - t23 * mrSges(6,1) + t53;
t51 = -m(6) * pkin(4) - m(7) * (pkin(5) * t24 + pkin(4)) - t24 * mrSges(6,1) - t16 * mrSges(7,1) + t23 * mrSges(6,2) + t15 * mrSges(7,2) - mrSges(4,1) - mrSges(5,1);
t28 = sin(qJ(1));
t49 = t27 * t28;
t31 = cos(qJ(1));
t48 = t27 * t31;
t47 = t28 * t30;
t46 = t30 * t31;
t22 = pkin(6) + r_base(3);
t44 = t27 * pkin(2) + t22;
t42 = t31 * pkin(1) + t28 * pkin(7) + r_base(1);
t39 = t28 * pkin(1) - pkin(7) * t31 + r_base(2);
t38 = pkin(2) * t46 + pkin(8) * t48 + t42;
t37 = -pkin(8) * t30 + t44;
t6 = t28 * t26 + t29 * t46;
t36 = t6 * pkin(3) + t38;
t35 = pkin(2) * t47 + pkin(8) * t49 + t39;
t4 = -t26 * t31 + t29 * t47;
t34 = t4 * pkin(3) + t35;
t5 = t26 * t46 - t28 * t29;
t3 = t26 * t47 + t29 * t31;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t37 - m(5) * (t37 + t58) + (-m(6) - m(7)) * (t44 + t58) + (-m(2) - m(3)) * t22 + (-mrSges(3,2) - m(6) * (-pkin(8) + qJ(5)) - m(7) * (-pkin(8) - t25) + t56) * t30 + (t51 * t29 - mrSges(3,1) + ((-m(7) * pkin(5) - mrSges(6,1)) * t23 + t53) * t26) * t27) * g(3) + (-m(3) * t39 - m(4) * t35 - m(7) * t34 - mrSges(1,2) + t60 * r_base(2) + t59 * (t3 * qJ(4) + t34) - t57 * t31 + t51 * t4 + t52 * t3 + t55 * t28 + t54 * t49) * g(2) + (-m(3) * t42 - m(4) * t38 - m(7) * t36 - mrSges(1,1) + t60 * r_base(1) + t59 * (t5 * qJ(4) + t36) + t51 * t6 + t52 * t5 + t55 * t31 + t57 * t28 + t54 * t48) * g(1);
U  = t1;
