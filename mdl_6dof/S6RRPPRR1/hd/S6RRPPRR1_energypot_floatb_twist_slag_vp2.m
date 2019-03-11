% Calculate potential energy for
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:52
% EndTime: 2019-03-09 08:45:52
% DurationCPUTime: 0.49s
% Computational Cost: add. (252->76), mult. (261->67), div. (0->0), fcn. (244->10), ass. (0->37)
t65 = -mrSges(4,1) - mrSges(5,1);
t64 = mrSges(4,2) - mrSges(5,3);
t29 = sin(qJ(5));
t56 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t63 = t56 * t29;
t25 = qJ(2) + pkin(10);
t21 = sin(t25);
t22 = cos(t25);
t30 = sin(qJ(2));
t34 = cos(qJ(2));
t33 = cos(qJ(5));
t5 = t21 * t29 + t22 * t33;
t51 = t21 * t33;
t28 = sin(qJ(6));
t32 = cos(qJ(6));
t55 = -m(7) * pkin(5) - t32 * mrSges(7,1) + t28 * mrSges(7,2) - mrSges(6,1);
t62 = -m(3) * pkin(1) - mrSges(3,1) * t34 + mrSges(3,2) * t30 + t64 * t21 + t65 * t22 + t5 * t55 - t56 * t51 - mrSges(2,1);
t61 = -m(2) - m(3);
t60 = -m(6) - m(7);
t35 = cos(qJ(1));
t48 = qJ(4) * t21;
t49 = t22 * t35;
t59 = pkin(3) * t49 + t35 * t48;
t58 = -m(1) + t61;
t53 = m(3) * pkin(7) - t28 * mrSges(7,1) - t32 * mrSges(7,2) - mrSges(2,2) + mrSges(5,2) + mrSges(3,3) + mrSges(4,3) - mrSges(6,3);
t31 = sin(qJ(1));
t50 = t22 * t31;
t26 = pkin(6) + r_base(3);
t20 = pkin(2) * t34 + pkin(1);
t47 = t35 * t20 + r_base(1);
t46 = t30 * pkin(2) + t26;
t27 = -qJ(3) - pkin(7);
t45 = t31 * t20 + t35 * t27 + r_base(2);
t42 = -t27 * t31 + t47;
t41 = pkin(3) * t50 + t31 * t48 + t45;
t39 = t21 * pkin(3) - qJ(4) * t22 + t46;
t1 = (-m(1) * r_base(3) - m(4) * t46 - m(5) * t39 - t30 * mrSges(3,1) - t34 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t55 * (-t22 * t29 + t51) + t60 * (t21 * pkin(4) + t39) + t56 * t5 + t61 * t26 - t64 * t22 + t65 * t21) * g(3) + (-m(4) * t45 - m(5) * t41 - mrSges(1,2) + t60 * (pkin(4) * t50 + t35 * pkin(8) + t41) + t50 * t63 + t58 * r_base(2) + t53 * t35 + t62 * t31) * g(2) + (-mrSges(1,1) - m(4) * t42 - m(5) * (t42 + t59) + t60 * (pkin(4) * t49 + t47 + t59) + t49 * t63 + t58 * r_base(1) + (t60 * (-pkin(8) - t27) - t53) * t31 + t62 * t35) * g(1);
U  = t1;
