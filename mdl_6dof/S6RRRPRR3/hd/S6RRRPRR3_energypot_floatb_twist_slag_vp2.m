% Calculate potential energy for
% S6RRRPRR3
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
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:10:47
% EndTime: 2019-03-09 18:10:48
% DurationCPUTime: 0.52s
% Computational Cost: add. (252->76), mult. (261->67), div. (0->0), fcn. (244->10), ass. (0->37)
t65 = -mrSges(4,1) - mrSges(5,1);
t64 = mrSges(4,2) - mrSges(5,3);
t28 = sin(qJ(5));
t56 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t63 = t56 * t28;
t26 = qJ(2) + qJ(3);
t21 = sin(t26);
t22 = cos(t26);
t29 = sin(qJ(2));
t33 = cos(qJ(2));
t32 = cos(qJ(5));
t5 = t21 * t28 + t22 * t32;
t51 = t21 * t32;
t27 = sin(qJ(6));
t31 = cos(qJ(6));
t55 = -m(7) * pkin(5) - t31 * mrSges(7,1) + t27 * mrSges(7,2) - mrSges(6,1);
t62 = -m(3) * pkin(1) - mrSges(3,1) * t33 + mrSges(3,2) * t29 + t64 * t21 + t65 * t22 + t5 * t55 - t56 * t51 - mrSges(2,1);
t61 = -m(2) - m(3);
t60 = -m(6) - m(7);
t34 = cos(qJ(1));
t48 = qJ(4) * t21;
t49 = t22 * t34;
t59 = pkin(3) * t49 + t34 * t48;
t58 = -m(1) + t61;
t53 = m(3) * pkin(7) - t27 * mrSges(7,1) - t31 * mrSges(7,2) - mrSges(2,2) + mrSges(5,2) + mrSges(3,3) + mrSges(4,3) - mrSges(6,3);
t30 = sin(qJ(1));
t50 = t22 * t30;
t25 = pkin(6) + r_base(3);
t19 = pkin(2) * t33 + pkin(1);
t47 = t34 * t19 + r_base(1);
t46 = t29 * pkin(2) + t25;
t35 = -pkin(8) - pkin(7);
t45 = t30 * t19 + t34 * t35 + r_base(2);
t42 = -t30 * t35 + t47;
t41 = pkin(3) * t50 + t30 * t48 + t45;
t39 = t21 * pkin(3) - qJ(4) * t22 + t46;
t1 = (-m(1) * r_base(3) - m(4) * t46 - m(5) * t39 - t29 * mrSges(3,1) - t33 * mrSges(3,2) - mrSges(1,3) - mrSges(2,3) + t55 * (-t22 * t28 + t51) + t60 * (t21 * pkin(4) + t39) + t56 * t5 + t61 * t25 - t64 * t22 + t65 * t21) * g(3) + (-m(4) * t45 - m(5) * t41 - mrSges(1,2) + t60 * (pkin(4) * t50 + t34 * pkin(9) + t41) + t50 * t63 + t58 * r_base(2) + t53 * t34 + t62 * t30) * g(2) + (-mrSges(1,1) - m(4) * t42 - m(5) * (t42 + t59) + t60 * (pkin(4) * t49 + t47 + t59) + t49 * t63 + t58 * r_base(1) + (t60 * (-pkin(9) - t35) - t53) * t30 + t62 * t34) * g(1);
U  = t1;
