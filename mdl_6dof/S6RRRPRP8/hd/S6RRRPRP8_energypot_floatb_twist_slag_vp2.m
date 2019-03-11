% Calculate potential energy for
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:09
% EndTime: 2019-03-09 17:14:10
% DurationCPUTime: 0.70s
% Computational Cost: add. (215->94), mult. (376->95), div. (0->0), fcn. (388->8), ass. (0->45)
t67 = -m(1) - m(2);
t66 = -m(5) - m(6);
t29 = sin(qJ(2));
t33 = cos(qJ(2));
t65 = -t33 * mrSges(3,1) + t29 * mrSges(3,2) - mrSges(2,1);
t64 = -mrSges(4,1) - mrSges(5,1);
t63 = -mrSges(6,1) - mrSges(7,1);
t62 = mrSges(2,2) - mrSges(3,3);
t61 = mrSges(4,2) - mrSges(5,3);
t60 = -mrSges(6,2) - mrSges(7,2);
t32 = cos(qJ(3));
t51 = t29 * t32;
t28 = sin(qJ(3));
t53 = t28 * t29;
t59 = pkin(3) * t51 + qJ(4) * t53;
t58 = mrSges(5,2) + mrSges(4,3) - mrSges(6,3) - mrSges(7,3);
t27 = sin(qJ(5));
t54 = pkin(5) * t27;
t57 = -m(7) * (qJ(4) + t54) + t61;
t31 = cos(qJ(5));
t20 = pkin(5) * t31 + pkin(4);
t56 = -m(6) * pkin(4) - m(7) * t20 + t64;
t26 = -qJ(6) - pkin(9);
t55 = m(6) * pkin(9) - m(7) * t26 - t58;
t30 = sin(qJ(1));
t52 = t29 * t30;
t34 = cos(qJ(1));
t50 = t29 * t34;
t49 = t30 * t33;
t48 = t33 * t34;
t25 = pkin(6) + r_base(3);
t47 = t29 * pkin(2) + t25;
t45 = t34 * pkin(1) + t30 * pkin(7) + r_base(1);
t43 = t30 * pkin(1) - pkin(7) * t34 + r_base(2);
t42 = t47 + t59;
t41 = pkin(2) * t48 + pkin(8) * t50 + t45;
t40 = -pkin(8) * t33 + t47;
t12 = t30 * t28 + t32 * t48;
t39 = t12 * pkin(3) + t41;
t38 = pkin(2) * t49 + pkin(8) * t52 + t43;
t10 = -t28 * t34 + t32 * t49;
t37 = t10 * pkin(3) + t38;
t11 = t28 * t48 - t30 * t32;
t9 = t28 * t49 + t32 * t34;
t1 = (-m(3) * t43 - m(4) * t38 - m(7) * t37 - mrSges(1,2) + t67 * r_base(2) + t57 * t9 + t66 * (t9 * qJ(4) + t37) - t62 * t34 + t65 * t30 + t63 * (t10 * t31 + t27 * t9) + t56 * t10 + t60 * (-t10 * t27 + t31 * t9) + t55 * t52) * g(2) + (-m(3) * t45 - m(4) * t41 - m(7) * t39 - mrSges(1,1) + t67 * r_base(1) + t63 * (t11 * t27 + t12 * t31) + t66 * (t11 * qJ(4) + t39) + t65 * t34 + t62 * t30 + t60 * (t11 * t31 - t12 * t27) + t56 * t12 + t57 * t11 + t55 * t50) * g(1) + (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - m(4) * t40 - m(5) * (t40 + t59) - m(6) * (pkin(4) * t51 + t42) - m(7) * (t20 * t51 + t53 * t54 + t42) + (-m(2) - m(3)) * t25 + (-mrSges(3,2) - m(6) * (-pkin(8) + pkin(9)) - m(7) * (-pkin(8) - t26) + t58) * t33 + (t63 * (t27 * t28 + t31 * t32) + t60 * (-t27 * t32 + t28 * t31) + t61 * t28 + t64 * t32 - mrSges(3,1)) * t29) * g(3);
U  = t1;
