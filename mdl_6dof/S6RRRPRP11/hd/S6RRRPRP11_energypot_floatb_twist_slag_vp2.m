% Calculate potential energy for
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:18
% EndTime: 2019-03-09 17:41:18
% DurationCPUTime: 0.67s
% Computational Cost: add. (310->103), mult. (626->112), div. (0->0), fcn. (727->10), ass. (0->52)
t75 = -m(1) - m(2);
t68 = pkin(4) + pkin(9);
t74 = -mrSges(6,1) - mrSges(7,1);
t73 = -mrSges(6,2) - mrSges(7,2);
t34 = sin(qJ(5));
t72 = -m(7) * (pkin(5) * t34 + qJ(4)) + mrSges(4,2) - mrSges(5,3);
t38 = cos(qJ(5));
t71 = m(6) * t68 + m(7) * (pkin(5) * t38 + t68) + mrSges(5,1) + mrSges(4,3);
t70 = mrSges(3,2) - t71;
t69 = -m(6) * pkin(10) + m(7) * (-qJ(6) - pkin(10)) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t40 = cos(qJ(1));
t37 = sin(qJ(1));
t59 = cos(pkin(6));
t53 = t37 * t59;
t20 = t40 * t36 + t39 * t53;
t67 = pkin(9) * t20;
t52 = t40 * t59;
t18 = t36 * t37 - t39 * t52;
t66 = t18 * pkin(9);
t64 = cos(qJ(3));
t32 = sin(pkin(6));
t63 = t32 * t36;
t62 = t32 * t37;
t61 = t32 * t39;
t60 = t32 * t40;
t58 = pkin(7) + r_base(3);
t57 = pkin(9) * t61;
t56 = t32 * t64;
t55 = t59 * pkin(8) + t58;
t51 = t40 * pkin(1) + pkin(8) * t62 + r_base(1);
t50 = pkin(2) * t63 + t55;
t21 = -t36 * t53 + t40 * t39;
t49 = t21 * pkin(2) + t51;
t35 = sin(qJ(3));
t17 = t59 * t35 + t36 * t56;
t48 = t17 * pkin(3) + t50;
t12 = t21 * t64 + t35 * t62;
t47 = t12 * pkin(3) + t49;
t46 = t37 * pkin(1) - pkin(8) * t60 + r_base(2);
t19 = t36 * t52 + t37 * t39;
t45 = t19 * pkin(2) + t46;
t10 = t19 * t64 - t35 * t60;
t44 = t10 * pkin(3) + t45;
t16 = t35 * t63 - t59 * t64;
t43 = t16 * qJ(4) + t48;
t11 = t21 * t35 - t37 * t56;
t42 = qJ(4) * t11 + t47;
t9 = t19 * t35 + t40 * t56;
t41 = t9 * qJ(4) + t44;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t58 - mrSges(2,3) - m(3) * t55 - t59 * mrSges(3,3) - (t36 * mrSges(3,1) + t39 * mrSges(3,2)) * t32 - m(4) * (t50 - t57) - m(5) * (t43 - t57) - m(6) * t43 - m(7) * t48 + t74 * (t16 * t34 - t38 * t61) + t73 * (t16 * t38 + t34 * t61) + t71 * t61 + t72 * t16 + t69 * t17) * g(3) + (-m(3) * t46 - mrSges(1,2) + mrSges(3,3) * t60 - m(7) * t44 - m(4) * (t45 + t66) - m(5) * (t41 + t66) - m(6) * t41 - t19 * mrSges(3,1) - t37 * mrSges(2,1) - t40 * mrSges(2,2) + t75 * r_base(2) + t72 * t9 + t74 * (t18 * t38 + t34 * t9) + t73 * (-t18 * t34 + t38 * t9) + t70 * t18 + t69 * t10) * g(2) + (-m(3) * t51 - mrSges(1,1) - m(7) * t47 - m(4) * (t49 + t67) - m(5) * (t42 + t67) - m(6) * t42 - mrSges(3,3) * t62 - t21 * mrSges(3,1) + t37 * mrSges(2,2) - t40 * mrSges(2,1) + t75 * r_base(1) + t74 * (t11 * t34 + t20 * t38) + t73 * (t11 * t38 - t20 * t34) + t72 * t11 + t70 * t20 + t69 * t12) * g(1);
U  = t1;
