% Calculate potential energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:22
% EndTime: 2019-03-09 16:21:22
% DurationCPUTime: 0.70s
% Computational Cost: add. (322->98), mult. (626->102), div. (0->0), fcn. (727->12), ass. (0->53)
t70 = -m(1) - m(2);
t65 = pkin(4) + pkin(9);
t69 = -m(6) * qJ(5) + m(7) * (-pkin(10) - qJ(5)) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t28 = pkin(11) + qJ(6);
t23 = sin(t28);
t24 = cos(t28);
t29 = sin(pkin(11));
t31 = cos(pkin(11));
t68 = -m(7) * (pkin(5) * t29 + qJ(4)) - t29 * mrSges(6,1) - t23 * mrSges(7,1) - t31 * mrSges(6,2) - t24 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t67 = m(6) * t65 + m(7) * (pkin(5) * t31 + t65) + t31 * mrSges(6,1) + t24 * mrSges(7,1) - t29 * mrSges(6,2) - t23 * mrSges(7,2) + mrSges(5,1) + mrSges(4,3);
t66 = mrSges(3,2) - t67;
t34 = sin(qJ(2));
t36 = cos(qJ(2));
t37 = cos(qJ(1));
t35 = sin(qJ(1));
t56 = cos(pkin(6));
t50 = t35 * t56;
t14 = t37 * t34 + t36 * t50;
t64 = pkin(9) * t14;
t49 = t37 * t56;
t12 = t34 * t35 - t36 * t49;
t63 = t12 * pkin(9);
t61 = cos(qJ(3));
t30 = sin(pkin(6));
t60 = t30 * t34;
t59 = t30 * t35;
t58 = t30 * t36;
t57 = t30 * t37;
t55 = pkin(7) + r_base(3);
t54 = pkin(9) * t58;
t53 = t30 * t61;
t52 = t56 * pkin(8) + t55;
t48 = t37 * pkin(1) + pkin(8) * t59 + r_base(1);
t47 = pkin(2) * t60 + t52;
t15 = -t34 * t50 + t37 * t36;
t46 = t15 * pkin(2) + t48;
t33 = sin(qJ(3));
t11 = t56 * t33 + t34 * t53;
t45 = t11 * pkin(3) + t47;
t6 = t15 * t61 + t33 * t59;
t44 = t6 * pkin(3) + t46;
t43 = t35 * pkin(1) - pkin(8) * t57 + r_base(2);
t13 = t34 * t49 + t35 * t36;
t42 = t13 * pkin(2) + t43;
t4 = t13 * t61 - t33 * t57;
t41 = t4 * pkin(3) + t42;
t5 = t15 * t33 - t35 * t53;
t40 = qJ(4) * t5 + t44;
t10 = t33 * t60 - t56 * t61;
t39 = qJ(4) * t10 + t45;
t3 = t13 * t33 + t37 * t53;
t38 = t3 * qJ(4) + t41;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t55 - mrSges(2,3) - m(3) * t52 - t56 * mrSges(3,3) - (t34 * mrSges(3,1) + t36 * mrSges(3,2)) * t30 - m(4) * (t47 - t54) - m(5) * (t39 - t54) - m(6) * t39 - m(7) * t45 + t67 * t58 + t68 * t10 + t69 * t11) * g(3) + (-m(3) * t43 - mrSges(1,2) + mrSges(3,3) * t57 - m(7) * t41 - m(4) * (t42 + t63) - m(5) * (t38 + t63) - m(6) * t38 - t13 * mrSges(3,1) - t35 * mrSges(2,1) - t37 * mrSges(2,2) + t70 * r_base(2) + t68 * t3 + t66 * t12 + t69 * t4) * g(2) + (-m(3) * t48 - mrSges(1,1) - m(7) * t44 - m(4) * (t46 + t64) - m(5) * (t40 + t64) - m(6) * t40 - mrSges(3,3) * t59 - t15 * mrSges(3,1) + t35 * mrSges(2,2) - t37 * mrSges(2,1) + t70 * r_base(1) + t68 * t5 + t66 * t14 + t69 * t6) * g(1);
U  = t1;
