% Calculate potential energy for
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR14_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:27
% EndTime: 2019-03-09 11:32:28
% DurationCPUTime: 0.68s
% Computational Cost: add. (270->100), mult. (528->100), div. (0->0), fcn. (592->10), ass. (0->55)
t74 = -m(1) - m(2);
t30 = cos(pkin(6));
t32 = sin(qJ(4));
t36 = cos(qJ(4));
t29 = sin(pkin(6));
t37 = cos(qJ(2));
t60 = t29 * t37;
t12 = t30 * t32 + t36 * t60;
t13 = t30 * t36 - t32 * t60;
t73 = t13 * pkin(4) + qJ(5) * t12;
t72 = mrSges(4,2) - mrSges(3,1);
t71 = mrSges(3,3) + mrSges(4,1);
t70 = -mrSges(4,3) + mrSges(3,2);
t69 = m(7) * pkin(10) + mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t31 = sin(qJ(6));
t35 = cos(qJ(6));
t68 = -t35 * mrSges(7,1) + t31 * mrSges(7,2) - mrSges(6,1) - mrSges(5,3);
t67 = t31 * mrSges(7,1) + t35 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t66 = -m(7) * (pkin(5) + pkin(9)) + t68 + t72;
t38 = cos(qJ(1));
t55 = t37 * t38;
t33 = sin(qJ(2));
t34 = sin(qJ(1));
t58 = t33 * t34;
t17 = -t30 * t58 + t55;
t64 = pkin(9) * t17;
t56 = t34 * t37;
t57 = t33 * t38;
t15 = t30 * t57 + t56;
t63 = t15 * pkin(9);
t62 = t29 * t33;
t61 = t29 * t34;
t59 = t29 * t38;
t54 = qJ(3) * t37;
t52 = pkin(7) + r_base(3);
t51 = pkin(8) * t59;
t50 = t34 * pkin(1) + r_base(2);
t49 = t30 * pkin(8) + t52;
t48 = t38 * pkin(1) + pkin(8) * t61 + r_base(1);
t47 = pkin(2) * t62 + t49;
t46 = t30 * pkin(3) + pkin(9) * t62 + t47;
t14 = -t30 * t55 + t58;
t45 = t15 * pkin(2) + t14 * qJ(3) + t50;
t16 = t30 * t56 + t57;
t44 = t17 * pkin(2) + qJ(3) * t16 + t48;
t43 = pkin(3) * t61 + t44;
t42 = -t54 * t29 + t46;
t41 = (-pkin(3) - pkin(8)) * t59 + t45;
t3 = -t16 * t36 + t32 * t61;
t4 = t16 * t32 + t36 * t61;
t40 = t4 * pkin(4) + qJ(5) * t3 + t43;
t5 = t14 * t36 + t32 * t59;
t6 = -t14 * t32 + t36 * t59;
t39 = -t6 * pkin(4) - t5 * qJ(5) + t41;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t52 - mrSges(2,3) - m(3) * t49 - m(4) * t47 - m(5) * t42 - m(6) * (t42 + t73) - m(7) * (t46 + t73) - t71 * t30 + (m(7) * t54 + (m(4) * qJ(3) - t70) * t37 + (-m(7) * pkin(5) + t72) * t33) * t29 + t68 * t62 - t67 * t12 - t69 * t13) * g(3) + (-m(5) * (t41 + t63) - m(6) * (t39 + t63) - m(7) * t39 - m(3) * (t50 - t51) - m(4) * (t45 - t51) - mrSges(1,2) - t34 * mrSges(2,1) - t38 * mrSges(2,2) + t74 * r_base(2) + t71 * t59 + t70 * t14 + t69 * t6 + t67 * t5 + t66 * t15) * g(2) + (-m(5) * (t43 + t64) - m(6) * (t40 + t64) - m(7) * t40 - m(4) * t44 - m(3) * t48 - mrSges(1,1) + t34 * mrSges(2,2) - t38 * mrSges(2,1) + t74 * r_base(1) - t71 * t61 + t70 * t16 - t69 * t4 - t67 * t3 + t66 * t17) * g(1);
U  = t1;
