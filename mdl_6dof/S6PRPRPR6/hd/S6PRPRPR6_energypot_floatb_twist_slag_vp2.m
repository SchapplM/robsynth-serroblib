% Calculate potential energy for
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:33
% EndTime: 2019-03-08 19:46:34
% DurationCPUTime: 0.77s
% Computational Cost: add. (293->108), mult. (550->119), div. (0->0), fcn. (622->12), ass. (0->50)
t71 = -m(1) - m(2);
t70 = -m(5) - m(6);
t69 = mrSges(4,2) - mrSges(3,1);
t68 = mrSges(3,3) + mrSges(4,1);
t67 = -mrSges(4,3) + mrSges(3,2);
t66 = m(6) * qJ(5) - m(7) * (-pkin(9) - qJ(5)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t27 = pkin(11) + qJ(6);
t21 = sin(t27);
t22 = cos(t27);
t28 = sin(pkin(11));
t31 = cos(pkin(11));
t63 = pkin(5) * t28;
t65 = -m(7) * (pkin(8) + t63) - t28 * mrSges(6,1) - t21 * mrSges(7,1) - t31 * mrSges(6,2) - t22 * mrSges(7,2) - mrSges(5,3) + t69;
t20 = pkin(5) * t31 + pkin(4);
t64 = -m(6) * pkin(4) - m(7) * t20 - t31 * mrSges(6,1) - t22 * mrSges(7,1) + t28 * mrSges(6,2) + t21 * mrSges(7,2) - mrSges(5,1);
t29 = sin(pkin(10));
t30 = sin(pkin(6));
t62 = t29 * t30;
t32 = cos(pkin(10));
t61 = t30 * t32;
t35 = sin(qJ(4));
t60 = t30 * t35;
t36 = sin(qJ(2));
t59 = t30 * t36;
t37 = cos(qJ(4));
t58 = t30 * t37;
t38 = cos(qJ(2));
t57 = t30 * t38;
t33 = cos(pkin(6));
t56 = t33 * t36;
t55 = t33 * t38;
t54 = qJ(3) * t38;
t53 = pkin(7) * t61;
t52 = t29 * pkin(1) + r_base(2);
t51 = qJ(1) + r_base(3);
t49 = t32 * pkin(1) + pkin(7) * t62 + r_base(1);
t48 = t33 * pkin(7) + t51;
t47 = pkin(2) * t59 + t48;
t8 = t29 * t36 - t32 * t55;
t9 = t29 * t38 + t32 * t56;
t46 = t9 * pkin(2) + qJ(3) * t8 + t52;
t45 = t33 * pkin(3) + pkin(8) * t59 + t47;
t10 = t29 * t55 + t32 * t36;
t11 = -t29 * t56 + t32 * t38;
t44 = t11 * pkin(2) + qJ(3) * t10 + t49;
t43 = pkin(3) * t62 + t44;
t42 = -t54 * t30 + t45;
t40 = (-pkin(3) - pkin(7)) * t61 + t46;
t13 = t33 * t37 - t35 * t57;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t51 - mrSges(2,3) - m(3) * t48 - m(4) * t47 - m(5) * t42 - t13 * mrSges(5,1) - mrSges(5,3) * t59 - m(6) * (t13 * pkin(4) + t42) - (t13 * t31 + t28 * t59) * mrSges(6,1) - (-t13 * t28 + t31 * t59) * mrSges(6,2) - m(7) * (t13 * t20 + t45) - (t13 * t22 + t21 * t59) * mrSges(7,1) - (-t13 * t21 + t22 * t59) * mrSges(7,2) - t68 * t33 + (m(7) * t54 + (m(4) * qJ(3) - t67) * t38 + (-m(7) * t63 + t69) * t36) * t30 - t66 * (t33 * t35 + t37 * t57)) * g(3) + (-mrSges(1,2) - m(7) * t40 - m(3) * (t52 - t53) - m(4) * (t46 - t53) - t29 * mrSges(2,1) - t32 * mrSges(2,2) + t71 * r_base(2) + t67 * t8 + t68 * t61 + t70 * (pkin(8) * t9 + t40) + t64 * (-t32 * t58 + t35 * t8) + t65 * t9 + t66 * (t32 * t60 + t37 * t8)) * g(2) + (-m(3) * t49 - m(4) * t44 - m(7) * t43 - t32 * mrSges(2,1) + t29 * mrSges(2,2) - mrSges(1,1) + t71 * r_base(1) - t68 * t62 + t70 * (pkin(8) * t11 + t43) + t67 * t10 + t64 * (t10 * t35 + t29 * t58) + t65 * t11 - t66 * (-t10 * t37 + t29 * t60)) * g(1);
U  = t1;
