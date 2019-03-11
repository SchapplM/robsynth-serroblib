% Calculate potential energy for
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:28:23
% EndTime: 2019-03-09 09:28:24
% DurationCPUTime: 0.91s
% Computational Cost: add. (249->102), mult. (473->104), div. (0->0), fcn. (516->10), ass. (0->50)
t73 = -m(1) - m(2);
t72 = pkin(9) - qJ(3);
t71 = -mrSges(3,1) + mrSges(4,2) - mrSges(5,3);
t70 = -mrSges(3,2) + mrSges(5,2) + mrSges(4,3);
t69 = -mrSges(3,3) - mrSges(4,1) - mrSges(5,1);
t68 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t67 = -m(7) * pkin(5) - t34 * mrSges(7,1) + t30 * mrSges(7,2) - mrSges(6,1);
t66 = m(6) * t72 + t30 * mrSges(7,1) + t34 * mrSges(7,2) + mrSges(6,3) - t70;
t65 = -pkin(3) - pkin(8);
t28 = sin(pkin(6));
t32 = sin(qJ(2));
t64 = t28 * t32;
t33 = sin(qJ(1));
t63 = t28 * t33;
t35 = cos(qJ(5));
t62 = t28 * t35;
t36 = cos(qJ(2));
t61 = t28 * t36;
t37 = cos(qJ(1));
t60 = t28 * t37;
t59 = t32 * t33;
t58 = t32 * t37;
t57 = t33 * t36;
t56 = t36 * t37;
t29 = cos(pkin(6));
t12 = t29 * t57 + t58;
t54 = qJ(3) * t12;
t10 = -t29 * t56 + t59;
t53 = t10 * qJ(3);
t52 = pkin(7) + r_base(3);
t51 = t33 * pkin(1) + r_base(2);
t50 = t29 * pkin(8) + t52;
t49 = t37 * pkin(1) + pkin(8) * t63 + r_base(1);
t48 = pkin(2) * t64 + t50;
t13 = -t29 * t59 + t56;
t47 = t13 * pkin(2) + t49;
t46 = (-pkin(4) + t65) * t60;
t44 = -pkin(8) * t60 + t51;
t11 = t29 * t58 + t57;
t6 = t11 * pkin(2);
t43 = t11 * qJ(4) + t51 + t6;
t42 = t29 * pkin(3) + qJ(4) * t64 + t48;
t41 = pkin(3) * t63 + qJ(4) * t13 + t47;
t40 = t43 + t53;
t38 = t29 * pkin(4) + t72 * t61 + t42;
t31 = sin(qJ(5));
t9 = t29 * t35 + t31 * t64;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t52 - mrSges(2,3) - m(3) * t50 - m(4) * t48 - m(5) * t42 - m(6) * t38 - t9 * mrSges(6,1) - mrSges(6,3) * t61 - m(7) * (pkin(5) * t9 + t38) - (t30 * t61 + t34 * t9) * mrSges(7,1) - (-t30 * t9 + t34 * t61) * mrSges(7,2) - t68 * (t29 * t31 - t32 * t62) + t69 * t29 + (((m(4) + m(5)) * qJ(3) + t70) * t36 + t71 * t32) * t28) * g(3) + (-m(3) * t44 - m(5) * t40 - (t43 + t46) * m(6) - m(4) * (t44 + t6 + t53) - mrSges(1,2) - t33 * mrSges(2,1) - t37 * mrSges(2,2) - (t40 + t46) * m(7) + t73 * r_base(2) + t68 * (t11 * t35 + t31 * t60) + (-m(5) * t65 - t69) * t60 + t67 * (t11 * t31 - t35 * t60) + t71 * t11 + (m(7) * pkin(9) + t66) * t10) * g(2) + (-m(3) * t49 - m(4) * (t47 + t54) - m(5) * (t41 + t54) - mrSges(1,1) + t33 * mrSges(2,2) - t37 * mrSges(2,1) + t73 * r_base(1) + (-m(6) - m(7)) * (pkin(4) * t63 + t41) - t68 * (-t13 * t35 + t31 * t63) + t69 * t63 + t67 * (t13 * t31 + t33 * t62) + t71 * t13 + (m(7) * t72 + t66) * t12) * g(1);
U  = t1;
