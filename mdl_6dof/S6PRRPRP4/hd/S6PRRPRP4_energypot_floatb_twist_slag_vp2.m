% Calculate potential energy for
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:05
% EndTime: 2019-03-08 21:40:05
% DurationCPUTime: 0.67s
% Computational Cost: add. (310->103), mult. (626->113), div. (0->0), fcn. (727->10), ass. (0->53)
t76 = -m(1) - m(2);
t69 = pkin(4) + pkin(8);
t75 = -mrSges(6,1) - mrSges(7,1);
t74 = -mrSges(6,2) - mrSges(7,2);
t36 = sin(qJ(5));
t73 = -m(7) * (pkin(5) * t36 + qJ(4)) + mrSges(4,2) - mrSges(5,3);
t39 = cos(qJ(5));
t72 = m(6) * t69 + m(7) * (pkin(5) * t39 + t69) + mrSges(5,1) + mrSges(4,3);
t71 = mrSges(3,2) - t72;
t70 = -m(6) * pkin(9) + m(7) * (-qJ(6) - pkin(9)) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t32 = sin(pkin(10));
t34 = cos(pkin(10));
t38 = sin(qJ(2));
t40 = cos(qJ(2));
t59 = cos(pkin(6));
t53 = t40 * t59;
t16 = t32 * t38 - t34 * t53;
t68 = pkin(8) * t16;
t18 = t32 * t53 + t34 * t38;
t67 = pkin(8) * t18;
t65 = cos(qJ(3));
t33 = sin(pkin(6));
t64 = t32 * t33;
t63 = t33 * t34;
t37 = sin(qJ(3));
t62 = t33 * t37;
t61 = t33 * t38;
t60 = t33 * t40;
t58 = pkin(8) * t60;
t57 = qJ(1) + r_base(3);
t56 = t33 * t65;
t54 = t38 * t59;
t52 = t34 * pkin(1) + pkin(7) * t64 + r_base(1);
t51 = t59 * pkin(7) + t57;
t19 = -t32 * t54 + t34 * t40;
t50 = t19 * pkin(2) + t52;
t49 = pkin(2) * t61 + t51;
t10 = t19 * t65 + t32 * t62;
t48 = t10 * pkin(3) + t50;
t21 = t59 * t37 + t38 * t56;
t47 = t21 * pkin(3) + t49;
t46 = t32 * pkin(1) - pkin(7) * t63 + r_base(2);
t17 = t32 * t40 + t34 * t54;
t45 = t17 * pkin(2) + t46;
t8 = t17 * t65 - t34 * t62;
t44 = t8 * pkin(3) + t45;
t9 = t19 * t37 - t32 * t56;
t43 = qJ(4) * t9 + t48;
t20 = t37 * t61 - t59 * t65;
t42 = t20 * qJ(4) + t47;
t7 = t17 * t37 + t34 * t56;
t41 = qJ(4) * t7 + t44;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t57 - mrSges(2,3) - m(3) * t51 - t59 * mrSges(3,3) - (t38 * mrSges(3,1) + t40 * mrSges(3,2)) * t33 - m(4) * (t49 - t58) - m(5) * (t42 - t58) - m(6) * t42 - m(7) * t47 + t72 * t60 + t73 * t20 + t75 * (t20 * t36 - t39 * t60) + t74 * (t20 * t39 + t36 * t60) + t70 * t21) * g(3) + (-m(3) * t46 - mrSges(1,2) + mrSges(3,3) * t63 - m(7) * t44 - m(4) * (t45 + t68) - m(5) * (t41 + t68) - m(6) * t41 - t17 * mrSges(3,1) - t32 * mrSges(2,1) - t34 * mrSges(2,2) + t76 * r_base(2) + t73 * t7 + t75 * (t16 * t39 + t36 * t7) + t74 * (-t16 * t36 + t39 * t7) + t71 * t16 + t70 * t8) * g(2) + (-m(3) * t52 - mrSges(1,1) - m(7) * t48 - m(4) * (t50 + t67) - m(5) * (t43 + t67) - m(6) * t43 - mrSges(3,3) * t64 - t19 * mrSges(3,1) + t32 * mrSges(2,2) - t34 * mrSges(2,1) + t76 * r_base(1) + t73 * t9 + t75 * (t18 * t39 + t36 * t9) + t74 * (-t18 * t36 + t39 * t9) + t71 * t18 + t70 * t10) * g(1);
U  = t1;
