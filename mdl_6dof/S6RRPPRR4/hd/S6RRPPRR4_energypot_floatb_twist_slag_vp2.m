% Calculate potential energy for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:04
% EndTime: 2019-03-09 09:01:05
% DurationCPUTime: 0.71s
% Computational Cost: add. (357->97), mult. (738->114), div. (0->0), fcn. (881->12), ass. (0->50)
t75 = -m(6) - m(7);
t74 = mrSges(4,2) - mrSges(5,3);
t73 = -m(3) - m(1) - m(2);
t72 = -mrSges(3,3) - mrSges(4,3) - mrSges(5,1);
t71 = -m(3) * pkin(1) - mrSges(2,1);
t70 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t69 = m(3) * pkin(8) - t72;
t37 = sin(qJ(6));
t41 = cos(qJ(6));
t68 = -t37 * mrSges(7,1) - t41 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t67 = -m(7) * pkin(5) - t41 * mrSges(7,1) + t37 * mrSges(7,2) - mrSges(6,1);
t34 = sin(pkin(6));
t39 = sin(qJ(2));
t66 = t34 * t39;
t40 = sin(qJ(1));
t65 = t34 * t40;
t44 = cos(qJ(1));
t64 = t34 * t44;
t35 = cos(pkin(11));
t43 = cos(qJ(2));
t63 = t35 * t43;
t62 = t39 * t44;
t61 = t40 * t39;
t60 = t40 * t43;
t59 = t43 * t44;
t58 = pkin(7) + r_base(3);
t36 = cos(pkin(6));
t57 = t36 * pkin(8) + t58;
t21 = pkin(2) * t36 * t39 + (-pkin(8) - qJ(3)) * t34;
t29 = pkin(2) * t43 + pkin(1);
t56 = t44 * t21 + t40 * t29 + r_base(2);
t33 = sin(pkin(11));
t55 = t33 * t43 + t35 * t39;
t23 = -t33 * t39 + t63;
t54 = -t21 * t40 + t44 * t29 + r_base(1);
t53 = pkin(2) * t66 + t36 * qJ(3) + t57;
t52 = t55 * t36;
t51 = t23 * t36;
t8 = -t40 * t55 + t44 * t51;
t9 = t40 * t23 + t44 * t52;
t50 = t9 * pkin(3) - t8 * qJ(4) + t56;
t10 = -t40 * t51 - t44 * t55;
t11 = t23 * t44 - t40 * t52;
t49 = t11 * pkin(3) - qJ(4) * t10 + t54;
t19 = t33 * t66 - t34 * t63;
t20 = t55 * t34;
t48 = t20 * pkin(3) + qJ(4) * t19 + t53;
t42 = cos(qJ(5));
t38 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t58 - mrSges(2,3) - m(3) * t57 - (t39 * mrSges(3,1) + t43 * mrSges(3,2)) * t34 - m(4) * t53 - m(5) * t48 + t75 * (t36 * pkin(4) + pkin(9) * t20 + t48) + t74 * t19 - t70 * (-t19 * t42 + t36 * t38) + t72 * t36 + t67 * (t19 * t38 + t36 * t42) + t68 * t20) * g(3) + (-(t36 * t59 - t61) * mrSges(3,2) - (t36 * t62 + t60) * mrSges(3,1) - m(5) * t50 - m(4) * t56 - mrSges(1,2) - t44 * mrSges(2,2) + t71 * t40 + t73 * r_base(2) - t74 * t8 + t75 * (-pkin(4) * t64 + t9 * pkin(9) + t50) + t70 * (t38 * t64 - t8 * t42) + t67 * (-t8 * t38 - t42 * t64) + t68 * t9 + t69 * t64) * g(2) + (-(-t36 * t61 + t59) * mrSges(3,1) - (-t36 * t60 - t62) * mrSges(3,2) - m(5) * t49 - m(4) * t54 - mrSges(1,1) + t40 * mrSges(2,2) + t71 * t44 + t73 * r_base(1) + t75 * (pkin(4) * t65 + pkin(9) * t11 + t49) - t74 * t10 - t70 * (t10 * t42 + t38 * t65) - t69 * t65 + t67 * (-t10 * t38 + t42 * t65) + t68 * t11) * g(1);
U  = t1;
