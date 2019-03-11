% Calculate potential energy for
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:37
% EndTime: 2019-03-08 20:40:38
% DurationCPUTime: 0.90s
% Computational Cost: add. (305->99), mult. (496->102), div. (0->0), fcn. (547->12), ass. (0->46)
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t75 = t39 * mrSges(5,1) - t36 * mrSges(5,2);
t72 = -m(1) - m(2);
t71 = -m(4) - m(5);
t70 = m(6) + m(7);
t69 = mrSges(4,1) + mrSges(3,3);
t68 = -m(5) * pkin(3) - t69;
t67 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t66 = t39 * mrSges(5,2) - mrSges(3,2) + mrSges(4,3);
t32 = sin(pkin(6));
t65 = t75 * t32 - mrSges(2,2);
t35 = sin(qJ(6));
t38 = cos(qJ(6));
t64 = -m(7) * pkin(5) - t38 * mrSges(7,1) + t35 * mrSges(7,2) - mrSges(6,1);
t63 = -m(5) * pkin(8) - t35 * mrSges(7,1) - t38 * mrSges(7,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t31 = sin(pkin(11));
t33 = cos(pkin(11));
t37 = sin(qJ(2));
t34 = cos(pkin(6));
t40 = cos(qJ(2));
t53 = t34 * t40;
t14 = t31 * t37 - t33 * t53;
t61 = t14 * t36;
t16 = t31 * t53 + t33 * t37;
t60 = t16 * t36;
t59 = t31 * t32;
t58 = t32 * t33;
t55 = t32 * t40;
t54 = t34 * t37;
t52 = pkin(7) * t58;
t51 = t31 * pkin(1) + r_base(2);
t50 = qJ(1) + r_base(3);
t49 = t33 * pkin(1) + pkin(7) * t59 + r_base(1);
t48 = t34 * pkin(7) + t50;
t47 = t32 * t37 * pkin(2) + t48;
t15 = t31 * t40 + t33 * t54;
t45 = t15 * pkin(2) + qJ(3) * t14 + t51;
t17 = -t31 * t54 + t33 * t40;
t44 = t17 * pkin(2) + qJ(3) * t16 + t49;
t41 = -pkin(9) - pkin(8);
t30 = qJ(4) + qJ(5);
t26 = cos(t30);
t25 = sin(t30);
t24 = pkin(4) * t39 + pkin(3);
t1 = (-m(1) * r_base(3) - m(2) * t50 - m(3) * t48 - mrSges(1,3) - mrSges(2,3) + t64 * (-t25 * t55 + t26 * t34) - t70 * (t34 * t24 + t47) - t67 * (t25 * t34 + t26 * t55) + t71 * t47 + (t68 - t75) * t34 + (((t70 - t71) * qJ(3) + t66 + (t70 * pkin(4) + mrSges(5,1)) * t36) * t40 + (t70 * t41 + t63) * t37) * t32) * g(3) + (-mrSges(1,2) - t31 * mrSges(2,1) - m(3) * (t51 - t52) - m(4) * (t45 - t52) - m(5) * t45 - t61 * mrSges(5,1) + t72 * r_base(2) + (-m(5) * (-pkin(3) - pkin(7)) + t69) * t58 - t70 * (-t15 * t41 + pkin(4) * t61 + (-pkin(7) - t24) * t58 + t45) + t65 * t33 + t67 * (t14 * t26 + t25 * t58) - t66 * t14 + t64 * (t14 * t25 - t26 * t58) + t63 * t15) * g(2) + (-m(3) * t49 - t33 * mrSges(2,1) - t60 * mrSges(5,1) - mrSges(1,1) + t72 * r_base(1) + t71 * t44 + t68 * t59 - t70 * (pkin(4) * t60 - t17 * t41 + t24 * t59 + t44) - t65 * t31 - t66 * t16 - t67 * (-t16 * t26 + t25 * t59) + t64 * (t16 * t25 + t26 * t59) + t63 * t17) * g(1);
U  = t1;
