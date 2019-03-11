% Calculate potential energy for
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR11_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:38:45
% EndTime: 2019-03-09 09:38:46
% DurationCPUTime: 0.83s
% Computational Cost: add. (305->98), mult. (496->99), div. (0->0), fcn. (547->12), ass. (0->46)
t71 = -m(1) - m(2);
t70 = -m(4) - m(5);
t69 = m(6) + m(7);
t68 = m(7) * pkin(10) - mrSges(6,2) + mrSges(7,3);
t33 = cos(pkin(11));
t67 = t33 * mrSges(5,2) - mrSges(3,2) + mrSges(4,3);
t31 = sin(pkin(11));
t66 = t33 * mrSges(5,1) - t31 * mrSges(5,2) + mrSges(4,1) + mrSges(3,3);
t36 = sin(qJ(6));
t39 = cos(qJ(6));
t65 = -m(7) * pkin(5) - t39 * mrSges(7,1) + t36 * mrSges(7,2) - mrSges(6,1);
t64 = -m(5) * pkin(3) - t66;
t63 = -m(5) * qJ(4) - t36 * mrSges(7,1) - t39 * mrSges(7,2) - mrSges(3,1) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t34 = cos(pkin(6));
t40 = cos(qJ(2));
t41 = cos(qJ(1));
t53 = t40 * t41;
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t56 = t37 * t38;
t14 = -t34 * t53 + t56;
t61 = t14 * t31;
t54 = t38 * t40;
t55 = t37 * t41;
t16 = t34 * t54 + t55;
t60 = t16 * t31;
t32 = sin(pkin(6));
t59 = t32 * t38;
t58 = t32 * t40;
t57 = t32 * t41;
t52 = pkin(7) + r_base(3);
t51 = pkin(8) * t57;
t50 = t38 * pkin(1) + r_base(2);
t49 = t34 * pkin(8) + t52;
t48 = t41 * pkin(1) + pkin(8) * t59 + r_base(1);
t47 = t32 * t37 * pkin(2) + t49;
t15 = t34 * t55 + t54;
t45 = t15 * pkin(2) + t14 * qJ(3) + t50;
t17 = -t34 * t56 + t53;
t44 = t17 * pkin(2) + qJ(3) * t16 + t48;
t35 = -pkin(9) - qJ(4);
t30 = pkin(11) + qJ(5);
t26 = cos(t30);
t25 = sin(t30);
t24 = pkin(4) * t33 + pkin(3);
t1 = (-m(1) * r_base(3) - m(2) * t52 - m(3) * t49 - mrSges(1,3) - mrSges(2,3) + t65 * (-t25 * t58 + t26 * t34) - t69 * (t34 * t24 + t47) - t68 * (t25 * t34 + t26 * t58) + t70 * t47 + t64 * t34 + (((t69 - t70) * qJ(3) + t67 + (t69 * pkin(4) + mrSges(5,1)) * t31) * t40 + (t69 * t35 + t63) * t37) * t32) * g(3) + (-t61 * mrSges(5,1) - m(5) * t45 - m(3) * (t50 - t51) - m(4) * (t45 - t51) - mrSges(1,2) - t38 * mrSges(2,1) - t41 * mrSges(2,2) + t71 * r_base(2) - t69 * (-t15 * t35 + pkin(4) * t61 + (-pkin(8) - t24) * t57 + t45) + t68 * (t14 * t26 + t25 * t57) + (-m(5) * (-pkin(3) - pkin(8)) + t66) * t57 - t67 * t14 + t65 * (t14 * t25 - t26 * t57) + t63 * t15) * g(2) + (-m(3) * t48 - t41 * mrSges(2,1) - t60 * mrSges(5,1) + t38 * mrSges(2,2) - mrSges(1,1) + t71 * r_base(1) + t70 * t44 - t69 * (pkin(4) * t60 - t17 * t35 + t24 * t59 + t44) + t64 * t59 - t67 * t16 - t68 * (-t16 * t26 + t25 * t59) + t65 * (t16 * t25 + t26 * t59) + t63 * t17) * g(1);
U  = t1;
