% Calculate potential energy for
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:29
% EndTime: 2019-03-09 16:03:30
% DurationCPUTime: 0.65s
% Computational Cost: add. (299->98), mult. (604->98), div. (0->0), fcn. (697->10), ass. (0->51)
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t70 = -t30 * mrSges(7,1) - t34 * mrSges(7,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t69 = -m(1) - m(2);
t31 = sin(qJ(3));
t32 = sin(qJ(2));
t29 = sin(pkin(6));
t61 = cos(qJ(3));
t52 = t29 * t61;
t54 = cos(pkin(6));
t14 = t31 * t54 + t32 * t52;
t35 = cos(qJ(2));
t58 = t29 * t35;
t68 = t14 * pkin(4) + qJ(5) * t58;
t67 = -m(7) * pkin(10) - mrSges(4,1) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t65 = -m(7) * (pkin(5) + qJ(4)) - t34 * mrSges(7,1) + t30 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t64 = mrSges(3,2) - t70 + (-m(6) - m(7)) * (pkin(9) - qJ(5));
t36 = cos(qJ(1));
t33 = sin(qJ(1));
t50 = t33 * t54;
t17 = t36 * t32 + t35 * t50;
t63 = pkin(9) * t17;
t49 = t36 * t54;
t15 = t32 * t33 - t35 * t49;
t62 = t15 * pkin(9);
t60 = t29 * t32;
t59 = t29 * t33;
t57 = t29 * t36;
t53 = pkin(7) + r_base(3);
t51 = t54 * pkin(8) + t53;
t48 = t36 * pkin(1) + pkin(8) * t59 + r_base(1);
t18 = -t32 * t50 + t36 * t35;
t47 = t18 * pkin(2) + t48;
t8 = t18 * t61 + t31 * t59;
t45 = t8 * pkin(3) + t47;
t44 = t33 * pkin(1) - pkin(8) * t57 + r_base(2);
t16 = t32 * t49 + t33 * t35;
t43 = t16 * pkin(2) + t44;
t42 = pkin(2) * t60 - pkin(9) * t58 + t51;
t6 = t16 * t61 - t31 * t57;
t41 = t6 * pkin(3) + t43;
t7 = t18 * t31 - t33 * t52;
t40 = qJ(4) * t7 + t45;
t39 = t14 * pkin(3) + t42;
t5 = t16 * t31 + t36 * t52;
t38 = t5 * qJ(4) + t41;
t13 = t31 * t60 - t54 * t61;
t37 = qJ(4) * t13 + t39;
t3 = t8 * pkin(4);
t1 = t6 * pkin(4);
t2 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t53 - mrSges(2,3) - m(3) * t51 - t54 * mrSges(3,3) - (t32 * mrSges(3,1) + t35 * mrSges(3,2)) * t29 - m(4) * t42 - m(5) * t37 - m(6) * (t37 + t68) - m(7) * (t39 + t68) + t70 * t58 + t65 * t13 + t67 * t14) * g(3) + (-m(4) * (t43 + t62) - m(5) * (t38 + t62) - m(7) * (t1 + t41) + mrSges(3,3) * t57 - m(3) * t44 - m(6) * (t1 + t38) - mrSges(1,2) - t16 * mrSges(3,1) - t33 * mrSges(2,1) - t36 * mrSges(2,2) + t69 * r_base(2) + t67 * t6 + t65 * t5 + t64 * t15) * g(2) + (-mrSges(3,3) * t59 - m(4) * (t47 + t63) - m(5) * (t40 + t63) - m(7) * (t3 + t45) - m(3) * t48 - m(6) * (t3 + t40) - mrSges(1,1) - t18 * mrSges(3,1) + t33 * mrSges(2,2) - t36 * mrSges(2,1) + t69 * r_base(1) + t67 * t8 + t65 * t7 + t64 * t17) * g(1);
U  = t2;
