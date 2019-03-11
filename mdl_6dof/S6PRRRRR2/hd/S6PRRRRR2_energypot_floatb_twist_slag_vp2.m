% Calculate potential energy for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:42:19
% EndTime: 2019-03-09 00:42:20
% DurationCPUTime: 0.96s
% Computational Cost: add. (376->115), mult. (559->130), div. (0->0), fcn. (633->14), ass. (0->46)
t32 = qJ(5) + qJ(6);
t25 = sin(t32);
t27 = cos(t32);
t38 = sin(qJ(5));
t41 = cos(qJ(5));
t45 = -pkin(9) - pkin(8);
t71 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t72 = -m(5) - m(6) - m(7);
t74 = -t25 * mrSges(7,1) - t41 * mrSges(6,2) - t27 * mrSges(7,2) - t72 * t45 - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t38 - t71;
t73 = -m(1) - m(2);
t69 = -m(6) * pkin(10) + m(7) * (-pkin(11) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t23 = pkin(5) * t41 + pkin(4);
t67 = -m(6) * pkin(4) - m(7) * t23 - t41 * mrSges(6,1) - t27 * mrSges(7,1) + t38 * mrSges(6,2) + t25 * mrSges(7,2) - mrSges(5,1);
t39 = sin(qJ(3));
t66 = pkin(3) * t39;
t34 = sin(pkin(12));
t35 = sin(pkin(6));
t63 = t34 * t35;
t36 = cos(pkin(12));
t62 = t35 * t36;
t61 = t35 * t39;
t40 = sin(qJ(2));
t60 = t35 * t40;
t42 = cos(qJ(3));
t59 = t35 * t42;
t43 = cos(qJ(2));
t58 = t35 * t43;
t37 = cos(pkin(6));
t57 = t37 * t40;
t56 = t37 * t43;
t55 = t34 * pkin(1) + r_base(2);
t54 = t34 * t61;
t53 = t38 * t58;
t52 = qJ(1) + r_base(3);
t51 = t36 * pkin(1) + pkin(7) * t63 + r_base(1);
t50 = t37 * pkin(7) + t52;
t49 = -pkin(7) * t62 + t55;
t24 = pkin(3) * t42 + pkin(2);
t47 = t24 * t60 + t37 * t66 + t45 * t58 + t50;
t33 = qJ(3) + qJ(4);
t28 = cos(t33);
t26 = sin(t33);
t14 = -t34 * t57 + t36 * t43;
t12 = t34 * t43 + t36 * t57;
t8 = t26 * t37 + t28 * t60;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t52 - mrSges(2,3) - m(5) * t47 - t8 * mrSges(5,1) + mrSges(5,3) * t58 - m(6) * (pkin(4) * t8 + t47) - (t41 * t8 - t53) * mrSges(6,1) - (-t38 * t8 - t41 * t58) * mrSges(6,2) - m(7) * (-pkin(5) * t53 + t8 * t23 + t47) - (-t25 * t58 + t27 * t8) * mrSges(7,1) - (-t25 * t8 - t27 * t58) * mrSges(7,2) + (-m(3) - m(4)) * t50 + (-t39 * mrSges(4,1) - t42 * mrSges(4,2) - mrSges(3,3)) * t37 + (t71 * t43 + (-m(4) * pkin(2) - t42 * mrSges(4,1) + t39 * mrSges(4,2) - mrSges(3,1)) * t40) * t35 + t69 * (t26 * t60 - t37 * t28)) * g(3) + (-mrSges(1,2) - m(4) * (pkin(2) * t12 + t49) - m(3) * t49 - (-t12 * t39 - t36 * t59) * mrSges(4,2) - (t12 * t42 - t36 * t61) * mrSges(4,1) + mrSges(3,3) * t62 - t12 * mrSges(3,1) - t34 * mrSges(2,1) - t36 * mrSges(2,2) + t73 * r_base(2) + t72 * (t12 * t24 + (-pkin(7) - t66) * t62 + t55) + t67 * (t12 * t28 - t26 * t62) + t69 * (t12 * t26 + t28 * t62) + t74 * (t34 * t40 - t36 * t56)) * g(2) + (-mrSges(3,3) * t63 - mrSges(1,1) - m(4) * (pkin(2) * t14 + t51) - m(3) * t51 - (t14 * t42 + t54) * mrSges(4,1) - (-t14 * t39 + t34 * t59) * mrSges(4,2) - t14 * mrSges(3,1) + t34 * mrSges(2,2) - t36 * mrSges(2,1) + t73 * r_base(1) + t72 * (pkin(3) * t54 + t14 * t24 + t51) + t69 * (t14 * t26 - t28 * t63) + t67 * (t14 * t28 + t26 * t63) + t74 * (t34 * t56 + t36 * t40)) * g(1);
U  = t1;
