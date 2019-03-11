% Calculate potential energy for
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:13:46
% EndTime: 2019-03-08 20:13:47
% DurationCPUTime: 0.69s
% Computational Cost: add. (281->108), mult. (550->122), div. (0->0), fcn. (622->10), ass. (0->51)
t74 = -m(1) - m(2);
t73 = -mrSges(3,1) + mrSges(4,2);
t72 = -mrSges(6,1) - mrSges(7,1);
t71 = -mrSges(6,2) - mrSges(7,2);
t70 = mrSges(3,3) + mrSges(4,1);
t69 = -mrSges(4,3) + mrSges(3,2);
t36 = sin(qJ(5));
t68 = -m(7) * (pkin(5) * t36 + pkin(8)) - mrSges(5,3) + t73;
t67 = m(6) * pkin(9) - m(7) * (-qJ(6) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t31 = sin(pkin(10));
t32 = sin(pkin(6));
t66 = t31 * t32;
t33 = cos(pkin(10));
t65 = t32 * t33;
t37 = sin(qJ(4));
t64 = t32 * t37;
t38 = sin(qJ(2));
t63 = t32 * t38;
t40 = cos(qJ(4));
t62 = t32 * t40;
t41 = cos(qJ(2));
t61 = t32 * t41;
t34 = cos(pkin(6));
t60 = t34 * t38;
t59 = t34 * t41;
t58 = t36 * t38;
t57 = qJ(3) * t41;
t56 = pkin(7) * t65;
t55 = t31 * pkin(1) + r_base(2);
t54 = qJ(1) + r_base(3);
t52 = t33 * pkin(1) + pkin(7) * t66 + r_base(1);
t51 = t34 * pkin(7) + t54;
t50 = pkin(2) * t63 + t51;
t14 = t31 * t38 - t33 * t59;
t15 = t31 * t41 + t33 * t60;
t49 = t15 * pkin(2) + qJ(3) * t14 + t55;
t48 = t34 * pkin(3) + pkin(8) * t63 + t50;
t16 = t31 * t59 + t33 * t38;
t17 = -t31 * t60 + t33 * t41;
t47 = t17 * pkin(2) + qJ(3) * t16 + t52;
t46 = pkin(3) * t66 + t47;
t45 = -t32 * t57 + t48;
t44 = pkin(8) * t17 + t46;
t43 = (-pkin(3) - pkin(7)) * t65 + t49;
t42 = pkin(8) * t15 + t43;
t39 = cos(qJ(5));
t26 = pkin(5) * t39 + pkin(4);
t19 = t34 * t40 - t37 * t61;
t8 = t14 * t37 - t33 * t62;
t6 = t16 * t37 + t31 * t62;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t54 - mrSges(2,3) - m(3) * t51 - m(4) * t50 - m(5) * t45 - t19 * mrSges(5,1) - mrSges(5,3) * t63 - m(6) * (t19 * pkin(4) + t45) - m(7) * (t19 * t26 + t48) + t71 * (-t19 * t36 + t39 * t63) - t70 * t34 + (-m(7) * (pkin(5) * t58 - t57) + (m(4) * qJ(3) - t69) * t41 + t73 * t38) * t32 + t72 * (t19 * t39 + t32 * t58) - t67 * (t34 * t37 + t40 * t61)) * g(3) + (-m(3) * (t55 - t56) - m(4) * (t49 - t56) - m(7) * (t26 * t8 + t43) - m(6) * (pkin(4) * t8 + t42) - m(5) * t42 - mrSges(1,2) - t8 * mrSges(5,1) - t31 * mrSges(2,1) - t33 * mrSges(2,2) + t74 * r_base(2) + t70 * t65 + t72 * (t15 * t36 + t39 * t8) + t71 * (t15 * t39 - t36 * t8) + t69 * t14 + t67 * (t14 * t40 + t33 * t64) + t68 * t15) * g(2) + (-m(7) * (t26 * t6 + t46) - m(6) * (pkin(4) * t6 + t44) - m(3) * t52 - m(4) * t47 - m(5) * t44 - mrSges(1,1) - t6 * mrSges(5,1) + t31 * mrSges(2,2) - t33 * mrSges(2,1) + t74 * r_base(1) - t70 * t66 + t72 * (t17 * t36 + t39 * t6) + t69 * t16 + t71 * (t17 * t39 - t36 * t6) - t67 * (-t16 * t40 + t31 * t64) + t68 * t17) * g(1);
U  = t1;
