% Calculate potential energy for
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:20
% EndTime: 2019-03-08 19:34:20
% DurationCPUTime: 0.72s
% Computational Cost: add. (393->98), mult. (832->107), div. (0->0), fcn. (1011->12), ass. (0->56)
t45 = cos(qJ(2));
t80 = t45 * mrSges(3,2);
t37 = sin(pkin(6));
t39 = cos(pkin(6));
t42 = sin(qJ(2));
t62 = t39 * t42;
t78 = -mrSges(3,3) - mrSges(4,3);
t79 = -t62 * mrSges(3,1) + (m(3) * pkin(7) - t78) * t37 - t39 * t80 - mrSges(2,2);
t35 = sin(pkin(11));
t59 = cos(pkin(11));
t24 = -t42 * t35 + t45 * t59;
t77 = -m(3) - m(2) - m(1);
t75 = -m(7) * pkin(9) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t40 = sin(qJ(6));
t43 = cos(qJ(6));
t73 = -t40 * mrSges(7,1) - t43 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t72 = -m(3) * pkin(1) - t45 * mrSges(3,1) + t42 * mrSges(3,2) - mrSges(2,1);
t71 = m(7) * (pkin(5) + pkin(8)) + t43 * mrSges(7,1) - t40 * mrSges(7,2) + mrSges(6,1) - mrSges(4,2) + mrSges(5,3);
t23 = -t45 * t35 - t42 * t59;
t36 = sin(pkin(10));
t38 = cos(pkin(10));
t48 = t39 * t24;
t9 = t36 * t23 + t38 * t48;
t70 = pkin(8) * t9;
t11 = t23 * t38 - t36 * t48;
t68 = pkin(8) * t11;
t19 = t24 * t37;
t67 = pkin(8) * t19;
t41 = sin(qJ(4));
t64 = t37 * t41;
t44 = cos(qJ(4));
t63 = t37 * t44;
t58 = qJ(1) + r_base(3);
t22 = pkin(2) * t62 + (-pkin(7) - qJ(3)) * t37;
t32 = pkin(2) * t45 + pkin(1);
t56 = t38 * t22 + t36 * t32 + r_base(2);
t55 = t39 * pkin(7) + t58;
t21 = t23 * t39;
t10 = -t21 * t38 + t24 * t36;
t54 = t10 * pkin(3) + t56;
t53 = -t22 * t36 + t38 * t32 + r_base(1);
t52 = t37 * t42 * pkin(2) + t39 * qJ(3) + t55;
t12 = t21 * t36 + t24 * t38;
t51 = t12 * pkin(3) + t53;
t20 = t23 * t37;
t50 = -t20 * pkin(3) + t52;
t3 = t10 * t41 + t38 * t63;
t4 = t10 * t44 - t38 * t64;
t49 = t4 * pkin(4) + qJ(5) * t3 + t54;
t5 = t12 * t41 - t36 * t63;
t6 = t12 * t44 + t36 * t64;
t47 = t6 * pkin(4) + qJ(5) * t5 + t51;
t14 = -t20 * t41 - t39 * t44;
t15 = -t20 * t44 + t39 * t41;
t46 = t15 * pkin(4) + qJ(5) * t14 + t50;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t58 - mrSges(2,3) - m(3) * t55 - (t42 * mrSges(3,1) + t80) * t37 - m(4) * t52 + t20 * mrSges(4,1) - m(5) * (t50 - t67) - m(6) * (t46 - t67) - m(7) * t46 + t78 * t39 + t73 * t14 + t71 * t19 + t75 * t15) * g(3) + (-m(7) * t49 - m(5) * (t54 - t70) - m(6) * (t49 - t70) - m(4) * t56 - mrSges(1,2) - t10 * mrSges(4,1) + t72 * t36 + t77 * r_base(2) + t73 * t3 + t71 * t9 + t75 * t4 + t79 * t38) * g(2) + (-m(5) * (t51 - t68) - m(6) * (t47 - t68) - m(7) * t47 - m(4) * t53 - mrSges(1,1) - t12 * mrSges(4,1) + t72 * t38 + t77 * r_base(1) + t75 * t6 + t73 * t5 + t71 * t11 - t79 * t36) * g(1);
U  = t1;
