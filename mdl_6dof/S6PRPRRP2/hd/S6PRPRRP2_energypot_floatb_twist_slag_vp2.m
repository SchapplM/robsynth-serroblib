% Calculate potential energy for
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:50
% EndTime: 2019-03-08 19:59:50
% DurationCPUTime: 0.71s
% Computational Cost: add. (439->99), mult. (952->118), div. (0->0), fcn. (1177->12), ass. (0->52)
t54 = cos(qJ(2));
t88 = t54 * mrSges(3,2);
t45 = sin(pkin(6));
t48 = cos(pkin(6));
t51 = sin(qJ(2));
t70 = t48 * t51;
t84 = -mrSges(3,3) - mrSges(4,3);
t87 = -t70 * mrSges(3,1) + t45 * (m(3) * pkin(7) - t84) - t48 * t88 - mrSges(2,2);
t86 = -m(6) - m(7);
t85 = mrSges(4,2) - mrSges(5,3);
t43 = sin(pkin(11));
t46 = cos(pkin(11));
t83 = t51 * t43 - t46 * t54;
t82 = -m(3) - m(1) - m(2);
t81 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t79 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t78 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t76 = -m(3) * pkin(1) - t54 * mrSges(3,1) + t51 * mrSges(3,2) - mrSges(2,1);
t50 = sin(qJ(4));
t73 = t45 * t50;
t53 = cos(qJ(4));
t72 = t45 * t53;
t67 = qJ(1) + r_base(3);
t31 = pkin(2) * t70 + (-pkin(7) - qJ(3)) * t45;
t40 = pkin(2) * t54 + pkin(1);
t44 = sin(pkin(10));
t47 = cos(pkin(10));
t66 = t47 * t31 + t44 * t40 + r_base(2);
t65 = t48 * pkin(7) + t67;
t64 = t43 * t54 + t51 * t46;
t63 = -t31 * t44 + t47 * t40 + r_base(1);
t62 = t45 * t51 * pkin(2) + t48 * qJ(3) + t65;
t61 = t83 * t48;
t17 = -t44 * t64 - t47 * t61;
t30 = t64 * t48;
t18 = t30 * t47 - t44 * t83;
t60 = t18 * pkin(3) - pkin(8) * t17 + t66;
t19 = t44 * t61 - t47 * t64;
t20 = -t30 * t44 - t47 * t83;
t59 = t20 * pkin(3) - pkin(8) * t19 + t63;
t28 = t83 * t45;
t29 = t64 * t45;
t58 = t29 * pkin(3) + pkin(8) * t28 + t62;
t52 = cos(qJ(5));
t49 = sin(qJ(5));
t23 = t29 * t53 + t48 * t50;
t22 = t29 * t50 - t48 * t53;
t12 = t20 * t53 + t44 * t73;
t11 = t20 * t50 - t44 * t72;
t10 = t18 * t53 - t47 * t73;
t9 = t18 * t50 + t47 * t72;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t67 - mrSges(2,3) - m(3) * t65 - (t51 * mrSges(3,1) + t88) * t45 - m(4) * t62 - t29 * mrSges(4,1) - m(5) * t58 - t23 * mrSges(5,1) + t86 * (t23 * pkin(4) + pkin(9) * t22 + t58) + t79 * (t23 * t52 + t28 * t49) + t78 * (t23 * t49 - t28 * t52) + t84 * t48 + t85 * t28 + t81 * t22) * g(3) + (-m(4) * t66 - m(5) * t60 - t18 * mrSges(4,1) - t10 * mrSges(5,1) - mrSges(1,2) + t76 * t44 + t82 * r_base(2) + t86 * (t10 * pkin(4) + pkin(9) * t9 + t60) + t79 * (t10 * t52 - t17 * t49) - t85 * t17 + t78 * (t10 * t49 + t17 * t52) + t81 * t9 + t87 * t47) * g(2) + (-m(4) * t63 - m(5) * t59 - t20 * mrSges(4,1) - t12 * mrSges(5,1) - mrSges(1,1) + t76 * t47 + t82 * r_base(1) + t86 * (t12 * pkin(4) + pkin(9) * t11 + t59) + t79 * (t12 * t52 - t19 * t49) + t78 * (t12 * t49 + t19 * t52) - t85 * t19 + t81 * t11 - t87 * t44) * g(1);
U  = t1;
