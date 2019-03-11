% Calculate Gravitation load on the joints for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRP4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:05:10
% EndTime: 2019-03-09 02:05:11
% DurationCPUTime: 0.62s
% Computational Cost: add. (293->71), mult. (593->91), div. (0->0), fcn. (670->8), ass. (0->40)
t69 = m(6) + m(7);
t59 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t58 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t71 = mrSges(7,2) + mrSges(6,3);
t27 = sin(qJ(5));
t29 = cos(qJ(5));
t70 = t69 * pkin(4) + t58 * t27 + t59 * t29;
t68 = mrSges(2,1) + mrSges(3,1);
t67 = mrSges(2,2) - mrSges(3,3);
t66 = mrSges(4,2) - mrSges(5,3);
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t65 = (-mrSges(5,2) + t71) * t30 + (-mrSges(5,1) - t70) * t28;
t64 = t69 * pkin(8);
t63 = m(4) + m(5) + t69;
t42 = t30 * mrSges(5,1) - t28 * mrSges(5,2);
t61 = t71 * t28 + mrSges(4,1) + t42;
t57 = g(3) * t28;
t56 = t28 * pkin(8);
t55 = cos(qJ(1));
t54 = sin(qJ(1));
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t18 = -t54 * t45 - t55 * t46;
t53 = t18 * t30;
t19 = t55 * t45 - t54 * t46;
t52 = t19 * t30;
t51 = t27 * t30;
t48 = t29 * t30;
t47 = t55 * pkin(1) + t54 * qJ(2);
t44 = t55 * pkin(2) + t47;
t43 = -t54 * pkin(1) + t55 * qJ(2);
t2 = t18 * t27 + t19 * t48;
t1 = -t18 * t29 + t19 * t51;
t38 = -t18 * pkin(3) + pkin(7) * t19 + t44;
t35 = -t54 * pkin(2) + t43;
t34 = t19 * pkin(3) + t18 * pkin(7) + t35;
t6 = -t18 * t48 + t19 * t27;
t5 = -t18 * t51 - t19 * t29;
t3 = [(-m(3) * t47 - m(4) * t44 - m(5) * t38 - t69 * (-pkin(4) * t53 - t18 * t56 + t38) - t59 * t6 - t68 * t55 + t67 * t54 - t58 * t5 + t66 * t19 + t61 * t18) * g(2) + (-m(3) * t43 - m(4) * t35 - m(5) * t34 + t67 * t55 + t68 * t54 - t69 * (pkin(4) * t52 + t19 * t56 + t34) - t59 * t2 + t66 * t18 - t58 * t1 - t61 * t19) * g(1) (-t54 * g(1) + t55 * g(2)) * (m(3) + t63) t63 * g(3) (t42 + t70 * t30 + (t64 + t71) * t28) * g(3) + (t65 * t19 + t52 * t64) * g(2) + (t65 * t18 + t53 * t64) * g(1) (-t59 * t27 + t58 * t29) * t57 + (-t1 * t59 + t58 * t2) * g(2) + (t59 * t5 - t58 * t6) * g(1) (-g(1) * t5 + g(2) * t1 + t27 * t57) * m(7)];
taug  = t3(:);
