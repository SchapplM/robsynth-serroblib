% Calculate Gravitation load on the joints for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:24
% EndTime: 2019-03-08 18:45:26
% DurationCPUTime: 0.66s
% Computational Cost: add. (757->74), mult. (1998->115), div. (0->0), fcn. (2518->16), ass. (0->55)
t85 = m(6) + m(7);
t81 = m(5) + t85;
t30 = pkin(13) + qJ(6);
t28 = sin(t30);
t29 = cos(t30);
t31 = sin(pkin(13));
t32 = cos(pkin(13));
t80 = mrSges(5,1) + m(7) * (pkin(5) * t32 + pkin(4)) + t29 * mrSges(7,1) - t28 * mrSges(7,2) + m(6) * pkin(4) + t32 * mrSges(6,1) - t31 * mrSges(6,2);
t79 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5)) - mrSges(7,3);
t34 = sin(qJ(4));
t36 = cos(qJ(4));
t86 = pkin(3) * t81 - t79 * t34 + t80 * t36 + mrSges(4,1);
t66 = sin(pkin(12));
t67 = sin(pkin(11));
t54 = t67 * t66;
t70 = cos(pkin(12));
t71 = cos(pkin(11));
t61 = t71 * t70;
t73 = cos(pkin(6));
t44 = -t61 * t73 + t54;
t68 = sin(pkin(7));
t69 = sin(pkin(6));
t58 = t69 * t68;
t72 = cos(pkin(7));
t84 = t44 * t72 + t71 * t58;
t55 = t67 * t70;
t59 = t71 * t66;
t45 = t55 * t73 + t59;
t57 = t69 * t67;
t83 = t45 * t72 - t68 * t57;
t82 = t70 * t72 * t69 + t73 * t68;
t76 = m(3) + m(4) + t81;
t75 = -t28 * mrSges(7,1) - t32 * mrSges(6,2) - t29 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t31 - t81 * pkin(9);
t74 = cos(qJ(3));
t60 = t71 * t69;
t56 = t69 * t66;
t43 = -t58 * t70 + t72 * t73;
t38 = t45 * t68 + t57 * t72;
t37 = t44 * t68 - t60 * t72;
t35 = sin(qJ(3));
t22 = -t54 * t73 + t61;
t21 = t59 * t73 + t55;
t17 = t82 * t35 + t74 * t56;
t16 = t35 * t56 - t82 * t74;
t12 = t17 * t36 + t34 * t43;
t11 = t17 * t34 - t36 * t43;
t10 = t22 * t74 - t83 * t35;
t9 = t22 * t35 + t83 * t74;
t8 = t21 * t74 - t84 * t35;
t7 = t21 * t35 + t84 * t74;
t4 = t10 * t36 + t34 * t38;
t3 = t10 * t34 - t36 * t38;
t2 = t34 * t37 + t8 * t36;
t1 = t34 * t8 - t36 * t37;
t5 = [(-m(2) - t76) * g(3) (-g(1) * t57 + g(2) * t60 - g(3) * t73) * t76 (t86 * t16 + t75 * t17) * g(3) + (t86 * t7 + t75 * t8) * g(2) + (t75 * t10 + t86 * t9) * g(1) (t80 * t11 + t79 * t12) * g(3) + (t80 * t1 + t79 * t2) * g(2) + (t80 * t3 + t79 * t4) * g(1), t85 * (-g(1) * t3 - g(2) * t1 - g(3) * t11) -g(1) * ((-t28 * t4 + t29 * t9) * mrSges(7,1) + (-t28 * t9 - t29 * t4) * mrSges(7,2)) - g(2) * ((-t2 * t28 + t29 * t7) * mrSges(7,1) + (-t2 * t29 - t28 * t7) * mrSges(7,2)) - g(3) * ((-t12 * t28 + t16 * t29) * mrSges(7,1) + (-t12 * t29 - t16 * t28) * mrSges(7,2))];
taug  = t5(:);
