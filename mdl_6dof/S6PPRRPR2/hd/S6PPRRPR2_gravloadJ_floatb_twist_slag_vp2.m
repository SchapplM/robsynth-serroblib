% Calculate Gravitation load on the joints for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:07
% EndTime: 2019-03-08 18:49:09
% DurationCPUTime: 0.67s
% Computational Cost: add. (702->89), mult. (1943->135), div. (0->0), fcn. (2447->14), ass. (0->63)
t96 = m(6) + m(7);
t101 = mrSges(5,1) - mrSges(6,2) + mrSges(7,3);
t36 = sin(qJ(6));
t39 = cos(qJ(6));
t100 = t36 * mrSges(7,1) + t39 * mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t97 = m(7) * pkin(10) + pkin(4) * t96 + t101;
t68 = sin(pkin(11));
t71 = cos(pkin(12));
t54 = t68 * t71;
t67 = sin(pkin(12));
t72 = cos(pkin(11));
t58 = t72 * t67;
t74 = cos(pkin(6));
t30 = t74 * t58 + t54;
t38 = sin(qJ(3));
t82 = cos(qJ(3));
t53 = t68 * t67;
t60 = t72 * t71;
t46 = -t74 * t60 + t53;
t69 = sin(pkin(7));
t70 = sin(pkin(6));
t57 = t70 * t69;
t73 = cos(pkin(7));
t92 = t46 * t73 + t72 * t57;
t14 = t30 * t38 + t92 * t82;
t37 = sin(qJ(4));
t75 = qJ(5) * t37;
t40 = cos(qJ(4));
t81 = t14 * t40;
t95 = -pkin(4) * t81 - t14 * t75;
t31 = -t74 * t53 + t60;
t47 = t74 * t54 + t58;
t56 = t70 * t68;
t91 = t47 * t73 - t69 * t56;
t16 = t31 * t38 + t91 * t82;
t80 = t16 * t40;
t94 = -pkin(4) * t80 - t16 * t75;
t55 = t70 * t67;
t90 = t71 * t73 * t70 + t74 * t69;
t25 = t38 * t55 - t90 * t82;
t79 = t25 * t40;
t93 = -pkin(4) * t79 - t25 * t75;
t88 = m(3) + m(4) + m(5) + t96;
t86 = -t96 * qJ(5) - t100;
t85 = t100 * t37 + t101 * t40 + mrSges(4,1);
t84 = -m(7) * (pkin(5) + pkin(9)) - t39 * mrSges(7,1) + t36 * mrSges(7,2) - mrSges(6,1) + mrSges(4,2) - mrSges(5,3);
t11 = t14 * pkin(3);
t15 = t30 * t82 - t92 * t38;
t66 = t15 * pkin(9) - t11;
t12 = t16 * pkin(3);
t17 = t31 * t82 - t91 * t38;
t65 = t17 * pkin(9) - t12;
t24 = t25 * pkin(3);
t26 = t90 * t38 + t82 * t55;
t64 = t26 * pkin(9) - t24;
t59 = t72 * t70;
t45 = -t71 * t57 + t74 * t73;
t42 = t47 * t69 + t73 * t56;
t41 = t46 * t69 - t73 * t59;
t18 = t26 * t37 - t45 * t40;
t5 = t17 * t37 - t42 * t40;
t3 = t15 * t37 - t41 * t40;
t1 = [(-m(2) - t88) * g(3) (-g(1) * t56 + g(2) * t59 - g(3) * t74) * t88 (-m(5) * t64 - m(6) * (t64 + t93) - m(7) * (-pkin(10) * t79 - t24 + t93) + t84 * t26 + t85 * t25) * g(3) + (-m(5) * t66 - m(6) * (t66 + t95) - m(7) * (-pkin(10) * t81 - t11 + t95) + t84 * t15 + t85 * t14) * g(2) + (-m(5) * t65 - m(6) * (t65 + t94) - m(7) * (-pkin(10) * t80 - t12 + t94) + t84 * t17 + t85 * t16) * g(1) (t86 * (t26 * t40 + t45 * t37) + t97 * t18) * g(3) + (t86 * (t15 * t40 + t41 * t37) + t97 * t3) * g(2) + (t86 * (t17 * t40 + t42 * t37) + t97 * t5) * g(1), t96 * (-g(1) * t5 - g(2) * t3 - g(3) * t18) -g(1) * ((-t16 * t36 + t5 * t39) * mrSges(7,1) + (-t16 * t39 - t5 * t36) * mrSges(7,2)) - g(2) * ((-t14 * t36 + t3 * t39) * mrSges(7,1) + (-t14 * t39 - t3 * t36) * mrSges(7,2)) - g(3) * ((t18 * t39 - t25 * t36) * mrSges(7,1) + (-t18 * t36 - t25 * t39) * mrSges(7,2))];
taug  = t1(:);
