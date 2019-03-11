% Calculate Gravitation load on the joints for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR8_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:41
% EndTime: 2019-03-09 05:22:43
% DurationCPUTime: 0.78s
% Computational Cost: add. (340->101), mult. (464->118), div. (0->0), fcn. (415->10), ass. (0->61)
t34 = -qJ(5) - pkin(8);
t95 = mrSges(4,2) - m(5) * pkin(8) + m(6) * t34 + m(7) * (-pkin(9) + t34) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t94 = t36 * mrSges(4,1) + t95 * t39;
t33 = qJ(4) + pkin(10);
t25 = cos(t33);
t38 = cos(qJ(4));
t29 = t38 * pkin(4);
t19 = pkin(5) * t25 + t29;
t93 = -m(5) * pkin(3) - m(6) * (t29 + pkin(3)) - m(7) * (pkin(3) + t19);
t83 = m(6) + m(7);
t87 = -m(4) - m(5) - t83;
t91 = -m(3) + t87;
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t81 = -g(1) * t37 + g(2) * t40;
t90 = m(6) * pkin(4);
t88 = -mrSges(2,1) + mrSges(3,2) - mrSges(4,3);
t24 = sin(t33);
t35 = sin(qJ(4));
t70 = t35 * pkin(4);
t18 = pkin(5) * t24 + t70;
t86 = m(6) * t70 + m(7) * t18;
t85 = t93 * t36 + mrSges(2,2) - mrSges(3,3) - t94;
t82 = -mrSges(5,1) - t90;
t26 = qJ(6) + t33;
t21 = sin(t26);
t22 = cos(t26);
t80 = t38 * mrSges(5,1) + t25 * mrSges(6,1) + t22 * mrSges(7,1) - t35 * mrSges(5,2) - t24 * mrSges(6,2) - t21 * mrSges(7,2) - t93;
t61 = t40 * t22;
t69 = t37 * t21;
t5 = -t36 * t69 + t61;
t62 = t40 * t21;
t68 = t37 * t22;
t6 = t36 * t68 + t62;
t75 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t7 = t36 * t62 + t68;
t8 = t36 * t61 - t69;
t74 = t7 * mrSges(7,1) + t8 * mrSges(7,2);
t71 = g(3) * t39;
t67 = t37 * t24;
t66 = t37 * t25;
t65 = t37 * t35;
t64 = t37 * t38;
t63 = t40 * t18;
t60 = t40 * t24;
t59 = t40 * t25;
t58 = t40 * t35;
t57 = t40 * t38;
t56 = t40 * pkin(1) + t37 * qJ(2);
t47 = -mrSges(7,1) * t21 - mrSges(7,2) * t22;
t15 = t36 * t58 + t64;
t13 = -t36 * t65 + t57;
t16 = t36 * t57 - t65;
t14 = t36 * t64 + t58;
t12 = t36 * t59 - t67;
t11 = t36 * t60 + t66;
t10 = t36 * t66 + t60;
t9 = -t36 * t67 + t59;
t1 = [(-t58 * t90 - m(3) * t56 - m(7) * t63 - t14 * mrSges(5,1) - t10 * mrSges(6,1) - t6 * mrSges(7,1) - t13 * mrSges(5,2) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + t87 * (t40 * pkin(7) + t56) + t88 * t40 + t85 * t37) * g(2) + (-t16 * mrSges(5,1) - t12 * mrSges(6,1) - t8 * mrSges(7,1) + t15 * mrSges(5,2) + t11 * mrSges(6,2) + t7 * mrSges(7,2) + (m(3) * pkin(1) + t87 * (-pkin(1) - pkin(7)) + t86 - t88) * t37 + (t91 * qJ(2) + t85) * t40) * g(1), -t81 * t91 (t80 * t36 + t94) * g(3) - t81 * ((-mrSges(4,1) - t80) * t39 + t95 * t36) (mrSges(5,1) * t35 + mrSges(6,1) * t24 + mrSges(5,2) * t38 + mrSges(6,2) * t25 - t47 + t86) * t71 + (-t16 * mrSges(5,2) - t11 * mrSges(6,1) - t12 * mrSges(6,2) - m(7) * (t37 * t19 + t36 * t63) - t74 + t82 * t15) * g(2) + (t14 * mrSges(5,2) - t9 * mrSges(6,1) + t10 * mrSges(6,2) - m(7) * (-t37 * t36 * t18 + t40 * t19) - t75 + t82 * t13) * g(1) (-t36 * g(3) - t81 * t39) * t83, -g(1) * t75 - g(2) * t74 - t47 * t71];
taug  = t1(:);
