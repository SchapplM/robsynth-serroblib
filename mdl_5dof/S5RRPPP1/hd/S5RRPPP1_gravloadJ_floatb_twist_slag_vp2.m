% Calculate Gravitation load on the joints for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPPP1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:21
% EndTime: 2019-12-31 19:23:23
% DurationCPUTime: 0.71s
% Computational Cost: add. (262->91), mult. (684->122), div. (0->0), fcn. (733->8), ass. (0->60)
t90 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1) + mrSges(4,3);
t89 = m(5) + m(6);
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t55 = t44 * mrSges(3,1) - t42 * mrSges(3,2);
t93 = -mrSges(2,1) - t55;
t92 = mrSges(2,2) - mrSges(3,3);
t45 = cos(qJ(1));
t36 = t45 * pkin(7);
t43 = sin(qJ(1));
t41 = cos(pkin(5));
t68 = qJ(3) * t41;
t39 = sin(pkin(5));
t69 = qJ(3) * t39;
t29 = t42 * t69;
t71 = t44 * pkin(2) + t29;
t91 = (-pkin(1) - t71) * t43 + t45 * t68 + t36;
t87 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2);
t85 = m(6) * qJ(5) + mrSges(4,1) - mrSges(5,2) + mrSges(6,3);
t82 = t39 * t44;
t86 = mrSges(3,2) * t44 + (m(6) * pkin(2) + mrSges(3,1)) * t42 - t90 * t82;
t83 = t39 * t43;
t81 = t39 * t45;
t38 = sin(pkin(8));
t80 = t42 * t38;
t40 = cos(pkin(8));
t79 = t42 * t40;
t78 = t42 * t43;
t77 = t42 * t45;
t76 = t43 * t41;
t75 = t44 * t38;
t74 = t44 * t40;
t73 = t44 * t45;
t70 = t45 * pkin(1) + t43 * pkin(7);
t67 = pkin(2) * t78;
t66 = pkin(2) * t77;
t65 = t38 * t78;
t62 = t41 * t74;
t61 = t38 * t77;
t60 = t40 * t73;
t59 = t44 * t69;
t51 = t41 * t75 + t79;
t10 = t51 * t43;
t23 = t43 * t59;
t9 = -t43 * t62 + t65;
t58 = -t10 * pkin(3) - t9 * qJ(4) + t23;
t11 = -t41 * t60 + t61;
t12 = t51 * t45;
t24 = t45 * t59;
t57 = -t12 * pkin(3) - t11 * qJ(4) + t24;
t56 = pkin(2) * t73 + t45 * t29 + t43 * t68 + t70;
t15 = t41 * t79 + t75;
t18 = t39 * t77 + t76;
t17 = t39 * t78 - t45 * t41;
t16 = -t41 * t80 + t74;
t6 = t38 * t83 - t41 * t61 + t60;
t5 = t15 * t45 - t40 * t83;
t4 = -t38 * t81 - t41 * t65 + t43 * t74;
t3 = t43 * t75 + (t42 * t76 + t81) * t40;
t1 = [(-m(3) * t70 - m(4) * t56 - t89 * (t6 * pkin(3) + t5 * qJ(4) + t56) + t93 * t45 + t92 * t43 - t85 * t6 + t87 * t5 - t90 * t18) * g(2) + (-m(3) * t36 - t91 * m(4) - t89 * (-t4 * pkin(3) - t3 * qJ(4) + t91) + t92 * t45 + (m(3) * pkin(1) - t93) * t43 + t85 * t4 - t87 * t3 + t90 * t17) * g(1), (-m(4) * t71 - t55 - t89 * (t16 * pkin(3) + t15 * qJ(4) + t71) - t90 * t39 * t42 - t85 * t16 + t87 * t15) * g(3) + (-m(4) * (t23 - t67) - m(5) * (t58 - t67) - m(6) * t58 - t87 * t9 + t85 * t10 + t86 * t43) * g(2) + (-m(4) * (t24 - t66) - m(5) * (t57 - t66) - m(6) * t57 + t85 * t12 - t87 * t11 + t86 * t45) * g(1), (m(4) + t89) * (-g(1) * t18 - g(2) * t17 + g(3) * t82), t89 * (-g(1) * t5 - g(2) * t3 - g(3) * (-t62 + t80)), (-g(1) * t6 - g(2) * t4 - g(3) * t51) * m(6)];
taug = t1(:);
