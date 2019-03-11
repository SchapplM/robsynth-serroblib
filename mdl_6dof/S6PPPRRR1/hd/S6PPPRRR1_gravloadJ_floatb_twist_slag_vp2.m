% Calculate Gravitation load on the joints for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:10
% EndTime: 2019-03-08 18:39:12
% DurationCPUTime: 0.59s
% Computational Cost: add. (1154->78), mult. (3283->136), div. (0->0), fcn. (4281->18), ass. (0->68)
t91 = m(6) + m(7);
t28 = sin(qJ(6));
t31 = cos(qJ(6));
t93 = m(7) * pkin(5) + t31 * mrSges(7,1) - t28 * mrSges(7,2) + mrSges(6,1);
t89 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t29 = sin(qJ(5));
t32 = cos(qJ(5));
t92 = pkin(4) * t91 - t89 * t29 + t93 * t32 + mrSges(5,1);
t69 = sin(pkin(13));
t70 = sin(pkin(12));
t57 = t70 * t69;
t75 = cos(pkin(13));
t76 = cos(pkin(12));
t65 = t76 * t75;
t79 = cos(pkin(6));
t48 = t79 * t65 - t57;
t78 = cos(pkin(7));
t44 = t48 * t78;
t58 = t70 * t75;
t63 = t76 * t69;
t49 = t79 * t63 + t58;
t72 = sin(pkin(7));
t73 = sin(pkin(6));
t62 = t73 * t72;
t74 = cos(pkin(14));
t53 = t74 * t62;
t68 = sin(pkin(14));
t35 = -t74 * t44 + t49 * t68 + t76 * t53;
t64 = t76 * t73;
t41 = t48 * t72 + t78 * t64;
t71 = sin(pkin(8));
t77 = cos(pkin(8));
t88 = t35 * t77 + t41 * t71;
t50 = -t79 * t58 - t63;
t45 = t50 * t78;
t51 = -t79 * t57 + t65;
t36 = -t74 * t45 + t51 * t68 - t70 * t53;
t61 = t73 * t70;
t42 = t50 * t72 - t78 * t61;
t87 = t36 * t77 + t42 * t71;
t54 = t78 * t75 * t73;
t60 = t73 * t69;
t40 = t68 * t60 + (-t72 * t79 - t54) * t74;
t47 = t75 * t62 - t79 * t78;
t86 = t40 * t77 + t47 * t71;
t83 = m(4) + m(5) + t91;
t82 = m(3) + t83;
t81 = -t28 * mrSges(7,1) - t31 * mrSges(7,2) - t91 * pkin(10) + mrSges(5,2) - mrSges(6,3);
t80 = cos(qJ(4));
t59 = t72 * t68;
t52 = t73 * t59;
t30 = sin(qJ(4));
t26 = t68 * t54 + t79 * t59 + t74 * t60;
t22 = t40 * t71 - t47 * t77;
t21 = t68 * t45 + t51 * t74 + t70 * t52;
t20 = t68 * t44 + t49 * t74 - t76 * t52;
t17 = t36 * t71 - t42 * t77;
t16 = t35 * t71 - t41 * t77;
t15 = t26 * t80 - t86 * t30;
t14 = t26 * t30 + t86 * t80;
t12 = t21 * t80 - t87 * t30;
t11 = t21 * t30 + t87 * t80;
t10 = t20 * t80 - t88 * t30;
t9 = t20 * t30 + t88 * t80;
t8 = t15 * t32 + t22 * t29;
t4 = t12 * t32 + t17 * t29;
t2 = t10 * t32 + t16 * t29;
t1 = [(-m(2) - t82) * g(3) (-t61 * g(1) + t64 * g(2) - t79 * g(3)) * t82 (t42 * g(1) + t41 * g(2) + t47 * g(3)) * t83 (t92 * t14 + t81 * t15) * g(3) + (t81 * t10 + t92 * t9) * g(2) + (t92 * t11 + t81 * t12) * g(1) (t89 * t8 - t93 * (-t15 * t29 + t22 * t32)) * g(3) + (t89 * t2 - t93 * (-t10 * t29 + t16 * t32)) * g(2) + (t89 * t4 - t93 * (-t12 * t29 + t17 * t32)) * g(1), -g(1) * ((t11 * t31 - t4 * t28) * mrSges(7,1) + (-t11 * t28 - t4 * t31) * mrSges(7,2)) - g(2) * ((-t2 * t28 + t9 * t31) * mrSges(7,1) + (-t2 * t31 - t9 * t28) * mrSges(7,2)) - g(3) * ((t14 * t31 - t8 * t28) * mrSges(7,1) + (-t14 * t28 - t8 * t31) * mrSges(7,2))];
taug  = t1(:);
