% Calculate Gravitation load on the joints for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:13
% EndTime: 2019-03-08 21:34:16
% DurationCPUTime: 0.89s
% Computational Cost: add. (587->94), mult. (1211->130), div. (0->0), fcn. (1420->12), ass. (0->58)
t115 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t53 = pkin(11) + qJ(5);
t51 = sin(t53);
t52 = cos(t53);
t76 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t121 = t115 * t51 + t76 * t52;
t54 = sin(pkin(11));
t56 = cos(pkin(11));
t119 = -m(5) * pkin(3) - t56 * mrSges(5,1) + t54 * mrSges(5,2) - mrSges(4,1);
t118 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(5,3);
t109 = -m(6) - m(7);
t50 = pkin(4) * t56 + pkin(3);
t117 = t109 * t50;
t116 = mrSges(6,3) + mrSges(7,2);
t114 = -m(4) + t109;
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t113 = t118 * t58 + t119 * t60 - mrSges(3,1);
t112 = -t54 * mrSges(5,1) - t56 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t100 = pkin(4) * t54;
t111 = -m(5) * pkin(8) + t109 * t100 + t115 * t52 - t76 * t51 + t112;
t57 = -pkin(9) - qJ(4);
t110 = -t113 + (t109 * t57 + t116) * t58 + (-t117 + t121) * t60;
t105 = m(5) - t109;
t103 = -t119 + t121;
t102 = -t116 + t118;
t55 = sin(pkin(6));
t59 = sin(qJ(2));
t95 = t55 * t59;
t61 = cos(qJ(2));
t94 = t55 * t61;
t90 = t60 * t61;
t87 = pkin(2) * t94 + pkin(8) * t95;
t86 = cos(pkin(6));
t85 = cos(pkin(10));
t84 = sin(pkin(10));
t83 = t58 * t94;
t82 = t51 * t94;
t78 = t55 * t85;
t77 = t55 * t84;
t71 = t86 * t85;
t70 = t86 * t84;
t38 = t86 * t58 + t60 * t95;
t37 = t58 * t95 - t86 * t60;
t36 = -t59 * t70 + t85 * t61;
t35 = t85 * t59 + t61 * t70;
t34 = t59 * t71 + t84 * t61;
t33 = t84 * t59 - t61 * t71;
t32 = t35 * pkin(2);
t31 = t33 * pkin(2);
t18 = t36 * t60 + t58 * t77;
t17 = t36 * t58 - t60 * t77;
t16 = t34 * t60 - t58 * t78;
t15 = t34 * t58 + t60 * t78;
t11 = t38 * t51 + t52 * t94;
t3 = t18 * t51 - t35 * t52;
t1 = t16 * t51 - t33 * t52;
t2 = [(-m(2) - m(3) - m(4) - t105) * g(3) ((-m(4) - m(5)) * t87 - t116 * t83 + t109 * (t95 * t100 - t57 * t83 + t87) - t115 * (-t52 * t95 + t60 * t82) + (t90 * t117 - t76 * (t51 * t59 + t52 * t90) + t113 * t61 + t112 * t59) * t55) * g(3) + (m(5) * t31 + t114 * (pkin(8) * t34 - t31) + t111 * t34 + t110 * t33) * g(2) + (m(5) * t32 + t114 * (pkin(8) * t36 - t32) + t111 * t36 + t110 * t35) * g(1) (t109 * (-t37 * t50 - t38 * t57) + t102 * t38 + t103 * t37) * g(3) + (t109 * (-t15 * t50 - t16 * t57) + t102 * t16 + t103 * t15) * g(2) + (t109 * (-t17 * t50 - t18 * t57) + t102 * t18 + t103 * t17) * g(1), t105 * (-g(1) * t17 - g(2) * t15 - g(3) * t37) (-t115 * (t38 * t52 - t82) + t76 * t11) * g(3) + (-t115 * (t16 * t52 + t33 * t51) + t76 * t1) * g(2) + (-t115 * (t18 * t52 + t35 * t51) + t76 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t11) * m(7)];
taug  = t2(:);
