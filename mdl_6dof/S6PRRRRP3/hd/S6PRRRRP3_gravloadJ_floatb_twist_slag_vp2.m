% Calculate Gravitation load on the joints for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:23
% EndTime: 2019-03-09 00:06:26
% DurationCPUTime: 1.05s
% Computational Cost: add. (643->111), mult. (1297->159), div. (0->0), fcn. (1517->12), ass. (0->58)
t53 = qJ(4) + qJ(5);
t50 = cos(t53);
t58 = cos(qJ(4));
t51 = t58 * pkin(4);
t39 = pkin(5) * t50 + t51;
t55 = sin(qJ(4));
t122 = -m(6) * (t51 + pkin(3)) - m(7) * (pkin(3) + t39) - mrSges(4,1) - m(5) * pkin(3) - t58 * mrSges(5,1) + t55 * mrSges(5,2);
t61 = -pkin(10) - pkin(9);
t108 = m(6) * t61 + m(7) * (-qJ(6) + t61) + mrSges(4,2) - mrSges(6,3) - mrSges(7,3) - m(5) * pkin(9) - mrSges(5,3);
t121 = mrSges(6,1) + mrSges(7,1);
t116 = -mrSges(6,2) - mrSges(7,2);
t49 = sin(t53);
t98 = pkin(4) * t55;
t38 = pkin(5) * t49 + t98;
t120 = -m(6) * t98 - m(7) * t38 - mrSges(5,1) * t55 - mrSges(5,2) * t58 + mrSges(3,2) - mrSges(4,3);
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t119 = t108 * t56 + t122 * t59 - mrSges(3,1);
t87 = cos(pkin(6));
t54 = sin(pkin(6));
t57 = sin(qJ(2));
t95 = t54 * t57;
t36 = t87 * t56 + t59 * t95;
t60 = cos(qJ(2));
t94 = t54 * t60;
t21 = -t36 * t49 - t50 * t94;
t115 = t116 * (-t36 * t50 + t49 * t94) - t121 * t21;
t114 = -m(6) * pkin(4) - mrSges(5,1);
t85 = sin(pkin(11));
t70 = t87 * t85;
t86 = cos(pkin(11));
t34 = -t57 * t70 + t86 * t60;
t78 = t54 * t85;
t26 = t34 * t59 + t56 * t78;
t33 = t86 * t57 + t60 * t70;
t11 = -t26 * t49 + t33 * t50;
t113 = t116 * (-t26 * t50 - t33 * t49) - t121 * t11;
t71 = t87 * t86;
t32 = t57 * t71 + t85 * t60;
t79 = t54 * t86;
t24 = t32 * t59 - t56 * t79;
t31 = t85 * t57 - t60 * t71;
t9 = -t24 * t49 + t31 * t50;
t112 = -t121 * t9 + t116 * (-t24 * t50 - t31 * t49);
t111 = -m(4) - m(6) - m(7);
t110 = -m(5) + t111;
t109 = t116 * t49 + t121 * t50 - t122;
t105 = -m(5) * pkin(8) + t120;
t103 = m(7) * pkin(5);
t97 = t49 * t59;
t96 = t50 * t59;
t91 = t59 * t60;
t35 = t56 * t95 - t87 * t59;
t30 = t33 * pkin(2);
t29 = t31 * pkin(2);
t25 = t34 * t56 - t59 * t78;
t23 = t32 * t56 + t59 * t79;
t1 = [(-m(2) - m(3) + t110) * g(3) (m(5) * t29 - t121 * (-t31 * t96 + t32 * t49) + t116 * (t31 * t97 + t32 * t50) + t111 * (t32 * pkin(8) - t29) + t105 * t32 - t119 * t31) * g(2) + (m(5) * t30 - t121 * (-t33 * t96 + t34 * t49) + t116 * (t33 * t97 + t34 * t50) + t111 * (t34 * pkin(8) - t30) + t105 * t34 - t119 * t33) * g(1) + (t110 * (pkin(2) * t94 + pkin(8) * t95) + (-t121 * (t49 * t57 + t50 * t91) + t116 * (-t49 * t91 + t50 * t57) + t119 * t60 + t120 * t57) * t54) * g(3) (t108 * t36 + t109 * t35) * g(3) + (t108 * t24 + t109 * t23) * g(2) + (t108 * t26 + t109 * t25) * g(1) (-(-t36 * t58 + t55 * t94) * mrSges(5,2) - m(7) * (-t36 * t38 - t39 * t94) + t114 * (-t36 * t55 - t58 * t94) + t115) * g(3) + (-(-t24 * t58 - t31 * t55) * mrSges(5,2) - m(7) * (-t24 * t38 + t31 * t39) + t114 * (-t24 * t55 + t31 * t58) + t112) * g(2) + (-(-t26 * t58 - t33 * t55) * mrSges(5,2) - m(7) * (-t26 * t38 + t33 * t39) + t114 * (-t26 * t55 + t33 * t58) + t113) * g(1) (-t21 * t103 + t115) * g(3) + (-t9 * t103 + t112) * g(2) + (-t11 * t103 + t113) * g(1) (-g(1) * t25 - g(2) * t23 - g(3) * t35) * m(7)];
taug  = t1(:);
