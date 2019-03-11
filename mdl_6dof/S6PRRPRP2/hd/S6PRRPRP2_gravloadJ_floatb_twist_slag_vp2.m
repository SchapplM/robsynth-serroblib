% Calculate Gravitation load on the joints for
% S6PRRPRP2
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:25
% EndTime: 2019-03-08 21:28:26
% DurationCPUTime: 1.01s
% Computational Cost: add. (614->104), mult. (1106->149), div. (0->0), fcn. (1279->12), ass. (0->60)
t55 = sin(qJ(5));
t74 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t127 = t74 * t55;
t76 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t113 = -m(6) - m(7);
t51 = qJ(3) + pkin(11);
t50 = cos(t51);
t126 = t50 * t76;
t125 = mrSges(6,3) + mrSges(7,2);
t124 = -m(4) * pkin(8) - t76 * t55 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t110 = m(5) - t113;
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t89 = cos(pkin(6));
t53 = sin(pkin(6));
t57 = sin(qJ(2));
t97 = t53 * t57;
t120 = -t56 * t97 + t89 * t59;
t49 = sin(t51);
t119 = -m(4) * pkin(2) - t59 * mrSges(4,1) - t50 * mrSges(5,1) + t56 * mrSges(4,2) + t49 * mrSges(5,2) - mrSges(3,1);
t60 = cos(qJ(2));
t52 = sin(pkin(10));
t81 = t52 * t89;
t88 = cos(pkin(10));
t37 = -t57 * t81 + t88 * t60;
t96 = t53 * t59;
t118 = -t37 * t56 + t52 * t96;
t58 = cos(qJ(5));
t117 = -t76 * t58 - mrSges(5,1) + t127;
t116 = mrSges(5,2) - t125;
t115 = -t74 * t58 + t124;
t107 = pkin(4) * t50;
t114 = t113 * (-pkin(9) * t49 - t107) + t58 * t126 - t50 * t127 - t119 + t125 * t49;
t69 = t89 * t88;
t35 = t52 * t60 + t57 * t69;
t80 = t53 * t88;
t65 = -t35 * t56 - t59 * t80;
t62 = t65 * pkin(3);
t98 = t52 * t53;
t95 = t53 * t60;
t92 = t58 * t60;
t85 = t49 * t95;
t84 = t55 * t95;
t75 = t118 * pkin(3);
t70 = t120 * pkin(3);
t54 = -qJ(4) - pkin(8);
t48 = pkin(3) * t59 + pkin(2);
t40 = t48 * t95;
t36 = t88 * t57 + t60 * t81;
t34 = t52 * t57 - t60 * t69;
t27 = t89 * t49 + t50 * t97;
t26 = -t49 * t97 + t89 * t50;
t15 = t27 * t55 + t53 * t92;
t14 = t37 * t50 + t49 * t98;
t13 = -t37 * t49 + t50 * t98;
t12 = t35 * t50 - t49 * t80;
t11 = -t35 * t49 - t50 * t80;
t3 = t14 * t55 - t36 * t58;
t1 = t12 * t55 - t34 * t58;
t2 = [(-m(2) - m(3) - m(4) - t110) * g(3) (-t110 * (-t34 * t48 - t35 * t54) + t115 * t35 + t114 * t34) * g(2) + (-t110 * (-t36 * t48 - t37 * t54) + t115 * t37 + t114 * t36) * g(1) + (-m(5) * t40 - t125 * t85 + t113 * (pkin(9) * t85 + t95 * t107 - t54 * t97 + t40) + t74 * (t50 * t84 - t58 * t97) + (-t92 * t126 + t119 * t60 + (m(5) * t54 + t124) * t57) * t53) * g(3) (-t120 * mrSges(4,1) - (-t89 * t56 - t57 * t96) * mrSges(4,2) - m(5) * t70 + t113 * (t26 * pkin(4) + pkin(9) * t27 + t70) + t116 * t27 + t117 * t26) * g(3) + (-t65 * mrSges(4,1) - (-t35 * t59 + t56 * t80) * mrSges(4,2) - m(5) * t62 + t116 * t12 + t117 * t11 + t113 * (t11 * pkin(4) + t12 * pkin(9) + t62)) * g(2) + (-t118 * mrSges(4,1) - (-t37 * t59 - t56 * t98) * mrSges(4,2) - m(5) * t75 + t113 * (t13 * pkin(4) + pkin(9) * t14 + t75) + t116 * t14 + t117 * t13) * g(1), t110 * (-g(1) * t36 - g(2) * t34 + g(3) * t95) (t74 * (t27 * t58 - t84) + t76 * t15) * g(3) + (t74 * (t12 * t58 + t34 * t55) + t76 * t1) * g(2) + (t74 * (t14 * t58 + t36 * t55) + t76 * t3) * g(1) (-g(1) * t3 - g(2) * t1 - g(3) * t15) * m(7)];
taug  = t2(:);
