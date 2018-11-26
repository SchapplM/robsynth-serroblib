% Calculate Gravitation load on the joints for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR6_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From joint_gravload_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:17:50
% EndTime: 2018-11-23 15:17:51
% DurationCPUTime: 1.14s
% Computational Cost: add. (3094->130), mult. (3066->187), div. (0->0), fcn. (2947->24), ass. (0->77)
t80 = sin(qJ(6));
t82 = cos(qJ(6));
t141 = m(7) * pkin(5) + t82 * mrSges(7,1) - t80 * mrSges(7,2) + mrSges(6,1);
t112 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t137 = m(6) + m(7);
t140 = m(5) + t137;
t119 = m(4) + t140;
t139 = pkin(2) * t119 + mrSges(3,1);
t75 = sin(pkin(13));
t77 = cos(pkin(13));
t99 = -m(5) * pkin(3) - t77 * mrSges(5,1) + t75 * mrSges(5,2) - mrSges(4,1);
t74 = pkin(13) + qJ(5);
t72 = sin(t74);
t73 = cos(t74);
t138 = -t112 * t72 + t141 * t73 - t99;
t93 = -m(5) * qJ(4) - t80 * mrSges(7,1) - t82 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t136 = sin(qJ(2));
t76 = sin(pkin(7));
t135 = t72 * t76;
t134 = t73 * t76;
t127 = cos(pkin(6));
t126 = cos(pkin(12));
t125 = sin(pkin(6));
t124 = sin(pkin(12));
t123 = pkin(6) - qJ(2);
t122 = pkin(6) + qJ(2);
t121 = pkin(7) - qJ(3);
t120 = pkin(7) + qJ(3);
t116 = cos(t123);
t115 = cos(t121);
t114 = sin(t123);
t113 = sin(t121);
t110 = cos(t122) / 0.2e1;
t109 = cos(t120) / 0.2e1;
t108 = sin(t122) / 0.2e1;
t107 = sin(t120) / 0.2e1;
t104 = t126 * t125;
t103 = t125 * t124;
t64 = t108 - t114 / 0.2e1;
t84 = cos(qJ(2));
t56 = -t124 * t64 + t126 * t84;
t54 = t124 * t84 + t126 * t64;
t97 = t116 / 0.2e1 + t110;
t96 = t115 / 0.2e1 + t109;
t95 = t108 + t114 / 0.2e1;
t94 = t107 + t113 / 0.2e1;
t92 = t94 * t125;
t90 = -mrSges(3,2) + (t77 * mrSges(5,2) + mrSges(4,3) + (pkin(4) * t137 + mrSges(5,1)) * t75 + t119 * pkin(9)) * t76;
t63 = t107 - t113 / 0.2e1;
t65 = t109 - t115 / 0.2e1;
t66 = t110 - t116 / 0.2e1;
t83 = cos(qJ(3));
t89 = -t127 * t65 + t63 * t95 - t66 * t83;
t88 = t124 * t97 + t126 * t136;
t87 = t124 * t136 - t126 * t97;
t86 = t104 * t65 + t54 * t83 - t63 * t87;
t85 = -t103 * t65 + t56 * t83 - t63 * t88;
t81 = sin(qJ(3));
t79 = -pkin(10) - qJ(4);
t78 = cos(pkin(7));
t71 = pkin(4) * t77 + pkin(3);
t53 = t127 * t78 - t76 * t95;
t41 = t103 * t78 + t76 * t88;
t40 = -t104 * t78 + t76 * t87;
t39 = t66 * t63 + t83 * t95;
t38 = -t66 * t96 + t81 * t95;
t33 = -t127 * t94 - t66 * t81 - t95 * t96;
t31 = -t56 * t63 - t83 * t88;
t30 = t56 * t96 - t81 * t88;
t29 = -t54 * t63 - t83 * t87;
t28 = t54 * t96 - t81 * t87;
t18 = -t124 * t92 + t56 * t81 + t88 * t96;
t15 = t126 * t92 + t54 * t81 + t87 * t96;
t10 = t53 * t72 + t73 * t89;
t4 = t41 * t72 + t73 * t85;
t2 = t40 * t72 + t73 * t86;
t1 = [(-m(2) - m(3) - t119) * g(3) (t112 * (t134 * t66 + t39 * t72) + t99 * t39 - t141 * (-t135 * t66 + t39 * t73) + t93 * t38 + t90 * t66 - t139 * t95 - t137 * (-t38 * t79 + t39 * t71)) * g(3) + (t112 * (-t134 * t54 + t29 * t72) - t141 * (t135 * t54 + t29 * t73) + t99 * t29 + t93 * t28 - t90 * t54 + t139 * t87 - t137 * (-t28 * t79 + t29 * t71)) * g(2) + (t112 * (-t134 * t56 + t31 * t72) - t141 * (t135 * t56 + t31 * t73) + t99 * t31 + t93 * t30 - t90 * t56 + t139 * t88 - t137 * (-t30 * t79 + t31 * t71)) * g(1) (-t137 * (-t33 * t71 - t79 * t89) + t93 * t89 + t138 * t33) * g(3) + (-t137 * (-t15 * t71 - t79 * t86) + t93 * t86 + t138 * t15) * g(2) + (-t137 * (-t18 * t71 - t79 * t85) + t93 * t85 + t138 * t18) * g(1), t140 * (-g(1) * t18 - g(2) * t15 - g(3) * t33) (-t141 * (t53 * t73 - t72 * t89) + t112 * t10) * g(3) + (t112 * t2 - t141 * (t40 * t73 - t72 * t86)) * g(2) + (t112 * t4 - t141 * (t41 * t73 - t72 * t85)) * g(1), -g(1) * ((t18 * t82 - t4 * t80) * mrSges(7,1) + (-t18 * t80 - t4 * t82) * mrSges(7,2)) - g(2) * ((t15 * t82 - t2 * t80) * mrSges(7,1) + (-t15 * t80 - t2 * t82) * mrSges(7,2)) - g(3) * ((-t10 * t80 + t33 * t82) * mrSges(7,1) + (-t10 * t82 - t33 * t80) * mrSges(7,2))];
taug  = t1(:);
