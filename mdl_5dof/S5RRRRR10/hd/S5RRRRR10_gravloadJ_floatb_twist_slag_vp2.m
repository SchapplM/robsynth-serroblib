% Calculate Gravitation load on the joints for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR10_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:31:59
% EndTime: 2019-12-31 22:32:03
% DurationCPUTime: 0.99s
% Computational Cost: add. (563->122), mult. (993->170), div. (0->0), fcn. (1143->12), ass. (0->62)
t125 = mrSges(5,2) - mrSges(6,3);
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t124 = t63 * mrSges(6,1) - t59 * mrSges(6,2) + mrSges(5,1);
t126 = m(6) * pkin(4) + t124;
t83 = -m(6) * pkin(10) + t125;
t75 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t117 = t59 * mrSges(6,1) + t63 * mrSges(6,2) - t75;
t58 = sin(pkin(5));
t61 = sin(qJ(2));
t107 = t58 * t61;
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t97 = cos(pkin(5));
t123 = -t60 * t107 + t97 * t64;
t105 = t58 * t64;
t109 = cos(qJ(1));
t65 = cos(qJ(2));
t62 = sin(qJ(1));
t85 = t62 * t97;
t39 = t109 * t65 - t61 * t85;
t19 = t62 * t105 - t39 * t60;
t57 = qJ(3) + qJ(4);
t54 = sin(t57);
t55 = cos(t57);
t114 = m(4) * pkin(2) + t64 * mrSges(4,1) - t60 * mrSges(4,2) + t126 * t55 - t83 * t54 + mrSges(3,1);
t121 = m(5) + m(6);
t30 = -t107 * t54 + t55 * t97;
t31 = t107 * t55 + t54 * t97;
t118 = -t124 * t30 + t125 * t31;
t106 = t58 * t62;
t17 = -t106 * t55 + t39 * t54;
t18 = t106 * t54 + t39 * t55;
t116 = t124 * t17 + t125 * t18;
t80 = t97 * t109;
t37 = t61 * t80 + t62 * t65;
t91 = t58 * t109;
t13 = -t37 * t54 - t55 * t91;
t14 = t37 * t55 - t54 * t91;
t115 = -t124 * t13 + t125 * t14;
t104 = t58 * t65;
t98 = t109 * pkin(1) + pkin(7) * t106;
t96 = t60 * t106;
t90 = -pkin(1) * t62 + pkin(7) * t91;
t89 = t13 * pkin(4) + t14 * pkin(10);
t88 = -t17 * pkin(4) + pkin(10) * t18;
t87 = t30 * pkin(4) + pkin(10) * t31;
t48 = t60 * t91;
t86 = -t37 * t64 + t48;
t82 = t19 * pkin(3);
t38 = t109 * t61 + t65 * t85;
t53 = pkin(3) * t64 + pkin(2);
t66 = -pkin(9) - pkin(8);
t81 = pkin(3) * t96 - t38 * t66 + t39 * t53 + t98;
t74 = t123 * pkin(3);
t69 = t37 * t60 + t64 * t91;
t68 = t69 * pkin(3);
t36 = t61 * t62 - t65 * t80;
t20 = t39 * t64 + t96;
t2 = t18 * t63 + t38 * t59;
t1 = -t18 * t59 + t38 * t63;
t3 = [(-t109 * mrSges(2,1) - m(3) * t98 - t39 * mrSges(3,1) - m(4) * (pkin(2) * t39 + t98) - t20 * mrSges(4,1) - t19 * mrSges(4,2) - m(5) * t81 - t18 * mrSges(5,1) - m(6) * (pkin(4) * t18 + t81) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + (-mrSges(3,3) * t58 + mrSges(2,2)) * t62 + t83 * t17 + t75 * t38) * g(2) + (t62 * mrSges(2,1) + t109 * mrSges(2,2) - m(3) * t90 + t37 * mrSges(3,1) - mrSges(3,3) * t91 - m(4) * (-pkin(2) * t37 + t90) - t86 * mrSges(4,1) - t69 * mrSges(4,2) + t83 * t13 + t126 * t14 + t117 * t36 + t121 * (-pkin(3) * t48 - t36 * t66 + t37 * t53 - t90)) * g(1), (-t121 * (-t36 * t53 - t37 * t66) - t117 * t37 + t114 * t36) * g(2) + (-t121 * (-t38 * t53 - t39 * t66) - t117 * t39 + t114 * t38) * g(1) + (-t121 * t53 * t104 + (-t114 * t65 + (t121 * t66 - t117) * t61) * t58) * g(3), (-t123 * mrSges(4,1) - (-t105 * t61 - t60 * t97) * mrSges(4,2) - m(5) * t74 - m(6) * (t74 + t87) + t118) * g(3) + (t69 * mrSges(4,1) - t86 * mrSges(4,2) + m(5) * t68 - m(6) * (-t68 + t89) + t115) * g(2) + (-t19 * mrSges(4,1) + t20 * mrSges(4,2) - m(5) * t82 - m(6) * (t82 + t88) + t116) * g(1), (-m(6) * t87 + t118) * g(3) + (-m(6) * t89 + t115) * g(2) + (-m(6) * t88 + t116) * g(1), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((-t14 * t59 + t36 * t63) * mrSges(6,1) + (-t14 * t63 - t36 * t59) * mrSges(6,2)) - g(3) * ((-t104 * t63 - t31 * t59) * mrSges(6,1) + (t104 * t59 - t31 * t63) * mrSges(6,2))];
taug = t3(:);
