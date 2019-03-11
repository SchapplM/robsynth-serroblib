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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

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
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:53
% EndTime: 2019-03-08 22:22:57
% DurationCPUTime: 1.07s
% Computational Cost: add. (988->126), mult. (2364->196), div. (0->0), fcn. (2947->16), ass. (0->70)
t65 = sin(qJ(6));
t67 = cos(qJ(6));
t127 = m(7) * pkin(5) + mrSges(7,1) * t67 - mrSges(7,2) * t65 + mrSges(6,1);
t98 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t120 = -m(6) - m(7);
t126 = m(5) - t120;
t103 = -m(4) - t126;
t125 = -pkin(2) * t103 + mrSges(3,1);
t61 = sin(pkin(13));
t63 = cos(pkin(13));
t76 = -m(5) * pkin(3) - t63 * mrSges(5,1) + t61 * mrSges(5,2) - mrSges(4,1);
t107 = cos(pkin(7));
t62 = sin(pkin(7));
t104 = sin(pkin(12));
t116 = sin(qJ(2));
t118 = cos(qJ(2));
t106 = cos(pkin(12));
t108 = cos(pkin(6));
t90 = t108 * t106;
t72 = t104 * t116 - t118 * t90;
t105 = sin(pkin(6));
t87 = t106 * t105;
t124 = t72 * t107 + t62 * t87;
t89 = t108 * t104;
t73 = t106 * t116 + t118 * t89;
t86 = t105 * t104;
t123 = t73 * t107 - t62 * t86;
t88 = t107 * t105;
t122 = t108 * t62 + t118 * t88;
t60 = pkin(13) + qJ(5);
t58 = sin(t60);
t59 = cos(t60);
t121 = t127 * t59 - t98 * t58 - t76;
t74 = -m(5) * qJ(4) - mrSges(7,1) * t65 - mrSges(7,2) * t67 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t117 = cos(qJ(3));
t115 = t58 * t62;
t114 = t59 * t62;
t93 = t105 * t116;
t85 = t62 * t93;
t94 = t118 * t105;
t109 = pkin(2) * t94 + pkin(9) * t85;
t66 = sin(qJ(3));
t100 = t66 * t107;
t95 = t107 * t117;
t82 = t61 * t85;
t78 = t116 * t88;
t68 = mrSges(3,2) + (-t63 * mrSges(5,2) - mrSges(4,3) + (t120 * pkin(4) - mrSges(5,1)) * t61 + t103 * pkin(9)) * t62;
t64 = -pkin(10) - qJ(4);
t57 = pkin(4) * t63 + pkin(3);
t47 = t106 * t118 - t116 * t89;
t46 = t104 * t118 + t116 * t90;
t45 = t108 * t107 - t62 * t94;
t42 = t117 * t94 - t66 * t78;
t41 = t117 * t78 + t66 * t94;
t34 = t107 * t86 + t73 * t62;
t33 = -t107 * t87 + t72 * t62;
t32 = t117 * t93 + t122 * t66;
t31 = -t122 * t117 + t66 * t93;
t26 = -t47 * t100 - t73 * t117;
t25 = t47 * t95 - t73 * t66;
t24 = -t46 * t100 - t72 * t117;
t23 = t46 * t95 - t72 * t66;
t18 = t47 * t117 - t123 * t66;
t17 = t123 * t117 + t47 * t66;
t16 = t46 * t117 - t124 * t66;
t15 = t124 * t117 + t46 * t66;
t12 = t32 * t59 + t45 * t58;
t4 = t18 * t59 + t34 * t58;
t2 = t16 * t59 + t33 * t58;
t1 = [(-m(2) - m(3) + t103) * g(3) (-mrSges(3,1) * t94 + mrSges(3,2) * t93 - m(4) * t109 - t42 * mrSges(4,1) - mrSges(4,3) * t85 - m(5) * (pkin(3) * t42 + t109) - (t42 * t63 + t82) * mrSges(5,1) - (-t42 * t61 + t63 * t85) * mrSges(5,2) + t98 * (t42 * t58 - t59 * t85) - t127 * (t42 * t59 + t58 * t85) + t74 * t41 + t120 * (pkin(4) * t82 - t41 * t64 + t42 * t57 + t109)) * g(3) + (t98 * (-t46 * t114 + t24 * t58) - t127 * (t46 * t115 + t24 * t59) + t76 * t24 + t74 * t23 + t68 * t46 + t125 * t72 + t120 * (-t23 * t64 + t24 * t57)) * g(2) + (t98 * (-t47 * t114 + t26 * t58) - t127 * (t47 * t115 + t26 * t59) + t76 * t26 + t74 * t25 + t68 * t47 + t125 * t73 + t120 * (-t25 * t64 + t26 * t57)) * g(1) (t120 * (-t31 * t57 - t32 * t64) + t74 * t32 + t121 * t31) * g(3) + (t120 * (-t15 * t57 - t16 * t64) + t74 * t16 + t121 * t15) * g(2) + (t120 * (-t17 * t57 - t18 * t64) + t74 * t18 + t121 * t17) * g(1), t126 * (-g(1) * t17 - g(2) * t15 - g(3) * t31) (t98 * t12 - t127 * (-t32 * t58 + t45 * t59)) * g(3) + (t98 * t2 - t127 * (-t16 * t58 + t33 * t59)) * g(2) + (t98 * t4 - t127 * (-t18 * t58 + t34 * t59)) * g(1), -g(1) * ((t17 * t67 - t4 * t65) * mrSges(7,1) + (-t17 * t65 - t4 * t67) * mrSges(7,2)) - g(2) * ((t15 * t67 - t2 * t65) * mrSges(7,1) + (-t15 * t65 - t2 * t67) * mrSges(7,2)) - g(3) * ((-t12 * t65 + t31 * t67) * mrSges(7,1) + (-t12 * t67 - t31 * t65) * mrSges(7,2))];
taug  = t1(:);
