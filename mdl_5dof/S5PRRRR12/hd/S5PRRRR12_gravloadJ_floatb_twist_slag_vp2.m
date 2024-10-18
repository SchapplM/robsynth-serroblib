% Calculate Gravitation load on the joints for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR12_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR12_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:20
% EndTime: 2024-09-28 18:07:22
% DurationCPUTime: 1.09s
% Computational Cost: add. (680->159), mult. (916->257), div. (0->0), fcn. (925->26), ass. (0->105)
t82 = sin(pkin(6));
t142 = pkin(10) * t82;
t81 = sin(pkin(11));
t110 = t81 * t142;
t84 = cos(pkin(11));
t86 = cos(pkin(5));
t48 = t84 * pkin(4) + t86 * t110;
t109 = t84 * t142;
t131 = t81 * t86;
t50 = pkin(4) * t131 - t109;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t13 = -t48 * t92 + t50 * t88;
t10 = -t84 * pkin(3) + t13;
t151 = t48 * t88 + t50 * t92;
t7 = -pkin(3) * t131 - t151;
t89 = sin(qJ(3));
t93 = cos(qJ(3));
t156 = t10 * t89 + t7 * t93;
t49 = -t81 * pkin(4) + t86 * t109;
t126 = t84 * t86;
t51 = pkin(4) * t126 + t110;
t14 = t49 * t92 - t51 * t88;
t11 = -pkin(3) * t81 + t14;
t102 = t49 * t88 + t51 * t92;
t8 = pkin(3) * t126 + t102;
t155 = t11 * t89 + t8 * t93;
t154 = t10 * t93 - t7 * t89;
t153 = t11 * t93 - t8 * t89;
t144 = m(4) * pkin(2);
t152 = -mrSges(3,1) - t144;
t85 = cos(pkin(6));
t91 = cos(qJ(5));
t124 = t85 * t91;
t87 = sin(qJ(5));
t125 = t85 * t87;
t83 = sin(pkin(5));
t130 = t82 * t83;
t80 = qJ(2) + qJ(3);
t75 = pkin(5) - t80;
t70 = -qJ(4) + t75;
t60 = sin(t70) / 0.2e1;
t74 = pkin(5) + t80;
t69 = qJ(4) + t74;
t61 = cos(t69) / 0.2e1;
t65 = sin(t69);
t66 = cos(t70);
t79 = qJ(4) + t80;
t71 = sin(t79);
t72 = cos(t79);
t150 = -t71 * mrSges(6,3) * t130 - ((-t71 * t125 + t72 * t91) * mrSges(6,1) + (-t71 * t124 - t72 * t87) * mrSges(6,2)) * t83 - (t65 / 0.2e1 + t60) * mrSges(5,1) - (t61 - t66 / 0.2e1) * mrSges(5,2);
t133 = t81 * t72;
t139 = mrSges(6,3) * t82;
t123 = t86 * t87;
t106 = t85 * t123;
t32 = t84 * t106 + t81 * t91;
t121 = t86 * t91;
t105 = t85 * t121;
t33 = t84 * t105 - t81 * t87;
t34 = -t84 * t121 + t81 * t125;
t35 = -t84 * t123 - t81 * t124;
t44 = t60 - t65 / 0.2e1;
t45 = t66 / 0.2e1 + t61;
t149 = -(-t32 * t71 - t34 * t72) * mrSges(6,1) - (-t33 * t71 + t35 * t72) * mrSges(6,2) - (t71 * t126 + t133) * t139 - (t84 * t45 - t81 * t71) * mrSges(5,1) - (t84 * t44 - t133) * mrSges(5,2);
t128 = t84 * t72;
t30 = t81 * t106 - t84 * t91;
t31 = -t81 * t105 - t84 * t87;
t36 = t81 * t121 + t84 * t125;
t37 = t81 * t123 - t84 * t124;
t148 = -(t30 * t71 - t36 * t72) * mrSges(6,1) - (-t31 * t71 + t37 * t72) * mrSges(6,2) - (-t71 * t131 + t128) * t139 - (-t81 * t45 - t84 * t71) * mrSges(5,1) - (-t81 * t44 - t128) * mrSges(5,2);
t62 = sin(t75) / 0.2e1;
t63 = cos(t74) / 0.2e1;
t67 = sin(t74);
t68 = cos(t75);
t147 = -(t67 / 0.2e1 + t62) * mrSges(4,1) - (t63 - t68 / 0.2e1) * mrSges(4,2) + t150;
t76 = sin(t80);
t132 = t81 * t76;
t53 = t62 - t67 / 0.2e1;
t54 = t68 / 0.2e1 + t63;
t77 = cos(t80);
t146 = -(t84 * t54 - t132) * mrSges(4,1) - (t84 * t53 - t81 * t77) * mrSges(4,2) + t149;
t127 = t84 * t76;
t145 = -(-t81 * t54 - t127) * mrSges(4,1) - (-t81 * t53 - t84 * t77) * mrSges(4,2) + t148;
t143 = m(5) * t83;
t58 = -t88 * pkin(4) + t92 * t142;
t134 = t58 * t89;
t94 = cos(qJ(2));
t129 = t83 * t94;
t90 = sin(qJ(2));
t122 = t86 * t90;
t120 = t86 * t94;
t119 = t89 * t90;
t108 = t87 * t130;
t107 = t91 * t130;
t57 = -pkin(4) * t92 - t88 * t142;
t56 = pkin(3) - t57;
t99 = -t56 * t93 - t134;
t98 = -pkin(3) * t119 + (t93 * pkin(3) + pkin(2)) * t94;
t95 = pkin(3) * (t93 * t94 - t119);
t59 = -t90 * pkin(2) - pkin(3) * t76;
t52 = t58 * t93;
t29 = t86 * t95;
t28 = t98 * t86;
t15 = (-t89 * t56 + t52) * t83 * t90;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), (-(mrSges(3,1) * t94 - mrSges(3,2) * t90) * t83 - t129 * t144 - t98 * t143 - m(6) * ((pkin(2) - t99) * t129 + t15) + t147) * g(3) + (-(-t84 * t122 - t81 * t94) * mrSges(3,2) - m(5) * (t84 * t28 + t81 * t59) - m(6) * (-(t81 * pkin(2) - t153) * t90 + (pkin(2) * t126 + t155) * t94) + t152 * (t84 * t120 - t81 * t90) + t146) * g(2) + (-(t81 * t122 - t84 * t94) * mrSges(3,2) - m(5) * (-t81 * t28 + t84 * t59) - m(6) * (-(t84 * pkin(2) - t154) * t90 + (-pkin(2) * t131 + t156) * t94) + t152 * (-t81 * t120 - t84 * t90) + t145) * g(1), (-t95 * t143 - m(6) * (-t99 * t129 + t15) + t147) * g(3) + (-m(5) * (-pkin(3) * t132 + t84 * t29) - m(6) * (t153 * t90 + t155 * t94) + t146) * g(2) + (-m(5) * (-pkin(3) * t127 - t81 * t29) - m(6) * (t154 * t90 + t156 * t94) + t145) * g(1), (-m(6) * ((t57 * t89 + t52) * t90 - (t57 * t93 - t134) * t94) * t83 + t150) * g(3) + (-m(6) * ((t102 * t93 + t14 * t89) * t94 + (-t102 * t89 + t14 * t93) * t90) + t149) * g(2) + (-m(6) * ((t13 * t89 - t151 * t93) * t94 + (t13 * t93 + t151 * t89) * t90) + t148) * g(1), -g(1) * ((t81 * t107 + t31 * t72 + t37 * t71) * mrSges(6,1) + (-t81 * t108 + t30 * t72 + t36 * t71) * mrSges(6,2)) - g(2) * ((-t84 * t107 + t33 * t72 + t35 * t71) * mrSges(6,1) + (t84 * t108 - t32 * t72 + t34 * t71) * mrSges(6,2)) - g(3) * ((mrSges(6,1) * t91 - mrSges(6,2) * t87) * t86 * t82 + ((t72 * t124 - t71 * t87) * mrSges(6,1) + (-t72 * t125 - t71 * t91) * mrSges(6,2)) * t83)];
taug = t1(:);
