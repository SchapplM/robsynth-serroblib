% Calculate joint inertia matrix for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:27:33
% EndTime: 2018-11-23 15:27:34
% DurationCPUTime: 1.45s
% Computational Cost: add. (1284->257), mult. (2684->349), div. (0->0), fcn. (2765->10), ass. (0->98)
t173 = Ifges(6,4) + Ifges(7,4);
t168 = Ifges(6,6) + Ifges(7,6);
t102 = sin(qJ(4));
t103 = sin(qJ(3));
t106 = cos(qJ(4));
t107 = cos(qJ(3));
t67 = t102 * t103 - t106 * t107;
t172 = t168 * t67;
t171 = Ifges(6,1) + Ifges(7,1);
t170 = Ifges(6,5) + Ifges(7,5);
t169 = Ifges(6,2) + Ifges(7,2);
t105 = cos(qJ(5));
t167 = t173 * t105;
t101 = sin(qJ(5));
t166 = t173 * t101;
t108 = cos(qJ(2));
t99 = sin(pkin(6));
t129 = t108 * t99;
t100 = cos(pkin(6));
t104 = sin(qJ(2));
t133 = t104 * t99;
t55 = t100 * t107 - t103 * t133;
t56 = t100 * t103 + t107 * t133;
t27 = t102 * t55 + t106 * t56;
t15 = -t101 * t27 - t105 * t129;
t16 = -t101 * t129 + t105 * t27;
t116 = -t101 * t15 + t105 * t16;
t73 = -mrSges(6,1) * t105 + mrSges(6,2) * t101;
t165 = -mrSges(5,1) + t73;
t163 = Ifges(6,3) + Ifges(7,3);
t68 = t102 * t107 + t103 * t106;
t162 = (-t169 * t101 + t167) * t68 + t172;
t161 = (t171 * t105 - t166) * t68 + t170 * t67;
t160 = t169 * t105 + t166;
t159 = t171 * t101 + t167;
t122 = t170 * t101 + t168 * t105;
t25 = t102 * t56 - t106 * t55;
t158 = t25 ^ 2;
t152 = -pkin(9) - pkin(8);
t127 = t152 * t103;
t80 = t152 * t107;
t45 = -t102 * t80 - t106 * t127;
t157 = t45 ^ 2;
t156 = 0.2e1 * t45;
t72 = -t105 * mrSges(7,1) + t101 * mrSges(7,2);
t155 = 0.2e1 * t72;
t98 = t107 ^ 2;
t154 = m(7) * pkin(5);
t149 = pkin(3) * t106;
t148 = pkin(10) * t105;
t87 = -pkin(3) * t107 - pkin(2);
t33 = pkin(4) * t67 - pkin(10) * t68 + t87;
t47 = t102 * t127 - t106 * t80;
t7 = t101 * t33 + t105 * t47;
t147 = t105 * t7;
t146 = t25 * t45;
t131 = t105 * t68;
t134 = t101 * t68;
t31 = mrSges(7,1) * t134 + mrSges(7,2) * t131;
t142 = t101 ^ 2 + t105 ^ 2;
t141 = t103 ^ 2 + t98;
t136 = t101 * mrSges(7,3);
t84 = pkin(3) * t102 + pkin(10);
t130 = t105 * t84;
t88 = t105 * qJ(6);
t128 = 0.2e1 * mrSges(7,3);
t86 = -pkin(5) * t105 - pkin(4);
t124 = t142 * t84;
t6 = -t101 * t47 + t105 * t33;
t123 = t170 * t131 + t163 * t67;
t121 = t159 * t101 + t160 * t105 + Ifges(5,3);
t2 = pkin(5) * t67 - t68 * t88 + t6;
t120 = -t6 * mrSges(6,3) - t2 * mrSges(7,3);
t119 = -t101 * t6 + t147;
t118 = 0.2e1 * mrSges(6,3) * t142;
t117 = mrSges(6,1) * t101 + mrSges(6,2) * t105;
t115 = -t103 * t55 + t107 * t56;
t114 = (mrSges(5,1) * t106 - mrSges(5,2) * t102) * pkin(3);
t113 = -t27 * mrSges(5,2) + (t72 + t165) * t25 + t116 * (mrSges(6,3) + mrSges(7,3));
t23 = pkin(5) * t134 + t45;
t4 = -qJ(6) * t134 + t7;
t112 = -t47 * mrSges(5,2) + mrSges(6,3) * t147 + Ifges(5,5) * t68 + t23 * t72 + t165 * t45 + t161 * t101 / 0.2e1 - t160 * t134 / 0.2e1 + t159 * t131 / 0.2e1 + (-Ifges(5,6) + t122 / 0.2e1) * t67 + (t4 * mrSges(7,3) + t162 / 0.2e1) * t105;
t94 = t99 ^ 2;
t85 = -pkin(4) - t149;
t82 = t94 * t108 ^ 2;
t75 = t88 + t148;
t74 = -mrSges(4,1) * t107 + mrSges(4,2) * t103;
t71 = (-qJ(6) - pkin(10)) * t101;
t70 = t86 - t149;
t58 = t88 + t130;
t57 = (-qJ(6) - t84) * t101;
t38 = mrSges(5,1) * t67 + mrSges(5,2) * t68;
t37 = mrSges(6,1) * t67 - mrSges(6,3) * t131;
t36 = mrSges(7,1) * t67 - mrSges(7,3) * t131;
t35 = -mrSges(6,2) * t67 - mrSges(6,3) * t134;
t34 = -mrSges(7,2) * t67 - mrSges(7,3) * t134;
t32 = t117 * t68;
t1 = [m(2) + m(5) * (t27 ^ 2 + t158 + t82) + m(4) * (t55 ^ 2 + t56 ^ 2 + t82) + m(3) * (t104 ^ 2 * t94 + t100 ^ 2 + t82) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t15 ^ 2 + t16 ^ 2 + t158); -t27 * t67 * mrSges(5,3) + (t34 + t35) * t16 + (t36 + t37) * t15 + t115 * mrSges(4,3) + (t68 * mrSges(5,3) + t31 + t32) * t25 + (-t104 * mrSges(3,2) + (mrSges(3,1) - t38 - t74) * t108) * t99 + m(6) * (t15 * t6 + t16 * t7 + t146) + m(7) * (t15 * t2 + t16 * t4 + t23 * t25) + m(5) * (-t129 * t87 + t27 * t47 + t146) + m(4) * (pkin(2) * t129 + pkin(8) * t115); Ifges(4,2) * t98 - 0.2e1 * pkin(2) * t74 + 0.2e1 * t2 * t36 + 0.2e1 * t23 * t31 + t32 * t156 + 0.2e1 * t4 * t34 + 0.2e1 * t7 * t35 + 0.2e1 * t6 * t37 + 0.2e1 * t87 * t38 + Ifges(3,3) + (Ifges(4,1) * t103 + 0.2e1 * Ifges(4,4) * t107) * t103 + 0.2e1 * t141 * pkin(8) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t47 + Ifges(5,2) * t67 + t123) * t67 + m(4) * (pkin(8) ^ 2 * t141 + pkin(2) ^ 2) + m(5) * (t47 ^ 2 + t87 ^ 2 + t157) + m(6) * (t6 ^ 2 + t7 ^ 2 + t157) + m(7) * (t2 ^ 2 + t23 ^ 2 + t4 ^ 2) + (mrSges(5,3) * t156 + Ifges(5,1) * t68 - 0.2e1 * Ifges(5,4) * t67 + t161 * t105 + (-t162 - t172) * t101) * t68; t55 * mrSges(4,1) - t56 * mrSges(4,2) + m(6) * (t116 * t84 + t25 * t85) + m(7) * (t15 * t57 + t16 * t58 + t25 * t70) + m(5) * (t102 * t27 - t106 * t25) * pkin(3) + t113; (-t103 * mrSges(4,1) - t107 * mrSges(4,2)) * pkin(8) + t112 + (m(5) * (t102 * t47 - t106 * t45) + (-t102 * t67 - t106 * t68) * mrSges(5,3)) * pkin(3) + Ifges(4,6) * t107 + Ifges(4,5) * t103 + t85 * t32 + t70 * t31 + t57 * t36 + t58 * t34 + t35 * t130 + m(7) * (t2 * t57 + t23 * t70 + t4 * t58) + (-t84 * t37 + t120) * t101 + m(6) * (t119 * t84 + t45 * t85); t70 * t155 + 0.2e1 * t85 * t73 + Ifges(4,3) + 0.2e1 * t114 + (-t57 * t101 + t58 * t105) * t128 + t84 * t118 + m(6) * (t142 * t84 ^ 2 + t85 ^ 2) + m(7) * (t57 ^ 2 + t58 ^ 2 + t70 ^ 2) + m(5) * (t102 ^ 2 + t106 ^ 2) * pkin(3) ^ 2 + t121; m(6) * (-pkin(4) * t25 + pkin(10) * t116) + m(7) * (t15 * t71 + t16 * t75 + t25 * t86) + t113; (-pkin(10) * t37 + t120) * t101 + m(7) * (t2 * t71 + t23 * t86 + t4 * t75) + t112 + t86 * t31 + t71 * t36 + t75 * t34 - pkin(4) * t32 + t35 * t148 + m(6) * (-pkin(4) * t45 + pkin(10) * t119); (t85 - pkin(4)) * t73 + (t70 + t86) * t72 + t114 + m(6) * (-pkin(4) * t85 + pkin(10) * t124) + m(7) * (t57 * t71 + t58 * t75 + t70 * t86) + ((t58 + t75) * t105 + (-t57 - t71) * t101) * mrSges(7,3) + (pkin(10) * t142 + t124) * mrSges(6,3) + t121; -0.2e1 * pkin(4) * t73 + t86 * t155 + (-t71 * t101 + t75 * t105) * t128 + pkin(10) * t118 + m(7) * (t71 ^ 2 + t75 ^ 2 + t86 ^ 2) + m(6) * (pkin(10) ^ 2 * t142 + pkin(4) ^ 2) + t121; (-mrSges(6,2) - mrSges(7,2)) * t16 + (mrSges(6,1) + mrSges(7,1) + t154) * t15; mrSges(6,1) * t6 + mrSges(7,1) * t2 - mrSges(6,2) * t7 - mrSges(7,2) * t4 - t168 * t134 + (m(7) * t2 + t36) * pkin(5) + t123; mrSges(7,1) * t57 - mrSges(7,2) * t58 - t117 * t84 + (m(7) * t57 - t136) * pkin(5) + t122; mrSges(7,1) * t71 - mrSges(7,2) * t75 - t117 * pkin(10) + (m(7) * t71 - t136) * pkin(5) + t122; (0.2e1 * mrSges(7,1) + t154) * pkin(5) + t163; m(7) * t25; m(7) * t23 + t31; m(7) * t70 + t72; m(7) * t86 + t72; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
