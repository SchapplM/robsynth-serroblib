% Calculate joint inertia matrix for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:32:40
% EndTime: 2018-11-23 15:32:41
% DurationCPUTime: 1.11s
% Computational Cost: add. (1971->290), mult. (4078->421), div. (0->0), fcn. (4495->12), ass. (0->115)
t117 = sin(qJ(5));
t122 = cos(qJ(5));
t92 = -mrSges(6,1) * t122 + mrSges(6,2) * t117;
t175 = t92 - mrSges(5,1);
t114 = sin(pkin(6));
t125 = cos(qJ(2));
t153 = t114 * t125;
t118 = sin(qJ(4));
t123 = cos(qJ(4));
t115 = cos(pkin(6));
t119 = sin(qJ(3));
t124 = cos(qJ(3));
t120 = sin(qJ(2));
t154 = t114 * t120;
t74 = t115 * t124 - t119 * t154;
t75 = t115 * t119 + t124 * t154;
t43 = t118 * t74 + t123 * t75;
t28 = -t117 * t43 - t122 * t153;
t29 = -t117 * t153 + t122 * t43;
t140 = -t117 * t28 + t122 * t29;
t105 = -pkin(3) * t124 - pkin(2);
t87 = t118 * t119 - t123 * t124;
t89 = t118 * t124 + t119 * t123;
t53 = pkin(4) * t87 - pkin(10) * t89 + t105;
t168 = -pkin(9) - pkin(8);
t145 = t168 * t119;
t98 = t168 * t124;
t69 = t118 * t145 - t123 * t98;
t21 = -t117 * t69 + t122 * t53;
t22 = t117 * t53 + t122 * t69;
t141 = -t117 * t21 + t122 * t22;
t155 = t122 * t89;
t174 = Ifges(6,5) * t155 + Ifges(6,3) * t87;
t116 = sin(qJ(6));
t152 = Ifges(6,5) * t117 + Ifges(6,6) * t122;
t121 = cos(qJ(6));
t86 = -t116 * t117 + t121 * t122;
t165 = t86 * mrSges(7,3);
t173 = t116 * pkin(5) * t165 + t152;
t41 = t118 * t75 - t123 * t74;
t38 = t41 ^ 2;
t66 = -t118 * t98 - t123 * t145;
t172 = t66 ^ 2;
t88 = t116 * t122 + t117 * t121;
t56 = -mrSges(7,1) * t86 + mrSges(7,2) * t88;
t171 = 0.2e1 * t56;
t170 = 0.2e1 * t66;
t113 = t124 ^ 2;
t167 = pkin(3) * t123;
t166 = t41 * t66;
t164 = t88 * mrSges(7,3);
t163 = Ifges(7,5) * t88 + Ifges(7,6) * t86;
t162 = Ifges(6,4) * t117;
t161 = Ifges(6,4) * t122;
t158 = t117 * t89;
t151 = t117 ^ 2 + t122 ^ 2;
t150 = t119 ^ 2 + t113;
t149 = 0.2e1 * mrSges(7,3);
t44 = t88 * t89;
t45 = t86 * t89;
t148 = Ifges(7,5) * t45 - Ifges(7,6) * t44 + Ifges(7,3) * t87;
t147 = t121 * t164;
t7 = -t116 * t29 + t121 * t28;
t8 = t116 * t28 + t121 * t29;
t146 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t104 = -pkin(5) * t122 - pkin(4);
t102 = pkin(3) * t118 + pkin(10);
t144 = t151 * t102;
t58 = Ifges(7,4) * t88 + Ifges(7,2) * t86;
t59 = Ifges(7,1) * t88 + Ifges(7,4) * t86;
t94 = Ifges(6,2) * t122 + t162;
t95 = Ifges(6,1) * t117 + t161;
t143 = t117 * t95 + t122 * t94 + t86 * t58 + t88 * t59 + Ifges(5,3);
t142 = mrSges(6,1) * t117 + mrSges(6,2) * t122;
t54 = -mrSges(6,2) * t87 - mrSges(6,3) * t158;
t55 = mrSges(6,1) * t87 - mrSges(6,3) * t155;
t139 = -t117 * t55 + t122 * t54;
t138 = -t119 * t74 + t124 * t75;
t82 = (-pkin(11) - t102) * t117;
t108 = t122 * pkin(11);
t83 = t102 * t122 + t108;
t50 = -t116 * t83 + t121 * t82;
t51 = t116 * t82 + t121 * t83;
t137 = t50 * mrSges(7,1) - t51 * mrSges(7,2) + t163;
t96 = (-pkin(11) - pkin(10)) * t117;
t97 = pkin(10) * t122 + t108;
t65 = -t116 * t97 + t121 * t96;
t68 = t116 * t96 + t121 * t97;
t136 = t65 * mrSges(7,1) - t68 * mrSges(7,2) + t163;
t135 = 0.2e1 * t151 * mrSges(6,3);
t11 = pkin(5) * t87 - pkin(11) * t155 + t21;
t14 = -pkin(11) * t158 + t22;
t3 = t11 * t121 - t116 * t14;
t4 = t11 * t116 + t121 * t14;
t134 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t148;
t133 = (mrSges(5,1) * t123 - mrSges(5,2) * t118) * pkin(3);
t132 = (mrSges(7,1) * t121 - mrSges(7,2) * t116) * pkin(5);
t131 = -t43 * mrSges(5,2) - t7 * t164 + t8 * t165 + t140 * mrSges(6,3) + (t56 + t175) * t41;
t12 = Ifges(7,4) * t45 - Ifges(7,2) * t44 + Ifges(7,6) * t87;
t13 = Ifges(7,1) * t45 - Ifges(7,4) * t44 + Ifges(7,5) * t87;
t31 = Ifges(6,6) * t87 + (-Ifges(6,2) * t117 + t161) * t89;
t32 = Ifges(6,5) * t87 + (Ifges(6,1) * t122 - t162) * t89;
t36 = pkin(5) * t158 + t66;
t130 = -t69 * mrSges(5,2) - t3 * t164 + t4 * t165 + t88 * t13 / 0.2e1 + t36 * t56 - t44 * t58 / 0.2e1 + t45 * t59 / 0.2e1 + t117 * t32 / 0.2e1 + t122 * t31 / 0.2e1 - t94 * t158 / 0.2e1 + t95 * t155 / 0.2e1 - Ifges(5,6) * t87 + Ifges(5,5) * t89 + t86 * t12 / 0.2e1 + t175 * t66 + (t163 + t152) * t87 / 0.2e1 + t141 * mrSges(6,3);
t109 = t114 ^ 2;
t103 = -pkin(4) - t167;
t100 = t109 * t125 ^ 2;
t93 = -mrSges(4,1) * t124 + mrSges(4,2) * t119;
t91 = t104 - t167;
t57 = mrSges(5,1) * t87 + mrSges(5,2) * t89;
t52 = t142 * t89;
t25 = mrSges(7,1) * t87 - mrSges(7,3) * t45;
t24 = -mrSges(7,2) * t87 - mrSges(7,3) * t44;
t15 = mrSges(7,1) * t44 + mrSges(7,2) * t45;
t1 = [m(2) + m(7) * (t7 ^ 2 + t8 ^ 2 + t38) + m(6) * (t28 ^ 2 + t29 ^ 2 + t38) + m(5) * (t43 ^ 2 + t100 + t38) + m(4) * (t74 ^ 2 + t75 ^ 2 + t100) + m(3) * (t109 * t120 ^ 2 + t115 ^ 2 + t100); -t43 * t87 * mrSges(5,3) + t8 * t24 + t7 * t25 + t28 * t55 + t29 * t54 + t138 * mrSges(4,3) + (t89 * mrSges(5,3) + t15 + t52) * t41 + (-t120 * mrSges(3,2) + (mrSges(3,1) - t57 - t93) * t125) * t114 + m(7) * (t3 * t7 + t36 * t41 + t4 * t8) + m(6) * (t21 * t28 + t22 * t29 + t166) + m(5) * (-t105 * t153 + t43 * t69 + t166) + m(4) * (pkin(2) * t153 + t138 * pkin(8)); Ifges(4,2) * t113 - 0.2e1 * pkin(2) * t93 + 0.2e1 * t105 * t57 - t44 * t12 + t45 * t13 + 0.2e1 * t36 * t15 + 0.2e1 * t21 * t55 + 0.2e1 * t22 * t54 + 0.2e1 * t4 * t24 + 0.2e1 * t3 * t25 + t52 * t170 + Ifges(3,3) + (Ifges(4,1) * t119 + 0.2e1 * Ifges(4,4) * t124) * t119 + 0.2e1 * t150 * pkin(8) * mrSges(4,3) + (mrSges(5,3) * t170 + Ifges(5,1) * t89 - t117 * t31 + t122 * t32) * t89 + (-0.2e1 * t69 * mrSges(5,3) + Ifges(5,2) * t87 + (-Ifges(6,6) * t117 - (2 * Ifges(5,4))) * t89 + t148 + t174) * t87 + m(4) * (t150 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t105 ^ 2 + t69 ^ 2 + t172) + m(6) * (t21 ^ 2 + t22 ^ 2 + t172) + m(7) * (t3 ^ 2 + t36 ^ 2 + t4 ^ 2); t74 * mrSges(4,1) - t75 * mrSges(4,2) + m(7) * (t41 * t91 + t50 * t7 + t51 * t8) + m(6) * (t140 * t102 + t103 * t41) + m(5) * (t118 * t43 - t123 * t41) * pkin(3) + t131; (m(5) * (t118 * t69 - t123 * t66) + (-t118 * t87 - t123 * t89) * mrSges(5,3)) * pkin(3) + m(6) * (t141 * t102 + t103 * t66) + t139 * t102 + m(7) * (t3 * t50 + t36 * t91 + t4 * t51) + (-t119 * mrSges(4,1) - t124 * mrSges(4,2)) * pkin(8) + Ifges(4,6) * t124 + Ifges(4,5) * t119 + t91 * t15 + t103 * t52 + t50 * t25 + t51 * t24 + t130; 0.2e1 * t103 * t92 + t91 * t171 + Ifges(4,3) + 0.2e1 * t133 + (-t50 * t88 + t51 * t86) * t149 + t102 * t135 + m(7) * (t50 ^ 2 + t51 ^ 2 + t91 ^ 2) + m(6) * (t151 * t102 ^ 2 + t103 ^ 2) + m(5) * (t118 ^ 2 + t123 ^ 2) * pkin(3) ^ 2 + t143; m(7) * (t104 * t41 + t65 * t7 + t68 * t8) + m(6) * (-pkin(4) * t41 + t140 * pkin(10)) + t131; m(6) * (-pkin(4) * t66 + t141 * pkin(10)) + t139 * pkin(10) + m(7) * (t104 * t36 + t3 * t65 + t4 * t68) + t104 * t15 - pkin(4) * t52 + t65 * t25 + t68 * t24 + t130; (t103 - pkin(4)) * t92 + (t91 + t104) * t56 + t133 + m(7) * (t104 * t91 + t50 * t65 + t51 * t68) + m(6) * (-pkin(4) * t103 + pkin(10) * t144) + ((-t50 - t65) * t88 + (t51 + t68) * t86) * mrSges(7,3) + (t151 * pkin(10) + t144) * mrSges(6,3) + t143; -0.2e1 * pkin(4) * t92 + t104 * t171 + (-t65 * t88 + t68 * t86) * t149 + pkin(10) * t135 + m(7) * (t104 ^ 2 + t65 ^ 2 + t68 ^ 2) + m(6) * (t151 * pkin(10) ^ 2 + pkin(4) ^ 2) + t143; t28 * mrSges(6,1) - t29 * mrSges(6,2) + m(7) * (t116 * t8 + t121 * t7) * pkin(5) + t146; -Ifges(6,6) * t158 + t21 * mrSges(6,1) - t22 * mrSges(6,2) + (m(7) * (t116 * t4 + t121 * t3) + t116 * t24 + t121 * t25) * pkin(5) + t134 + t174; -t142 * t102 + (m(7) * (t116 * t51 + t121 * t50) - t147) * pkin(5) + t137 + t173; -t142 * pkin(10) + (m(7) * (t116 * t68 + t121 * t65) - t147) * pkin(5) + t136 + t173; Ifges(6,3) + Ifges(7,3) + m(7) * (t116 ^ 2 + t121 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t132; t146; t134; t137; t136; Ifges(7,3) + t132; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
