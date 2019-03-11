% Calculate joint inertia matrix for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:50
% EndTime: 2019-03-09 20:53:54
% DurationCPUTime: 1.49s
% Computational Cost: add. (1518->293), mult. (2852->372), div. (0->0), fcn. (2727->6), ass. (0->110)
t195 = Ifges(5,1) + Ifges(7,3);
t194 = Ifges(5,5) + Ifges(7,5);
t192 = Ifges(6,5) + Ifges(7,4);
t193 = Ifges(7,2) + Ifges(6,3);
t122 = cos(qJ(4));
t119 = sin(qJ(4));
t160 = Ifges(7,6) * t119;
t164 = Ifges(5,4) * t119;
t120 = sin(qJ(3));
t121 = sin(qJ(2));
t123 = cos(qJ(3));
t124 = cos(qJ(2));
t76 = t120 * t121 - t123 * t124;
t77 = t120 * t124 + t121 * t123;
t191 = (t195 * t122 + t160 - t164) * t77 + t194 * t76;
t159 = Ifges(7,6) * t122;
t161 = Ifges(6,6) * t122;
t190 = (t193 * t119 + t159 - t161) * t77 + t192 * t76;
t162 = Ifges(6,6) * t119;
t189 = -t193 * t122 + t160 - t162;
t163 = Ifges(5,4) * t122;
t188 = t195 * t119 - t159 + t163;
t187 = t119 ^ 2 + t122 ^ 2;
t101 = -pkin(2) * t124 - pkin(1);
t32 = pkin(3) * t76 - pkin(9) * t77 + t101;
t176 = -pkin(8) - pkin(7);
t147 = t176 * t121;
t93 = t176 * t124;
t50 = t120 * t147 - t123 * t93;
t39 = t119 * t50;
t12 = t122 * t32 - t39;
t13 = t119 * t32 + t122 * t50;
t137 = -t119 * t12 + t122 * t13;
t186 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t185 = (Ifges(5,6) - t192) * t122 + t194 * t119;
t184 = (mrSges(5,3) + mrSges(6,1)) * t187;
t48 = -t120 * t93 - t123 * t147;
t183 = t48 ^ 2;
t182 = 2 * mrSges(7,1);
t181 = 0.2e1 * t48;
t82 = mrSges(6,2) * t122 - mrSges(6,3) * t119;
t180 = 0.2e1 * t82;
t84 = -mrSges(7,2) * t119 - mrSges(7,3) * t122;
t179 = 0.2e1 * t84;
t178 = 0.2e1 * t101;
t172 = Ifges(6,4) * t76;
t171 = Ifges(5,6) * t76;
t170 = pkin(2) * t123;
t169 = pkin(9) * t119;
t7 = -qJ(5) * t76 - t13;
t168 = t122 * t7;
t118 = -pkin(4) - qJ(6);
t154 = t122 * t77;
t38 = mrSges(6,1) * t154 + t76 * mrSges(6,2);
t99 = pkin(2) * t120 + pkin(9);
t167 = t187 * pkin(9) * t99;
t165 = t187 * t99 ^ 2;
t105 = t119 * mrSges(6,1);
t104 = t119 * mrSges(7,1);
t157 = t119 * t77;
t156 = t119 * t99;
t106 = t122 * mrSges(7,1);
t153 = qJ(5) * t122;
t152 = t104 + t105;
t150 = t187 * pkin(9) ^ 2;
t148 = t121 ^ 2 + t124 ^ 2;
t35 = mrSges(7,1) * t154 - t76 * mrSges(7,3);
t144 = -qJ(5) * t119 - pkin(3);
t89 = Ifges(5,2) * t122 + t164;
t143 = t188 * t119 + t122 * t89 + Ifges(4,3);
t37 = -mrSges(7,1) * t157 + t76 * mrSges(7,2);
t141 = pkin(4) * t157 + t48;
t8 = -pkin(4) * t76 - t12;
t140 = t119 * t8 - t168;
t139 = t119 * mrSges(5,1) + t122 * mrSges(5,2);
t138 = -t119 * mrSges(6,2) - t122 * mrSges(6,3);
t79 = -pkin(4) * t122 + t144;
t136 = t194 * t154 + t192 * t157 + t186 * t76;
t135 = (mrSges(4,1) * t123 - mrSges(4,2) * t120) * pkin(2);
t62 = t118 * t122 + t144;
t33 = -mrSges(5,2) * t76 - mrSges(5,3) * t157;
t34 = mrSges(5,1) * t76 - mrSges(5,3) * t154;
t36 = mrSges(6,1) * t157 - mrSges(6,3) * t76;
t134 = (t33 - t36) * t122 + (-t34 + t38) * t119;
t132 = 0.2e1 * t184;
t131 = t118 * t104 + (-pkin(4) * mrSges(6,1) - Ifges(6,4)) * t119 + (mrSges(7,1) + mrSges(6,1)) * t153 + t185;
t130 = m(6) * (-pkin(4) * t119 + t153) - t138 - t139;
t14 = -t77 * t153 + t141;
t2 = t39 + t118 * t76 + (pkin(5) * t77 - t32) * t122;
t21 = t171 + (-Ifges(5,2) * t119 + t163) * t77;
t26 = t172 + (-Ifges(6,2) * t122 + t162) * t77;
t4 = -pkin(5) * t157 - t7;
t83 = -mrSges(5,1) * t122 + mrSges(5,2) * t119;
t88 = -Ifges(6,2) * t119 - t161;
t9 = (qJ(6) * t119 - t153) * t77 + t141;
t129 = t4 * t106 + t2 * t104 + t8 * t105 - t119 * t26 / 0.2e1 + t14 * t82 + t9 * t84 - Ifges(4,6) * t76 + Ifges(4,5) * t77 - t50 * mrSges(4,2) - mrSges(6,1) * t168 + (t83 - mrSges(4,1)) * t48 + t191 * t119 / 0.2e1 + t137 * mrSges(5,3) + (-Ifges(6,4) * t119 + t185) * t76 / 0.2e1 + (-t89 / 0.2e1 + t189 / 0.2e1) * t157 + (-t88 / 0.2e1 + t188 / 0.2e1) * t154 + (t21 / 0.2e1 - t190 / 0.2e1) * t122;
t125 = qJ(5) ^ 2;
t113 = t122 * pkin(5);
t111 = t119 * pkin(5);
t100 = -pkin(3) - t170;
t92 = pkin(9) * t122 + t113;
t91 = t111 + t169;
t72 = t122 * t99 + t113;
t71 = t111 + t156;
t63 = t79 - t170;
t54 = t62 - t170;
t31 = t139 * t77;
t30 = t138 * t77;
t29 = (-mrSges(7,2) * t122 + mrSges(7,3) * t119) * t77;
t1 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t124 + mrSges(3,2) * t121) + t124 * (Ifges(3,4) * t121 + Ifges(3,2) * t124) + t121 * (Ifges(3,1) * t121 + Ifges(3,4) * t124) + 0.2e1 * t8 * t38 + t31 * t181 + 0.2e1 * t9 * t29 + 0.2e1 * t14 * t30 + 0.2e1 * t13 * t33 + 0.2e1 * t12 * t34 + 0.2e1 * t2 * t35 + 0.2e1 * t7 * t36 + 0.2e1 * t4 * t37 + Ifges(2,3) + 0.2e1 * t148 * pkin(7) * mrSges(3,3) + (mrSges(4,1) * t178 - 0.2e1 * t50 * mrSges(4,3) + Ifges(4,2) * t76 + t136) * t76 + m(3) * (t148 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t101 ^ 2 + t50 ^ 2 + t183) + m(5) * (t12 ^ 2 + t13 ^ 2 + t183) + m(6) * (t14 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t2 ^ 2 + t4 ^ 2 + t9 ^ 2) + (mrSges(4,2) * t178 + mrSges(4,3) * t181 + Ifges(4,1) * t77 - 0.2e1 * Ifges(4,4) * t76 + (-t26 - t172 + t191) * t122 + (-t21 - t171 + t190) * t119) * t77; t129 + (-t121 * mrSges(3,1) - t124 * mrSges(3,2)) * pkin(7) + m(6) * (t14 * t63 + t140 * t99) + m(5) * (t100 * t48 + t137 * t99) + Ifges(3,6) * t124 + Ifges(3,5) * t121 + t100 * t31 + t72 * t37 + t54 * t29 + t63 * t30 + t71 * t35 + t134 * t99 + m(7) * (t2 * t71 + t4 * t72 + t54 * t9) + (m(4) * (t120 * t50 - t123 * t48) + (-t120 * t76 - t123 * t77) * mrSges(4,3)) * pkin(2); 0.2e1 * t100 * t83 + t54 * t179 + t63 * t180 + Ifges(3,3) + (t71 * t182 - t88) * t119 + 0.2e1 * t135 + (t72 * t182 - t189) * t122 + t132 * t99 + m(6) * (t63 ^ 2 + t165) + m(7) * (t54 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(5) * (t100 ^ 2 + t165) + m(4) * (t120 ^ 2 + t123 ^ 2) * pkin(2) ^ 2 + t143; t129 + t134 * pkin(9) + m(7) * (t2 * t91 + t4 * t92 + t62 * t9) + m(6) * (t140 * pkin(9) + t14 * t79) + m(5) * (-pkin(3) * t48 + t137 * pkin(9)) + t79 * t30 + t91 * t35 + t92 * t37 + t62 * t29 - pkin(3) * t31; -t119 * t88 + (t62 + t54) * t84 + (t100 - pkin(3)) * t83 + (t79 + t63) * t82 - t189 * t122 + t135 + m(7) * (t54 * t62 + t71 * t91 + t72 * t92) + m(5) * (-pkin(3) * t100 + t167) + m(6) * (t63 * t79 + t167) + ((t72 + t92) * t122 + (t71 + t91) * t119) * mrSges(7,1) + t143 + (pkin(9) + t99) * t184; -0.2e1 * pkin(3) * t83 + t62 * t179 + t79 * t180 + (t91 * t182 - t88) * t119 + (t92 * t182 - t189) * t122 + m(6) * (t79 ^ 2 + t150) + m(7) * (t62 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (pkin(3) ^ 2 + t150) + t132 * pkin(9) + t143; t136 + (-t36 + t37) * qJ(5) + m(7) * (qJ(5) * t4 + t118 * t2) + m(6) * (-pkin(4) * t8 - qJ(5) * t7) + (-Ifges(6,4) * t122 - Ifges(5,6) * t119) * t77 + t118 * t35 - pkin(4) * t38 + t4 * mrSges(7,2) - t7 * mrSges(6,3) + t8 * mrSges(6,2) + t12 * mrSges(5,1) - t13 * mrSges(5,2) - t2 * mrSges(7,3); m(7) * (qJ(5) * t72 + t118 * t71) + t72 * mrSges(7,2) - t71 * mrSges(7,3) + t130 * t99 + t131; m(7) * (qJ(5) * t92 + t118 * t91) - t91 * mrSges(7,3) + t92 * mrSges(7,2) + t130 * pkin(9) + t131; -0.2e1 * pkin(4) * mrSges(6,2) - 0.2e1 * t118 * mrSges(7,3) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJ(5) + m(6) * (pkin(4) ^ 2 + t125) + m(7) * (t118 ^ 2 + t125) + t186; m(6) * t8 + m(7) * t2 + t35 + t38; m(6) * t156 + m(7) * t71 + t152; m(6) * t169 + m(7) * t91 + t152; -m(6) * pkin(4) + m(7) * t118 + mrSges(6,2) - mrSges(7,3); m(6) + m(7); m(7) * t4 + t37; m(7) * t72 + t106; m(7) * t92 + t106; m(7) * qJ(5) + mrSges(7,2); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
