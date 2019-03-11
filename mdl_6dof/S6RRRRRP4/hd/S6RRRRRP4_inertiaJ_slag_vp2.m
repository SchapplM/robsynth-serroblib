% Calculate joint inertia matrix for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:12:08
% EndTime: 2019-03-10 01:12:12
% DurationCPUTime: 1.74s
% Computational Cost: add. (2720->317), mult. (5128->430), div. (0->0), fcn. (5495->8), ass. (0->122)
t219 = Ifges(6,1) + Ifges(7,1);
t218 = -Ifges(6,4) + Ifges(7,5);
t217 = Ifges(7,4) + Ifges(6,5);
t216 = -Ifges(6,6) + Ifges(7,6);
t215 = -mrSges(6,1) - mrSges(7,1);
t214 = -mrSges(6,2) + mrSges(7,3);
t138 = sin(qJ(4));
t142 = cos(qJ(4));
t137 = sin(qJ(5));
t141 = cos(qJ(5));
t109 = t137 * t138 - t141 * t142;
t111 = t137 * t142 + t138 * t141;
t164 = t216 * t109 + t217 * t111;
t213 = Ifges(5,5) * t138 + Ifges(5,6) * t142 + t164;
t209 = (mrSges(6,3) + mrSges(7,2));
t212 = 2 * t209;
t189 = Ifges(7,2) + Ifges(6,3);
t139 = sin(qJ(3));
t140 = sin(qJ(2));
t143 = cos(qJ(3));
t144 = cos(qJ(2));
t110 = t139 * t140 - t143 * t144;
t112 = t139 * t144 + t140 * t143;
t55 = t111 * t112;
t56 = t109 * t112;
t211 = t217 * t110 + t218 * t55 - t219 * t56;
t210 = t218 * t109 + t219 * t111;
t129 = -pkin(2) * t144 - pkin(1);
t68 = pkin(3) * t110 - pkin(9) * t112 + t129;
t198 = -pkin(8) - pkin(7);
t120 = t198 * t144;
t169 = t198 * t140;
t89 = -t143 * t120 + t139 * t169;
t31 = -t138 * t89 + t142 * t68;
t32 = t138 * t68 + t142 * t89;
t162 = -t138 * t31 + t142 * t32;
t177 = t112 * t142;
t208 = Ifges(5,5) * t177 + Ifges(5,3) * t110;
t125 = pkin(2) * t139 + pkin(9);
t132 = t142 * pkin(10);
t106 = t125 * t142 + t132;
t166 = (-pkin(10) - t125) * t138;
t63 = t106 * t137 - t141 * t166;
t65 = t141 * t106 + t137 * t166;
t207 = t214 * t65 + t215 * t63;
t119 = pkin(9) * t142 + t132;
t168 = (-pkin(10) - pkin(9)) * t138;
t84 = t119 * t137 - t141 * t168;
t88 = t141 * t119 + t137 * t168;
t206 = t214 * t88 + t215 * t84;
t86 = -t120 * t139 - t143 * t169;
t205 = t86 ^ 2;
t204 = 2 * mrSges(7,3);
t71 = mrSges(7,1) * t109 - mrSges(7,3) * t111;
t203 = 0.2e1 * t71;
t72 = mrSges(6,1) * t109 + mrSges(6,2) * t111;
t202 = 0.2e1 * t72;
t201 = 0.2e1 * t86;
t200 = 0.2e1 * t129;
t195 = pkin(2) * t143;
t194 = pkin(4) * t137;
t16 = pkin(4) * t110 - pkin(10) * t177 + t31;
t178 = t112 * t138;
t22 = -pkin(10) * t178 + t32;
t8 = t137 * t16 + t141 * t22;
t35 = mrSges(6,1) * t110 + mrSges(6,3) * t56;
t36 = -t110 * mrSges(7,1) - t56 * mrSges(7,2);
t188 = t36 - t35;
t34 = -mrSges(6,2) * t110 - mrSges(6,3) * t55;
t37 = -mrSges(7,2) * t55 + mrSges(7,3) * t110;
t187 = t37 + t34;
t186 = mrSges(6,3) * t109;
t185 = Ifges(5,4) * t138;
t184 = Ifges(5,4) * t142;
t183 = t109 * mrSges(7,2);
t96 = t111 * mrSges(7,2);
t182 = t111 * mrSges(6,3);
t174 = t138 ^ 2 + t142 ^ 2;
t173 = t140 ^ 2 + t144 ^ 2;
t172 = t63 ^ 2 + t65 ^ 2;
t171 = t84 ^ 2 + t88 ^ 2;
t170 = t141 * t182;
t128 = -pkin(4) * t142 - pkin(3);
t167 = t63 * t84 + t88 * t65;
t165 = t174 * t125;
t48 = pkin(4) * t178 + t86;
t163 = mrSges(5,1) * t138 + mrSges(5,2) * t142;
t7 = -t137 * t22 + t141 * t16;
t69 = -mrSges(5,2) * t110 - mrSges(5,3) * t178;
t70 = mrSges(5,1) * t110 - mrSges(5,3) * t177;
t161 = -t138 * t70 + t142 * t69;
t160 = 0.2e1 * mrSges(5,3) * t174;
t159 = t189 * t110 + t216 * t55 - t217 * t56;
t158 = (t143 * mrSges(4,1) - t139 * mrSges(4,2)) * pkin(2);
t157 = (mrSges(6,1) * t141 - mrSges(6,2) * t137) * pkin(4);
t117 = Ifges(5,2) * t142 + t185;
t118 = Ifges(5,1) * t138 + t184;
t73 = Ifges(7,5) * t111 + Ifges(7,3) * t109;
t74 = Ifges(6,4) * t111 - Ifges(6,2) * t109;
t156 = t142 * t117 + t138 * t118 + Ifges(4,3) + t210 * t111 + (t73 - t74) * t109;
t67 = pkin(5) * t109 - qJ(6) * t111 + t128;
t152 = (-pkin(5) * t111 - qJ(6) * t109) * mrSges(7,2) + t164;
t122 = qJ(6) + t194;
t126 = -pkin(4) * t141 - pkin(5);
t151 = -t122 * t183 + t126 * t96 - t186 * t194 + t213;
t3 = qJ(6) * t110 + t8;
t5 = -pkin(5) * t110 - t7;
t150 = t7 * mrSges(6,1) - t5 * mrSges(7,1) - t8 * mrSges(6,2) + t3 * mrSges(7,3) + t159;
t10 = pkin(5) * t55 + qJ(6) * t56 + t48;
t116 = -mrSges(5,1) * t142 + mrSges(5,2) * t138;
t17 = -Ifges(7,5) * t56 + Ifges(7,6) * t110 + Ifges(7,3) * t55;
t18 = -Ifges(6,4) * t56 - Ifges(6,2) * t55 + Ifges(6,6) * t110;
t40 = Ifges(5,6) * t110 + (-Ifges(5,2) * t138 + t184) * t112;
t41 = Ifges(5,5) * t110 + (Ifges(5,1) * t142 - t185) * t112;
t149 = t142 * t40 / 0.2e1 + t138 * t41 / 0.2e1 + Ifges(4,5) * t112 - t89 * mrSges(4,2) + t10 * t71 + t48 * t72 - t7 * t182 - t3 * t183 + t118 * t177 / 0.2e1 - t117 * t178 / 0.2e1 + t5 * t96 - t8 * t186 + (t116 - mrSges(4,1)) * t86 + (t73 / 0.2e1 - t74 / 0.2e1) * t55 - t210 * t56 / 0.2e1 + t211 * t111 / 0.2e1 + (-t18 / 0.2e1 + t17 / 0.2e1) * t109 + t162 * mrSges(5,3) + (-Ifges(4,6) + t213 / 0.2e1) * t110;
t127 = -pkin(3) - t195;
t115 = t128 - t195;
t66 = t163 * t112;
t54 = t67 - t195;
t24 = mrSges(6,1) * t55 - mrSges(6,2) * t56;
t23 = mrSges(7,1) * t55 + mrSges(7,3) * t56;
t1 = [t140 * (Ifges(3,1) * t140 + Ifges(3,4) * t144) + t144 * (Ifges(3,4) * t140 + Ifges(3,2) * t144) - 0.2e1 * pkin(1) * (-mrSges(3,1) * t144 + mrSges(3,2) * t140) + t66 * t201 + 0.2e1 * t32 * t69 + 0.2e1 * t31 * t70 + 0.2e1 * t3 * t37 + 0.2e1 * t48 * t24 + 0.2e1 * t8 * t34 + 0.2e1 * t7 * t35 + 0.2e1 * t5 * t36 + 0.2e1 * t10 * t23 + Ifges(2,3) - t211 * t56 + (t17 - t18) * t55 + 0.2e1 * t173 * pkin(7) * mrSges(3,3) + (mrSges(4,2) * t200 + mrSges(4,3) * t201 + Ifges(4,1) * t112 - t138 * t40 + t142 * t41) * t112 + (mrSges(4,1) * t200 - 0.2e1 * t89 * mrSges(4,3) + Ifges(4,2) * t110 + (-Ifges(5,6) * t138 - (2 * Ifges(4,4))) * t112 + t159 + t208) * t110 + m(3) * (pkin(7) ^ 2 * t173 + pkin(1) ^ 2) + m(4) * (t129 ^ 2 + t89 ^ 2 + t205) + m(5) * (t31 ^ 2 + t32 ^ 2 + t205) + m(6) * (t48 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t10 ^ 2 + t3 ^ 2 + t5 ^ 2); (-t140 * mrSges(3,1) - t144 * mrSges(3,2)) * pkin(7) + t149 + Ifges(3,5) * t140 + Ifges(3,6) * t144 + t127 * t66 + t115 * t24 + t54 * t23 + m(5) * (t125 * t162 + t127 * t86) + t187 * t65 + t188 * t63 + t161 * t125 + m(6) * (t115 * t48 - t63 * t7 + t65 * t8) + m(7) * (t10 * t54 + t3 * t65 + t5 * t63) + (m(4) * (t139 * t89 - t143 * t86) + (-t110 * t139 - t112 * t143) * mrSges(4,3)) * pkin(2); t156 + t125 * t160 + 0.2e1 * t127 * t116 + t115 * t202 + t54 * t203 + m(4) * (t139 ^ 2 + t143 ^ 2) * pkin(2) ^ 2 + m(5) * (t125 ^ 2 * t174 + t127 ^ 2) + 0.2e1 * t158 + Ifges(3,3) + m(7) * (t54 ^ 2 + t172) + m(6) * (t115 ^ 2 + t172) + (-t65 * t109 + t63 * t111) * t212; t187 * t88 + t188 * t84 + t161 * pkin(9) + m(6) * (t128 * t48 - t7 * t84 + t8 * t88) + m(7) * (t10 * t67 + t3 * t88 + t5 * t84) + t149 + t128 * t24 - pkin(3) * t66 + t67 * t23 + m(5) * (-pkin(3) * t86 + pkin(9) * t162); t156 + m(6) * (t115 * t128 + t167) + m(7) * (t54 * t67 + t167) + (t128 + t115) * t72 + (t67 + t54) * t71 + m(5) * (-pkin(3) * t127 + pkin(9) * t165) + (pkin(9) * t174 + t165) * mrSges(5,3) + t158 + (t127 - pkin(3)) * t116 + t209 * ((t63 + t84) * t111 + (-t65 - t88) * t109); -0.2e1 * pkin(3) * t116 + t128 * t202 + t67 * t203 + pkin(9) * t160 + m(7) * (t67 ^ 2 + t171) + m(6) * (t128 ^ 2 + t171) + m(5) * (pkin(9) ^ 2 * t174 + pkin(3) ^ 2) + t156 + (-t88 * t109 + t84 * t111) * t212; t150 + t126 * t36 + t122 * t37 + t31 * mrSges(5,1) - t32 * mrSges(5,2) + (m(6) * (t137 * t8 + t141 * t7) + t137 * t34 + t141 * t35) * pkin(4) - Ifges(5,6) * t178 + m(7) * (t122 * t3 + t126 * t5) + t208; m(7) * (t122 * t65 + t126 * t63) - t163 * t125 + t151 + (m(6) * (t137 * t65 - t141 * t63) - t170) * pkin(4) + t207; m(7) * (t122 * t88 + t126 * t84) - t163 * pkin(9) + t151 + (m(6) * (t137 * t88 - t141 * t84) - t170) * pkin(4) + t206; -0.2e1 * t126 * mrSges(7,1) + t122 * t204 + Ifges(5,3) + 0.2e1 * t157 + m(7) * (t122 ^ 2 + t126 ^ 2) + m(6) * (t137 ^ 2 + t141 ^ 2) * pkin(4) ^ 2 + t189; -pkin(5) * t36 + m(7) * (-pkin(5) * t5 + qJ(6) * t3) + qJ(6) * t37 + t150; m(7) * (-pkin(5) * t63 + qJ(6) * t65) + t152 + t207; m(7) * (-pkin(5) * t84 + qJ(6) * t88) + t152 + t206; m(7) * (-pkin(5) * t126 + qJ(6) * t122) + t157 + (qJ(6) + t122) * mrSges(7,3) + (pkin(5) - t126) * mrSges(7,1) + t189; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t204 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t189; m(7) * t5 + t36; m(7) * t63 + t96; m(7) * t84 + t96; m(7) * t126 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
