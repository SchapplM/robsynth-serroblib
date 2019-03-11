% Calculate time derivative of joint inertia matrix for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:21:16
% EndTime: 2019-03-09 15:21:22
% DurationCPUTime: 3.34s
% Computational Cost: add. (9064->370), mult. (20055->562), div. (0->0), fcn. (19729->10), ass. (0->168)
t163 = cos(pkin(11));
t160 = t163 ^ 2;
t161 = sin(pkin(11));
t193 = t161 ^ 2 + t160;
t248 = t193 * mrSges(6,3);
t247 = -mrSges(6,1) * t163 + mrSges(6,2) * t161 - mrSges(5,1);
t165 = sin(qJ(3));
t166 = sin(qJ(2));
t168 = cos(qJ(3));
t169 = cos(qJ(2));
t141 = t165 * t169 + t168 * t166;
t162 = sin(pkin(10));
t205 = cos(pkin(10));
t182 = t205 * t165;
t221 = pkin(2) * qJD(3);
t124 = (t162 * t168 + t182) * t221;
t246 = 0.2e1 * t124;
t230 = -pkin(8) - pkin(7);
t148 = t230 * t166;
t142 = t165 * t148;
t149 = t230 * t169;
t118 = -t168 * t149 + t142;
t164 = sin(qJ(6));
t167 = cos(qJ(6));
t174 = t161 * t164 - t163 * t167;
t131 = t174 * qJD(6);
t244 = qJD(2) + qJD(3);
t139 = t161 * t167 + t163 * t164;
t110 = mrSges(7,1) * t174 + mrSges(7,2) * t139;
t243 = t110 + t247;
t229 = pkin(3) * t162;
t153 = qJ(5) + t229;
t133 = (-pkin(9) - t153) * t161;
t158 = t163 * pkin(9);
t198 = t153 * t163;
t134 = t158 + t198;
t100 = t133 * t164 + t134 * t167;
t132 = t139 * qJD(6);
t99 = t133 * t167 - t134 * t164;
t91 = -qJD(5) * t174 + qJD(6) * t99;
t92 = -qJD(5) * t139 - qJD(6) * t100;
t242 = -t100 * t132 + t99 * t131 - t92 * t139 - t174 * t91;
t191 = qJD(3) * t168;
t186 = pkin(2) * t191;
t192 = qJD(3) * t165;
t187 = pkin(2) * t192;
t125 = -t162 * t187 + t205 * t186;
t121 = qJD(5) + t125;
t155 = pkin(2) * t168 + pkin(3);
t129 = pkin(2) * t182 + t162 * t155;
t122 = qJ(5) + t129;
t115 = (-pkin(9) - t122) * t161;
t116 = t122 * t163 + t158;
t86 = t115 * t167 - t116 * t164;
t58 = qJD(6) * t86 - t121 * t174;
t87 = t115 * t164 + t116 * t167;
t59 = -qJD(6) * t87 - t121 * t139;
t241 = t86 * t131 - t87 * t132 - t59 * t139 - t174 * t58;
t102 = t132 * mrSges(7,1) - t131 * mrSges(7,2);
t103 = -Ifges(7,4) * t131 - Ifges(7,2) * t132;
t104 = -Ifges(7,1) * t131 - Ifges(7,4) * t132;
t195 = t168 * t169;
t140 = -t165 * t166 + t195;
t106 = -t140 * t205 + t141 * t162;
t107 = t162 * t140 + t141 * t205;
t201 = t107 * t163;
t156 = -pkin(2) * t169 - pkin(1);
t120 = -t140 * pkin(3) + t156;
t66 = t106 * pkin(4) - t107 * qJ(5) + t120;
t117 = t168 * t148 + t149 * t165;
t172 = -qJ(4) * t141 + t117;
t98 = qJ(4) * t140 + t118;
t68 = t162 * t172 + t205 * t98;
t46 = -t161 * t68 + t163 * t66;
t23 = pkin(5) * t106 - pkin(9) * t201 + t46;
t202 = t107 * t161;
t47 = t161 * t66 + t163 * t68;
t30 = -pkin(9) * t202 + t47;
t11 = -t164 * t30 + t167 * t23;
t111 = Ifges(7,4) * t139 - Ifges(7,2) * t174;
t112 = Ifges(7,1) * t139 - Ifges(7,4) * t174;
t113 = t244 * t140;
t114 = t244 * t141;
t12 = t164 * t23 + t167 * t30;
t194 = -Ifges(7,5) * t131 - Ifges(7,6) * t132;
t171 = (t195 * t230 - t142) * qJD(2);
t173 = -t113 * qJ(4) - t141 * qJD(4);
t89 = t141 * qJD(2) * t230 + t148 * t191 + t149 * t192;
t61 = -qJ(4) * t114 + qJD(4) * t140 + t89;
t36 = t205 * t61 + (-t148 * t192 + t149 * t191 + t171 + t173) * t162;
t101 = qJD(2) * t166 * pkin(2) + pkin(3) * t114;
t84 = t113 * t162 + t114 * t205;
t85 = t113 * t205 - t162 * t114;
t43 = pkin(4) * t84 - qJ(5) * t85 - qJD(5) * t107 + t101;
t15 = -t161 * t36 + t163 * t43;
t213 = t163 * t85;
t4 = pkin(5) * t84 - pkin(9) * t213 + t15;
t16 = t161 * t43 + t163 * t36;
t214 = t161 * t85;
t9 = -pkin(9) * t214 + t16;
t2 = qJD(6) * t11 + t164 * t4 + t167 * t9;
t90 = -qJD(3) * t118 + t171;
t35 = t162 * t61 - t205 * (t173 + t90);
t20 = pkin(5) * t214 + t35;
t215 = t16 * t163;
t223 = Ifges(6,4) * t163;
t224 = Ifges(6,4) * t161;
t27 = -t107 * t132 - t174 * t85;
t28 = t107 * t131 - t139 * t85;
t3 = -qJD(6) * t12 - t164 * t9 + t167 * t4;
t222 = Ifges(6,2) * t161;
t41 = Ifges(6,6) * t84 + (-t222 + t223) * t85;
t42 = Ifges(6,5) * t84 + (Ifges(6,1) * t163 - t224) * t85;
t72 = t139 * t107;
t73 = t174 * t107;
t44 = -Ifges(7,4) * t73 - Ifges(7,2) * t72 + Ifges(7,6) * t106;
t45 = -Ifges(7,1) * t73 - Ifges(7,4) * t72 + Ifges(7,5) * t106;
t67 = t162 * t98 - t205 * t172;
t53 = pkin(5) * t202 + t67;
t7 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t84;
t8 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t84;
t240 = t247 * t35 + mrSges(6,3) * t215 + t163 * t41 / 0.2e1 + t161 * t42 / 0.2e1 + t139 * t8 / 0.2e1 - t132 * t44 / 0.2e1 - t131 * t45 / 0.2e1 - Ifges(4,6) * t114 + t20 * t110 + t28 * t111 / 0.2e1 + t27 * t112 / 0.2e1 + Ifges(4,5) * t113 + t53 * t102 - t72 * t103 / 0.2e1 - t73 * t104 / 0.2e1 - t89 * mrSges(4,2) + t90 * mrSges(4,1) - Ifges(5,6) * t84 + Ifges(5,5) * t85 - t36 * mrSges(5,2) + t106 * t194 / 0.2e1 + (Ifges(6,1) * t161 + t223) * t213 / 0.2e1 - (Ifges(6,2) * t163 + t224) * t214 / 0.2e1 - t174 * t7 / 0.2e1 + (Ifges(6,5) * t161 + Ifges(7,5) * t139 + Ifges(6,6) * t163 - Ifges(7,6) * t174) * t84 / 0.2e1 + (t11 * t131 - t12 * t132 - t3 * t139 - t174 * t2) * mrSges(7,3);
t239 = 2 * m(5);
t238 = 2 * m(6);
t237 = 2 * m(7);
t236 = 0.2e1 * t35;
t235 = 0.2e1 * t67;
t234 = 0.2e1 * t101;
t233 = 0.2e1 * t156;
t232 = m(5) * pkin(3);
t228 = t163 * pkin(5);
t226 = t35 * t67;
t225 = t84 * mrSges(5,3);
t50 = mrSges(6,1) * t214 + mrSges(6,2) * t213;
t218 = t124 * t67;
t217 = t125 * mrSges(5,2);
t216 = t15 * t161;
t128 = -t162 * t165 * pkin(2) + t155 * t205;
t123 = -pkin(4) - t128;
t119 = t123 - t228;
t200 = t119 * t102;
t184 = t205 * pkin(3);
t154 = -t184 - pkin(4);
t145 = t154 - t228;
t199 = t145 * t102;
t190 = 0.2e1 * mrSges(7,3);
t189 = 0.2e1 * t169;
t188 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t84;
t183 = t84 * mrSges(5,1) + t85 * mrSges(5,2);
t13 = -t28 * mrSges(7,1) + t27 * mrSges(7,2);
t181 = -t103 * t174 + t139 * t104 - t132 * t111 - t131 * t112;
t180 = t193 * t122;
t179 = t193 * t153;
t178 = Ifges(6,5) * t163 - Ifges(6,6) * t161;
t177 = t215 - t216;
t176 = -t161 * t46 + t163 * t47;
t175 = 0.2e1 * t248;
t76 = mrSges(6,1) * t106 - mrSges(6,3) * t201;
t75 = -mrSges(6,2) * t106 - mrSges(6,3) * t202;
t74 = (mrSges(6,1) * t161 + mrSges(6,2) * t163) * t107;
t57 = mrSges(7,1) * t106 + mrSges(7,3) * t73;
t56 = -mrSges(7,2) * t106 - mrSges(7,3) * t72;
t52 = mrSges(6,1) * t84 - mrSges(6,3) * t213;
t51 = -mrSges(6,2) * t84 - mrSges(6,3) * t214;
t48 = mrSges(7,1) * t72 - mrSges(7,2) * t73;
t19 = -mrSges(7,2) * t84 + mrSges(7,3) * t28;
t18 = mrSges(7,1) * t84 - mrSges(7,3) * t27;
t1 = [t84 * (-Ifges(7,5) * t73 - Ifges(7,6) * t72) + t50 * t235 + t74 * t236 + (t11 * t3 + t12 * t2 + t20 * t53) * t237 + (t15 * t46 + t16 * t47 + t226) * t238 + (t101 * t120 + t36 * t68 + t226) * t239 + 0.2e1 * m(4) * (t117 * t90 + t118 * t89) + (mrSges(5,3) * t235 + (Ifges(6,1) * t160 + (2 * Ifges(5,1)) + (t222 - 0.2e1 * t223) * t161) * t107) * t85 + 0.2e1 * t120 * t183 + 0.2e1 * t16 * t75 + 0.2e1 * t15 * t76 - t72 * t7 - t73 * t8 - 0.2e1 * t140 * Ifges(4,2) * t114 + 0.2e1 * t2 * t56 + 0.2e1 * t3 * t57 + 0.2e1 * t47 * t51 + 0.2e1 * t46 * t52 + 0.2e1 * t53 * t13 + t28 * t44 + t27 * t45 + 0.2e1 * t20 * t48 + (0.2e1 * (-Ifges(5,4) + t178) * t85 + mrSges(5,1) * t234 - 0.2e1 * t36 * mrSges(5,3) + ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t84 + t188) * t106 + 0.2e1 * t12 * t19 + 0.2e1 * t11 * t18 + 0.2e1 * t113 * t141 * Ifges(4,1) - 0.2e1 * t68 * t225 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t169) * t189 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t140 + mrSges(4,2) * t141) + m(4) * pkin(2) * t233 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t166 + (-Ifges(3,2) + Ifges(3,1)) * t189) * t166) * qJD(2) + (mrSges(5,2) * t234 + mrSges(5,3) * t236 - t161 * t41 + t163 * t42 + (-0.2e1 * Ifges(5,4) + t178) * t84) * t107 + 0.2e1 * (t113 * t140 - t114 * t141) * Ifges(4,4) + 0.2e1 * (-t113 * t117 - t114 * t118 + t140 * t89 - t141 * t90) * mrSges(4,3) + (t114 * mrSges(4,1) + t113 * mrSges(4,2)) * t233; (-t125 * t106 + t124 * t107 - t128 * t85 - t129 * t84) * mrSges(5,3) + m(7) * (t11 * t59 + t119 * t20 + t12 * t58 + t124 * t53 + t2 * t87 + t3 * t86) + m(5) * (t125 * t68 - t128 * t35 + t129 * t36 + t218) + (m(4) * (t165 * t89 + t168 * t90 + (-t117 * t165 + t118 * t168) * qJD(3)) + (-t168 * t113 - t165 * t114 + (t140 * t168 + t141 * t165) * qJD(3)) * mrSges(4,3)) * pkin(2) + m(6) * (t121 * t176 + t122 * t177 + t123 * t35 + t218) + (Ifges(3,5) * t169 - Ifges(3,6) * t166 + (-mrSges(3,1) * t169 + mrSges(3,2) * t166) * pkin(7)) * qJD(2) + (t121 * t75 + t122 * t51) * t163 + (-t15 * mrSges(6,3) - t121 * t76 - t122 * t52) * t161 + (t74 + t48) * t124 + t123 * t50 + t119 * t13 + t86 * t18 + t87 * t19 + t58 * t56 + t59 * t57 + t240; -0.2e1 * t217 + 0.2e1 * t200 + t121 * t175 + (t119 * t124 + t58 * t87 + t59 * t86) * t237 + (-t124 * t128 + t125 * t129) * t239 + (t121 * t180 + t123 * t124) * t238 + t241 * t190 + t181 + 0.2e1 * (-mrSges(4,1) * t165 - mrSges(4,2) * t168) * t221 + t243 * t246; -t225 * t229 - t85 * mrSges(5,3) * t184 + (t162 * t36 - t205 * t35) * t232 + t51 * t198 - t161 * t153 * t52 - mrSges(6,3) * t216 + m(6) * (t153 * t177 + t154 * t35) + t154 * t50 + t145 * t13 + t99 * t18 + t100 * t19 + t91 * t56 + t92 * t57 + m(7) * (t100 * t2 + t11 * t92 + t12 * t91 + t145 * t20 + t3 * t99) + (m(6) * t176 - t161 * t76 + t163 * t75) * qJD(5) + t240; -mrSges(4,1) * t187 - mrSges(4,2) * t186 + t181 + t125 * t162 * t232 + m(7) * (t100 * t58 + t59 * t99 + t86 * t92 + t87 * t91) + m(6) * (qJD(5) * t180 + t121 * t179) + t199 - t217 + t200 + (m(6) * t154 + m(7) * t145 - t205 * t232 + t243) * t124 + (t241 + t242) * mrSges(7,3) + (qJD(5) + t121) * t248; 0.2e1 * t199 + (t100 * t91 + t92 * t99) * t237 + (t179 * t238 + t175) * qJD(5) + t242 * t190 + t181; -t131 * t56 - t132 * t57 - t174 * t18 + t139 * t19 + t161 * t51 + t163 * t52 + m(7) * (-t11 * t132 - t12 * t131 + t139 * t2 - t174 * t3) + m(6) * (t15 * t163 + t16 * t161) + m(5) * t101 + t183; m(7) * (-t131 * t87 - t132 * t86 + t139 * t58 - t174 * t59); m(7) * (-t100 * t131 - t132 * t99 + t139 * t91 - t174 * t92); (-t131 * t139 + t132 * t174) * t237; m(6) * t35 + m(7) * t20 + t13 + t50; (m(7) / 0.2e1 + m(6) / 0.2e1) * t246 + t102; t102; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t188; mrSges(7,1) * t59 - mrSges(7,2) * t58 + t194; mrSges(7,1) * t92 - mrSges(7,2) * t91 + t194; -t102; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
