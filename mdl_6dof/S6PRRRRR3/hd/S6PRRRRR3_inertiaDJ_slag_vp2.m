% Calculate time derivative of joint inertia matrix for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:48:12
% EndTime: 2019-03-09 00:48:22
% DurationCPUTime: 4.78s
% Computational Cost: add. (7980->547), mult. (20102->821), div. (0->0), fcn. (18924->12), ass. (0->217)
t185 = sin(qJ(4));
t186 = sin(qJ(3));
t190 = cos(qJ(4));
t226 = qJD(4) * t190;
t191 = cos(qJ(3));
t229 = qJD(3) * t191;
t196 = t185 * t229 + t186 * t226;
t184 = sin(qJ(5));
t189 = cos(qJ(5));
t203 = t184 * t185 - t189 * t190;
t129 = t203 * t186;
t252 = -pkin(10) - pkin(9);
t164 = t252 * t185;
t165 = t252 * t190;
t117 = t184 * t164 - t189 * t165;
t212 = t190 * t229;
t230 = qJD(3) * t186;
t259 = -Ifges(5,5) * t212 - Ifges(5,3) * t230;
t160 = -pkin(3) * t191 - pkin(9) * t186 - pkin(2);
t233 = t190 * t191;
t171 = pkin(8) * t233;
t124 = t185 * t160 + t171;
t258 = qJD(4) + qJD(5);
t146 = t184 * t190 + t185 * t189;
t107 = t258 * t146;
t72 = -t107 * t186 - t203 * t229;
t73 = t129 * t258 - t146 * t229;
t257 = -Ifges(6,5) * t72 - Ifges(6,6) * t73 - Ifges(6,3) * t230;
t256 = 2 * m(5);
t255 = 2 * m(6);
t254 = 2 * m(7);
t253 = 0.2e1 * pkin(8);
t246 = Ifges(5,4) * t185;
t162 = Ifges(5,2) * t190 + t246;
t251 = -t162 / 0.2e1;
t250 = -t185 / 0.2e1;
t249 = pkin(8) * t185;
t173 = pkin(4) * t189 + pkin(5);
t188 = cos(qJ(6));
t222 = qJD(6) * t188;
t183 = sin(qJ(6));
t223 = qJD(6) * t183;
t237 = t183 * t184;
t94 = t173 * t222 + (-t184 * t223 + (t188 * t189 - t237) * qJD(5)) * pkin(4);
t248 = t94 * mrSges(7,2);
t245 = Ifges(5,4) * t190;
t244 = Ifges(5,6) * t185;
t243 = t191 * Ifges(5,6);
t161 = -mrSges(5,1) * t190 + mrSges(5,2) * t185;
t242 = -mrSges(4,1) + t161;
t144 = t190 * t160;
t234 = t186 * t190;
t101 = -pkin(10) * t234 + t144 + (-pkin(4) - t249) * t191;
t235 = t185 * t186;
t115 = -pkin(10) * t235 + t124;
t62 = t184 * t101 + t189 * t115;
t182 = cos(pkin(6));
t181 = sin(pkin(6));
t187 = sin(qJ(2));
t239 = t181 * t187;
t133 = t182 * t186 + t191 * t239;
t192 = cos(qJ(2));
t231 = qJD(2) * t192;
t216 = t181 * t231;
t111 = qJD(3) * t133 + t186 * t216;
t241 = t111 * t186;
t132 = -t182 * t191 + t186 * t239;
t112 = -qJD(3) * t132 + t191 * t216;
t240 = t112 * t191;
t81 = t132 * t111;
t238 = t181 * t192;
t236 = t184 * t188;
t158 = (pkin(3) * t186 - pkin(9) * t191) * qJD(3);
t232 = t190 * t158 + t230 * t249;
t159 = pkin(4) * t235 + t186 * pkin(8);
t228 = qJD(4) * t185;
t227 = qJD(4) * t186;
t225 = qJD(5) * t184;
t224 = qJD(5) * t189;
t128 = t146 * t186;
t84 = -t128 * t188 + t129 * t183;
t30 = qJD(6) * t84 + t183 * t73 + t188 * t72;
t85 = -t128 * t183 - t129 * t188;
t31 = -qJD(6) * t85 - t183 * t72 + t188 * t73;
t221 = -Ifges(7,5) * t30 - Ifges(7,6) * t31 - Ifges(7,3) * t230;
t220 = pkin(4) * t228;
t179 = pkin(8) * t229;
t200 = -t133 * t190 + t185 * t238;
t217 = qJD(2) * t239;
t51 = qJD(4) * t200 - t112 * t185 + t190 * t217;
t113 = -t133 * t185 - t190 * t238;
t52 = qJD(4) * t113 + t112 * t190 + t185 * t217;
t64 = t113 * t189 + t184 * t200;
t16 = qJD(5) * t64 + t184 * t51 + t189 * t52;
t65 = t113 * t184 - t189 * t200;
t17 = -qJD(5) * t65 - t184 * t52 + t189 * t51;
t32 = -t183 * t65 + t188 * t64;
t5 = qJD(6) * t32 + t16 * t188 + t17 * t183;
t33 = t183 * t64 + t188 * t65;
t6 = -qJD(6) * t33 - t16 * t183 + t17 * t188;
t219 = t6 * mrSges(7,1) - t5 * mrSges(7,2);
t122 = pkin(4) * t196 + t179;
t174 = -pkin(4) * t190 - pkin(3);
t218 = qJD(4) * t252;
t215 = t185 * t227;
t95 = -t173 * t223 + (-t184 * t222 + (-t183 * t189 - t236) * qJD(5)) * pkin(4);
t92 = t95 * mrSges(7,1);
t211 = t92 - t248;
t210 = (2 * Ifges(4,4)) + t244;
t61 = t189 * t101 - t115 * t184;
t123 = -t191 * t249 + t144;
t79 = t185 * t158 + t160 * t226 + (-t190 * t230 - t191 * t228) * pkin(8);
t209 = -t123 * qJD(4) + t79;
t116 = t189 * t164 + t165 * t184;
t208 = mrSges(5,1) * t185 + mrSges(5,2) * t190;
t207 = Ifges(5,1) * t190 - t246;
t163 = Ifges(5,1) * t185 + t245;
t206 = -Ifges(5,2) * t185 + t245;
t205 = Ifges(5,5) * t185 + Ifges(5,6) * t190;
t45 = -pkin(5) * t191 + pkin(11) * t129 + t61;
t49 = -pkin(11) * t128 + t62;
t24 = -t183 * t49 + t188 * t45;
t25 = t183 * t45 + t188 * t49;
t89 = -pkin(11) * t146 + t116;
t90 = -pkin(11) * t203 + t117;
t47 = -t183 * t90 + t188 * t89;
t48 = t183 * t89 + t188 * t90;
t156 = t185 * t218;
t157 = t190 * t218;
t75 = t189 * t156 + t184 * t157 + t164 * t224 + t165 * t225;
t43 = -pkin(11) * t107 + t75;
t106 = t258 * t203;
t76 = -qJD(5) * t117 - t156 * t184 + t189 * t157;
t44 = pkin(11) * t106 + t76;
t11 = qJD(6) * t47 + t183 * t44 + t188 * t43;
t12 = -qJD(6) * t48 - t183 * t43 + t188 * t44;
t100 = t146 * t188 - t183 * t203;
t40 = -qJD(6) * t100 + t106 * t183 - t107 * t188;
t37 = Ifges(7,6) * t40;
t99 = -t146 * t183 - t188 * t203;
t39 = qJD(6) * t99 - t106 * t188 - t107 * t183;
t38 = Ifges(7,5) * t39;
t204 = t12 * mrSges(7,1) - t11 * mrSges(7,2) + t37 + t38;
t58 = (pkin(4) * t186 - pkin(10) * t233) * qJD(3) + (-t171 + (pkin(10) * t186 - t160) * t185) * qJD(4) + t232;
t69 = -pkin(10) * t196 + t79;
t20 = -qJD(5) * t62 - t184 * t69 + t189 * t58;
t13 = pkin(5) * t230 - pkin(11) * t72 + t20;
t19 = t101 * t224 - t115 * t225 + t184 * t58 + t189 * t69;
t14 = pkin(11) * t73 + t19;
t2 = qJD(6) * t24 + t13 * t183 + t14 * t188;
t3 = -qJD(6) * t25 + t13 * t188 - t14 * t183;
t202 = t3 * mrSges(7,1) - t2 * mrSges(7,2) - t221;
t201 = t17 * mrSges(6,1) - t16 * mrSges(6,2) + t219;
t199 = t132 * t229 + t241;
t198 = (-mrSges(6,1) * t184 - mrSges(6,2) * t189) * qJD(5) * pkin(4);
t197 = t212 - t215;
t102 = Ifges(6,6) * t107;
t103 = Ifges(6,5) * t106;
t195 = t76 * mrSges(6,1) - t75 * mrSges(6,2) - t102 - t103 + t204;
t194 = t20 * mrSges(6,1) - t19 * mrSges(6,2) + t202 - t257;
t193 = -t51 * t185 + t52 * t190 + (-t113 * t190 + t185 * t200) * qJD(4);
t178 = Ifges(5,5) * t226;
t155 = -mrSges(5,1) * t191 - mrSges(5,3) * t234;
t154 = mrSges(5,2) * t191 - mrSges(5,3) * t235;
t153 = t207 * qJD(4);
t152 = t206 * qJD(4);
t151 = (mrSges(4,1) * t186 + mrSges(4,2) * t191) * qJD(3);
t150 = t208 * qJD(4);
t141 = (-mrSges(7,1) * t183 - mrSges(7,2) * t188) * qJD(6) * pkin(5);
t139 = t208 * t186;
t131 = pkin(4) * t236 + t173 * t183;
t130 = -pkin(4) * t237 + t173 * t188;
t127 = -Ifges(5,5) * t191 + t186 * t207;
t126 = t186 * t206 - t243;
t125 = pkin(5) * t203 + t174;
t121 = -mrSges(5,2) * t230 - mrSges(5,3) * t196;
t120 = mrSges(5,1) * t230 - mrSges(5,3) * t197;
t119 = -mrSges(6,1) * t191 + mrSges(6,3) * t129;
t118 = mrSges(6,2) * t191 - mrSges(6,3) * t128;
t110 = Ifges(6,1) * t146 - Ifges(6,4) * t203;
t109 = Ifges(6,4) * t146 - Ifges(6,2) * t203;
t108 = mrSges(6,1) * t203 + mrSges(6,2) * t146;
t104 = pkin(5) * t128 + t159;
t98 = mrSges(5,1) * t196 + mrSges(5,2) * t197;
t91 = pkin(5) * t107 + t220;
t88 = mrSges(6,1) * t128 - mrSges(6,2) * t129;
t87 = -t163 * t227 + (Ifges(5,5) * t186 + t191 * t207) * qJD(3);
t86 = -t162 * t227 + (Ifges(5,6) * t186 + t191 * t206) * qJD(3);
t83 = -Ifges(6,1) * t129 - Ifges(6,4) * t128 - Ifges(6,5) * t191;
t82 = -Ifges(6,4) * t129 - Ifges(6,2) * t128 - Ifges(6,6) * t191;
t80 = -qJD(4) * t124 + t232;
t78 = -mrSges(7,1) * t191 - mrSges(7,3) * t85;
t77 = mrSges(7,2) * t191 + mrSges(7,3) * t84;
t68 = -Ifges(6,1) * t106 - Ifges(6,4) * t107;
t67 = -Ifges(6,4) * t106 - Ifges(6,2) * t107;
t66 = mrSges(6,1) * t107 - mrSges(6,2) * t106;
t60 = -mrSges(6,2) * t230 + mrSges(6,3) * t73;
t59 = mrSges(6,1) * t230 - mrSges(6,3) * t72;
t57 = Ifges(7,1) * t100 + Ifges(7,4) * t99;
t56 = Ifges(7,4) * t100 + Ifges(7,2) * t99;
t55 = -mrSges(7,1) * t99 + mrSges(7,2) * t100;
t50 = -pkin(5) * t73 + t122;
t46 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t42 = Ifges(7,1) * t85 + Ifges(7,4) * t84 - Ifges(7,5) * t191;
t41 = Ifges(7,4) * t85 + Ifges(7,2) * t84 - Ifges(7,6) * t191;
t36 = -mrSges(6,1) * t73 + mrSges(6,2) * t72;
t35 = Ifges(6,1) * t72 + Ifges(6,4) * t73 + Ifges(6,5) * t230;
t34 = Ifges(6,4) * t72 + Ifges(6,2) * t73 + Ifges(6,6) * t230;
t27 = -mrSges(7,2) * t230 + mrSges(7,3) * t31;
t26 = mrSges(7,1) * t230 - mrSges(7,3) * t30;
t23 = Ifges(7,1) * t39 + Ifges(7,4) * t40;
t22 = Ifges(7,4) * t39 + Ifges(7,2) * t40;
t21 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t9 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t8 = Ifges(7,1) * t30 + Ifges(7,4) * t31 + Ifges(7,5) * t230;
t7 = Ifges(7,4) * t30 + Ifges(7,2) * t31 + Ifges(7,6) * t230;
t1 = [0.2e1 * m(7) * (t32 * t6 + t33 * t5 + t81) + 0.2e1 * m(6) * (t16 * t65 + t17 * t64 + t81) + 0.2e1 * m(5) * (t113 * t51 - t200 * t52 + t81) + 0.2e1 * m(4) * (-t181 ^ 2 * t187 * t231 + t133 * t112 + t81); t113 * t120 - t200 * t121 + t16 * t118 + t17 * t119 + t52 * t154 + t51 * t155 + t32 * t26 + t33 * t27 + t5 * t77 + t64 * t59 + t6 * t78 + t65 * t60 + (t98 + t36 + t9) * t132 + (t139 + t88 + t46) * t111 + (-t192 * t151 + (-t192 * mrSges(3,2) + (-mrSges(4,1) * t191 + mrSges(4,2) * t186 - mrSges(3,1)) * t187) * qJD(2)) * t181 + (t241 + t240 + (t132 * t191 - t133 * t186) * qJD(3)) * mrSges(4,3) + m(5) * (t113 * t80 + t123 * t51 + t124 * t52 - t200 * t79) - m(4) * pkin(2) * t217 + m(6) * (t111 * t159 + t122 * t132 + t16 * t62 + t17 * t61 + t19 * t65 + t20 * t64) + m(7) * (t104 * t111 + t132 * t50 + t2 * t33 + t24 * t6 + t25 * t5 + t3 * t32) + (m(5) * t199 / 0.2e1 + m(4) * (-t133 * t230 + t199 + t240) / 0.2e1) * t253; ((-t185 * t126 + t190 * t127 + t139 * t253 + t191 * t210) * qJD(3) + t221 + t257 + t259) * t191 + (t104 * t50 + t2 * t25 + t24 * t3) * t254 + (t122 * t159 + t19 * t62 + t20 * t61) * t255 + (t123 * t80 + t124 * t79) * t256 + (t98 * t253 - t185 * t86 + t190 * t87 + (-t190 * t126 - t185 * t127 + t191 * t205) * qJD(4) + (-Ifges(6,5) * t129 - Ifges(6,6) * t128 + Ifges(7,5) * t85 + Ifges(7,6) * t84 + (Ifges(5,5) * t190 - t210) * t186 + (pkin(8) ^ 2 * t256 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(5,3) - Ifges(6,3) - Ifges(7,3)) * t191) * qJD(3)) * t186 + 0.2e1 * t79 * t154 + 0.2e1 * t80 * t155 + 0.2e1 * t159 * t36 - 0.2e1 * pkin(2) * t151 - t128 * t34 - t129 * t35 + 0.2e1 * t19 * t118 + 0.2e1 * t20 * t119 + 0.2e1 * t122 * t88 + 0.2e1 * t123 * t120 + 0.2e1 * t124 * t121 + 0.2e1 * t104 * t9 + t84 * t7 + t85 * t8 + 0.2e1 * t2 * t77 + 0.2e1 * t3 * t78 + t73 * t82 + t72 * t83 + 0.2e1 * t61 * t59 + 0.2e1 * t62 * t60 + 0.2e1 * t50 * t46 + t31 * t41 + t30 * t42 + 0.2e1 * t24 * t26 + 0.2e1 * t25 * t27; -t112 * mrSges(4,2) + (t150 + t21 + t66) * t132 + (t108 + t55 + t242) * t111 + m(7) * (t11 * t33 + t111 * t125 + t12 * t32 + t132 * t91 + t47 * t6 + t48 * t5) + m(6) * (t111 * t174 + t116 * t17 + t117 * t16 + t132 * t220 + t64 * t76 + t65 * t75) + (-t100 * t6 - t32 * t39 + t33 * t40 + t5 * t99) * mrSges(7,3) + (t106 * t64 - t107 * t65 - t146 * t17 - t16 * t203) * mrSges(6,3) + t193 * mrSges(5,3) + (-pkin(3) * t111 + pkin(9) * t193) * m(5); (t106 * t61 - t107 * t62 - t146 * t20 - t19 * t203) * mrSges(6,3) + (t190 * t153 / 0.2e1 + t152 * t250 - Ifges(4,6) * qJD(3) + (t163 * t250 + t190 * t251) * qJD(4) + (qJD(3) * mrSges(4,2) + t150) * pkin(8) + (Ifges(6,5) * t146 + Ifges(7,5) * t100 - Ifges(6,6) * t203 + Ifges(7,6) * t99 + t205) * qJD(3) / 0.2e1) * t186 - t203 * t34 / 0.2e1 + m(7) * (t104 * t91 + t11 * t25 + t12 * t24 + t125 * t50 + t2 * t48 + t3 * t47) + (t103 / 0.2e1 + t102 / 0.2e1 - t38 / 0.2e1 - t37 / 0.2e1 - t178 / 0.2e1 + (pkin(8) * t242 + Ifges(4,5)) * qJD(3)) * t191 + m(5) * (-pkin(3) * t179 + (-t124 * t228 - t80 * t185) * pkin(9)) + t174 * t36 + t159 * t66 + t146 * t35 / 0.2e1 - t128 * t67 / 0.2e1 - t129 * t68 / 0.2e1 + t116 * t59 + t117 * t60 + t75 * t118 + t76 * t119 + t122 * t108 + t125 * t9 + t104 * t21 - t106 * t83 / 0.2e1 - t107 * t82 / 0.2e1 + t73 * t109 / 0.2e1 + t72 * t110 / 0.2e1 + (-t100 * t3 + t2 * t99 - t24 * t39 + t25 * t40) * mrSges(7,3) - pkin(3) * t98 + t99 * t7 / 0.2e1 + t100 * t8 / 0.2e1 + t84 * t22 / 0.2e1 + t85 * t23 / 0.2e1 + t91 * t46 + t11 * t77 + t12 * t78 + t50 * t55 + t31 * t56 / 0.2e1 + t30 * t57 / 0.2e1 + t47 * t26 + t48 * t27 + t40 * t41 / 0.2e1 + t39 * t42 / 0.2e1 + m(6) * (t116 * t20 + t117 * t19 + t122 * t174 + t159 * t220 + t61 * t76 + t62 * t75) + (qJD(4) * t127 / 0.2e1 + t163 * t229 / 0.2e1 + t86 / 0.2e1 + t209 * mrSges(5,3) + (m(5) * t209 - qJD(4) * t155 + t121) * pkin(9)) * t190 + (-pkin(9) * t120 - t80 * mrSges(5,3) + t229 * t251 + t87 / 0.2e1 + (-t126 / 0.2e1 + t243 / 0.2e1 + pkin(4) * t88 - t124 * mrSges(5,3) - pkin(9) * t154) * qJD(4)) * t185; -0.2e1 * pkin(3) * t150 + t100 * t23 - t106 * t110 - t107 * t109 + 0.2e1 * t125 * t21 - t203 * t67 + t146 * t68 + t190 * t152 + t185 * t153 + 0.2e1 * t174 * t66 + t99 * t22 + t39 * t57 + t40 * t56 + 0.2e1 * t91 * t55 + (t190 * t163 + (0.2e1 * pkin(4) * t108 - t162) * t185) * qJD(4) + (t116 * t76 + t117 * t75 + t174 * t220) * t255 + (t11 * t48 + t12 * t47 + t125 * t91) * t254 + 0.2e1 * (-t100 * t12 + t11 * t99 - t39 * t47 + t40 * t48) * mrSges(7,3) + 0.2e1 * (t106 * t116 - t107 * t117 - t146 * t76 - t203 * t75) * mrSges(6,3); m(7) * (t130 * t6 + t131 * t5 + t32 * t95 + t33 * t94) - t52 * mrSges(5,2) + t51 * mrSges(5,1) + m(6) * (t16 * t184 + t17 * t189 + (-t184 * t64 + t189 * t65) * qJD(5)) * pkin(4) + t201; t194 + m(7) * (t130 * t3 + t131 * t2 + t24 * t95 + t25 * t94) + (m(6) * (t184 * t19 + t189 * t20 + t224 * t62 - t225 * t61) + t189 * t59 + t184 * t60 + t118 * t224 - t119 * t225) * pkin(4) - t196 * Ifges(5,6) + t130 * t26 + t131 * t27 + t94 * t77 + t95 * t78 - t79 * mrSges(5,2) + t80 * mrSges(5,1) - Ifges(5,5) * t215 - t259; m(7) * (t11 * t131 + t12 * t130 + t47 * t95 + t48 * t94) + t178 + (pkin(9) * t161 - t244) * qJD(4) + (-t100 * t95 - t130 * t39 + t131 * t40 + t94 * t99) * mrSges(7,3) + (m(6) * (t184 * t75 + t189 * t76 + (-t116 * t184 + t117 * t189) * qJD(5)) + (t189 * t106 - t184 * t107 + (t146 * t184 - t189 * t203) * qJD(5)) * mrSges(6,3)) * pkin(4) + t195; (t130 * t95 + t131 * t94) * t254 - 0.2e1 * t248 + 0.2e1 * t92 + 0.2e1 * t198; m(7) * (t183 * t5 + t188 * t6 + (-t183 * t32 + t188 * t33) * qJD(6)) * pkin(5) + t201; (m(7) * (t183 * t2 + t188 * t3 + t222 * t25 - t223 * t24) + t77 * t222 + t183 * t27 - t78 * t223 + t188 * t26) * pkin(5) + t194; (m(7) * (t11 * t183 + t12 * t188 + (-t183 * t47 + t188 * t48) * qJD(6)) + (t183 * t40 - t188 * t39 + (t100 * t183 + t188 * t99) * qJD(6)) * mrSges(7,3)) * pkin(5) + t195; t198 + (m(7) * (-t130 * t223 + t131 * t222 + t183 * t94 + t188 * t95) - mrSges(7,2) * t222 - mrSges(7,1) * t223) * pkin(5) + t211; 0.2e1 * t141; t219; t202; t204; t211; t141; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
