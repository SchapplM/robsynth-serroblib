% Calculate time derivative of joint inertia matrix for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:36:05
% EndTime: 2019-03-08 22:36:14
% DurationCPUTime: 4.06s
% Computational Cost: add. (4110->522), mult. (11538->775), div. (0->0), fcn. (10778->12), ass. (0->223)
t237 = -mrSges(4,1) + mrSges(5,2);
t269 = mrSges(5,1) + mrSges(4,3);
t152 = sin(qJ(6));
t156 = cos(qJ(6));
t157 = cos(qJ(5));
t197 = qJD(6) * t157;
t153 = sin(qJ(5));
t202 = qJD(5) * t153;
t163 = t152 * t202 - t156 * t197;
t268 = t152 / 0.2e1;
t245 = t156 / 0.2e1;
t125 = t153 * mrSges(6,1) + t157 * mrSges(6,2);
t267 = mrSges(5,3) + t125;
t148 = sin(pkin(7));
t265 = 0.2e1 * (Ifges(5,6) + Ifges(4,4)) * t148;
t208 = t152 ^ 2 + t156 ^ 2;
t154 = sin(qJ(3));
t218 = t148 * t154;
t138 = pkin(9) * t218;
t150 = cos(pkin(7));
t158 = cos(qJ(3));
t242 = pkin(2) * t158;
t190 = -pkin(3) - t242;
t63 = pkin(4) * t218 + t138 + (-pkin(10) + t190) * t150;
t160 = -pkin(3) - pkin(10);
t223 = qJ(4) * t154;
t77 = (t160 * t158 - pkin(2) - t223) * t148;
t234 = t153 * t63 + t157 * t77;
t264 = qJD(5) * t234;
t206 = qJD(3) * t148;
t188 = t154 * t206;
t135 = pkin(3) * t188;
t204 = qJD(4) * t154;
t59 = t135 + (-t204 + (pkin(10) * t154 - qJ(4) * t158) * qJD(3)) * t148;
t215 = t150 * t154;
t141 = pkin(2) * t215;
t217 = t148 * t158;
t252 = pkin(4) + pkin(9);
t78 = (t252 * t217 + t141) * qJD(3);
t14 = -t153 * t59 + t157 * t78 - t264;
t155 = sin(qJ(2));
t149 = sin(pkin(6));
t207 = qJD(2) * t149;
t189 = t155 * t207;
t180 = t148 * t189;
t151 = cos(pkin(6));
t191 = t149 * t154 * t155;
t159 = cos(qJ(2));
t216 = t149 * t159;
t192 = t150 * t216;
t64 = -t151 * t217 - t158 * t192 + t191;
t97 = -t148 * t216 + t150 * t151;
t42 = t153 * t97 - t64 * t157;
t222 = qJD(5) * t42;
t212 = t155 * t158;
t213 = t154 * t159;
t165 = t150 * t213 + t212;
t40 = t151 * t188 + (t165 * qJD(3) + (t150 * t212 + t213) * qJD(2)) * t149;
t16 = t153 * t40 + t157 * t180 - t222;
t43 = t153 * t64 + t157 * t97;
t221 = qJD(5) * t43;
t15 = t153 * t180 - t40 * t157 + t221;
t227 = t15 * t157;
t263 = qJD(5) * (t153 * t42 + t157 * t43) + t16 * t153 - t227;
t124 = -mrSges(7,1) * t156 + mrSges(7,2) * t152;
t262 = -m(7) * pkin(5) - mrSges(6,1) + t124;
t261 = 0.2e1 * m(7);
t260 = 2 * mrSges(5,1);
t259 = -2 * mrSges(4,3);
t258 = 0.2e1 * t160;
t257 = m(6) / 0.2e1;
t256 = m(7) / 0.2e1;
t193 = t153 * t217;
t99 = t150 * t157 - t193;
t73 = -t152 * t99 + t156 * t218;
t255 = t73 / 0.2e1;
t74 = t152 * t218 + t156 * t99;
t254 = t74 / 0.2e1;
t230 = Ifges(7,4) * t152;
t176 = Ifges(7,1) * t156 - t230;
t92 = Ifges(7,5) * t153 + t176 * t157;
t253 = t92 / 0.2e1;
t198 = qJD(6) * t156;
t143 = Ifges(7,5) * t198;
t199 = qJD(6) * t152;
t251 = -Ifges(7,6) * t199 / 0.2e1 + t143 / 0.2e1;
t250 = Ifges(7,5) * t268 + Ifges(7,6) * t245;
t127 = Ifges(7,2) * t156 + t230;
t249 = t127 / 0.2e1;
t248 = -t152 / 0.2e1;
t246 = -t156 / 0.2e1;
t205 = qJD(3) * t158;
t187 = t148 * t205;
t98 = t150 * t153 + t157 * t217;
t71 = -qJD(5) * t98 + t153 * t188;
t50 = mrSges(6,1) * t187 - mrSges(6,3) * t71;
t28 = qJD(6) * t73 + t152 * t187 + t156 * t71;
t29 = -qJD(6) * t74 - t152 * t71 + t156 * t187;
t8 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t244 = -t8 + t50;
t243 = m(7) * t153;
t241 = pkin(5) * t153;
t240 = pkin(11) * t157;
t239 = t15 * t42;
t102 = pkin(9) * t217 + t141;
t96 = t102 * qJD(3);
t238 = t64 * t96;
t41 = -t189 * t215 - qJD(3) * t191 + (t159 * t207 + (t148 * t151 + t192) * qJD(3)) * t158;
t65 = t149 * t165 + t151 * t218;
t22 = t65 * t41;
t38 = -mrSges(7,1) * t73 + mrSges(7,2) * t74;
t81 = mrSges(6,1) * t218 - mrSges(6,3) * t99;
t235 = t38 - t81;
t233 = mrSges(7,3) * t157;
t232 = Ifges(6,4) * t153;
t231 = Ifges(6,4) * t157;
t229 = Ifges(7,4) * t156;
t228 = Ifges(7,6) * t152;
t24 = Ifges(7,4) * t74 + Ifges(7,2) * t73 + Ifges(7,6) * t98;
t226 = t152 * t24;
t25 = Ifges(7,1) * t74 + Ifges(7,4) * t73 + Ifges(7,5) * t98;
t225 = t156 * t25;
t108 = -mrSges(5,1) * t217 - mrSges(5,3) * t150;
t56 = mrSges(6,1) * t98 + mrSges(6,2) * t99;
t224 = -t108 + t56;
t120 = qJ(4) - t240 + t241;
t214 = t153 * t160;
t85 = t156 * t120 - t152 * t214;
t220 = qJD(6) * t85;
t86 = t152 * t120 + t156 * t214;
t219 = qJD(6) * t86;
t211 = t157 * t160;
t210 = t237 * t150 + t269 * t218;
t209 = Ifges(4,5) * t187 + Ifges(5,5) * t188;
t203 = qJD(5) * t152;
t201 = qJD(5) * t156;
t200 = qJD(5) * t157;
t72 = -qJD(5) * t193 + t150 * t200 - t157 * t188;
t5 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t72;
t195 = t150 * t242;
t194 = Ifges(6,5) * t71 - Ifges(6,6) * t72 + Ifges(6,3) * t187;
t87 = -t150 * qJ(4) - t102;
t186 = t160 * t200;
t104 = qJD(4) + (pkin(5) * t157 + pkin(11) * t153) * qJD(5);
t44 = t152 * t104 + t156 * t186 + t220;
t182 = t44 - t220;
t45 = t156 * t104 - t152 * t186 - t219;
t181 = -t45 - t219;
t76 = pkin(4) * t217 - t87;
t31 = pkin(11) * t218 + t234;
t37 = pkin(5) * t98 - pkin(11) * t99 + t76;
t11 = -t152 * t31 + t156 * t37;
t136 = qJD(3) * t195;
t144 = t150 * qJD(4);
t62 = -t252 * t188 + t136 + t144;
t21 = pkin(5) * t72 - pkin(11) * t71 + t62;
t13 = t153 * t78 + t157 * t59 + t63 * t200 - t77 * t202;
t9 = pkin(11) * t187 + t13;
t1 = qJD(6) * t11 + t152 * t21 + t156 * t9;
t12 = t152 * t37 + t156 * t31;
t2 = -qJD(6) * t12 - t152 * t9 + t156 * t21;
t179 = t1 * t156 - t152 * t2;
t20 = t152 * t65 + t156 * t43;
t3 = -qJD(6) * t20 - t152 * t16 + t156 * t41;
t19 = -t152 * t43 + t156 * t65;
t4 = qJD(6) * t19 + t152 * t41 + t156 * t16;
t178 = -t3 * t152 + t4 * t156;
t177 = mrSges(7,1) * t152 + mrSges(7,2) * t156;
t129 = Ifges(7,1) * t152 + t229;
t175 = -Ifges(7,2) * t152 + t229;
t174 = -pkin(3) * t158 - t223;
t17 = mrSges(7,1) * t72 - mrSges(7,3) * t28;
t18 = -mrSges(7,2) * t72 + mrSges(7,3) * t29;
t172 = -t152 * t17 + t156 * t18;
t34 = -t153 * t77 + t157 * t63;
t171 = t153 * t34 - t157 * t234;
t166 = t41 * qJ(4) + t65 * qJD(4);
t95 = -pkin(9) * t188 + t136;
t164 = t152 * t197 + t153 * t201;
t161 = -t11 * t198 - t12 * t199 + t179;
t52 = -Ifges(7,5) * t164 + t163 * Ifges(7,6) + Ifges(7,3) * t200;
t130 = Ifges(6,1) * t157 - t232;
t128 = -Ifges(6,2) * t153 + t231;
t118 = mrSges(7,1) * t153 - t156 * t233;
t117 = -mrSges(7,2) * t153 - t152 * t233;
t116 = (-Ifges(6,1) * t153 - t231) * qJD(5);
t115 = t176 * qJD(6);
t114 = (-Ifges(6,2) * t157 - t232) * qJD(5);
t113 = t175 * qJD(6);
t111 = (mrSges(6,1) * t157 - mrSges(6,2) * t153) * qJD(5);
t110 = t177 * qJD(6);
t107 = -mrSges(4,2) * t150 + mrSges(4,3) * t217;
t103 = t177 * t157;
t101 = -t138 + t195;
t100 = (mrSges(5,2) * t158 - mrSges(5,3) * t154) * t148;
t94 = (mrSges(4,1) * t154 + mrSges(4,2) * t158) * t206;
t93 = (-mrSges(5,2) * t154 - mrSges(5,3) * t158) * t206;
t91 = Ifges(7,6) * t153 + t175 * t157;
t90 = Ifges(7,3) * t153 + (Ifges(7,5) * t156 - t228) * t157;
t89 = t190 * t150 + t138;
t88 = (-pkin(2) + t174) * t148;
t84 = -t144 - t95;
t83 = -mrSges(7,2) * t200 + mrSges(7,3) * t163;
t82 = mrSges(7,1) * t200 + mrSges(7,3) * t164;
t80 = -mrSges(6,2) * t218 - mrSges(6,3) * t98;
t79 = t135 + (-qJ(4) * t205 - t204) * t148;
t60 = -mrSges(7,1) * t163 - mrSges(7,2) * t164;
t54 = -t129 * t197 + (Ifges(7,5) * t157 - t176 * t153) * qJD(5);
t53 = -t127 * t197 + (Ifges(7,6) * t157 - t175 * t153) * qJD(5);
t51 = -mrSges(6,2) * t187 - mrSges(6,3) * t72;
t49 = Ifges(6,1) * t99 - Ifges(6,4) * t98 + Ifges(6,5) * t218;
t48 = Ifges(6,4) * t99 - Ifges(6,2) * t98 + Ifges(6,6) * t218;
t47 = mrSges(7,1) * t98 - mrSges(7,3) * t74;
t46 = -mrSges(7,2) * t98 + mrSges(7,3) * t73;
t36 = mrSges(6,1) * t72 + mrSges(6,2) * t71;
t33 = Ifges(6,1) * t71 - Ifges(6,4) * t72 + Ifges(6,5) * t187;
t32 = Ifges(6,4) * t71 - Ifges(6,2) * t72 + Ifges(6,6) * t187;
t30 = -pkin(5) * t218 - t34;
t23 = Ifges(7,5) * t74 + Ifges(7,6) * t73 + Ifges(7,3) * t98;
t10 = -pkin(5) * t187 - t14;
t7 = Ifges(7,1) * t28 + Ifges(7,4) * t29 + Ifges(7,5) * t72;
t6 = Ifges(7,4) * t28 + Ifges(7,2) * t29 + Ifges(7,6) * t72;
t26 = [0.2e1 * m(7) * (t19 * t3 + t20 * t4 + t239) + 0.2e1 * m(6) * (t16 * t43 + t22 + t239) + 0.2e1 * (m(5) + m(4)) * (t97 * t180 + t40 * t64 + t22); t16 * t80 + t19 * t17 + t20 * t18 + t3 * t47 + t65 * t36 + t4 * t46 + t43 * t51 + (t93 + t94) * t97 - t244 * t42 + t210 * t40 + t235 * t15 + (-mrSges(3,1) * t155 - mrSges(3,2) * t159) * t207 + (t107 + t224) * t41 + ((t100 + (-mrSges(4,1) * t158 + mrSges(4,2) * t154) * t148) * t189 + t269 * qJD(3) * (-t154 * t65 + t158 * t64)) * t148 + m(5) * (t88 * t180 + t40 * t89 - t41 * t87 - t65 * t84 + t79 * t97 + t238) + m(4) * (-pkin(2) * t148 ^ 2 * t189 - t101 * t40 + t102 * t41 + t65 * t95 + t238) + m(6) * (t13 * t43 - t14 * t42 - t15 * t34 + t16 * t234 + t41 * t76 + t62 * t65) + m(7) * (t1 * t20 + t10 * t42 + t11 * t3 + t12 * t4 + t15 * t30 + t19 * t2); t209 * t150 + (t23 - t48) * t72 + (t154 * t194 - 0.2e1 * pkin(2) * t94 + ((t87 * t260 + t102 * t259 - t154 * t265 + (Ifges(5,5) - (2 * Ifges(4,6))) * t150) * t154 + (Ifges(6,5) * t99 - Ifges(6,6) * t98 + t89 * t260 + t101 * t259 + t158 * t265 + (-(2 * Ifges(5,4)) + Ifges(4,5)) * t150 + ((2 * Ifges(4,1)) - (2 * Ifges(4,2)) + (2 * Ifges(5,2)) - (2 * Ifges(5,3)) + Ifges(6,3)) * t218) * t158) * qJD(3)) * t148 + 0.2e1 * t210 * t96 + 0.2e1 * t234 * t51 + 0.2e1 * m(6) * (t13 * t234 + t14 * t34 + t62 * t76) + (t1 * t12 + t10 * t30 + t11 * t2) * t261 + 0.2e1 * t95 * t107 + 0.2e1 * t84 * t108 + t99 * t33 + 0.2e1 * t79 * t100 + 0.2e1 * t88 * t93 + 0.2e1 * t13 * t80 + 0.2e1 * t14 * t81 + t73 * t6 + t74 * t7 + 0.2e1 * t76 * t36 + t71 * t49 + 0.2e1 * t62 * t56 + 0.2e1 * t34 * t50 + 0.2e1 * t1 * t46 + 0.2e1 * t2 * t47 + 0.2e1 * t10 * t38 + t29 * t24 + 0.2e1 * t30 * t8 + t28 * t25 + 0.2e1 * t11 * t17 + 0.2e1 * t12 * t18 + 0.2e1 * m(5) * (t79 * t88 + t84 * t87 + t89 * t96) + 0.2e1 * m(4) * (-t101 * t96 + t102 * t95) + (t5 - t32) * t98; t15 * t103 + t65 * t111 + t4 * t117 + t3 * t118 + t19 * t82 + t20 * t83 + t42 * t60 + t237 * t40 + (-mrSges(4,2) + t267) * t41 + m(7) * (t45 * t19 + t44 * t20 + t85 * t3 + t86 * t4) + m(6) * t166 + m(5) * (-pkin(3) * t40 + t166) + ((t42 * t202 - t227) * t256 + t263 * t257) * t258 - t263 * mrSges(6,3); m(6) * (qJ(4) * t62 + qJD(4) * t76 + t13 * t214 + t14 * t211) + m(7) * (t86 * t1 - t10 * t211 + t45 * t11 + t44 * t12 + t85 * t2) + (-t32 / 0.2e1 + t5 / 0.2e1 + t160 * t51 - t13 * mrSges(6,3)) * t153 + t224 * qJD(4) + (t153 * t226 / 0.2e1 + (-Ifges(6,5) * t153 - Ifges(6,6) * t157) * t218 / 0.2e1 + t157 * t23 / 0.2e1 - t157 * t48 / 0.2e1 + t171 * mrSges(6,3) + (-m(6) * t171 + t235 * t153 + t157 * t80 + t30 * t243) * t160 - (t225 + t49) * t153 / 0.2e1) * qJD(5) + t209 + t28 * t253 + t54 * t254 + t53 * t255 + m(5) * (-pkin(3) * t96 - qJ(4) * t84 - qJD(4) * t87) + (-t114 / 0.2e1 + t52 / 0.2e1) * t98 + (-t128 / 0.2e1 + t90 / 0.2e1) * t72 + (t158 * (Ifges(6,5) * t157 - Ifges(6,6) * t153) / 0.2e1 - t158 * Ifges(5,4) - t154 * Ifges(4,6) + t174 * mrSges(5,1)) * t206 + t237 * t96 + (t33 / 0.2e1 + t7 * t245 + t6 * t248 - t14 * mrSges(6,3) + t244 * t160 + (t24 * t246 + t25 * t248) * qJD(6)) * t157 + t62 * t125 + t71 * t130 / 0.2e1 + t1 * t117 + t2 * t118 + t76 * t111 + t99 * t116 / 0.2e1 + t10 * t103 + t29 * t91 / 0.2e1 - t95 * mrSges(4,2) + t11 * t82 + t12 * t83 - t84 * mrSges(5,3) + t85 * t17 + t86 * t18 + t30 * t60 + t44 * t46 + t45 * t47 + qJ(4) * t36; 0.2e1 * t44 * t117 + 0.2e1 * t86 * t83 + (t86 * t44 + t85 * t45) * t261 + 0.2e1 * t45 * t118 + 0.2e1 * t85 * t82 + 0.2e1 * qJ(4) * t111 + 0.2e1 * ((m(5) + m(6)) * qJ(4) + t267) * qJD(4) + (-t114 + t52 + (t103 * t258 + t152 * t91 - t156 * t92 - t130) * qJD(5)) * t153 + (-t152 * t53 + t156 * t54 - 0.2e1 * t160 * t60 + t116 + (-t152 * t92 - t156 * t91) * qJD(6) + (-0.2e1 * t160 ^ 2 * t243 - t128 + t90) * qJD(5)) * t157; m(5) * t40 + 0.2e1 * ((-t19 * t203 + t20 * t201 - t15) * t256 + (-t15 + t221) * t257) * t157 + 0.2e1 * ((-t19 * t198 - t20 * t199 + t178 + t222) * t256 + (t16 + t222) * t257) * t153; m(5) * t96 + mrSges(5,1) * t187 + ((-t152 * t47 + t156 * t46 + t80) * qJD(5) + m(7) * (-t11 * t203 + t12 * t201 - t10) + m(6) * (t14 + t264) + t244) * t157 + (t51 + (-t152 * t46 - t156 * t47) * qJD(6) + t235 * qJD(5) + m(7) * (qJD(5) * t30 + t161) + m(6) * (-qJD(5) * t34 + t13) + t172) * t153; (-t60 + (m(7) * (-t152 * t85 + t156 * t86) + t156 * t117 - t152 * t118) * qJD(5)) * t157 + (m(7) * (-t152 * t45 + t156 * t44 - t85 * t198 - t86 * t199 - 0.2e1 * t186) - t117 * t199 + t156 * t83 - t118 * t198 - t152 * t82 + qJD(5) * t103) * t153; 0.2e1 * (-0.1e1 + t208) * t200 * t243; -t16 * mrSges(6,2) + t42 * t110 + (m(7) * pkin(11) + mrSges(7,3)) * ((-t152 * t20 - t156 * t19) * qJD(6) + t178) + t262 * t15; t6 * t245 + t7 * t268 + t72 * t250 + t29 * t249 + t28 * t129 / 0.2e1 + t10 * t124 + t30 * t110 + t98 * t251 + t113 * t255 + t115 * t254 - t13 * mrSges(6,2) + t14 * mrSges(6,1) + (t225 / 0.2e1 - t226 / 0.2e1) * qJD(6) + (-m(7) * t10 - t8) * pkin(5) + ((-t11 * t156 - t12 * t152) * qJD(6) + t179) * mrSges(7,3) + (m(7) * t161 - t47 * t198 - t46 * t199 + t172) * pkin(11) + t194; -pkin(5) * t60 + (t251 + (t160 * t262 - Ifges(6,5)) * qJD(5)) * t153 + (-t129 * t202 / 0.2e1 + t53 / 0.2e1 + qJD(6) * t253 + t182 * mrSges(7,3) + (m(7) * t182 - qJD(6) * t118 + t83) * pkin(11)) * t156 + (t202 * t249 + t54 / 0.2e1 - qJD(6) * t91 / 0.2e1 + t181 * mrSges(7,3) + (m(7) * t181 - qJD(6) * t117 - t82) * pkin(11)) * t152 + (-t160 * t110 + t115 * t245 + t113 * t248 + (t127 * t246 + t129 * t248) * qJD(6) + (-t160 * mrSges(6,2) - Ifges(6,6) + t250) * qJD(5)) * t157; -t157 * t110 + (t153 * t124 + m(7) * (t208 * t240 - t241) + t208 * t233 - t125) * qJD(5); -0.2e1 * pkin(5) * t110 + t113 * t156 + t115 * t152 + (-t127 * t152 + t129 * t156) * qJD(6); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t45 - mrSges(7,2) * t44 + t52; (t153 * t199 - t156 * t200) * mrSges(7,2) + (-t152 * t200 - t153 * t198) * mrSges(7,1); t143 + (t124 * pkin(11) - t228) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t26(1) t26(2) t26(4) t26(7) t26(11) t26(16); t26(2) t26(3) t26(5) t26(8) t26(12) t26(17); t26(4) t26(5) t26(6) t26(9) t26(13) t26(18); t26(7) t26(8) t26(9) t26(10) t26(14) t26(19); t26(11) t26(12) t26(13) t26(14) t26(15) t26(20); t26(16) t26(17) t26(18) t26(19) t26(20) t26(21);];
Mq  = res;
