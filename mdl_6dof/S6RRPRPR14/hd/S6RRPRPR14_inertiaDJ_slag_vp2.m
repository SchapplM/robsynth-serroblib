% Calculate time derivative of joint inertia matrix for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:36
% EndTime: 2019-03-09 11:32:48
% DurationCPUTime: 4.59s
% Computational Cost: add. (4140->539), mult. (10435->754), div. (0->0), fcn. (9049->8), ass. (0->225)
t275 = Ifges(6,1) + Ifges(5,3);
t271 = Ifges(5,5) - Ifges(6,4);
t270 = -Ifges(5,6) + Ifges(6,5);
t168 = sin(qJ(6));
t169 = sin(qJ(4));
t171 = cos(qJ(6));
t212 = qJD(6) * t171;
t172 = cos(qJ(4));
t215 = qJD(4) * t172;
t178 = t168 * t215 + t169 * t212;
t274 = 2 * qJ(3);
t273 = m(6) + m(7);
t252 = t168 / 0.2e1;
t272 = -t171 / 0.2e1;
t130 = mrSges(7,1) * t168 + mrSges(7,2) * t171;
t269 = mrSges(6,3) + t130;
t268 = -t169 * pkin(4) + qJ(5) * t172;
t166 = sin(pkin(6));
t170 = sin(qJ(2));
t217 = qJD(4) * t169;
t232 = t166 * t170;
t206 = qJD(2) * t232;
t146 = pkin(2) * t206;
t173 = cos(qJ(2));
t219 = qJD(3) * t170;
t63 = t146 + (-t219 + (pkin(9) * t170 - qJ(3) * t173) * qJD(2)) * t166;
t151 = pkin(8) * t232;
t167 = cos(pkin(6));
t249 = pkin(1) * t173;
t207 = -pkin(2) - t249;
t69 = pkin(3) * t232 + t151 + (-pkin(9) + t207) * t167;
t175 = -pkin(2) - pkin(9);
t234 = qJ(3) * t170;
t83 = (t173 * t175 - pkin(1) - t234) * t166;
t154 = t167 * t170 * pkin(1);
t231 = t166 * t173;
t258 = pkin(3) + pkin(8);
t84 = (t231 * t258 + t154) * qJD(2);
t182 = -t169 * t84 - t172 * t63 - t69 * t215 + t217 * t83;
t220 = qJD(2) * t173;
t11 = -t166 * (qJ(5) * t220 + qJD(5) * t170) + t182;
t205 = t166 * t220;
t208 = t169 * t231;
t79 = -qJD(4) * t208 + t167 * t215 - t172 * t206;
t103 = t167 * t169 + t172 * t231;
t81 = t103 * t168 + t171 * t232;
t29 = -qJD(6) * t81 - t168 * t205 + t171 * t79;
t78 = qJD(4) * t103 - t169 * t206;
t17 = mrSges(7,2) * t78 + mrSges(7,3) * t29;
t80 = t103 * t171 - t168 * t232;
t30 = qJD(6) * t80 + t168 * t79 + t171 * t205;
t18 = -mrSges(7,1) * t78 - mrSges(7,3) * t30;
t104 = t167 * t172 - t208;
t43 = -mrSges(7,2) * t104 + mrSges(7,3) * t80;
t44 = mrSges(7,1) * t104 - mrSges(7,3) * t81;
t186 = t168 * t44 - t171 * t43;
t267 = -t186 * qJD(6) + t168 * t17 + t171 * t18;
t266 = 2 * m(7);
t265 = -0.2e1 * pkin(1);
t264 = 2 * mrSges(4,1);
t263 = -2 * mrSges(3,3);
t189 = -pkin(2) * t173 - t234;
t95 = (-pkin(1) + t189) * t166;
t262 = -0.2e1 * t95;
t261 = t29 / 0.2e1;
t260 = t80 / 0.2e1;
t259 = t81 / 0.2e1;
t257 = pkin(4) + pkin(10);
t190 = Ifges(7,5) * t168 + Ifges(7,6) * t171;
t256 = -t190 * qJD(6) / 0.2e1;
t255 = Ifges(7,5) * t272 + Ifges(7,6) * t252;
t239 = Ifges(7,4) * t168;
t138 = Ifges(7,1) * t171 - t239;
t254 = t138 / 0.2e1;
t253 = -t168 / 0.2e1;
t251 = t171 / 0.2e1;
t248 = Ifges(4,6) + Ifges(3,4);
t247 = pkin(5) - t175;
t49 = mrSges(6,1) * t79 - mrSges(6,3) * t205;
t52 = -mrSges(5,2) * t205 - mrSges(5,3) * t79;
t246 = -t49 + t52;
t50 = -t78 * mrSges(6,1) + mrSges(6,2) * t205;
t51 = mrSges(5,1) * t205 + mrSges(5,3) * t78;
t245 = -t50 + t51;
t38 = t169 * t69 + t172 * t83;
t86 = mrSges(6,1) * t103 - mrSges(6,3) * t232;
t88 = -mrSges(5,2) * t232 - mrSges(5,3) * t103;
t244 = -t86 + t88;
t87 = mrSges(6,1) * t104 + mrSges(6,2) * t232;
t89 = mrSges(5,1) * t232 - mrSges(5,3) * t104;
t243 = t87 - t89;
t242 = mrSges(7,3) * t169;
t241 = Ifges(5,4) * t169;
t240 = Ifges(5,4) * t172;
t238 = Ifges(7,4) * t171;
t237 = Ifges(6,6) * t169;
t236 = Ifges(6,6) * t172;
t209 = t167 * t249;
t147 = qJD(2) * t209;
t101 = -pkin(8) * t206 + t147;
t235 = t101 * mrSges(3,2);
t230 = t168 * t138;
t229 = t168 * t257;
t228 = t169 * t175;
t136 = -Ifges(7,2) * t168 + t238;
t227 = t171 * t136;
t226 = t171 * t257;
t225 = t172 * t175;
t224 = Ifges(3,5) * t205 + Ifges(4,5) * t206;
t106 = pkin(8) * t231 + t154;
t222 = Ifges(6,4) * t217 + Ifges(6,5) * t215;
t221 = t168 ^ 2 + t171 ^ 2;
t218 = qJD(4) * t168;
t216 = qJD(4) * t171;
t214 = qJD(6) * t168;
t213 = qJD(6) * t169;
t211 = qJD(6) * t172;
t210 = 0.2e1 * t173;
t6 = Ifges(7,5) * t30 + Ifges(7,6) * t29 - Ifges(7,3) * t78;
t94 = -t167 * qJ(3) - t106;
t202 = t171 * t215;
t201 = -t232 / 0.2e1;
t200 = pkin(4) * t215 + qJ(5) * t217 + qJD(3);
t199 = m(6) * t175 - mrSges(6,1);
t198 = m(7) * t221;
t102 = t106 * qJD(2);
t197 = (mrSges(4,2) - mrSges(3,1)) * t102;
t37 = -t169 * t83 + t172 * t69;
t196 = qJD(4) * t247;
t82 = pkin(3) * t231 - t94;
t160 = t167 * qJD(3);
t68 = -t206 * t258 + t147 + t160;
t176 = qJ(5) * t78 - qJD(5) * t104 + t68;
t13 = t257 * t79 + t176;
t19 = pkin(5) * t104 - t232 * t257 - t37;
t180 = -qJ(5) * t104 + t82;
t21 = t103 * t257 + t180;
t3 = -t168 * t21 + t171 * t19;
t15 = -t169 * t63 + t172 * t84 - t83 * t215 - t69 * t217;
t5 = -pkin(5) * t78 - t205 * t257 - t15;
t1 = qJD(6) * t3 + t13 * t171 + t168 * t5;
t4 = t168 * t19 + t171 * t21;
t2 = -qJD(6) * t4 - t13 * t168 + t171 * t5;
t195 = t1 * t168 + t171 * t2;
t194 = t168 * t3 - t171 * t4;
t124 = qJ(3) - t268;
t132 = t169 * mrSges(5,1) + t172 * mrSges(5,2);
t193 = mrSges(7,1) * t171 - mrSges(7,2) * t168;
t131 = -t169 * mrSges(6,2) - t172 * mrSges(6,3);
t192 = Ifges(7,1) * t168 + t238;
t191 = Ifges(7,2) * t171 + t239;
t109 = t169 * t196;
t107 = pkin(10) * t169 + t124;
t126 = t247 * t172;
t65 = -t107 * t168 + t126 * t171;
t90 = (pkin(10) * qJD(4) - qJD(5)) * t172 + t200;
t22 = qJD(6) * t65 - t109 * t168 + t171 * t90;
t66 = t107 * t171 + t126 * t168;
t23 = -qJD(6) * t66 - t109 * t171 - t168 * t90;
t187 = t168 * t22 + t171 * t23;
t185 = -t168 * t65 + t171 * t66;
t31 = -qJ(5) * t232 - t38;
t184 = t205 * t275 + t270 * t79 - t271 * t78;
t25 = Ifges(7,4) * t81 + Ifges(7,2) * t80 + Ifges(7,6) * t104;
t26 = Ifges(7,1) * t81 + Ifges(7,4) * t80 + Ifges(7,5) * t104;
t183 = t25 * t272 + t26 * t253;
t179 = -t168 * t213 + t202;
t53 = Ifges(7,6) * t202 + (-Ifges(7,6) * t214 - Ifges(7,3) * qJD(4)) * t169 + t178 * Ifges(7,5);
t122 = mrSges(7,1) * t172 - t168 * t242;
t123 = -mrSges(7,2) * t172 + t171 * t242;
t91 = mrSges(7,2) * t217 + mrSges(7,3) * t179;
t92 = -mrSges(7,1) * t217 - mrSges(7,3) * t178;
t177 = -t122 * t214 + t123 * t212 + t168 * t91 + t171 * t92;
t139 = Ifges(5,1) * t172 - t241;
t137 = -Ifges(5,2) * t169 + t240;
t134 = -Ifges(6,2) * t172 + t237;
t133 = Ifges(6,3) * t169 - t236;
t125 = t247 * t169;
t121 = (-Ifges(5,1) * t169 - t240) * qJD(4);
t120 = t192 * qJD(6);
t119 = (-Ifges(5,2) * t172 - t241) * qJD(4);
t118 = t191 * qJD(6);
t116 = (Ifges(6,2) * t169 + t236) * qJD(4);
t115 = (Ifges(6,3) * t172 + t237) * qJD(4);
t114 = (mrSges(5,1) * t172 - mrSges(5,2) * t169) * qJD(4);
t113 = (-mrSges(6,2) * t172 + mrSges(6,3) * t169) * qJD(4);
t112 = t193 * qJD(6);
t111 = -mrSges(4,1) * t231 - mrSges(4,3) * t167;
t110 = t172 * t196;
t108 = t193 * t169;
t105 = -t151 + t209;
t100 = -qJD(5) * t172 + t200;
t99 = Ifges(7,5) * t172 + t169 * t192;
t98 = Ifges(7,6) * t172 + t169 * t191;
t97 = Ifges(7,3) * t172 + t169 * t190;
t96 = t167 * t207 + t151;
t93 = -t101 - t160;
t85 = t146 + (-qJ(3) * t220 - t219) * t166;
t64 = -mrSges(7,1) * t179 + mrSges(7,2) * t178;
t59 = -mrSges(6,2) * t103 - mrSges(6,3) * t104;
t58 = mrSges(5,1) * t103 + mrSges(5,2) * t104;
t55 = t138 * t213 + (-Ifges(7,5) * t169 + t172 * t192) * qJD(4);
t54 = t136 * t213 + (-Ifges(7,6) * t169 + t172 * t191) * qJD(4);
t48 = Ifges(5,1) * t104 - Ifges(5,4) * t103 + Ifges(5,5) * t232;
t47 = Ifges(5,4) * t104 - Ifges(5,2) * t103 + Ifges(5,6) * t232;
t46 = Ifges(6,4) * t232 - Ifges(6,2) * t104 + Ifges(6,6) * t103;
t45 = Ifges(6,5) * t232 - Ifges(6,6) * t104 + Ifges(6,3) * t103;
t42 = -mrSges(7,1) * t80 + mrSges(7,2) * t81;
t41 = mrSges(5,1) * t79 - mrSges(5,2) * t78;
t40 = -mrSges(6,2) * t79 + mrSges(6,3) * t78;
t39 = pkin(4) * t103 + t180;
t36 = -Ifges(5,1) * t78 - Ifges(5,4) * t79 + Ifges(5,5) * t205;
t35 = -Ifges(5,4) * t78 - Ifges(5,2) * t79 + Ifges(5,6) * t205;
t34 = Ifges(6,4) * t205 + Ifges(6,2) * t78 + Ifges(6,6) * t79;
t33 = Ifges(6,5) * t205 + Ifges(6,6) * t78 + Ifges(6,3) * t79;
t32 = -pkin(4) * t232 - t37;
t24 = Ifges(7,5) * t81 + Ifges(7,6) * t80 + Ifges(7,3) * t104;
t20 = -pkin(5) * t103 - t31;
t16 = pkin(4) * t79 + t176;
t12 = -pkin(4) * t205 - t15;
t10 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t9 = -pkin(5) * t79 - t11;
t8 = Ifges(7,1) * t30 + Ifges(7,4) * t29 - Ifges(7,5) * t78;
t7 = Ifges(7,4) * t30 + Ifges(7,2) * t29 - Ifges(7,6) * t78;
t14 = [0.2e1 * m(6) * (t11 * t31 + t12 * t32 + t16 * t39) + 0.2e1 * m(4) * (t102 * t96 + t85 * t95 + t93 * t94) + 0.2e1 * m(3) * (t101 * t106 - t102 * t105) + (-t34 + t36 + t6) * t104 + (t33 - t35) * t103 + (0.2e1 * t197 + t224 - 0.2e1 * t235) * t167 + (t45 - t47) * t79 + ((t85 * mrSges(4,2) + t101 * mrSges(3,3)) * t210 + (-0.2e1 * t85 * mrSges(4,3) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * t102 + t184) * t170 + ((t94 * t264 + mrSges(4,2) * t262 + t106 * t263 + (Ifges(4,5) - (2 * Ifges(3,6))) * t167 + (mrSges(3,1) * t265 - 0.2e1 * t170 * t248) * t166) * t170 + (t96 * t264 + t105 * t263 + mrSges(4,3) * t262 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t167 + (mrSges(3,2) * t265 + t210 * t248) * t166 + t271 * t104 + t270 * t103 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + t275) * t232) * t173) * qJD(2)) * t166 + 0.2e1 * t93 * t111 + t80 * t7 + t81 * t8 + 0.2e1 * t82 * t41 + 0.2e1 * t11 * t86 + 0.2e1 * t12 * t87 + 0.2e1 * t15 * t89 + 0.2e1 * t16 * t59 + 0.2e1 * t68 * t58 + 0.2e1 * t31 * t49 + 0.2e1 * t32 * t50 + 0.2e1 * t37 * t51 + 0.2e1 * t38 * t52 + 0.2e1 * t39 * t40 + 0.2e1 * t9 * t42 + 0.2e1 * t1 * t43 + 0.2e1 * t2 * t44 + t29 * t25 + t30 * t26 + 0.2e1 * t4 * t17 + 0.2e1 * t3 * t18 + 0.2e1 * t20 * t10 + 0.2e1 * m(5) * (t15 * t37 - t182 * t38 + t68 * t82) - 0.2e1 * t182 * t88 + (t46 - t48 - t24) * t78 + (t1 * t4 + t2 * t3 + t20 * t9) * t266; t197 + (t133 / 0.2e1 - t137 / 0.2e1) * t79 + (t134 / 0.2e1 - t139 / 0.2e1 - t97 / 0.2e1) * t78 + t224 + (t170 * t222 / 0.2e1 + (t189 * mrSges(4,1) - t173 * Ifges(4,4) - t170 * Ifges(3,6) + (t270 * t169 + t271 * t172) * t173 / 0.2e1) * qJD(2)) * t166 - t235 + (-t116 / 0.2e1 + t121 / 0.2e1 + t53 / 0.2e1) * t104 + (t115 / 0.2e1 - t119 / 0.2e1) * t103 + m(4) * (-pkin(2) * t102 - qJ(3) * t93 - qJD(3) * t94) + t2 * t122 + t1 * t123 + t124 * t40 - t125 * t10 + t16 * t131 + t68 * t132 - t9 * t108 - t110 * t42 + t39 * t113 + t82 * t114 + t30 * t99 / 0.2e1 + t100 * t59 + t4 * t91 + t3 * t92 - t93 * mrSges(4,3) + t20 * t64 + t65 * t18 + t66 * t17 + qJ(3) * t41 + t22 * t43 + t23 * t44 + (-t111 + t58) * qJD(3) + m(7) * (t1 * t66 - t110 * t20 - t125 * t9 + t2 * t65 + t22 * t4 + t23 * t3) + m(5) * (qJ(3) * t68 + qJD(3) * t82 + t15 * t225 - t182 * t228) + (t33 / 0.2e1 - t35 / 0.2e1 + t7 * t251 + t8 * t252 + t11 * mrSges(6,1) + t182 * mrSges(5,3) + t246 * t175 + (t25 * t253 + t251 * t26) * qJD(6)) * t169 + t55 * t259 + t54 * t260 + t98 * t261 + m(6) * (t100 * t39 - t11 * t228 - t12 * t225 + t124 * t16) + ((-t32 * mrSges(6,1) + t37 * mrSges(5,3) + Ifges(5,5) * t201 - t24 / 0.2e1 + t46 / 0.2e1 - t48 / 0.2e1) * t169 + (-t38 * mrSges(5,3) + t31 * mrSges(6,1) + Ifges(5,6) * t201 + t45 / 0.2e1 - t47 / 0.2e1 - t183) * t172 + (t244 * t172 + t243 * t169 + m(6) * (t32 * t169 - t31 * t172) + m(5) * (-t37 * t169 + t38 * t172)) * t175) * qJD(4) + (-t34 / 0.2e1 + t36 / 0.2e1 + t6 / 0.2e1 + t12 * mrSges(6,1) - t15 * mrSges(5,3) + t245 * t175) * t172; 0.2e1 * t23 * t122 + 0.2e1 * t22 * t123 + 0.2e1 * t124 * t113 - 0.2e1 * t125 * t64 + 0.2e1 * t110 * t108 + t114 * t274 + 0.2e1 * t66 * t91 + 0.2e1 * t65 * t92 + (t110 * t125 + t22 * t66 + t23 * t65) * t266 + (-t116 + t121 + t53 + (t168 * t99 + t171 * t98 + t133 - t137) * qJD(4)) * t172 + (t168 * t55 + t171 * t54 + t115 - t119 + (-t168 * t98 + t171 * t99) * qJD(6) + (t134 - t139 - t97) * qJD(4)) * t169 + 0.2e1 * (m(6) * t124 + t131) * t100 + (0.2e1 * mrSges(4,3) + 0.2e1 * t132 + (m(4) + m(5)) * t274) * qJD(3); m(4) * t102 + mrSges(4,1) * t205 + (t10 + (t168 * t43 + t171 * t44 + t243) * qJD(4) + m(7) * (t216 * t3 + t218 * t4 + t9) + m(6) * (qJD(4) * t32 - t11) + m(5) * (-qJD(4) * t37 - t182) + t246) * t169 + ((t42 + t244) * qJD(4) + m(7) * (qJD(4) * t20 - t212 * t4 + t214 * t3 - t195) + m(6) * (-qJD(4) * t31 - t12) + m(5) * (qJD(4) * t38 + t15) + t245 - t267) * t172; (m(7) * (t216 * t65 + t218 * t66 - t110) + t122 * t216 + t123 * t218 + t64) * t169 + (m(7) * (-qJD(4) * t125 - t212 * t66 + t214 * t65 - t187) - qJD(4) * t108 - t177) * t172; (0.1e1 - t221) * t169 * t215 * t266; t184 + (t194 * mrSges(7,3) - (-m(7) * t194 - t186) * t257 + t183) * qJD(6) + t78 * t255 + t136 * t261 + t30 * t254 + t104 * t256 - t118 * t260 - t120 * t259 + t9 * t130 + t20 * t112 - pkin(4) * t50 - t11 * mrSges(6,3) + t12 * mrSges(6,2) + t182 * mrSges(5,2) + t15 * mrSges(5,1) + m(6) * (-pkin(4) * t12 - qJ(5) * t11 - qJD(5) * t31) + m(7) * (qJ(5) * t9 + qJD(5) * t20 - t1 * t229 - t2 * t226) + (-t49 + t10) * qJ(5) + (-t86 + t42) * qJD(5) + (-t7 / 0.2e1 - t257 * t17 - t1 * mrSges(7,3)) * t168 + (t8 / 0.2e1 - t257 * t18 - t2 * mrSges(7,3)) * t171; t172 * t256 + t55 * t251 + t54 * t253 - t125 * t112 - t110 * t130 - qJD(5) * t108 + qJ(5) * t64 + m(7) * (-qJ(5) * t110 - qJD(5) * t125 - t22 * t229 - t226 * t23) - t92 * t226 - t91 * t229 - t187 * mrSges(7,3) + (t199 * qJD(5) - t118 * t251 - t120 * t252) * t169 + ((t169 * t254 - t66 * mrSges(7,3) - t98 / 0.2e1) * t171 + (-t169 * t136 / 0.2e1 + t65 * mrSges(7,3) - t99 / 0.2e1) * t168 - (m(7) * t185 - t168 * t122 + t171 * t123) * t257) * qJD(6) + ((pkin(4) * mrSges(6,1) - Ifges(5,5) + t255) * t169 + (-qJ(5) * mrSges(6,1) + t227 / 0.2e1 + t230 / 0.2e1 - Ifges(5,6)) * t172 + (m(6) * t268 - t131 - t132) * t175) * qJD(4) + t222; t169 * t112 + ((-mrSges(5,2) + t269) * t172 + (-m(6) * pkin(4) - mrSges(7,3) * t221 - t198 * t257 - mrSges(5,1) + mrSges(6,2)) * t169) * qJD(4) + t273 * (qJ(5) * t215 + qJD(5) * t169); 0.2e1 * qJ(5) * t112 + t118 * t168 - t120 * t171 + (-t227 - t230) * qJD(6) + 0.2e1 * (qJ(5) * t273 + t269) * qJD(5); m(7) * (-qJD(6) * t194 + t195) + m(6) * t12 + t50 + t267; m(7) * (qJD(6) * t185 + t187) + t199 * t217 + t177; (m(6) + t198) * t217; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t6; mrSges(7,1) * t23 - mrSges(7,2) * t22 + t53; (-t168 * t217 + t171 * t211) * mrSges(7,2) + (t168 * t211 + t169 * t216) * mrSges(7,1); ((mrSges(7,2) * t257 - Ifges(7,6)) * t171 + (mrSges(7,1) * t257 - Ifges(7,5)) * t168) * qJD(6); -t130 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
