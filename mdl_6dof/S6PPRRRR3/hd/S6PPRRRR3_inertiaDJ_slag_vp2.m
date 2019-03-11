% Calculate time derivative of joint inertia matrix for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_inertiaDJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:35
% EndTime: 2019-03-08 19:07:44
% DurationCPUTime: 4.20s
% Computational Cost: add. (6042->521), mult. (18844->802), div. (0->0), fcn. (20028->16), ass. (0->219)
t157 = sin(qJ(6));
t161 = cos(qJ(6));
t158 = sin(qJ(5));
t200 = qJD(6) * t158;
t162 = cos(qJ(5));
t202 = qJD(5) * t162;
t168 = -t157 * t200 + t161 * t202;
t263 = t157 / 0.2e1;
t242 = t161 / 0.2e1;
t149 = sin(pkin(14));
t151 = sin(pkin(7));
t152 = sin(pkin(6));
t156 = cos(pkin(6));
t160 = sin(qJ(3));
t164 = cos(qJ(3));
t153 = cos(pkin(14));
t155 = cos(pkin(7));
t213 = t153 * t155;
t74 = (t149 * t164 + t160 * t213) * t152 + t156 * t151 * t160;
t262 = qJD(3) * t74;
t150 = sin(pkin(8));
t206 = qJD(3) * t160;
t192 = t151 * t206;
t185 = t150 * t192;
t159 = sin(qJ(4));
t217 = t150 * t159;
t141 = pkin(10) * t217;
t154 = cos(pkin(8));
t163 = cos(qJ(4));
t211 = t154 * t163;
t112 = pkin(3) * t211 - t141;
t207 = t163 * t164;
t210 = t159 * t160;
t261 = t154 * t207 - t210;
t260 = -m(6) * pkin(4) - mrSges(6,1) * t162 + mrSges(6,2) * t158 - mrSges(5,1);
t259 = m(7) * pkin(12) + mrSges(7,3);
t205 = qJD(4) * t150;
t104 = (pkin(4) * t159 - pkin(11) * t163) * t205;
t105 = t112 * qJD(4);
t212 = t154 * t159;
t215 = t150 * t163;
t113 = pkin(3) * t212 + pkin(10) * t215;
t98 = pkin(11) * t154 + t113;
t99 = (-pkin(4) * t163 - pkin(11) * t159 - pkin(3)) * t150;
t232 = t158 * t99 + t162 * t98;
t32 = -qJD(5) * t232 + t104 * t162 - t105 * t158;
t130 = -mrSges(7,1) * t161 + mrSges(7,2) * t157;
t258 = -m(7) * pkin(5) - mrSges(6,1) + t130;
t257 = 0.2e1 * m(7);
t256 = 0.2e1 * pkin(11);
t255 = -2 * mrSges(5,3);
t254 = m(6) / 0.2e1;
t253 = m(7) / 0.2e1;
t216 = t150 * t162;
t110 = t154 * t158 + t159 * t216;
t171 = -t110 * t161 + t157 * t215;
t204 = qJD(4) * t159;
t191 = t150 * t204;
t109 = -t162 * t154 + t158 * t217;
t190 = t163 * t205;
t80 = -t109 * qJD(5) + t162 * t190;
t40 = t171 * qJD(6) - t157 * t80 + t161 * t191;
t251 = t40 / 0.2e1;
t82 = -t110 * t157 - t161 * t215;
t250 = t82 / 0.2e1;
t249 = -t171 / 0.2e1;
t107 = -t151 * t152 * t153 + t155 * t156;
t214 = t151 * t164;
t73 = t156 * t214 + (-t149 * t160 + t164 * t213) * t152;
t177 = t107 * t150 + t154 * t73;
t28 = t177 * t159 + t163 * t74;
t48 = t107 * t154 - t150 * t73;
t16 = t158 * t28 - t48 * t162;
t222 = t159 * t74;
t70 = t73 * qJD(3);
t15 = -t262 * t212 + t163 * t70 + (t177 * t163 - t222) * qJD(4);
t17 = t158 * t48 + t162 * t28;
t3 = qJD(5) * t17 + t15 * t158 - t216 * t262;
t248 = t16 * t3;
t228 = Ifges(7,4) * t157;
t180 = Ifges(7,1) * t161 - t228;
t102 = -Ifges(7,5) * t162 + t180 * t158;
t247 = t102 / 0.2e1;
t246 = Ifges(7,5) * t263 + Ifges(7,6) * t242;
t227 = Ifges(7,4) * t161;
t135 = Ifges(7,1) * t157 + t227;
t245 = t135 / 0.2e1;
t244 = -t157 / 0.2e1;
t243 = -t161 / 0.2e1;
t241 = pkin(3) * t150 ^ 2;
t240 = pkin(11) * t162;
t14 = t28 * qJD(4) + t159 * t70 + t211 * t262;
t27 = -t107 * t215 - t73 * t211 + t222;
t239 = t14 * t27;
t108 = -t150 * t214 + t154 * t155;
t208 = t160 * t163;
t209 = t159 * t164;
t169 = t154 * t209 + t208;
t76 = t169 * t151 + t155 * t217;
t176 = t162 * t108 - t158 * t76;
t46 = t155 * t190 + (t261 * qJD(4) + (-t154 * t210 + t207) * qJD(3)) * t151;
t50 = t108 * t158 + t162 * t76;
t20 = t50 * qJD(5) + t158 * t46 - t162 * t185;
t238 = t20 * t176;
t237 = t3 * t158;
t223 = t150 * t262;
t4 = -qJD(5) * t16 + t15 * t162 + t158 * t223;
t236 = t4 * t162;
t47 = t155 * t191 + (t169 * qJD(4) + (t154 * t208 + t209) * qJD(3)) * t151;
t75 = -t151 * t261 - t155 * t215;
t235 = t47 * t75;
t39 = t82 * qJD(6) + t157 * t191 + t161 * t80;
t18 = -mrSges(7,1) * t40 + mrSges(7,2) * t39;
t63 = mrSges(6,1) * t191 - mrSges(6,3) * t80;
t234 = t18 - t63;
t44 = -mrSges(7,1) * t82 - mrSges(7,2) * t171;
t85 = -mrSges(6,1) * t215 - mrSges(6,3) * t110;
t233 = t44 - t85;
t231 = mrSges(7,3) * t158;
t230 = Ifges(6,4) * t158;
t229 = Ifges(6,4) * t162;
t226 = Ifges(7,6) * t157;
t106 = t113 * qJD(4);
t225 = t106 * t27;
t224 = t106 * t75;
t19 = t176 * qJD(5) + t158 * t185 + t162 * t46;
t221 = t19 * t162;
t220 = t20 * t158;
t218 = -mrSges(5,1) * t154 + mrSges(6,1) * t109 + mrSges(6,2) * t110 + mrSges(5,3) * t217;
t203 = qJD(5) * t158;
t201 = qJD(6) * t157;
t199 = qJD(6) * t161;
t198 = qJD(6) * t162;
t81 = t110 * qJD(5) + t158 * t190;
t11 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t81;
t196 = Ifges(6,5) * t80 - Ifges(6,6) * t81 + Ifges(6,3) * t191;
t195 = Ifges(6,6) * t215;
t127 = (pkin(5) * t158 - pkin(12) * t162) * qJD(5);
t129 = -pkin(5) * t162 - pkin(12) * t158 - pkin(4);
t56 = t129 * t199 + t127 * t157 + (-t157 * t198 - t161 * t203) * pkin(11);
t94 = t129 * t161 - t157 * t240;
t187 = -qJD(6) * t94 + t56;
t57 = -t129 * t201 + t127 * t161 + (t157 * t203 - t161 * t198) * pkin(11);
t95 = t129 * t157 + t161 * t240;
t186 = -qJD(6) * t95 - t57;
t184 = t16 * t20 - t176 * t3;
t183 = t14 * t75 + t27 * t47;
t97 = t141 + (-pkin(3) * t163 - pkin(4)) * t154;
t51 = pkin(5) * t109 - pkin(12) * t110 + t97;
t53 = -pkin(12) * t215 + t232;
t21 = -t157 * t53 + t161 * t51;
t31 = t158 * t104 + t162 * t105 + t99 * t202 - t98 * t203;
t25 = pkin(12) * t191 + t31;
t36 = pkin(5) * t81 - pkin(12) * t80 + t106;
t5 = t21 * qJD(6) + t157 * t36 + t161 * t25;
t22 = t157 * t51 + t161 * t53;
t6 = -t22 * qJD(6) - t157 * t25 + t161 * t36;
t182 = -t6 * t157 + t5 * t161;
t181 = mrSges(7,1) * t157 + mrSges(7,2) * t161;
t179 = -Ifges(7,2) * t157 + t227;
t133 = Ifges(7,2) * t161 + t228;
t10 = t157 * t27 + t161 * t17;
t9 = -t157 * t17 + t161 * t27;
t30 = t157 * t75 + t161 * t50;
t29 = -t157 * t50 + t161 * t75;
t58 = -t158 * t98 + t162 * t99;
t174 = t16 * t202 + t237;
t34 = -Ifges(7,4) * t171 + Ifges(7,2) * t82 + Ifges(7,6) * t109;
t35 = -Ifges(7,1) * t171 + Ifges(7,4) * t82 + Ifges(7,5) * t109;
t173 = t35 * t242 + t34 * t244;
t172 = -t176 * t202 + t220;
t167 = t157 * t202 + t158 * t199;
t65 = Ifges(7,5) * t168 - t167 * Ifges(7,6) + Ifges(7,3) * t203;
t147 = Ifges(6,5) * t202;
t146 = Ifges(7,5) * t199;
t139 = Ifges(5,5) * t190;
t136 = Ifges(6,1) * t158 + t229;
t134 = Ifges(6,2) * t162 + t230;
t126 = -mrSges(7,1) * t162 - t161 * t231;
t125 = mrSges(7,2) * t162 - t157 * t231;
t124 = (Ifges(6,1) * t162 - t230) * qJD(5);
t123 = t180 * qJD(6);
t122 = (-Ifges(6,2) * t158 + t229) * qJD(5);
t121 = t179 * qJD(6);
t120 = -Ifges(7,6) * t201 + t146;
t119 = (mrSges(6,1) * t158 + mrSges(6,2) * t162) * qJD(5);
t118 = t181 * qJD(6);
t117 = -mrSges(5,2) * t154 + mrSges(5,3) * t215;
t114 = t181 * t158;
t111 = (-mrSges(5,1) * t163 + mrSges(5,2) * t159) * t150;
t103 = (mrSges(5,1) * t159 + mrSges(5,2) * t163) * t205;
t101 = -Ifges(7,6) * t162 + t179 * t158;
t100 = -Ifges(7,3) * t162 + (Ifges(7,5) * t161 - t226) * t158;
t88 = -mrSges(7,2) * t203 - t167 * mrSges(7,3);
t87 = mrSges(7,1) * t203 - t168 * mrSges(7,3);
t84 = mrSges(6,2) * t215 - mrSges(6,3) * t109;
t72 = t167 * mrSges(7,1) + t168 * mrSges(7,2);
t67 = -t135 * t200 + (Ifges(7,5) * t158 + t180 * t162) * qJD(5);
t66 = -t133 * t200 + (Ifges(7,6) * t158 + t179 * t162) * qJD(5);
t64 = -mrSges(6,2) * t191 - mrSges(6,3) * t81;
t61 = Ifges(6,1) * t110 - Ifges(6,4) * t109 - Ifges(6,5) * t215;
t60 = Ifges(6,4) * t110 - Ifges(6,2) * t109 - t195;
t55 = mrSges(7,1) * t109 + mrSges(7,3) * t171;
t54 = -mrSges(7,2) * t109 + mrSges(7,3) * t82;
t52 = pkin(5) * t215 - t58;
t43 = mrSges(6,1) * t81 + mrSges(6,2) * t80;
t42 = Ifges(6,1) * t80 - Ifges(6,4) * t81 + Ifges(6,5) * t191;
t41 = Ifges(6,4) * t80 - Ifges(6,2) * t81 + Ifges(6,6) * t191;
t33 = -Ifges(7,5) * t171 + Ifges(7,6) * t82 + Ifges(7,3) * t109;
t26 = -pkin(5) * t191 - t32;
t24 = -mrSges(7,2) * t81 + mrSges(7,3) * t40;
t23 = mrSges(7,1) * t81 - mrSges(7,3) * t39;
t13 = Ifges(7,1) * t39 + Ifges(7,4) * t40 + Ifges(7,5) * t81;
t12 = Ifges(7,4) * t39 + Ifges(7,2) * t40 + Ifges(7,6) * t81;
t8 = -t30 * qJD(6) - t157 * t19 + t161 * t47;
t7 = t29 * qJD(6) + t157 * t47 + t161 * t19;
t2 = qJD(6) * t9 + t14 * t157 + t161 * t4;
t1 = -qJD(6) * t10 + t14 * t161 - t157 * t4;
t37 = [0.2e1 * m(7) * (t1 * t9 + t10 * t2 + t248) + 0.2e1 * m(6) * (t17 * t4 + t239 + t248) + 0.2e1 * m(5) * (t15 * t28 + t223 * t48 + t239) + 0.2e1 * m(4) * (-t262 * t73 + t70 * t74); m(7) * (t1 * t29 + t10 * t7 + t2 * t30 + t8 * t9 + t184) + m(6) * (t17 * t19 + t4 * t50 + t183 + t184) + m(4) * (t160 * t70 - t73 * t206) * t151 + (t108 * t223 + t15 * t76 + t185 * t48 + t28 * t46 + t183) * m(5); 0.2e1 * m(7) * (t29 * t8 + t30 * t7 - t238) + 0.2e1 * m(6) * (t19 * t50 + t235 - t238) + 0.2e1 * m(5) * (t108 * t185 + t46 * t76 + t235); -t262 * mrSges(4,1) - t70 * mrSges(4,2) + t1 * t55 + t10 * t24 + t48 * t103 + t15 * t117 + t17 * t64 + t2 * t54 + t9 * t23 + t27 * t43 + t4 * t84 + t233 * t3 + t234 * t16 + t218 * t14 + (t262 * t111 + (-t159 * t28 + t163 * t27) * qJD(4) * mrSges(5,3)) * t150 + m(7) * (t1 * t21 + t10 * t5 + t16 * t26 + t2 * t22 + t3 * t52 + t6 * t9) + m(6) * (t14 * t97 - t16 * t32 + t17 * t31 + t232 * t4 - t3 * t58 + t225) + m(5) * (t105 * t28 - t112 * t14 + t113 * t15 - t241 * t262 + t225); t108 * t103 + t46 * t117 + t19 * t84 + t29 * t23 + t30 * t24 + t75 * t43 + t50 * t64 + t7 * t54 + t8 * t55 - t234 * t176 + t218 * t47 + t233 * t20 + (-t159 * t76 + t163 * t75) * mrSges(5,3) * t205 + (-mrSges(4,2) * t164 + (t111 * t150 - mrSges(4,1)) * t160) * t151 * qJD(3) + m(7) * (-t176 * t26 + t20 * t52 + t21 * t8 + t22 * t7 + t29 * t6 + t30 * t5) + m(6) * (t176 * t32 + t19 * t232 - t20 * t58 + t31 * t50 + t47 * t97 + t224) + m(5) * (t105 * t76 - t112 * t47 + t113 * t46 - t192 * t241 + t224); 0.2e1 * t218 * t106 + 0.2e1 * m(6) * (t106 * t97 + t232 * t31 + t32 * t58) + 0.2e1 * t232 * t64 + (t21 * t6 + t22 * t5 + t26 * t52) * t257 - t171 * t13 + (t11 - t41) * t109 + 0.2e1 * m(5) * (t105 * t113 - t106 * t112) + (t33 - t60) * t81 + t154 * t139 + (-t163 * t196 - 0.2e1 * pkin(3) * t103 + ((0.2e1 * Ifges(5,4) * t215 + Ifges(5,5) * t154 + t112 * t255) * t163 + (-0.2e1 * Ifges(5,4) * t217 + t113 * t255 + Ifges(6,5) * t110 - 0.2e1 * Ifges(5,6) * t154 - Ifges(6,6) * t109 + ((2 * Ifges(5,1)) - (2 * Ifges(5,2)) - Ifges(6,3)) * t215) * t159) * qJD(4)) * t150 + 0.2e1 * t21 * t23 + 0.2e1 * t22 * t24 + t39 * t35 + t40 * t34 + 0.2e1 * t26 * t44 + 0.2e1 * t52 * t18 + 0.2e1 * t5 * t54 + 0.2e1 * t6 * t55 + 0.2e1 * t58 * t63 + t80 * t61 + t82 * t12 + 0.2e1 * t31 * t84 + 0.2e1 * t32 * t85 + 0.2e1 * t97 * t43 + t110 * t42 + 0.2e1 * t105 * t117; -t15 * mrSges(5,2) + t1 * t126 + t10 * t88 + t3 * t114 + t27 * t119 + t2 * t125 + t16 * t72 + t9 * t87 + m(7) * (t1 * t94 + t10 * t56 + t2 * t95 + t57 * t9) + (t174 * t253 + (-t17 * t203 + t174 + t236) * t254) * t256 + (t237 + t236 + (-t158 * t17 + t16 * t162) * qJD(5)) * mrSges(6,3) + t260 * t14; -t46 * mrSges(5,2) + t20 * t114 + t75 * t119 + t7 * t125 + t8 * t126 + t29 * t87 + t30 * t88 - t176 * t72 + m(7) * (t29 * t57 + t30 * t56 + t7 * t95 + t8 * t94) + (t172 * t253 + (-t203 * t50 + t172 + t221) * t254) * t256 + (t220 + t221 + (-t158 * t50 - t162 * t176) * qJD(5)) * mrSges(6,3) + t260 * t47; m(7) * (t21 * t57 + t22 * t56 + t5 * t95 + t6 * t94) + ((t61 / 0.2e1 - t58 * mrSges(6,3) + t173) * t162 + (t195 / 0.2e1 - t60 / 0.2e1 + t33 / 0.2e1 - t232 * mrSges(6,3)) * t158) * qJD(5) + (t162 * t64 + t234 * t158 + (-t158 * t84 + t162 * t233) * qJD(5) + m(6) * (-t32 * t158 + t31 * t162 - t202 * t58 - t203 * t232) + m(7) * (t158 * t26 + t202 * t52)) * pkin(11) + t39 * t247 + t67 * t249 + t66 * t250 + t101 * t251 + (t100 / 0.2e1 - t134 / 0.2e1) * t81 + (t65 / 0.2e1 - t122 / 0.2e1) * t109 + t139 + t260 * t106 + (-t163 * t147 / 0.2e1 + (Ifges(6,5) * t158 / 0.2e1 + Ifges(6,6) * t162 / 0.2e1 - Ifges(5,6)) * t204) * t150 + (t13 * t242 - t32 * mrSges(6,3) + t12 * t244 + t42 / 0.2e1 + (t243 * t34 + t244 * t35) * qJD(6)) * t158 - pkin(4) * t43 + t56 * t54 + t57 * t55 + t52 * t72 + t21 * t87 + t22 * t88 + t94 * t23 + t95 * t24 - t105 * mrSges(5,2) + t26 * t114 + t97 * t119 + t110 * t124 / 0.2e1 + t5 * t125 + t6 * t126 + t80 * t136 / 0.2e1 + (t31 * mrSges(6,3) - t11 / 0.2e1 + t41 / 0.2e1) * t162; 0.2e1 * t57 * t126 + 0.2e1 * t94 * t87 + 0.2e1 * t56 * t125 + 0.2e1 * t95 * t88 + (t56 * t95 + t57 * t94) * t257 - 0.2e1 * pkin(4) * t119 + (t122 - t65 + (-t101 * t157 + t102 * t161 + t114 * t256 + t136) * qJD(5)) * t162 + (t72 * t256 - t157 * t66 + t161 * t67 + t124 + (-t101 * t161 - t102 * t157) * qJD(6) + (pkin(11) ^ 2 * t162 * t257 + t100 - t134) * qJD(5)) * t158; -t4 * mrSges(6,2) + t16 * t118 + t259 * (-t1 * t157 + t2 * t161 + (-t10 * t157 - t161 * t9) * qJD(6)) + t258 * t3; -t19 * mrSges(6,2) - t176 * t118 + t259 * (-t8 * t157 + t7 * t161 + (-t157 * t30 - t161 * t29) * qJD(6)) + t258 * t20; -t31 * mrSges(6,2) + t32 * mrSges(6,1) + t52 * t118 + t109 * t120 / 0.2e1 + t121 * t250 + t123 * t249 + t26 * t130 + t81 * t246 + t133 * t251 + t39 * t245 + t13 * t263 + t12 * t242 + t173 * qJD(6) + (-m(7) * t26 - t18) * pkin(5) + ((-t157 * t22 - t161 * t21) * qJD(6) + t182) * mrSges(7,3) + (m(7) * (-t199 * t21 - t201 * t22 + t182) + t161 * t24 - t157 * t23 - t55 * t199 - t54 * t201) * pkin(12) + t196; -pkin(5) * t72 + t147 + (-t120 / 0.2e1 + t258 * qJD(5) * pkin(11)) * t162 + (qJD(6) * t247 + t202 * t245 + t66 / 0.2e1 + t187 * mrSges(7,3) + (m(7) * t187 - qJD(6) * t126 + t88) * pkin(12)) * t161 + (-qJD(6) * t101 / 0.2e1 - t133 * t202 / 0.2e1 + t67 / 0.2e1 + t186 * mrSges(7,3) + (m(7) * t186 - qJD(6) * t125 - t87) * pkin(12)) * t157 + (t123 * t242 + t121 * t244 + pkin(11) * t118 + (t133 * t243 + t135 * t244) * qJD(6) + (pkin(11) * mrSges(6,2) - Ifges(6,6) + t246) * qJD(5)) * t158; -0.2e1 * pkin(5) * t118 + t121 * t161 + t123 * t157 + (-t133 * t157 + t135 * t161) * qJD(6); mrSges(7,1) * t1 - mrSges(7,2) * t2; mrSges(7,1) * t8 - mrSges(7,2) * t7; mrSges(7,1) * t6 - mrSges(7,2) * t5 + t11; mrSges(7,1) * t57 - mrSges(7,2) * t56 + t65; t146 + (pkin(12) * t130 - t226) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t37(1) t37(2) t37(4) t37(7) t37(11) t37(16); t37(2) t37(3) t37(5) t37(8) t37(12) t37(17); t37(4) t37(5) t37(6) t37(9) t37(13) t37(18); t37(7) t37(8) t37(9) t37(10) t37(14) t37(19); t37(11) t37(12) t37(13) t37(14) t37(15) t37(20); t37(16) t37(17) t37(18) t37(19) t37(20) t37(21);];
Mq  = res;
