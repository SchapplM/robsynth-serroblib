% Calculate time derivative of joint inertia matrix for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:22:01
% EndTime: 2019-03-09 12:22:11
% DurationCPUTime: 4.36s
% Computational Cost: add. (7028->469), mult. (16692->675), div. (0->0), fcn. (15537->8), ass. (0->191)
t259 = Ifges(7,4) + Ifges(6,5);
t257 = Ifges(7,6) - Ifges(6,6);
t260 = -mrSges(6,1) - mrSges(7,1);
t258 = -Ifges(7,2) - Ifges(6,3);
t184 = sin(pkin(10));
t185 = cos(pkin(10));
t187 = sin(qJ(4));
t189 = cos(qJ(4));
t159 = t184 * t189 + t185 * t187;
t188 = sin(qJ(2));
t140 = t159 * t188;
t204 = t184 * t187 - t185 * t189;
t141 = t204 * t188;
t186 = sin(qJ(5));
t242 = cos(qJ(5));
t101 = -t140 * t186 - t141 * t242;
t190 = cos(qJ(2));
t88 = -mrSges(6,1) * t190 - mrSges(6,3) * t101;
t89 = mrSges(7,1) * t190 + mrSges(7,2) * t101;
t256 = -t88 + t89;
t163 = -pkin(2) * t190 - qJ(3) * t188 - pkin(1);
t155 = t185 * t163;
t229 = t185 * t188;
t119 = -pkin(8) * t229 + t155 + (-pkin(7) * t184 - pkin(3)) * t190;
t228 = t185 * t190;
t134 = pkin(7) * t228 + t163 * t184;
t231 = t184 * t188;
t126 = -pkin(8) * t231 + t134;
t82 = t119 * t187 + t126 * t189;
t148 = t204 * qJD(4);
t149 = t159 * qJD(4);
t199 = -t159 * t186 - t204 * t242;
t73 = qJD(5) * t199 - t148 * t242 - t149 * t186;
t118 = t159 * t242 - t186 * t204;
t74 = qJD(5) * t118 - t186 * t148 + t149 * t242;
t255 = t257 * t74 + t259 * t73;
t225 = qJD(2) * t190;
t107 = -t149 * t188 - t204 * t225;
t108 = t148 * t188 - t159 * t225;
t226 = qJD(2) * t188;
t254 = -Ifges(5,5) * t107 - Ifges(5,6) * t108 - Ifges(5,3) * t226;
t146 = -t188 * qJD(3) + (pkin(2) * t188 - qJ(3) * t190) * qJD(2);
t216 = pkin(7) * t226;
t124 = t146 * t185 + t184 * t216;
t102 = (pkin(3) * t188 - pkin(8) * t228) * qJD(2) + t124;
t137 = t184 * t146;
t230 = t184 * t190;
t110 = t137 + (-pkin(7) * t229 - pkin(8) * t230) * qJD(2);
t29 = -qJD(4) * t82 + t102 * t189 - t110 * t187;
t20 = pkin(4) * t226 - pkin(9) * t107 + t29;
t222 = qJD(4) * t189;
t223 = qJD(4) * t187;
t28 = t102 * t187 + t110 * t189 + t119 * t222 - t126 * t223;
t23 = pkin(9) * t108 + t28;
t81 = t119 * t189 - t126 * t187;
t56 = -pkin(4) * t190 + pkin(9) * t141 + t81;
t65 = -pkin(9) * t140 + t82;
t237 = t186 * t56 + t242 * t65;
t6 = -qJD(5) * t237 - t186 * t23 + t20 * t242;
t253 = 2 * m(4);
t252 = 2 * m(5);
t251 = 2 * m(6);
t250 = 2 * m(7);
t249 = -0.2e1 * pkin(1);
t248 = 0.2e1 * pkin(7);
t247 = m(6) * pkin(4);
t246 = -t204 / 0.2e1;
t245 = t159 / 0.2e1;
t244 = t185 / 0.2e1;
t241 = pkin(4) * t149;
t240 = pkin(4) * t186;
t239 = t73 * mrSges(7,2);
t238 = pkin(8) + qJ(3);
t234 = mrSges(4,2) * t185;
t233 = Ifges(4,4) * t184;
t232 = Ifges(4,4) * t185;
t227 = -Ifges(5,5) * t148 - Ifges(5,6) * t149;
t164 = t238 * t185;
t212 = t238 * t184;
t128 = t189 * t164 - t187 * t212;
t214 = t184 * t225;
t142 = mrSges(4,1) * t214 + t225 * t234;
t179 = pkin(7) * t225;
t152 = pkin(3) * t214 + t179;
t162 = pkin(3) * t231 + pkin(7) * t188;
t224 = qJD(3) * t184;
t221 = qJD(5) * t186;
t220 = t185 * qJD(3);
t219 = t189 * qJD(3);
t218 = t242 * pkin(4);
t217 = pkin(4) * t221;
t174 = -pkin(3) * t185 - pkin(2);
t106 = -pkin(9) * t204 + t128;
t95 = -t187 * t220 - t164 * t222 + (t223 * t238 - t219) * t184;
t191 = t148 * pkin(9) + t95;
t208 = t189 * t212;
t127 = -t187 * t164 - t208;
t198 = -t159 * pkin(9) + t127;
t193 = t242 * t198;
t94 = -qJD(4) * t208 + t185 * t219 + (-qJD(4) * t164 - t224) * t187;
t84 = -pkin(9) * t149 + t94;
t17 = qJD(5) * t193 - t106 * t221 + t186 * t191 + t242 * t84;
t195 = t186 * t198;
t213 = qJD(5) * t242;
t18 = qJD(5) * t195 + t106 * t213 + t186 * t84 - t191 * t242;
t63 = t186 * t106 - t193;
t64 = t106 * t242 + t195;
t215 = t17 * t64 + t18 * t63;
t200 = -t140 * t242 + t141 * t186;
t47 = qJD(5) * t200 + t107 * t242 + t108 * t186;
t48 = qJD(5) * t101 + t107 * t186 - t108 * t242;
t14 = t48 * mrSges(6,1) + mrSges(6,2) * t47;
t31 = t74 * mrSges(6,1) + mrSges(6,2) * t73;
t13 = mrSges(7,1) * t48 - t47 * mrSges(7,3);
t30 = mrSges(7,1) * t74 - t73 * mrSges(7,3);
t66 = -t108 * mrSges(5,1) + mrSges(5,2) * t107;
t210 = t118 * t217;
t209 = pkin(4) * t213;
t85 = -pkin(4) * t108 + t152;
t121 = pkin(4) * t140 + t162;
t38 = -mrSges(7,1) * t226 + mrSges(7,2) * t47;
t207 = Ifges(4,1) * t185 - t233;
t206 = -t184 * Ifges(4,2) + t232;
t205 = -Ifges(4,5) * t185 + Ifges(4,6) * t184;
t136 = pkin(4) * t204 + t174;
t203 = t226 * t258 - t257 * t48 - t259 * t47;
t26 = -t186 * t65 + t242 * t56;
t5 = t186 * t20 + t213 * t56 - t221 * t65 + t23 * t242;
t197 = t255 + t260 * t18 + (-mrSges(6,2) + mrSges(7,3)) * t17;
t2 = qJ(6) * t226 - qJD(6) * t190 + t5;
t3 = -pkin(5) * t226 - t6;
t196 = mrSges(6,1) * t6 - t3 * mrSges(7,1) - t5 * mrSges(6,2) + mrSges(7,3) * t2 - t203;
t169 = t209 + qJD(6);
t192 = -mrSges(6,2) * t209 + t169 * mrSges(7,3) + t217 * t260;
t183 = qJD(6) * mrSges(7,3);
t175 = -t218 - pkin(5);
t173 = qJ(6) + t240;
t161 = -mrSges(4,1) * t190 - mrSges(4,3) * t229;
t160 = mrSges(4,2) * t190 - mrSges(4,3) * t231;
t151 = (mrSges(4,1) * t188 - mrSges(4,3) * t228) * qJD(2);
t150 = (-mrSges(4,2) * t188 - mrSges(4,3) * t230) * qJD(2);
t143 = t148 * mrSges(5,2);
t133 = -pkin(7) * t230 + t155;
t132 = (t188 * Ifges(4,5) + t190 * t207) * qJD(2);
t131 = (t188 * Ifges(4,6) + t190 * t206) * qJD(2);
t130 = -mrSges(5,1) * t190 + mrSges(5,3) * t141;
t129 = mrSges(5,2) * t190 - mrSges(5,3) * t140;
t125 = -t185 * t216 + t137;
t123 = Ifges(5,1) * t159 - Ifges(5,4) * t204;
t122 = Ifges(5,4) * t159 - Ifges(5,2) * t204;
t116 = -Ifges(5,1) * t148 - Ifges(5,4) * t149;
t115 = -Ifges(5,4) * t148 - Ifges(5,2) * t149;
t114 = t149 * mrSges(5,1) - t143;
t98 = -Ifges(5,1) * t141 - Ifges(5,4) * t140 - Ifges(5,5) * t190;
t97 = -Ifges(5,4) * t141 - Ifges(5,2) * t140 - Ifges(5,6) * t190;
t91 = -mrSges(5,2) * t226 + mrSges(5,3) * t108;
t90 = mrSges(5,1) * t226 - mrSges(5,3) * t107;
t87 = mrSges(6,2) * t190 + mrSges(6,3) * t200;
t86 = mrSges(7,2) * t200 - mrSges(7,3) * t190;
t80 = Ifges(6,1) * t118 + Ifges(6,4) * t199;
t79 = Ifges(7,1) * t118 - Ifges(7,5) * t199;
t78 = Ifges(6,4) * t118 + Ifges(6,2) * t199;
t77 = Ifges(7,5) * t118 - Ifges(7,3) * t199;
t76 = -mrSges(6,1) * t199 + mrSges(6,2) * t118;
t75 = -mrSges(7,1) * t199 - mrSges(7,3) * t118;
t62 = -pkin(5) * t199 - qJ(6) * t118 + t136;
t60 = -mrSges(6,1) * t200 + mrSges(6,2) * t101;
t59 = -mrSges(7,1) * t200 - mrSges(7,3) * t101;
t58 = Ifges(5,1) * t107 + Ifges(5,4) * t108 + Ifges(5,5) * t226;
t57 = Ifges(5,4) * t107 + Ifges(5,2) * t108 + Ifges(5,6) * t226;
t54 = Ifges(6,1) * t101 + Ifges(6,4) * t200 - Ifges(6,5) * t190;
t53 = Ifges(7,1) * t101 - Ifges(7,4) * t190 - Ifges(7,5) * t200;
t52 = Ifges(6,4) * t101 + Ifges(6,2) * t200 - Ifges(6,6) * t190;
t51 = Ifges(7,5) * t101 - Ifges(7,6) * t190 - Ifges(7,3) * t200;
t49 = -pkin(5) * t200 - qJ(6) * t101 + t121;
t39 = -mrSges(6,2) * t226 - mrSges(6,3) * t48;
t37 = mrSges(6,1) * t226 - mrSges(6,3) * t47;
t36 = -mrSges(7,2) * t48 + mrSges(7,3) * t226;
t35 = Ifges(6,1) * t73 - Ifges(6,4) * t74;
t34 = Ifges(7,1) * t73 + Ifges(7,5) * t74;
t33 = Ifges(6,4) * t73 - Ifges(6,2) * t74;
t32 = Ifges(7,5) * t73 + Ifges(7,3) * t74;
t25 = pkin(5) * t190 - t26;
t24 = -qJ(6) * t190 + t237;
t21 = pkin(5) * t74 - qJ(6) * t73 - qJD(6) * t118 + t241;
t12 = Ifges(6,1) * t47 - Ifges(6,4) * t48 + Ifges(6,5) * t226;
t11 = Ifges(7,1) * t47 + Ifges(7,4) * t226 + Ifges(7,5) * t48;
t10 = Ifges(6,4) * t47 - Ifges(6,2) * t48 + Ifges(6,6) * t226;
t9 = Ifges(7,5) * t47 + Ifges(7,6) * t226 + Ifges(7,3) * t48;
t7 = pkin(5) * t48 - qJ(6) * t47 - qJD(6) * t101 + t85;
t1 = [0.2e1 * t152 * (mrSges(5,1) * t140 - mrSges(5,2) * t141) + (t121 * t85 + t237 * t5 + t26 * t6) * t251 + 0.2e1 * t237 * t39 - (t9 - t10) * t200 + (t203 + t254) * t190 + (t2 * t24 + t25 * t3 + t49 * t7) * t250 + (t152 * t162 + t28 * t82 + t29 * t81) * t252 + (t124 * t133 + t125 * t134) * t253 + (t11 + t12) * t101 + ((mrSges(3,2) * t249 + 0.2e1 * (Ifges(3,4) + t205) * t190) * t190 + (-Ifges(5,5) * t141 - Ifges(5,6) * t140 + mrSges(3,1) * t249 + (-0.2e1 * Ifges(3,4) - t205) * t188 + t259 * t101 - t257 * t200 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - Ifges(5,3) + pkin(7) ^ 2 * t253 + t185 * t207 - t184 * t206 + (mrSges(4,1) * t184 + t234) * t248 + t258) * t190) * t188) * qJD(2) + (t51 - t52) * t48 + (t53 + t54) * t47 + 0.2e1 * t24 * t36 + 0.2e1 * t26 * t37 + 0.2e1 * t25 * t38 + (-t131 * t184 + t132 * t185 + t142 * t248) * t188 + 0.2e1 * t49 * t13 + 0.2e1 * t7 * t59 + 0.2e1 * t85 * t60 + 0.2e1 * t2 * t86 + 0.2e1 * t5 * t87 + 0.2e1 * t6 * t88 + 0.2e1 * t3 * t89 + 0.2e1 * t81 * t90 + 0.2e1 * t82 * t91 + t107 * t98 + t108 * t97 + 0.2e1 * t121 * t14 + 0.2e1 * t28 * t129 + 0.2e1 * t29 * t130 - t140 * t57 - t141 * t58 + 0.2e1 * t134 * t150 + 0.2e1 * t133 * t151 + 0.2e1 * t125 * t160 + 0.2e1 * t124 * t161 + 0.2e1 * t162 * t66; m(6) * (t121 * t241 + t136 * t85 + t17 * t237 - t18 * t26 + t5 * t64 - t6 * t63) + m(4) * (-t133 * t224 + t134 * t220 + (-t124 * t184 + t125 * t185) * qJ(3)) - (t32 / 0.2e1 - t33 / 0.2e1) * t200 + (-t118 * t6 + t199 * t5 - t237 * t74 - t26 * t73) * mrSges(6,3) + (t118 * t3 + t199 * t2 - t24 * t74 + t25 * t73) * mrSges(7,2) + ((-Ifges(3,6) + Ifges(4,5) * t184 / 0.2e1 + Ifges(4,6) * t244 + Ifges(5,5) * t245 + Ifges(5,6) * t246 + pkin(7) * mrSges(3,2) + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t118 - (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t199) * t188 + (Ifges(3,5) + (Ifges(4,1) * t184 + t232) * t244 - t184 * (Ifges(4,2) * t185 + t233) / 0.2e1 + (-m(4) * pkin(2) - mrSges(4,1) * t185 + mrSges(4,2) * t184 - mrSges(3,1)) * pkin(7)) * t190) * qJD(2) - (t9 / 0.2e1 - t10 / 0.2e1) * t199 + (t148 * t81 - t149 * t82 - t159 * t29 - t204 * t28) * mrSges(5,3) + t152 * (mrSges(5,1) * t204 + mrSges(5,2) * t159) + (t36 + t39) * t64 + (t38 - t37) * t63 - (-pkin(4) * t60 + t97 / 0.2e1) * t149 + (t34 / 0.2e1 + t35 / 0.2e1) * t101 + t256 * t18 - (t227 + t255) * t190 / 0.2e1 + t58 * t245 + t57 * t246 + (t86 + t87) * t17 + (qJD(3) * t160 + qJ(3) * t150 + t125 * mrSges(4,3) + t131 / 0.2e1) * t185 + (-qJD(3) * t161 - qJ(3) * t151 - t124 * mrSges(4,3) + t132 / 0.2e1) * t184 + (t11 / 0.2e1 + t12 / 0.2e1) * t118 + m(5) * (t127 * t29 + t128 * t28 + t152 * t174 + t81 * t95 + t82 * t94) + m(7) * (t17 * t24 + t18 * t25 + t2 * t64 + t21 * t49 + t3 * t63 + t62 * t7) + (t77 / 0.2e1 - t78 / 0.2e1) * t48 + (t79 / 0.2e1 + t80 / 0.2e1) * t47 + (t53 / 0.2e1 + t54 / 0.2e1) * t73 + (t51 / 0.2e1 - t52 / 0.2e1) * t74 + t49 * t30 + t21 * t59 + t62 * t13 + t7 * t75 + t85 * t76 + t121 * t31 + t108 * t122 / 0.2e1 + t107 * t123 / 0.2e1 + t127 * t90 + t128 * t91 + t94 * t129 + t95 * t130 + t136 * t14 - t140 * t115 / 0.2e1 - t141 * t116 / 0.2e1 - pkin(2) * t142 - t148 * t98 / 0.2e1 + t162 * t114 + t174 * t66; 0.2e1 * t174 * t114 - t204 * t115 + t159 * t116 - t148 * t123 + 0.2e1 * t136 * t31 + 0.2e1 * t21 * t75 + 0.2e1 * t62 * t30 + (t77 - t78) * t74 + (t79 + t80) * t73 - (-0.2e1 * pkin(4) * t76 + t122) * t149 + (t34 + t35) * t118 - (t32 - t33) * t199 + (t136 * t241 + t215) * t251 + (t127 * t95 + t128 * t94) * t252 + (t21 * t62 + t215) * t250 + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (t118 * t18 + t17 * t199 + t63 * t73 - t64 * t74) + 0.2e1 * (t127 * t148 - t128 * t149 - t159 * t95 - t204 * t94) * mrSges(5,3) + (qJ(3) * t253 + 0.2e1 * mrSges(4,3)) * qJD(3) * (t184 ^ 2 + t185 ^ 2); m(4) * t179 + m(5) * t152 + m(6) * t85 + m(7) * t7 + t13 + t14 + t142 + t66; m(7) * t21 - t143 - (-mrSges(5,1) - t247) * t149 + t30 + t31; 0; t196 + t87 * t209 + (t242 * t6 + t186 * t5 + (-t186 * t26 + t237 * t242) * qJD(5)) * t247 + m(7) * (t169 * t24 + t173 * t2 + t175 * t3) + t39 * t240 + t37 * t218 - t28 * mrSges(5,2) + t29 * mrSges(5,1) + t169 * t86 + t173 * t36 + t175 * t38 + (m(7) * t25 + t256) * t217 - t254; t197 + (-t242 * t18 + t17 * t186 + (t186 * t63 + t242 * t64) * qJD(5)) * t247 + m(7) * (t169 * t64 + t17 * t173 + t175 * t18 + t217 * t63) + t175 * t239 - t94 * mrSges(5,2) + t95 * mrSges(5,1) + t227 + (t169 * t199 - t173 * t74 + t210) * mrSges(7,2) + (t199 * t209 - t218 * t73 - t240 * t74 + t210) * mrSges(6,3); 0; 0.2e1 * m(7) * (t169 * t173 + t175 * t217) + 0.2e1 * t192; t196 + m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t24) + qJ(6) * t36 - pkin(5) * t38 + qJD(6) * t86; m(7) * (-pkin(5) * t18 + qJ(6) * t17 + qJD(6) * t64) + (-pkin(5) * t73 - qJ(6) * t74 + qJD(6) * t199) * mrSges(7,2) + t197; 0; m(7) * (-pkin(5) * t217 + qJ(6) * t169 + qJD(6) * t173) + t183 + t192; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t183; m(7) * t3 + t38; m(7) * t18 + t239; 0; m(7) * t217; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
