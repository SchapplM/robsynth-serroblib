% Calculate time derivative of joint inertia matrix for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:21
% EndTime: 2019-03-09 09:54:28
% DurationCPUTime: 3.45s
% Computational Cost: add. (2827->425), mult. (7025->582), div. (0->0), fcn. (5895->6), ass. (0->165)
t226 = Ifges(7,4) - Ifges(5,6) + Ifges(6,5);
t225 = Ifges(7,5) + Ifges(5,5) - Ifges(6,4);
t227 = -Ifges(6,1) - Ifges(7,1) - Ifges(5,3);
t224 = 2 * mrSges(6,1) + 2 * mrSges(5,3);
t160 = sin(pkin(9));
t161 = cos(pkin(9));
t212 = cos(qJ(4));
t184 = qJD(4) * t212;
t163 = sin(qJ(4));
t190 = qJD(4) * t163;
t223 = -t160 * t190 + t161 * t184;
t131 = t212 * t160 + t163 * t161;
t122 = t131 * qJD(4);
t222 = t226 * t122 + t225 * t223;
t208 = pkin(8) + qJ(3);
t183 = t208 * t160;
t129 = t212 * t183;
t138 = t208 * t161;
t182 = t212 * qJD(3);
t192 = qJD(3) * t160;
t38 = (qJD(4) * t138 + t192) * t163 + qJD(4) * t129 - t161 * t182;
t221 = 2 * m(4);
t220 = 2 * m(5);
t219 = 2 * m(6);
t218 = 2 * m(7);
t217 = -2 * pkin(1);
t216 = 2 * pkin(7);
t215 = -2 * mrSges(7,1);
t214 = t161 / 0.2e1;
t210 = -mrSges(6,2) + mrSges(5,1);
t209 = pkin(4) + qJ(6);
t164 = sin(qJ(2));
t165 = cos(qJ(2));
t135 = -pkin(2) * t165 - t164 * qJ(3) - pkin(1);
t128 = t161 * t135;
t201 = t161 * t164;
t80 = -pkin(8) * t201 + t128 + (-pkin(7) * t160 - pkin(3)) * t165;
t200 = t161 * t165;
t105 = pkin(7) * t200 + t160 * t135;
t203 = t160 * t164;
t93 = -pkin(8) * t203 + t105;
t25 = t163 * t80 + t212 * t93;
t207 = Ifges(4,4) * t160;
t206 = Ifges(4,4) * t161;
t193 = qJD(2) * t165;
t155 = pkin(7) * t193;
t187 = t160 * t193;
t125 = pkin(3) * t187 + t155;
t205 = t125 * mrSges(5,1);
t204 = t125 * mrSges(5,2);
t202 = t160 * t165;
t186 = t161 * t193;
t110 = mrSges(4,1) * t187 + mrSges(4,2) * t186;
t120 = -t164 * qJD(3) + (pkin(2) * t164 - qJ(3) * t165) * qJD(2);
t194 = qJD(2) * t164;
t189 = pkin(7) * t194;
t91 = t161 * t120 + t160 * t189;
t134 = pkin(3) * t203 + t164 * pkin(7);
t195 = (qJ(5) * qJD(5));
t191 = qJD(3) * t161;
t149 = -pkin(3) * t161 - pkin(2);
t188 = t212 * t161;
t59 = t164 * t122 + t163 * t187 - t212 * t186;
t60 = t131 * t193 + t223 * t164;
t19 = t59 * mrSges(7,2) + t60 * mrSges(7,3);
t34 = -t59 * mrSges(6,1) + mrSges(6,2) * t194;
t33 = -t60 * mrSges(7,1) + mrSges(7,2) * t194;
t71 = -mrSges(7,2) * t223 + t122 * mrSges(7,3);
t94 = t138 * t163 + t129;
t21 = qJ(5) * t165 - t25;
t24 = -t163 * t93 + t212 * t80;
t31 = -t59 * mrSges(7,1) - mrSges(7,3) * t194;
t39 = t138 * t184 + t163 * t191 + (-t208 * t190 + t182) * t160;
t95 = t212 * t138 - t163 * t183;
t177 = -t38 * t95 + t39 * t94;
t109 = -t163 * t203 + t164 * t188;
t176 = -qJ(5) * t109 + t134;
t175 = t161 * Ifges(4,1) - t207;
t174 = -t160 * Ifges(4,2) + t206;
t173 = -Ifges(4,5) * t161 + Ifges(4,6) * t160;
t172 = -qJ(5) * t223 - qJD(5) * t131;
t130 = t160 * t163 - t188;
t171 = -qJ(5) * t122 - qJD(5) * t130;
t22 = t165 * pkin(4) - t24;
t169 = -qJ(5) * t131 + t149;
t46 = (pkin(3) * t164 - pkin(8) * t200) * qJD(2) + t91;
t106 = t160 * t120;
t65 = t106 + (-pkin(7) * t201 - pkin(8) * t202) * qJD(2);
t168 = t163 * t65 + t93 * t184 + t80 * t190 - t212 * t46;
t6 = t163 * t46 + t80 * t184 - t93 * t190 + t212 * t65;
t167 = qJ(5) * t59 - qJD(5) * t109 + t125;
t166 = t227 * t194 + t225 * t59 - t226 * t60;
t4 = -qJ(5) * t194 + qJD(5) * t165 - t6;
t133 = -mrSges(4,1) * t165 - mrSges(4,3) * t201;
t132 = mrSges(4,2) * t165 - mrSges(4,3) * t203;
t124 = (mrSges(4,1) * t164 - mrSges(4,3) * t200) * qJD(2);
t123 = (-mrSges(4,2) * t164 - mrSges(4,3) * t202) * qJD(2);
t112 = t223 * mrSges(5,2);
t111 = t223 * mrSges(6,3);
t108 = t131 * t164;
t104 = -pkin(7) * t202 + t128;
t103 = (t164 * Ifges(4,5) + t175 * t165) * qJD(2);
t102 = (t164 * Ifges(4,6) + t174 * t165) * qJD(2);
t101 = -mrSges(5,1) * t165 - t109 * mrSges(5,3);
t100 = mrSges(5,2) * t165 - t108 * mrSges(5,3);
t99 = t109 * mrSges(6,1) - mrSges(6,2) * t165;
t98 = -t108 * mrSges(7,1) - mrSges(7,2) * t165;
t97 = t108 * mrSges(6,1) + mrSges(6,3) * t165;
t96 = t109 * mrSges(7,1) + mrSges(7,3) * t165;
t92 = -t161 * t189 + t106;
t90 = Ifges(5,1) * t131 - Ifges(5,4) * t130;
t89 = Ifges(5,4) * t131 - Ifges(5,2) * t130;
t88 = -Ifges(6,2) * t131 + Ifges(6,6) * t130;
t87 = Ifges(7,2) * t130 + Ifges(7,6) * t131;
t86 = -Ifges(6,6) * t131 + Ifges(6,3) * t130;
t85 = Ifges(7,6) * t130 + Ifges(7,3) * t131;
t84 = -mrSges(6,2) * t130 - mrSges(6,3) * t131;
t83 = -mrSges(7,2) * t131 + mrSges(7,3) * t130;
t78 = pkin(4) * t130 + t169;
t77 = Ifges(5,1) * t223 - Ifges(5,4) * t122;
t76 = Ifges(5,4) * t223 - Ifges(5,2) * t122;
t75 = -Ifges(6,2) * t223 + Ifges(6,6) * t122;
t74 = Ifges(7,2) * t122 + Ifges(7,6) * t223;
t73 = -Ifges(6,6) * t223 + Ifges(6,3) * t122;
t72 = Ifges(7,6) * t122 + Ifges(7,3) * t223;
t70 = t122 * mrSges(5,1) + t112;
t69 = -t122 * mrSges(6,2) - t111;
t64 = -mrSges(6,2) * t108 - mrSges(6,3) * t109;
t63 = -mrSges(7,2) * t109 + mrSges(7,3) * t108;
t58 = -t130 * pkin(5) + t95;
t57 = pkin(5) * t131 + t94;
t48 = t59 * mrSges(5,2);
t47 = t59 * mrSges(6,3);
t45 = Ifges(5,1) * t109 - Ifges(5,4) * t108 - Ifges(5,5) * t165;
t44 = Ifges(5,4) * t109 - Ifges(5,2) * t108 - Ifges(5,6) * t165;
t43 = -Ifges(6,4) * t165 - Ifges(6,2) * t109 + Ifges(6,6) * t108;
t42 = -Ifges(7,4) * t165 + Ifges(7,2) * t108 + Ifges(7,6) * t109;
t41 = -Ifges(6,5) * t165 - Ifges(6,6) * t109 + Ifges(6,3) * t108;
t40 = -Ifges(7,5) * t165 + Ifges(7,6) * t108 + Ifges(7,3) * t109;
t36 = -mrSges(5,2) * t194 - mrSges(5,3) * t60;
t35 = mrSges(5,1) * t194 + mrSges(5,3) * t59;
t32 = mrSges(6,1) * t60 - mrSges(6,3) * t194;
t30 = t209 * t130 + t169;
t29 = pkin(4) * t122 + t172;
t28 = pkin(4) * t108 + t176;
t27 = pkin(5) * t223 + t39;
t26 = -pkin(5) * t122 - t38;
t23 = t209 * t108 + t176;
t20 = qJD(6) * t130 + t209 * t122 + t172;
t18 = t60 * mrSges(5,1) - t48;
t17 = -t60 * mrSges(6,2) + t47;
t16 = -Ifges(5,1) * t59 - Ifges(5,4) * t60 + Ifges(5,5) * t194;
t15 = -Ifges(5,4) * t59 - Ifges(5,2) * t60 + Ifges(5,6) * t194;
t14 = Ifges(6,4) * t194 + Ifges(6,2) * t59 + Ifges(6,6) * t60;
t13 = Ifges(7,4) * t194 + Ifges(7,2) * t60 - Ifges(7,6) * t59;
t12 = Ifges(6,5) * t194 + Ifges(6,6) * t59 + Ifges(6,3) * t60;
t11 = Ifges(7,5) * t194 + Ifges(7,6) * t60 - Ifges(7,3) * t59;
t10 = -t108 * pkin(5) - t21;
t9 = t109 * pkin(5) + t165 * qJ(6) + t22;
t8 = pkin(4) * t60 + t167;
t5 = -pkin(4) * t194 + t168;
t3 = qJD(6) * t108 + t209 * t60 + t167;
t2 = -t60 * pkin(5) - t4;
t1 = -t59 * pkin(5) + t165 * qJD(6) - t209 * t194 + t168;
t7 = [(t11 - t14 + t16 + 0.2e1 * t204) * t109 + (t12 + t13 - t15 + 0.2e1 * t205) * t108 + t166 * t165 + (t1 * t9 + t10 * t2 + t23 * t3) * t218 + (t21 * t4 + t22 * t5 + t28 * t8) * t219 + (t104 * t91 + t105 * t92) * t221 + (((mrSges(3,2) * t217) + 0.2e1 * (Ifges(3,4) + t173) * t165) * t165 + ((mrSges(3,1) * t217) + (-0.2e1 * Ifges(3,4) - t173) * t164 + t225 * t109 + t226 * t108 + (t161 * t175 - t160 * t174 + (mrSges(4,1) * t160 + mrSges(4,2) * t161) * t216 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (pkin(7) ^ 2 * t221) - (2 * Ifges(4,3)) + t227) * t165) * t164) * qJD(2) + (t125 * t134 - t168 * t24 + t25 * t6) * t220 - 0.2e1 * t168 * t101 + (-t160 * t102 + t161 * t103 + t110 * t216) * t164 + (t43 - t40 - t45) * t59 + 0.2e1 * t23 * t19 + 0.2e1 * t28 * t17 + 0.2e1 * t9 * t31 + 0.2e1 * t21 * t32 + 0.2e1 * t10 * t33 + 0.2e1 * t22 * t34 + 0.2e1 * t24 * t35 + 0.2e1 * t25 * t36 + 0.2e1 * t3 * t63 + 0.2e1 * t8 * t64 + 0.2e1 * t1 * t96 + 0.2e1 * t4 * t97 + 0.2e1 * t2 * t98 + 0.2e1 * t5 * t99 + 0.2e1 * t6 * t100 + (t41 + t42 - t44) * t60 + 0.2e1 * t105 * t123 + 0.2e1 * t104 * t124 + 0.2e1 * t92 * t132 + 0.2e1 * t91 * t133 + 0.2e1 * t134 * t18; (((pkin(7) * mrSges(3,2)) + Ifges(4,5) * t160 / 0.2e1 + Ifges(4,6) * t214 - Ifges(3,6) + (Ifges(7,5) / 0.2e1 + Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t131 + (Ifges(7,4) / 0.2e1 - Ifges(5,6) / 0.2e1 + Ifges(6,5) / 0.2e1) * t130) * t164 + ((Ifges(4,1) * t160 + t206) * t214 - t160 * (Ifges(4,2) * t161 + t207) / 0.2e1 + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(4,1) * t161 + mrSges(4,2) * t160 - mrSges(3,1)) * pkin(7)) * t165) * qJD(2) + (-t2 * mrSges(7,1) + t4 * mrSges(6,1) - t6 * mrSges(5,3) + t12 / 0.2e1 + t13 / 0.2e1 - t15 / 0.2e1 + t205) * t130 + m(4) * (-t104 * t192 + t105 * t191 + (-t91 * t160 + t92 * t161) * qJ(3)) - (-t22 * mrSges(6,1) - t9 * mrSges(7,1) + t24 * mrSges(5,3) + t43 / 0.2e1 - t40 / 0.2e1 - t45 / 0.2e1) * t223 + (t72 / 0.2e1 + t77 / 0.2e1 - t75 / 0.2e1) * t109 + (-t10 * mrSges(7,1) + t21 * mrSges(6,1) - t25 * mrSges(5,3) + t41 / 0.2e1 + t42 / 0.2e1 - t44 / 0.2e1) * t122 + m(7) * (t1 * t57 + t10 * t26 + t2 * t58 + t20 * t23 + t27 * t9 + t3 * t30) + m(6) * (t21 * t38 + t22 * t39 + t28 * t29 - t4 * t95 + t5 * t94 + t78 * t8) + (t1 * mrSges(7,1) + t168 * mrSges(5,3) + t5 * mrSges(6,1) + t11 / 0.2e1 + t16 / 0.2e1 + t204 - t14 / 0.2e1) * t131 + m(5) * (t125 * t149 + t168 * t94 - t24 * t39 - t25 * t38 + t6 * t95) - t222 * t165 / 0.2e1 + (t86 / 0.2e1 + t87 / 0.2e1 - t89 / 0.2e1) * t60 + (t36 - t32) * t95 + (t34 - t35) * t94 + (t99 - t101) * t39 + (t97 - t100) * t38 + (-t85 / 0.2e1 + t88 / 0.2e1 - t90 / 0.2e1) * t59 + t30 * t19 + t57 * t31 + t58 * t33 + t20 * t63 + t29 * t64 + t28 * t69 + t23 * t71 + t78 * t17 + t3 * t83 + t8 * t84 + t27 * t96 + t26 * t98 - pkin(2) * t110 + t134 * t70 + t149 * t18 + (t73 / 0.2e1 + t74 / 0.2e1 - t76 / 0.2e1) * t108 + (-t91 * mrSges(4,3) - qJ(3) * t124 - qJD(3) * t133 + t103 / 0.2e1) * t160 + (t92 * mrSges(4,3) + qJ(3) * t123 + qJD(3) * t132 + t102 / 0.2e1) * t161; 0.2e1 * t149 * t70 + 0.2e1 * t20 * t83 + 0.2e1 * t29 * t84 + 0.2e1 * t30 * t71 + 0.2e1 * t78 * t69 + (t20 * t30 + t26 * t58 + t27 * t57) * t218 + t177 * t220 + (t29 * t78 + t177) * t219 + (0.2e1 * t27 * mrSges(7,1) + t224 * t39 + t72 - t75 + t77) * t131 + (t26 * t215 + t224 * t38 + t73 + t74 - t76) * t130 + (t58 * t215 - t224 * t95 + t86 + t87 - t89) * t122 - (t57 * t215 - t224 * t94 - t85 + t88 - t90) * t223 + (qJ(3) * t221 + 0.2e1 * mrSges(4,3)) * (t160 ^ 2 + t161 ^ 2) * qJD(3); m(4) * t155 + m(5) * t125 + m(6) * t8 + m(7) * t3 + t210 * t60 + t110 + t19 + t47 - t48; m(6) * t29 + m(7) * t20 + t210 * t122 - t111 + t112 + t71; 0; -t166 + (-t97 + t98) * qJD(5) + (-t32 + t33) * qJ(5) + m(7) * (qJ(5) * t2 + qJD(5) * t10 - qJD(6) * t9 - t1 * t209) + m(6) * (-pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t21) - t1 * mrSges(7,3) + t2 * mrSges(7,2) - t4 * mrSges(6,3) + t5 * mrSges(6,2) - t6 * mrSges(5,2) - t168 * mrSges(5,1) - pkin(4) * t34 - qJD(6) * t96 - t209 * t31; t26 * mrSges(7,2) - t27 * mrSges(7,3) - t210 * t39 + (mrSges(5,2) - mrSges(6,3)) * t38 + m(6) * (-pkin(4) * t39 - qJ(5) * t38 + qJD(5) * t95) + m(7) * (qJ(5) * t26 + qJD(5) * t58 - qJD(6) * t57 - t209 * t27) + (-pkin(4) * t223 + t171) * mrSges(6,1) + (-qJD(6) * t131 - t209 * t223 + t171) * mrSges(7,1) + t222; 0; (2 * m(6) * t195) + 0.2e1 * qJD(6) * mrSges(7,3) + 0.2e1 * m(7) * (qJD(6) * t209 + t195) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJD(5); m(6) * t5 + m(7) * t1 + t31 + t34; -(-mrSges(7,1) - mrSges(6,1)) * t223 + m(6) * t39 + m(7) * t27; 0; -m(7) * qJD(6); 0; m(7) * t2 + t33; m(7) * t26 - t122 * mrSges(7,1); 0; m(7) * qJD(5); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
