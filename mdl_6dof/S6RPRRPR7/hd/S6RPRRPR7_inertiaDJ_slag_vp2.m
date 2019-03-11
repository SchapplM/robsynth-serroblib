% Calculate time derivative of joint inertia matrix for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:20:00
% EndTime: 2019-03-09 05:20:04
% DurationCPUTime: 2.70s
% Computational Cost: add. (5193->293), mult. (10358->446), div. (0->0), fcn. (9833->8), ass. (0->145)
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t164 = t114 ^ 2 + t117 ^ 2;
t113 = sin(pkin(10));
t171 = cos(pkin(10));
t115 = sin(qJ(4));
t116 = sin(qJ(3));
t118 = cos(qJ(4));
t119 = cos(qJ(3));
t89 = -t115 * t119 - t118 * t116;
t165 = t118 * t119;
t90 = -t115 * t116 + t165;
t126 = t113 * t90 - t171 * t89;
t162 = qJD(4) * t115;
t163 = qJD(3) * t116;
t212 = qJD(3) + qJD(4);
t68 = -t115 * t163 - t116 * t162 + t212 * t165;
t69 = t212 * t89;
t127 = t113 * t69 + t171 * t68;
t195 = t127 * t126;
t43 = t113 * t68 - t171 * t69;
t64 = t113 * t89 + t171 * t90;
t197 = t43 * t64;
t224 = t195 - t197;
t167 = t113 * t115;
t180 = pkin(3) * qJD(4);
t76 = (t118 * t171 - t167) * t180;
t193 = t126 * t76;
t150 = t171 * t115;
t75 = (t113 * t118 + t150) * t180;
t194 = t64 * t75;
t107 = pkin(3) * t118 + pkin(4);
t77 = -pkin(3) * t167 + t107 * t171;
t78 = pkin(3) * t150 + t113 * t107;
t223 = -t127 * t78 + t43 * t77 - t193 + t194;
t188 = mrSges(7,1) * t117;
t94 = mrSges(7,2) * t114 - t188;
t217 = -mrSges(6,1) + t94;
t203 = pkin(1) + pkin(7);
t159 = pkin(8) + t203;
t146 = qJD(4) * t159;
t147 = t159 * t119;
t222 = qJD(3) * t147 + t119 * t146;
t221 = t116 * t146 + t159 * t163;
t160 = qJD(6) * t117;
t130 = -t114 * t43 + t64 * t160;
t220 = qJD(3) * t119;
t219 = t164 * mrSges(7,3) - mrSges(6,2);
t215 = t164 * t76;
t161 = qJD(6) * t114;
t156 = t64 * t161;
t173 = t117 * t43;
t129 = t156 + t173;
t15 = mrSges(7,1) * t127 + mrSges(7,3) * t129;
t16 = -mrSges(7,2) * t127 - mrSges(7,3) * t130;
t214 = -t114 * t15 + t117 * t16;
t213 = -t116 * mrSges(4,1) - mrSges(4,2) * t119;
t211 = -t115 * t68 - t118 * t69 + (t115 * t90 + t118 * t89) * qJD(4);
t210 = 2 * m(6);
t209 = 2 * m(7);
t49 = t222 * t115 + t221 * t118;
t120 = -t69 * qJ(5) - t90 * qJD(5) + t49;
t48 = t221 * t115 - t222 * t118;
t21 = -t68 * qJ(5) + t89 * qJD(5) + t48;
t10 = t113 * t21 - t120 * t171;
t208 = 0.2e1 * t10;
t148 = t116 * t159;
t70 = t115 * t148 - t118 * t147;
t122 = -t90 * qJ(5) + t70;
t71 = -t115 * t147 - t118 * t148;
t57 = qJ(5) * t89 + t71;
t25 = t113 * t57 - t122 * t171;
t207 = 0.2e1 * t25;
t98 = pkin(3) * t220 + qJD(2);
t58 = pkin(4) * t68 + t98;
t206 = 0.2e1 * t58;
t205 = m(6) * pkin(4);
t202 = pkin(4) * t113;
t201 = t10 * t25;
t11 = t113 * t120 + t171 * t21;
t26 = t113 * t122 + t171 * t57;
t104 = t116 * pkin(3) + qJ(2);
t72 = -pkin(4) * t89 + t104;
t27 = pkin(5) * t126 - pkin(9) * t64 + t72;
t14 = t114 * t27 + t117 * t26;
t169 = qJD(6) * t14;
t17 = pkin(5) * t127 + pkin(9) * t43 + t58;
t3 = -t11 * t114 + t117 * t17 - t169;
t200 = t114 * t3;
t13 = -t114 * t26 + t117 * t27;
t170 = qJD(6) * t13;
t2 = t11 * t117 + t114 * t17 + t170;
t199 = t117 * t2;
t198 = t25 * t75;
t196 = t127 * mrSges(6,3);
t192 = t68 * t89;
t191 = t69 * t90;
t73 = -pkin(5) - t77;
t138 = mrSges(7,1) * t114 + mrSges(7,2) * t117;
t91 = t138 * qJD(6);
t190 = t73 * t91;
t189 = -Ifges(7,5) * t173 + Ifges(7,3) * t127;
t186 = mrSges(5,2) * t118;
t183 = Ifges(7,4) * t114;
t182 = Ifges(7,4) * t117;
t181 = Ifges(7,6) * t114;
t154 = t171 * pkin(4);
t103 = -t154 - pkin(5);
t179 = t103 * t91;
t176 = t114 * t64;
t172 = t117 * t64;
t74 = pkin(9) + t78;
t168 = qJD(6) * t74;
t158 = 2 * mrSges(5,3);
t157 = mrSges(7,3) * t170;
t153 = mrSges(6,1) * t127 - mrSges(6,2) * t43;
t151 = -(2 * Ifges(6,4)) - t181;
t102 = pkin(9) + t202;
t149 = t102 * t164;
t145 = -t10 * t64 + t25 * t43;
t143 = -t191 + t192;
t137 = Ifges(7,1) * t117 - t183;
t136 = -Ifges(7,2) * t114 + t182;
t135 = Ifges(7,5) * t114 + Ifges(7,6) * t117;
t134 = -t114 * t13 + t117 * t14;
t33 = -mrSges(7,2) * t126 - mrSges(7,3) * t176;
t34 = mrSges(7,1) * t126 - mrSges(7,3) * t172;
t133 = -t114 * t34 + t117 * t33;
t92 = t136 * qJD(6);
t93 = t137 * qJD(6);
t95 = Ifges(7,2) * t117 + t183;
t96 = Ifges(7,1) * t114 + t182;
t128 = t114 * t93 + t117 * t92 + t96 * t160 - t161 * t95;
t125 = t48 * t89 - t49 * t90 - t68 * t71 - t69 * t70;
t124 = t69 * mrSges(5,1) - t68 * mrSges(5,2) + t219 * t127 + t217 * t43 - t64 * t91;
t123 = -t200 + t199 + (-t114 * t14 - t117 * t13) * qJD(6);
t108 = Ifges(7,5) * t160;
t22 = Ifges(7,6) * t126 + t136 * t64;
t23 = Ifges(7,5) * t126 + t137 * t64;
t6 = -Ifges(7,4) * t129 - Ifges(7,2) * t130 + Ifges(7,6) * t127;
t7 = -Ifges(7,1) * t129 - Ifges(7,4) * t130 + Ifges(7,5) * t127;
t121 = -t48 * mrSges(5,2) - t11 * mrSges(6,2) + mrSges(7,3) * t199 + t23 * t160 / 0.2e1 + t25 * t91 - t96 * t173 / 0.2e1 + t114 * t7 / 0.2e1 - Ifges(6,5) * t43 + t49 * mrSges(5,1) + t117 * t6 / 0.2e1 - t92 * t176 / 0.2e1 + t93 * t172 / 0.2e1 + t126 * (-Ifges(7,6) * t161 + t108) / 0.2e1 - Ifges(5,6) * t68 + Ifges(5,5) * t69 - t130 * t95 / 0.2e1 - (t64 * t96 + t22) * t161 / 0.2e1 + (t135 / 0.2e1 - Ifges(6,6)) * t127 + t217 * t10;
t32 = t138 * t64;
t12 = mrSges(7,1) * t130 - mrSges(7,2) * t129;
t1 = [0.2e1 * t98 * (-mrSges(5,1) * t89 + mrSges(5,2) * t90) + 0.2e1 * t104 * (mrSges(5,1) * t68 + mrSges(5,2) * t69) - 0.2e1 * Ifges(5,2) * t192 + 0.2e1 * Ifges(5,1) * t191 + 0.2e1 * t72 * t153 + t32 * t208 + 0.2e1 * t2 * t33 + 0.2e1 * t3 * t34 + t12 * t207 + 0.2e1 * t13 * t15 + 0.2e1 * t14 * t16 - 0.2e1 * t26 * t196 - (mrSges(6,3) * t207 - t114 * t22 + t117 * t23) * t43 + 0.2e1 * m(5) * (t104 * t98 + t48 * t71 + t49 * t70) + (t11 * t26 + t58 * t72 + t201) * t210 + (t13 * t3 + t14 * t2 + t201) * t209 + 0.2e1 * (-t68 * t90 + t69 * t89) * Ifges(5,4) + t125 * t158 + (mrSges(6,1) * t206 - 0.2e1 * t11 * mrSges(6,3) - t151 * t43 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t127 + t189) * t126 + (mrSges(6,2) * t206 + mrSges(6,3) * t208 - 0.2e1 * Ifges(6,1) * t43 - t114 * t6 + t117 * t7 + (Ifges(7,5) * t117 + t151) * t127 + (-t114 * t23 - t117 * t22 - t126 * t135) * qJD(6)) * t64 + 0.2e1 * (qJ(2) * (mrSges(4,1) * t119 - mrSges(4,2) * t116) + (t116 ^ 2 - t119 ^ 2) * Ifges(4,4)) * qJD(3) + 0.2e1 * (mrSges(3,3) + (m(3) + m(4)) * qJ(2) - t213) * qJD(2) + 0.2e1 * (-Ifges(4,1) + Ifges(4,2)) * t116 * t220; -t64 * t12 + t43 * t32 + t133 * t127 + t143 * t158 + ((-t114 * t33 - t117 * t34) * qJD(6) + t214) * t126 + m(7) * (t123 * t126 + t127 * t134 + t145) + m(6) * (t11 * t126 + t127 * t26 + t145) - m(5) * t125 - 0.2e1 * t224 * mrSges(6,3); 0.2e1 * m(7) * (t164 * t195 - t197) + 0.2e1 * m(6) * t224 - 0.2e1 * m(5) * t143; t121 + m(6) * (-t10 * t77 + t11 * t78 + t26 * t76 + t198) + ((mrSges(4,2) * t203 - Ifges(4,6)) * t119 + (mrSges(4,1) * t203 - Ifges(4,5)) * t116) * qJD(3) + m(7) * (t10 * t73 + t198) + t223 * mrSges(6,3) + t73 * t12 + t75 * t32 + (-t157 - t34 * t168 + m(7) * (-t13 * t168 + t14 * t76 + t2 * t74) + t76 * t33 + t74 * t16) * t117 + (m(5) * (t115 * t48 + t118 * t49 + (-t115 * t70 + t118 * t71) * qJD(4)) + t211 * mrSges(5,3)) * pkin(3) + (-t33 * t168 + m(7) * (-t13 * t76 - t14 * t168 - t3 * t74) - t76 * t34 - t74 * t15 + (-t3 - t169) * mrSges(7,3)) * t114; t213 * qJD(3) + m(7) * (t43 * t73 - t194 + t164 * (t127 * t74 + t193)) - m(6) * t223 - m(5) * t211 * pkin(3) + t124; 0.2e1 * t190 - 0.2e1 * t76 * mrSges(6,2) + (-t75 * t77 + t76 * t78) * t210 + (t74 * t215 + t73 * t75) * t209 + t128 + 0.2e1 * t217 * t75 + 0.2e1 * (-mrSges(5,1) * t115 - t186) * t180 + 0.2e1 * mrSges(7,3) * t215; t121 + (-t10 * t171 + t11 * t113) * t205 - t117 * t157 - t196 * t202 + t43 * mrSges(6,3) * t154 + (m(7) * t10 + t12) * t103 + (-t14 * t161 - t200) * mrSges(7,3) + (m(7) * t123 - t160 * t34 - t161 * t33 + t214) * t102; m(7) * (t103 * t43 + t127 * t149) + (t113 * t127 - t171 * t43) * t205 + t124; -pkin(3) * mrSges(5,1) * t162 - t180 * t186 + t128 + t179 + t190 + (m(7) * t103 - t171 * t205 + t217) * t75 + (m(7) * t149 + t113 * t205 + t219) * t76; t128 + 0.2e1 * t179; t114 * t16 + t117 * t15 + t133 * qJD(6) + m(7) * (qJD(6) * t134 + t114 * t2 + t117 * t3) + m(6) * t58 + t153; 0; 0; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t156 - Ifges(7,6) * t130 + t189; (-t117 * t127 + t126 * t161) * mrSges(7,2) + (-t114 * t127 - t126 * t160) * mrSges(7,1); t108 - t138 * t76 + (-t74 * t188 + (mrSges(7,2) * t74 - Ifges(7,6)) * t114) * qJD(6); t108 + (t102 * t94 - t181) * qJD(6); -t91; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
