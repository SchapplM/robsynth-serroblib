% Calculate time derivative of joint inertia matrix for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:41
% EndTime: 2019-03-09 04:07:47
% DurationCPUTime: 2.65s
% Computational Cost: add. (4822->418), mult. (10772->632), div. (0->0), fcn. (9913->8), ass. (0->170)
t156 = (-pkin(1) - pkin(7));
t211 = 2 * t156;
t148 = sin(pkin(10));
t149 = cos(pkin(10));
t151 = sin(qJ(5));
t154 = cos(qJ(5));
t131 = t148 * t154 + t149 * t151;
t118 = t131 * qJD(5);
t155 = cos(qJ(3));
t159 = t148 * t151 - t149 * t154;
t152 = sin(qJ(3));
t178 = qJD(3) * t152;
t74 = -t118 * t155 + t159 * t178;
t208 = qJD(5) * t159;
t76 = t131 * t178 + t155 * t208;
t38 = -t76 * mrSges(6,1) + t74 * mrSges(6,2);
t150 = sin(qJ(6));
t153 = cos(qJ(6));
t108 = t131 * t155;
t110 = t159 * t155;
t65 = -t108 * t153 + t110 * t150;
t26 = qJD(6) * t65 + t150 * t76 + t153 * t74;
t67 = -t108 * t150 - t110 * t153;
t28 = -qJD(6) * t67 - t150 * t74 + t153 * t76;
t6 = -t28 * mrSges(7,1) + t26 * mrSges(7,2);
t210 = -t38 - t6;
t193 = pkin(8) + qJ(4);
t135 = t193 * t148;
t137 = t193 * t149;
t97 = -t151 * t135 + t154 * t137;
t111 = -qJD(4) * t155 + qJD(2) + (pkin(3) * t155 + qJ(4) * t152) * qJD(3);
t105 = t149 * t111;
t176 = qJD(3) * t156;
t167 = t155 * t176;
t90 = -t148 * t167 + t105;
t91 = t148 * t111 + t149 * t167;
t209 = -t90 * t148 + t91 * t149;
t177 = qJD(3) * t155;
t207 = Ifges(6,5) * t74 + Ifges(6,6) * t76 + Ifges(6,3) * t177;
t206 = 2 * m(5);
t205 = 2 * m(6);
t204 = 2 * m(7);
t203 = 2 * mrSges(4,1);
t147 = t149 ^ 2;
t202 = m(7) * pkin(5);
t87 = -t131 * t150 - t153 * t159;
t201 = t87 / 0.2e1;
t88 = t131 * t153 - t150 * t159;
t200 = t88 / 0.2e1;
t199 = -t159 / 0.2e1;
t198 = t131 / 0.2e1;
t197 = t148 / 0.2e1;
t195 = pkin(3) * t152;
t194 = pkin(5) * t118;
t42 = qJD(6) * t87 - t118 * t150 - t153 * t208;
t43 = -qJD(6) * t88 - t118 * t153 + t150 * t208;
t192 = Ifges(7,5) * t42 + Ifges(7,6) * t43;
t134 = -qJ(4) * t155 + qJ(2) + t195;
t124 = t149 * t134;
t165 = -t148 * t156 + pkin(4);
t182 = t149 * t155;
t85 = -pkin(8) * t182 + t152 * t165 + t124;
t181 = t152 * t156;
t103 = t148 * t134 + t149 * t181;
t184 = t148 * t155;
t92 = -pkin(8) * t184 + t103;
t45 = t151 * t85 + t154 * t92;
t191 = Ifges(5,4) * t148;
t190 = Ifges(5,4) * t149;
t189 = t148 * Ifges(5,2);
t188 = t151 * t92;
t185 = -mrSges(5,1) * t149 + mrSges(5,2) * t148 - mrSges(4,1);
t183 = t149 * t152;
t180 = -Ifges(6,5) * t208 - Ifges(6,6) * t118;
t179 = t148 ^ 2 + t147;
t175 = qJD(4) * t148;
t174 = qJD(4) * t149;
t173 = qJD(5) * t154;
t172 = qJD(6) * t150;
t171 = qJD(6) * t153;
t170 = Ifges(7,5) * t26 + Ifges(7,6) * t28 + Ifges(7,3) * t177;
t143 = -pkin(4) * t149 - pkin(3);
t169 = t148 * t178;
t168 = t152 * t177;
t16 = -mrSges(7,1) * t43 + t42 * mrSges(7,2);
t107 = t131 * t152;
t109 = t159 * t152;
t64 = -t107 * t153 + t109 * t150;
t73 = -qJD(3) * t110 - t118 * t152;
t75 = -qJD(3) * t108 + t152 * t208;
t25 = qJD(6) * t64 + t150 * t75 + t153 * t73;
t66 = -t107 * t150 - t109 * t153;
t27 = -qJD(6) * t66 - t150 * t73 + t153 * t75;
t166 = t27 * mrSges(7,1) - t25 * mrSges(7,2);
t44 = t154 * t85 - t188;
t164 = t179 * mrSges(5,3);
t163 = t179 * qJ(4);
t96 = -t154 * t135 - t137 * t151;
t126 = pkin(4) * t184 - t155 * t156;
t71 = -pkin(9) * t131 + t96;
t72 = -pkin(9) * t159 + t97;
t36 = -t150 * t72 + t153 * t71;
t59 = -t135 * t173 + t154 * t174 + (-qJD(5) * t137 - t175) * t151;
t49 = -pkin(9) * t118 + t59;
t60 = -t131 * qJD(4) - t97 * qJD(5);
t50 = pkin(9) * t208 + t60;
t8 = qJD(6) * t36 + t150 * t50 + t153 * t49;
t37 = t150 * t71 + t153 * t72;
t9 = -qJD(6) * t37 - t150 * t49 + t153 * t50;
t162 = t9 * mrSges(7,1) - t8 * mrSges(7,2) + t192;
t161 = mrSges(5,1) * t148 + mrSges(5,2) * t149;
t160 = -Ifges(5,5) * t149 + Ifges(5,6) * t148;
t29 = pkin(5) * t152 + pkin(9) * t110 + t44;
t34 = -pkin(9) * t108 + t45;
t12 = -t150 * t34 + t153 * t29;
t13 = t150 * t29 + t153 * t34;
t116 = -pkin(4) * t169 + t152 * t176;
t58 = t105 + (pkin(8) * t183 + t155 * t165) * qJD(3);
t77 = pkin(8) * t169 + t91;
t15 = -qJD(5) * t45 - t151 * t77 + t154 * t58;
t10 = pkin(5) * t177 - pkin(9) * t74 + t15;
t14 = -qJD(5) * t188 + t151 * t58 + t154 * t77 + t85 * t173;
t11 = pkin(9) * t76 + t14;
t2 = qJD(6) * t12 + t10 * t150 + t11 * t153;
t3 = -qJD(6) * t13 + t10 * t153 - t11 * t150;
t158 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t170;
t133 = mrSges(5,1) * t152 - mrSges(5,3) * t182;
t132 = -mrSges(5,2) * t152 - mrSges(5,3) * t184;
t125 = (-mrSges(7,1) * t150 - mrSges(7,2) * t153) * qJD(6) * pkin(5);
t121 = (mrSges(5,1) * t155 + mrSges(5,3) * t183) * qJD(3);
t120 = (mrSges(5,3) * t148 * t152 - mrSges(5,2) * t155) * qJD(3);
t119 = t161 * t155;
t113 = t208 * mrSges(6,2);
t112 = t161 * t178;
t106 = pkin(5) * t159 + t143;
t102 = -t148 * t181 + t124;
t101 = (t155 * Ifges(5,5) + (-t149 * Ifges(5,1) + t191) * t152) * qJD(3);
t100 = (t155 * Ifges(5,6) + (t189 - t190) * t152) * qJD(3);
t99 = mrSges(6,1) * t152 + mrSges(6,3) * t110;
t98 = -mrSges(6,2) * t152 - mrSges(6,3) * t108;
t95 = Ifges(6,1) * t131 - Ifges(6,4) * t159;
t94 = Ifges(6,4) * t131 - Ifges(6,2) * t159;
t93 = mrSges(6,1) * t159 + mrSges(6,2) * t131;
t89 = t108 * pkin(5) + t126;
t84 = -Ifges(6,1) * t208 - Ifges(6,4) * t118;
t83 = -Ifges(6,4) * t208 - Ifges(6,2) * t118;
t82 = mrSges(6,1) * t118 - t113;
t78 = mrSges(6,1) * t108 - mrSges(6,2) * t110;
t62 = -Ifges(6,1) * t110 - Ifges(6,4) * t108 + Ifges(6,5) * t152;
t61 = -Ifges(6,4) * t110 - Ifges(6,2) * t108 + Ifges(6,6) * t152;
t57 = -mrSges(6,2) * t177 + mrSges(6,3) * t76;
t56 = mrSges(6,1) * t177 - mrSges(6,3) * t74;
t53 = mrSges(7,1) * t152 - mrSges(7,3) * t67;
t52 = -mrSges(7,2) * t152 + mrSges(7,3) * t65;
t51 = -pkin(5) * t76 + t116;
t48 = Ifges(7,1) * t88 + Ifges(7,4) * t87;
t47 = Ifges(7,4) * t88 + Ifges(7,2) * t87;
t46 = -mrSges(7,1) * t87 + mrSges(7,2) * t88;
t35 = -mrSges(7,1) * t65 + mrSges(7,2) * t67;
t33 = Ifges(6,1) * t74 + Ifges(6,4) * t76 + Ifges(6,5) * t177;
t32 = Ifges(6,4) * t74 + Ifges(6,2) * t76 + Ifges(6,6) * t177;
t31 = Ifges(7,1) * t67 + Ifges(7,4) * t65 + Ifges(7,5) * t152;
t30 = Ifges(7,4) * t67 + Ifges(7,2) * t65 + Ifges(7,6) * t152;
t20 = -mrSges(7,2) * t177 + mrSges(7,3) * t28;
t19 = mrSges(7,1) * t177 - mrSges(7,3) * t26;
t18 = Ifges(7,1) * t42 + Ifges(7,4) * t43;
t17 = Ifges(7,4) * t42 + Ifges(7,2) * t43;
t5 = Ifges(7,1) * t26 + Ifges(7,4) * t28 + Ifges(7,5) * t177;
t4 = Ifges(7,4) * t26 + Ifges(7,2) * t28 + Ifges(7,6) * t177;
t1 = [((qJ(2) * t203 - Ifges(6,5) * t110 + Ifges(7,5) * t67 - Ifges(6,6) * t108 + Ifges(7,6) * t65 + (-(2 * Ifges(4,4)) - t160) * t155) * t155 + (t119 * t211 - 0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * (Ifges(4,4) + t160) * t152 + (-Ifges(5,1) * t147 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + (2 * Ifges(5,3)) + Ifges(7,3) + Ifges(6,3) - (2 * m(5) * t156 ^ 2) + (-t189 + 0.2e1 * t190) * t148) * t155) * t152) * qJD(3) + (-t148 * t100 + t149 * t101 + t112 * t211) * t155 + (t152 * t203 + 0.2e1 * t155 * mrSges(4,2) + (2 * mrSges(3,3)) + 0.2e1 * (m(3) + m(4)) * qJ(2)) * qJD(2) + 0.2e1 * t91 * t132 + 0.2e1 * t90 * t133 + 0.2e1 * t103 * t120 + 0.2e1 * t102 * t121 + 0.2e1 * t126 * t38 - t108 * t32 - t110 * t33 + 0.2e1 * t116 * t78 + 0.2e1 * t14 * t98 + 0.2e1 * t15 * t99 + 0.2e1 * t89 * t6 + t74 * t62 + t76 * t61 + t65 * t4 + t67 * t5 + 0.2e1 * t44 * t56 + 0.2e1 * t45 * t57 + 0.2e1 * t51 * t35 + 0.2e1 * t2 * t52 + 0.2e1 * t3 * t53 + t28 * t30 + t26 * t31 + 0.2e1 * t13 * t20 + 0.2e1 * t12 * t19 + (t170 + t207) * t152 + (t12 * t3 + t13 * t2 + t51 * t89) * t204 + (t116 * t126 + t14 * t45 + t15 * t44) * t205 + (t102 * t90 + t103 * t91) * t206; -t107 * t56 - t109 * t57 + t64 * t19 + t66 * t20 + t25 * t52 + t27 * t53 + t73 * t98 + t75 * t99 + (t149 * t120 - t148 * t121) * t152 + (t112 + t210) * t155 + ((t132 * t149 - t133 * t148) * t155 + (t119 + t35 + t78) * t152) * qJD(3) + m(7) * (t12 * t27 + t13 * t25 - t155 * t51 + t178 * t89 + t2 * t66 + t3 * t64) + m(6) * (-t107 * t15 - t109 * t14 - t116 * t155 + t126 * t178 + t44 * t75 + t45 * t73) + m(5) * (t209 * t152 + (-t102 * t148 + t103 * t149 - 0.2e1 * t181) * t177); 0.2e1 * m(7) * (t25 * t66 + t27 * t64 - t168) + 0.2e1 * m(6) * (-t107 * t75 - t109 * t73 - t168) + 0.2e1 * m(5) * (-0.1e1 + t179) * t168; -t208 * t62 / 0.2e1 + (-t118 * t45 - t131 * t15 - t14 * t159 + t208 * t44) * mrSges(6,3) + (-t12 * t42 + t13 * t43 + t2 * t87 - t3 * t88) * mrSges(7,3) + m(6) * (t116 * t143 + t14 * t97 + t15 * t96 + t44 * t60 + t45 * t59) + (t100 / 0.2e1 + t91 * mrSges(5,3) + qJ(4) * t120 + qJD(4) * t132) * t149 + (t101 / 0.2e1 - t90 * mrSges(5,3) - qJ(4) * t121 - qJD(4) * t133) * t148 - (t61 / 0.2e1 - pkin(5) * t35) * t118 + m(7) * (t106 * t51 + t12 * t9 + t13 * t8 + t194 * t89 + t2 * t37 + t3 * t36) + (t192 + t180) * t152 / 0.2e1 + ((-t156 * mrSges(4,2) - Ifges(4,6) + Ifges(7,5) * t200 + Ifges(7,6) * t201 + Ifges(6,5) * t198 + Ifges(6,6) * t199 + Ifges(5,5) * t197 + Ifges(5,6) * t149 / 0.2e1) * t155 + (-t149 * (Ifges(5,1) * t148 + t190) / 0.2e1 + (Ifges(5,2) * t149 + t191) * t197 - Ifges(4,5) + (-m(5) * pkin(3) + t185) * t156) * t152) * qJD(3) + t143 * t38 + t126 * t82 - t108 * t83 / 0.2e1 - t110 * t84 / 0.2e1 + pkin(3) * t112 + t116 * t93 + t96 * t56 + t97 * t57 + t59 * t98 + t60 * t99 + t106 * t6 + t89 * t16 + t76 * t94 / 0.2e1 + t74 * t95 / 0.2e1 + t65 * t17 / 0.2e1 + t67 * t18 / 0.2e1 + t28 * t47 / 0.2e1 + t26 * t48 / 0.2e1 + t51 * t46 + t8 * t52 + t9 * t53 + t42 * t31 / 0.2e1 + t43 * t30 / 0.2e1 + t36 * t19 + t37 * t20 + m(5) * (t209 * qJ(4) - t102 * t175 + t103 * t174) + t33 * t198 + t32 * t199 + t5 * t200 + t4 * t201; (-t16 - t82) * t155 + ((-mrSges(4,2) + t164) * t155 + (t46 + t93 + t185) * t152) * qJD(3) + m(7) * (t106 * t178 - t155 * t194 + t25 * t37 + t27 * t36 + t64 * t9 + t66 * t8) + m(6) * (-t107 * t60 - t109 * t59 + t143 * t178 + t73 * t97 + t75 * t96) + m(5) * (t179 * t152 * qJD(4) + (t155 * t163 - t195) * qJD(3)) + (t25 * t87 - t27 * t88 - t42 * t64 + t43 * t66) * mrSges(7,3) + (-t107 * t208 + t109 * t118 - t131 * t75 - t159 * t73) * mrSges(6,3); 0.2e1 * t106 * t16 - t208 * t95 - t159 * t83 + t131 * t84 + 0.2e1 * t143 * t82 + t87 * t17 + t88 * t18 + t42 * t48 + t43 * t47 - (-0.2e1 * pkin(5) * t46 + t94) * t118 + (t106 * t194 + t36 * t9 + t37 * t8) * t204 + (t59 * t97 + t60 * t96) * t205 + 0.2e1 * (-t36 * t42 + t37 * t43 + t8 * t87 - t88 * t9) * mrSges(7,3) + 0.2e1 * (-t118 * t97 - t131 * t60 - t159 * t59 + t208 * t96) * mrSges(6,3) + (t163 * t206 + 0.2e1 * t164) * qJD(4); m(6) * t116 + m(7) * t51 + ((m(5) * t156) - t161) * t178 - t210; (m(5) + m(6) + m(7)) * t178; -t113 - (-mrSges(6,1) - t202) * t118 + t16; 0; t15 * mrSges(6,1) - t14 * mrSges(6,2) + (m(7) * (-t12 * t172 + t13 * t171 + t150 * t2 + t153 * t3) + t52 * t171 + t150 * t20 - t53 * t172 + t153 * t19) * pkin(5) + t158 + t207; t75 * mrSges(6,1) - t73 * mrSges(6,2) + (t150 * t25 + t153 * t27 + (-t150 * t64 + t153 * t66) * qJD(6)) * t202 + t166; t60 * mrSges(6,1) - t59 * mrSges(6,2) + (m(7) * (t150 * t8 + t153 * t9 + (-t150 * t36 + t153 * t37) * qJD(6)) + (t150 * t43 - t153 * t42 + (t150 * t88 + t153 * t87) * qJD(6)) * mrSges(7,3)) * pkin(5) + t162 + t180; 0; 0.2e1 * t125; t158; t166; t162; 0; t125; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
