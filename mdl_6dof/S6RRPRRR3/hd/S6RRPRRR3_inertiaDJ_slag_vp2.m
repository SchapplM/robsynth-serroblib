% Calculate time derivative of joint inertia matrix for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:53
% EndTime: 2019-03-09 13:22:05
% DurationCPUTime: 5.05s
% Computational Cost: add. (12624->503), mult. (27313->752), div. (0->0), fcn. (27372->10), ass. (0->199)
t176 = sin(pkin(11));
t177 = cos(pkin(11));
t181 = sin(qJ(2));
t185 = cos(qJ(2));
t157 = t176 * t181 - t177 * t185;
t153 = t157 * qJD(2);
t158 = t176 * t185 + t177 * t181;
t180 = sin(qJ(4));
t184 = cos(qJ(4));
t215 = qJD(4) * t184;
t192 = -t153 * t180 + t158 * t215;
t250 = -mrSges(5,1) * t184 + mrSges(5,2) * t180;
t179 = sin(qJ(5));
t183 = cos(qJ(5));
t194 = t179 * t180 - t183 * t184;
t105 = t194 * t158;
t174 = -pkin(2) * t185 - pkin(1);
t114 = pkin(3) * t157 - pkin(8) * t158 + t174;
t234 = -qJ(3) - pkin(7);
t166 = t234 * t181;
t167 = t234 * t185;
t129 = t166 * t176 - t167 * t177;
t122 = t184 * t129;
t81 = t114 * t180 + t122;
t152 = t158 * qJD(2);
t218 = t184 * t153;
t249 = -Ifges(5,5) * t218 + Ifges(5,3) * t152;
t171 = pkin(2) * t176 + pkin(8);
t235 = pkin(9) + t171;
t154 = t235 * t180;
t155 = t235 * t184;
t111 = -t154 * t179 + t155 * t183;
t198 = mrSges(5,1) * t180 + mrSges(5,2) * t184;
t162 = t198 * qJD(4);
t248 = qJD(4) + qJD(5);
t161 = t179 * t184 + t180 * t183;
t124 = t248 * t161;
t48 = -t124 * t158 + t153 * t194;
t49 = t105 * t248 + t161 * t153;
t247 = Ifges(6,5) * t48 + Ifges(6,6) * t49 + Ifges(6,3) * t152;
t246 = 2 * m(6);
t245 = 2 * m(7);
t244 = -2 * mrSges(4,3);
t128 = -t166 * t177 - t167 * t176;
t242 = 0.2e1 * t128;
t241 = 0.2e1 * t174;
t240 = m(6) * pkin(4);
t237 = -t158 / 0.2e1;
t232 = Ifges(5,4) * t180;
t168 = Ifges(5,2) * t184 + t232;
t236 = -t168 / 0.2e1;
t178 = sin(qJ(6));
t182 = cos(qJ(6));
t117 = -t161 * t178 - t182 * t194;
t123 = t248 * t194;
t62 = qJD(6) * t117 - t123 * t182 - t124 * t178;
t118 = t161 * t182 - t178 * t194;
t63 = -qJD(6) * t118 + t123 * t178 - t124 * t182;
t233 = Ifges(7,5) * t62 + Ifges(7,6) * t63;
t222 = t158 * t184;
t80 = t114 * t184 - t129 * t180;
t64 = pkin(4) * t157 - pkin(9) * t222 + t80;
t223 = t158 * t180;
t70 = -pkin(9) * t223 + t81;
t38 = t179 * t64 + t183 * t70;
t231 = Ifges(5,4) * t184;
t230 = Ifges(5,6) * t180;
t173 = pkin(4) * t183 + pkin(5);
t211 = qJD(6) * t182;
t212 = qJD(6) * t178;
t221 = t178 * t179;
t112 = t173 * t211 + (-t179 * t212 + (t182 * t183 - t221) * qJD(5)) * pkin(4);
t229 = t112 * mrSges(7,2);
t228 = t152 * Ifges(5,5);
t227 = t152 * Ifges(5,6);
t226 = t157 * Ifges(5,6);
t203 = qJD(2) * t234;
t149 = qJD(3) * t185 + t181 * t203;
t188 = -qJD(3) * t181 + t185 * t203;
t100 = t149 * t176 - t177 * t188;
t225 = t100 * t128;
t220 = t179 * t182;
t217 = -Ifges(6,5) * t123 - Ifges(6,6) * t124;
t216 = qJD(4) * t180;
t214 = qJD(5) * t179;
t213 = qJD(5) * t183;
t104 = t161 * t158;
t71 = -t104 * t182 + t105 * t178;
t18 = qJD(6) * t71 + t178 * t49 + t182 * t48;
t72 = -t104 * t178 - t105 * t182;
t19 = -qJD(6) * t72 - t178 * t48 + t182 * t49;
t209 = Ifges(7,5) * t18 + Ifges(7,6) * t19 + Ifges(7,3) * t152;
t208 = pkin(2) * qJD(2) * t181;
t207 = pkin(4) * t216;
t172 = -pkin(2) * t177 - pkin(3);
t205 = t158 * t216;
t32 = -mrSges(7,1) * t63 + t62 * mrSges(7,2);
t204 = -(2 * Ifges(4,4)) - t230;
t37 = -t179 * t70 + t183 * t64;
t202 = qJD(4) * t235;
t201 = 0.2e1 * t208;
t85 = mrSges(6,1) * t124 - t123 * mrSges(6,2);
t113 = -t173 * t212 + (-t179 * t211 + (-t178 * t183 - t220) * qJD(5)) * pkin(4);
t109 = t113 * mrSges(7,1);
t200 = t109 - t229;
t101 = t177 * t149 + t176 * t188;
t102 = pkin(3) * t152 + pkin(8) * t153 + t208;
t199 = -t101 * t180 + t102 * t184;
t110 = -t154 * t183 - t155 * t179;
t99 = pkin(4) * t223 + t128;
t197 = Ifges(5,1) * t184 - t232;
t196 = -Ifges(5,2) * t180 + t231;
t26 = pkin(5) * t157 + pkin(10) * t105 + t37;
t27 = -pkin(10) * t104 + t38;
t12 = -t178 * t27 + t182 * t26;
t13 = t178 * t26 + t182 * t27;
t94 = -pkin(10) * t161 + t110;
t95 = -pkin(10) * t194 + t111;
t52 = -t178 * t95 + t182 * t94;
t53 = t178 * t94 + t182 * t95;
t147 = t180 * t202;
t148 = t184 * t202;
t78 = -t147 * t183 - t148 * t179 - t154 * t213 - t155 * t214;
t50 = -pkin(10) * t124 + t78;
t79 = -qJD(5) * t111 + t147 * t179 - t148 * t183;
t51 = pkin(10) * t123 + t79;
t21 = qJD(6) * t52 + t178 * t51 + t182 * t50;
t22 = -qJD(6) * t53 - t178 * t50 + t182 * t51;
t195 = mrSges(7,1) * t22 - t21 * mrSges(7,2) + t233;
t165 = -pkin(4) * t184 + t172;
t76 = pkin(4) * t192 + t100;
t30 = pkin(9) * t218 + pkin(4) * t152 + (-t122 + (pkin(9) * t158 - t114) * t180) * qJD(4) + t199;
t42 = t101 * t184 + t102 * t180 + t114 * t215 - t129 * t216;
t36 = -pkin(9) * t192 + t42;
t11 = -qJD(5) * t38 - t179 * t36 + t183 * t30;
t4 = pkin(5) * t152 - pkin(10) * t48 + t11;
t10 = t179 * t30 + t183 * t36 + t213 * t64 - t214 * t70;
t5 = pkin(10) * t49 + t10;
t2 = qJD(6) * t12 + t178 * t4 + t182 * t5;
t3 = -qJD(6) * t13 - t178 * t5 + t182 * t4;
t193 = mrSges(7,1) * t3 - t2 * mrSges(7,2) + t209;
t191 = t205 + t218;
t190 = -t32 - t85;
t189 = (-mrSges(6,1) * t179 - mrSges(6,2) * t183) * qJD(5) * pkin(4);
t187 = mrSges(6,1) * t79 - t78 * mrSges(6,2) + t195 + t217;
t186 = mrSges(6,1) * t11 - t10 * mrSges(6,2) + t193 + t247;
t175 = Ifges(5,5) * t215;
t169 = Ifges(5,1) * t180 + t231;
t164 = t197 * qJD(4);
t163 = t196 * qJD(4);
t156 = (-mrSges(7,1) * t178 - mrSges(7,2) * t182) * qJD(6) * pkin(5);
t151 = pkin(4) * t220 + t173 * t178;
t150 = -pkin(4) * t221 + t173 * t182;
t139 = t153 * mrSges(4,2);
t130 = pkin(5) * t194 + t165;
t127 = Ifges(6,1) * t161 - Ifges(6,4) * t194;
t126 = Ifges(6,4) * t161 - Ifges(6,2) * t194;
t125 = mrSges(6,1) * t194 + mrSges(6,2) * t161;
t116 = mrSges(5,1) * t157 - mrSges(5,3) * t222;
t115 = -mrSges(5,2) * t157 - mrSges(5,3) * t223;
t108 = pkin(5) * t124 + t207;
t93 = Ifges(5,5) * t157 + t158 * t197;
t92 = t158 * t196 + t226;
t91 = mrSges(6,1) * t157 + mrSges(6,3) * t105;
t90 = -mrSges(6,2) * t157 - mrSges(6,3) * t104;
t89 = -mrSges(5,2) * t152 - mrSges(5,3) * t192;
t88 = mrSges(5,1) * t152 + mrSges(5,3) * t191;
t87 = -Ifges(6,1) * t123 - Ifges(6,4) * t124;
t86 = -Ifges(6,4) * t123 - Ifges(6,2) * t124;
t84 = Ifges(7,1) * t118 + Ifges(7,4) * t117;
t83 = Ifges(7,4) * t118 + Ifges(7,2) * t117;
t82 = -mrSges(7,1) * t117 + mrSges(7,2) * t118;
t77 = mrSges(5,1) * t192 - mrSges(5,2) * t191;
t74 = mrSges(6,1) * t104 - mrSges(6,2) * t105;
t73 = pkin(5) * t104 + t99;
t68 = -Ifges(6,1) * t105 - Ifges(6,4) * t104 + Ifges(6,5) * t157;
t67 = -Ifges(6,4) * t105 - Ifges(6,2) * t104 + Ifges(6,6) * t157;
t66 = -Ifges(5,1) * t191 - Ifges(5,4) * t192 + t228;
t65 = -Ifges(5,4) * t191 - Ifges(5,2) * t192 + t227;
t61 = mrSges(7,1) * t157 - mrSges(7,3) * t72;
t60 = -mrSges(7,2) * t157 + mrSges(7,3) * t71;
t45 = -mrSges(6,2) * t152 + mrSges(6,3) * t49;
t44 = mrSges(6,1) * t152 - mrSges(6,3) * t48;
t43 = -qJD(4) * t81 + t199;
t41 = -mrSges(7,1) * t71 + mrSges(7,2) * t72;
t40 = Ifges(7,1) * t72 + Ifges(7,4) * t71 + Ifges(7,5) * t157;
t39 = Ifges(7,4) * t72 + Ifges(7,2) * t71 + Ifges(7,6) * t157;
t34 = Ifges(7,1) * t62 + Ifges(7,4) * t63;
t33 = Ifges(7,4) * t62 + Ifges(7,2) * t63;
t31 = -pkin(5) * t49 + t76;
t25 = -mrSges(6,1) * t49 + mrSges(6,2) * t48;
t24 = Ifges(6,1) * t48 + Ifges(6,4) * t49 + t152 * Ifges(6,5);
t23 = Ifges(6,4) * t48 + Ifges(6,2) * t49 + t152 * Ifges(6,6);
t15 = -mrSges(7,2) * t152 + mrSges(7,3) * t19;
t14 = mrSges(7,1) * t152 - mrSges(7,3) * t18;
t8 = -mrSges(7,1) * t19 + mrSges(7,2) * t18;
t7 = Ifges(7,1) * t18 + Ifges(7,4) * t19 + t152 * Ifges(7,5);
t6 = Ifges(7,4) * t18 + Ifges(7,2) * t19 + t152 * Ifges(7,6);
t1 = [0.2e1 * m(4) * (t101 * t129 + t174 * t208 + t225) + 0.2e1 * m(5) * (t42 * t81 + t43 * t80 + t225) - (mrSges(4,3) * t242 - t180 * t92 + t184 * t93) * t153 + (mrSges(4,1) * t201 + t101 * t244 - t204 * t153 + ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3) + Ifges(7,3)) * t152 + t209 + t247 + t249) * t157 + (t12 * t3 + t13 * t2 + t31 * t73) * t245 + (t10 * t38 + t11 * t37 + t76 * t99) * t246 - t139 * t241 + t77 * t242 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t181 + mrSges(3,2) * t185) + (-Ifges(3,2) + Ifges(3,1)) * t181 * t185 + (-t181 ^ 2 + t185 ^ 2) * Ifges(3,4)) * qJD(2) + (mrSges(4,2) * t201 - 0.2e1 * Ifges(4,1) * t153 - t180 * t65 + t184 * t66 + (Ifges(5,5) * t184 + t204) * t152 + (-t184 * t92 - t180 * t93 + t157 * (-Ifges(5,5) * t180 - Ifges(5,6) * t184)) * qJD(4) + 0.2e1 * (t198 + mrSges(4,3)) * t100) * t158 + (mrSges(4,1) * t241 - Ifges(6,5) * t105 + Ifges(7,5) * t72 - Ifges(6,6) * t104 + Ifges(7,6) * t71 + t129 * t244) * t152 + 0.2e1 * t42 * t115 + 0.2e1 * t43 * t116 - t104 * t23 - t105 * t24 + 0.2e1 * t99 * t25 + 0.2e1 * t80 * t88 + 0.2e1 * t81 * t89 + 0.2e1 * t10 * t90 + 0.2e1 * t11 * t91 + t71 * t6 + t72 * t7 + 0.2e1 * t73 * t8 + 0.2e1 * t76 * t74 + 0.2e1 * t2 * t60 + 0.2e1 * t3 * t61 + t49 * t67 + t48 * t68 + 0.2e1 * t31 * t41 + 0.2e1 * t37 * t44 + 0.2e1 * t38 * t45 + t19 * t39 + t18 * t40 + 0.2e1 * t12 * t14 + 0.2e1 * t13 * t15; (Ifges(3,5) * t185 - Ifges(3,6) * t181 + (-mrSges(3,1) * t185 + mrSges(3,2) * t181) * pkin(7)) * qJD(2) + ((-t152 * t176 + t153 * t177) * mrSges(4,3) + m(4) * (-t100 * t177 + t101 * t176)) * pkin(2) + m(6) * (t10 * t111 + t11 * t110 + t165 * t76 + t37 * t79 + t38 * t78) + (t117 * t2 - t118 * t3 - t12 * t62 + t13 * t63) * mrSges(7,3) + (t233 + t217 + t175) * t157 / 0.2e1 + (Ifges(6,5) * t161 + Ifges(7,5) * t118 - Ifges(6,6) * t194 + Ifges(7,6) * t117) * t152 / 0.2e1 - t194 * t23 / 0.2e1 + m(7) * (t108 * t73 + t12 * t22 + t13 * t21 + t130 * t31 + t2 * t53 + t3 * t52) + (m(5) * t172 - mrSges(4,1) + t250) * t100 + (t227 / 0.2e1 + t65 / 0.2e1 + t42 * mrSges(5,3) + t158 * t164 / 0.2e1 - t153 * t169 / 0.2e1 + (t158 * t236 - t80 * mrSges(5,3) + t93 / 0.2e1) * qJD(4) + (-qJD(4) * t116 + t89 + m(5) * (-qJD(4) * t80 + t42)) * t171) * t184 + (t228 / 0.2e1 + t66 / 0.2e1 - t43 * mrSges(5,3) + t163 * t237 - t153 * t236 + (t169 * t237 + pkin(4) * t74 - t81 * mrSges(5,3) - t226 / 0.2e1 - t92 / 0.2e1 + t99 * t240) * qJD(4) + (-m(5) * t43 - t88 + (-m(5) * t81 - t115) * qJD(4)) * t171) * t180 + t128 * t162 + t165 * t25 + t172 * t77 + t161 * t24 / 0.2e1 - Ifges(4,5) * t153 - Ifges(4,6) * t152 + t130 * t8 - t123 * t68 / 0.2e1 - t124 * t67 / 0.2e1 + t76 * t125 + t49 * t126 / 0.2e1 + t48 * t127 / 0.2e1 + t117 * t6 / 0.2e1 + t118 * t7 / 0.2e1 - t104 * t86 / 0.2e1 - t105 * t87 / 0.2e1 + t108 * t41 + t110 * t44 + t111 * t45 + t99 * t85 - t101 * mrSges(4,2) + t78 * t90 + t79 * t91 + t31 * t82 + t19 * t83 / 0.2e1 + t18 * t84 / 0.2e1 + t71 * t33 / 0.2e1 + t72 * t34 / 0.2e1 + t73 * t32 + (-t10 * t194 - t11 * t161 + t123 * t37 - t124 * t38) * mrSges(6,3) + t21 * t60 + t22 * t61 + t62 * t40 / 0.2e1 + t63 * t39 / 0.2e1 + t52 * t14 + t53 * t15; 0.2e1 * t108 * t82 + t117 * t33 + t118 * t34 - t123 * t127 - t124 * t126 + 0.2e1 * t130 * t32 - t194 * t86 + t161 * t87 + 0.2e1 * t172 * t162 + t184 * t163 + t180 * t164 + 0.2e1 * t165 * t85 + t62 * t84 + t63 * t83 + (t184 * t169 + (0.2e1 * pkin(4) * t125 - t168) * t180) * qJD(4) + (t108 * t130 + t21 * t53 + t22 * t52) * t245 + (t110 * t79 + t111 * t78 + t165 * t207) * t246 + 0.2e1 * (t117 * t21 - t118 * t22 - t52 * t62 + t53 * t63) * mrSges(7,3) + 0.2e1 * (t110 * t123 - t111 * t124 - t161 * t79 - t194 * t78) * mrSges(6,3); m(4) * t208 + t152 * mrSges(4,1) + t117 * t14 + t118 * t15 - t123 * t90 - t124 * t91 - t194 * t44 + t161 * t45 + t180 * t89 + t184 * t88 + t62 * t60 + t63 * t61 - t139 + (t115 * t184 - t116 * t180) * qJD(4) + m(7) * (t117 * t3 + t118 * t2 + t12 * t63 + t13 * t62) + m(6) * (t10 * t161 - t11 * t194 - t123 * t38 - t124 * t37) + m(5) * (t180 * t42 + t184 * t43 + (-t180 * t80 + t184 * t81) * qJD(4)); m(7) * (t117 * t22 + t118 * t21 + t52 * t63 + t53 * t62) + m(6) * (-t110 * t124 - t111 * t123 + t161 * t78 - t194 * t79); 0.2e1 * m(6) * (-t123 * t161 + t124 * t194) + 0.2e1 * m(7) * (t117 * t63 + t118 * t62); -Ifges(5,5) * t205 + t186 + m(7) * (t112 * t13 + t113 * t12 + t150 * t3 + t151 * t2) - t192 * Ifges(5,6) + (t90 * t213 - t91 * t214 + t179 * t45 + t183 * t44 + m(6) * (t10 * t179 + t11 * t183 + t213 * t38 - t214 * t37)) * pkin(4) + t150 * t14 + t151 * t15 + t113 * t61 + t112 * t60 - t42 * mrSges(5,2) + t43 * mrSges(5,1) + t249; m(7) * (t112 * t53 + t113 * t52 + t150 * t22 + t151 * t21) + t175 + (t171 * t250 - t230) * qJD(4) + (t112 * t117 - t113 * t118 - t150 * t62 + t151 * t63) * mrSges(7,3) + (m(6) * (t179 * t78 + t183 * t79 + (-t110 * t179 + t111 * t183) * qJD(5)) + (t183 * t123 - t179 * t124 + (t161 * t179 - t183 * t194) * qJD(5)) * mrSges(6,3)) * pkin(4) + t187; m(7) * (t112 * t118 + t113 * t117 + t150 * t63 + t151 * t62) - t162 + (-t123 * t179 - t124 * t183 + (t161 * t183 + t179 * t194) * qJD(5)) * t240 + t190; (t112 * t151 + t113 * t150) * t245 - 0.2e1 * t229 + 0.2e1 * t109 + 0.2e1 * t189; (m(7) * (-t12 * t212 + t13 * t211 + t178 * t2 + t182 * t3) + t60 * t211 + t178 * t15 - t61 * t212 + t182 * t14) * pkin(5) + t186; (m(7) * (t178 * t21 + t182 * t22 + (-t178 * t52 + t182 * t53) * qJD(6)) + (t178 * t63 - t182 * t62 + (t117 * t182 + t118 * t178) * qJD(6)) * mrSges(7,3)) * pkin(5) + t187; m(7) * (t178 * t62 + t182 * t63 + (-t117 * t178 + t118 * t182) * qJD(6)) * pkin(5) + t190; t189 + (m(7) * (t112 * t178 + t113 * t182 - t150 * t212 + t151 * t211) - mrSges(7,2) * t211 - mrSges(7,1) * t212) * pkin(5) + t200; 0.2e1 * t156; t193; t195; -t32; t200; t156; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
