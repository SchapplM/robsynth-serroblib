% Calculate time derivative of joint inertia matrix for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR11_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:51
% EndTime: 2019-03-09 11:11:59
% DurationCPUTime: 3.74s
% Computational Cost: add. (5490->439), mult. (11857->646), div. (0->0), fcn. (10294->8), ass. (0->198)
t177 = sin(pkin(10));
t178 = cos(pkin(10));
t183 = cos(qJ(4));
t219 = qJD(4) * t183;
t180 = sin(qJ(4));
t220 = qJD(4) * t180;
t129 = t177 * t220 - t178 * t219;
t193 = t177 * t183 + t178 * t180;
t130 = t193 * qJD(4);
t179 = sin(qJ(6));
t182 = cos(qJ(6));
t192 = t177 * t180 - t178 * t183;
t256 = -t179 * t193 - t182 * t192;
t47 = qJD(6) * t256 - t129 * t182 - t130 * t179;
t90 = t179 * t192 - t182 * t193;
t265 = t47 * t90;
t184 = cos(qJ(2));
t120 = t193 * t184;
t181 = sin(qJ(2));
t244 = pkin(3) + pkin(7);
t157 = t244 * t181;
t143 = t183 * t157;
t185 = -pkin(2) - pkin(8);
t231 = qJ(3) * t181;
t138 = t184 * t185 - pkin(1) - t231;
t209 = qJ(5) * t184 - t138;
t76 = pkin(4) * t181 + t180 * t209 + t143;
t225 = t183 * t184;
t94 = t183 * t138 + t180 * t157;
t83 = -qJ(5) * t225 + t94;
t40 = -t177 * t83 + t178 * t76;
t26 = pkin(5) * t181 + pkin(9) * t120 + t40;
t119 = t192 * t184;
t41 = t177 * t76 + t178 * t83;
t27 = pkin(9) * t119 + t41;
t12 = -t179 * t27 + t182 * t26;
t13 = t179 * t26 + t182 * t27;
t186 = qJD(6) * t90 + t129 * t179 - t182 * t130;
t222 = qJD(2) * t181;
t207 = pkin(2) * t222 - qJD(3) * t181;
t107 = (pkin(8) * t181 - qJ(3) * t184) * qJD(2) + t207;
t221 = qJD(2) * t184;
t150 = t244 * t221;
t133 = t183 * t150;
t30 = pkin(4) * t221 + t133 + t209 * t219 + (-qJ(5) * t222 - qJD(4) * t157 + qJD(5) * t184 - t107) * t180;
t218 = qJD(4) * t184;
t214 = t180 * t218;
t191 = t183 * t222 + t214;
t217 = t183 * qJD(5);
t51 = t183 * t107 - t138 * t220 + t180 * t150 + t157 * t219;
t32 = qJ(5) * t191 - t184 * t217 + t51;
t14 = -t177 * t32 + t178 * t30;
t78 = qJD(4) * t119 + t193 * t222;
t7 = pkin(5) * t221 - pkin(9) * t78 + t14;
t15 = t177 * t30 + t178 * t32;
t77 = -t192 * t222 + t193 * t218;
t8 = pkin(9) * t77 + t15;
t2 = qJD(6) * t12 + t179 * t7 + t182 * t8;
t3 = -qJD(6) * t13 - t179 * t8 + t182 * t7;
t268 = t12 * t186 + t13 * t47 - t2 * t90 + t256 * t3;
t224 = qJ(5) - t185;
t151 = t224 * t180;
t152 = t224 * t183;
t98 = t151 * t177 - t178 * t152;
t69 = pkin(9) * t192 + t98;
t99 = -t178 * t151 - t177 * t152;
t70 = -pkin(9) * t193 + t99;
t35 = -t179 * t70 + t182 * t69;
t117 = t220 * t224 - t217;
t118 = -qJD(4) * t152 - t180 * qJD(5);
t67 = t178 * t117 - t118 * t177;
t57 = pkin(9) * t130 + t67;
t68 = t177 * t117 + t178 * t118;
t58 = pkin(9) * t129 + t68;
t10 = qJD(6) * t35 + t179 * t57 + t182 * t58;
t36 = t179 * t69 + t182 * t70;
t11 = -qJD(6) * t36 - t179 * t58 + t182 * t57;
t267 = -t10 * t90 + t11 * t256 + t186 * t35 + t36 * t47;
t167 = pkin(4) * t178 + pkin(5);
t238 = pkin(4) * t177;
t125 = t167 * t182 - t179 * t238;
t112 = t125 * qJD(6);
t126 = t167 * t179 + t182 * t238;
t113 = t126 * qJD(6);
t266 = -t112 * t90 - t113 * t256 + t125 * t186 + t126 * t47;
t264 = Ifges(5,3) + Ifges(6,3);
t260 = 0.2e1 * qJ(3);
t259 = m(4) * pkin(7);
t258 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t257 = t186 * t256;
t52 = -qJD(4) * t94 - t107 * t180 + t133;
t255 = t51 * t180 + t52 * t183;
t251 = 2 * m(6);
t250 = 2 * m(7);
t249 = -0.2e1 * pkin(1);
t127 = -qJ(3) * t221 + t207;
t248 = 0.2e1 * t127;
t200 = -pkin(2) * t184 - t231;
t153 = -pkin(1) + t200;
t247 = -0.2e1 * t153;
t246 = t90 / 0.2e1;
t245 = t256 / 0.2e1;
t242 = -t193 / 0.2e1;
t241 = -t192 / 0.2e1;
t239 = -t184 / 0.2e1;
t175 = t184 * pkin(7);
t237 = Ifges(7,5) * t186 - Ifges(7,6) * t47;
t236 = Ifges(5,4) * t180;
t235 = Ifges(5,4) * t183;
t234 = t181 * Ifges(5,5);
t203 = Ifges(5,1) * t180 + t235;
t116 = -t184 * t203 + t234;
t229 = t180 * t116;
t156 = Ifges(5,1) * t183 - t236;
t228 = t180 * t156;
t202 = Ifges(5,2) * t183 + t236;
t115 = t181 * Ifges(5,6) - t184 * t202;
t227 = t183 * t115;
t155 = -Ifges(5,2) * t180 + t235;
t226 = t183 * t155;
t168 = t180 * pkin(4) + qJ(3);
t223 = -Ifges(6,5) * t130 + Ifges(6,6) * t129;
t158 = t184 * pkin(3) + t175;
t162 = pkin(4) * t219 + qJD(3);
t216 = 2 * mrSges(6,3);
t71 = t119 * t182 + t120 * t179;
t24 = qJD(6) * t71 + t179 * t77 + t182 * t78;
t72 = t119 * t179 - t120 * t182;
t25 = -qJD(6) * t72 - t179 * t78 + t182 * t77;
t215 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t221;
t131 = pkin(4) * t225 + t158;
t213 = t183 * t218;
t212 = t180 * t222;
t42 = -t77 * mrSges(6,1) + t78 * mrSges(6,2);
t6 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t16 = mrSges(7,1) * t47 + mrSges(7,2) * t186;
t210 = mrSges(7,1) * t186 - t47 * mrSges(7,2);
t84 = -t129 * mrSges(6,1) - t130 * mrSges(6,2);
t206 = t11 * mrSges(7,1) - t10 * mrSges(7,2) + t237;
t205 = mrSges(5,1) * t183 - mrSges(5,2) * t180;
t154 = mrSges(5,1) * t180 + mrSges(5,2) * t183;
t204 = -t113 * mrSges(7,1) - t112 * mrSges(7,2);
t201 = -Ifges(5,5) * t180 - Ifges(5,6) * t183;
t93 = -t138 * t180 + t143;
t199 = t93 * t180 - t94 * t183;
t198 = -t129 * t193 + t130 * t192;
t197 = t129 * t177 + t130 * t178;
t147 = mrSges(5,3) * t180 * t184 + mrSges(5,1) * t181;
t148 = -mrSges(5,2) * t181 - mrSges(5,3) * t225;
t194 = -t180 * t147 + t183 * t148;
t190 = t212 - t213;
t189 = Ifges(5,5) * t212 + Ifges(6,5) * t78 + Ifges(5,6) * t191 + Ifges(6,6) * t77 + t221 * t264 + t215;
t188 = t129 * t41 + t130 * t40 + t14 * t192 - t15 * t193;
t187 = t129 * t99 + t130 * t98 + t192 * t67 - t193 * t68;
t100 = -pkin(4) * t214 + (-pkin(4) * t183 - t244) * t222;
t149 = t244 * t222;
t146 = t203 * qJD(4);
t145 = t202 * qJD(4);
t144 = t205 * qJD(4);
t134 = t205 * t184;
t111 = pkin(5) * t193 + t168;
t106 = mrSges(5,1) * t221 - mrSges(5,3) * t190;
t105 = -mrSges(5,2) * t221 + mrSges(5,3) * t191;
t104 = mrSges(6,1) * t181 + mrSges(6,3) * t120;
t103 = -mrSges(6,2) * t181 + mrSges(6,3) * t119;
t102 = -pkin(5) * t129 + t162;
t97 = -Ifges(6,1) * t192 - Ifges(6,4) * t193;
t96 = -Ifges(6,4) * t192 - Ifges(6,2) * t193;
t95 = mrSges(6,1) * t193 - mrSges(6,2) * t192;
t88 = -mrSges(5,1) * t191 + mrSges(5,2) * t190;
t87 = -pkin(5) * t119 + t131;
t86 = -Ifges(6,1) * t130 + Ifges(6,4) * t129;
t85 = -Ifges(6,4) * t130 + Ifges(6,2) * t129;
t82 = -t156 * t218 + (t184 * Ifges(5,5) + t181 * t203) * qJD(2);
t81 = -t155 * t218 + (t184 * Ifges(5,6) + t181 * t202) * qJD(2);
t80 = -mrSges(6,1) * t119 - mrSges(6,2) * t120;
t66 = -Ifges(6,1) * t120 + Ifges(6,4) * t119 + Ifges(6,5) * t181;
t65 = -Ifges(6,4) * t120 + Ifges(6,2) * t119 + Ifges(6,6) * t181;
t62 = mrSges(6,1) * t221 - mrSges(6,3) * t78;
t61 = -mrSges(6,2) * t221 + mrSges(6,3) * t77;
t60 = mrSges(7,1) * t181 - mrSges(7,3) * t72;
t59 = -mrSges(7,2) * t181 + mrSges(7,3) * t71;
t56 = -pkin(5) * t77 + t100;
t55 = Ifges(7,1) * t256 + Ifges(7,4) * t90;
t54 = Ifges(7,4) * t256 + Ifges(7,2) * t90;
t53 = -mrSges(7,1) * t90 + mrSges(7,2) * t256;
t39 = -mrSges(7,1) * t71 + mrSges(7,2) * t72;
t38 = Ifges(6,1) * t78 + Ifges(6,4) * t77 + Ifges(6,5) * t221;
t37 = Ifges(6,4) * t78 + Ifges(6,2) * t77 + Ifges(6,6) * t221;
t34 = Ifges(7,1) * t72 + Ifges(7,4) * t71 + Ifges(7,5) * t181;
t33 = Ifges(7,4) * t72 + Ifges(7,2) * t71 + Ifges(7,6) * t181;
t20 = -mrSges(7,2) * t221 + mrSges(7,3) * t25;
t19 = mrSges(7,1) * t221 - mrSges(7,3) * t24;
t18 = Ifges(7,1) * t186 - Ifges(7,4) * t47;
t17 = Ifges(7,4) * t186 - Ifges(7,2) * t47;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t25 + Ifges(7,5) * t221;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t25 + Ifges(7,6) * t221;
t1 = [(t12 * t3 + t13 * t2 + t56 * t87) * t250 + (t100 * t131 + t14 * t40 + t15 * t41) * t251 + (-0.2e1 * t127 * mrSges(4,3) + t189) * t181 + 0.2e1 * m(5) * (-t149 * t158 + t51 * t94 + t52 * t93) + m(4) * t153 * t248 + 0.2e1 * t158 * t88 + 0.2e1 * t52 * t147 + 0.2e1 * t51 * t148 - 0.2e1 * t149 * t134 + 0.2e1 * t131 * t42 + t119 * t37 - t120 * t38 + 0.2e1 * t100 * t80 + 0.2e1 * t15 * t103 + 0.2e1 * t14 * t104 + 0.2e1 * t94 * t105 + 0.2e1 * t93 * t106 + 0.2e1 * t87 * t6 + t77 * t65 + t78 * t66 + t71 * t4 + t72 * t5 + 0.2e1 * t2 * t59 + 0.2e1 * t3 * t60 + 0.2e1 * t41 * t61 + 0.2e1 * t40 * t62 + 0.2e1 * t56 * t39 + t25 * t33 + t24 * t34 + 0.2e1 * t12 * t19 + 0.2e1 * t13 * t20 + ((mrSges(3,1) * t249 + mrSges(4,2) * t247 + t227 + t229 + 0.2e1 * (-Ifges(3,4) - Ifges(4,6)) * t181) * t181 + (mrSges(3,2) * t249 + mrSges(4,3) * t247 - Ifges(6,5) * t120 + Ifges(7,5) * t72 + Ifges(6,6) * t119 + Ifges(7,6) * t71 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + t201) * t184 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(7,3) + t264) * t181) * t184) * qJD(2) + (mrSges(4,2) * t248 - t180 * t82 - t183 * t81 + (t115 * t180 + (-t116 - t234) * t183) * qJD(4)) * t184; (t200 * t259 + (-pkin(2) * mrSges(4,1) + Ifges(5,5) * t183 / 0.2e1 - Ifges(5,6) * t180 / 0.2e1 - Ifges(4,4) + Ifges(3,5) + Ifges(7,5) * t245 + Ifges(7,6) * t246 + Ifges(6,5) * t241 + Ifges(6,6) * t242 + (-mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t184 + (t226 / 0.2e1 + t228 / 0.2e1 - qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + (mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t181) * qJD(2) + (qJD(4) * t201 + t223 + t237) * t181 / 0.2e1 + t5 * t245 + t4 * t246 + t38 * t241 + t37 * t242 + t188 * mrSges(6,3) + (-t227 / 0.2e1 - t229 / 0.2e1 + (-t183 * t156 / 0.2e1 + t180 * t155 / 0.2e1) * t184 + t199 * mrSges(5,3) + (-m(5) * t199 + t194) * t185) * qJD(4) + (t82 / 0.2e1 + t185 * t106 - t145 * t239 - t52 * mrSges(5,3)) * t183 + (-t81 / 0.2e1 + t185 * t105 - t146 * t239 - t51 * mrSges(5,3)) * t180 + t186 * t34 / 0.2e1 - t47 * t33 / 0.2e1 - t268 * mrSges(7,3) + t168 * t42 - t149 * t154 + t158 * t144 + t162 * t80 - t130 * t66 / 0.2e1 + t131 * t84 + t129 * t65 / 0.2e1 + t119 * t85 / 0.2e1 - t120 * t86 / 0.2e1 + t111 * t6 + t77 * t96 / 0.2e1 + t78 * t97 / 0.2e1 + t98 * t62 + t99 * t61 + t100 * t95 + t102 * t39 + t68 * t103 + t67 * t104 + t87 * t16 + qJ(3) * t88 + t71 * t17 / 0.2e1 + t72 * t18 / 0.2e1 + t10 * t59 + t11 * t60 + t25 * t54 / 0.2e1 + t24 * t55 / 0.2e1 + t56 * t53 + t35 * t19 + t36 * t20 + (m(4) * t175 + m(5) * t158 + t184 * mrSges(4,1) + t134) * qJD(3) + m(6) * (t100 * t168 + t131 * t162 + t14 * t98 + t15 * t99 + t40 * t67 + t41 * t68) + m(7) * (t10 * t13 + t102 * t87 + t11 * t12 + t111 * t56 + t2 * t36 + t3 * t35) + m(5) * (-qJ(3) * t149 + t185 * t255); t144 * t260 + 0.2e1 * t102 * t53 + 0.2e1 * t111 * t16 + t129 * t96 - t130 * t97 - t193 * t85 - t192 * t86 + t180 * t145 - t183 * t146 + 0.2e1 * t162 * t95 + 0.2e1 * t168 * t84 + t90 * t17 + t256 * t18 + t186 * t55 - t47 * t54 + (-t226 - t228) * qJD(4) + (t162 * t168 + t67 * t98 + t68 * t99) * t251 + (t10 * t36 + t102 * t111 + t11 * t35) * t250 + t187 * t216 - 0.2e1 * t267 * mrSges(7,3) + (0.2e1 * mrSges(4,3) + 0.2e1 * t154 + (m(4) + m(5)) * t260) * qJD(3); -t129 * t103 - t130 * t104 + t180 * t105 + t183 * t106 + t193 * t61 - t192 * t62 + t256 * t19 - t90 * t20 + t47 * t59 + t186 * t60 + t194 * qJD(4) + (mrSges(4,1) + t259) * t221 + m(7) * t268 - m(6) * t188 + m(5) * (-t199 * qJD(4) + t255); -t198 * t216 + m(7) * t267 - m(6) * t187 + (-0.2e1 * t257 + 0.2e1 * t265) * mrSges(7,3); 0.2e1 * m(6) * t198 + 0.2e1 * m(7) * (t257 - t265); t189 + t125 * t19 + t126 * t20 + t112 * t59 - t113 * t60 - t51 * mrSges(5,2) + t52 * mrSges(5,1) - t15 * mrSges(6,2) + t14 * mrSges(6,1) - Ifges(5,5) * t213 + (t177 * t61 + m(6) * (t14 * t178 + t15 * t177) + t178 * t62) * pkin(4) + m(7) * (t112 * t13 - t113 * t12 + t125 * t3 + t126 * t2) + t258; m(7) * (t10 * t126 + t11 * t125 + t112 * t36 - t113 * t35) - t68 * mrSges(6,2) + t67 * mrSges(6,1) + ((-mrSges(5,2) * t185 - Ifges(5,6)) * t183 + (-mrSges(5,1) * t185 - Ifges(5,5)) * t180) * qJD(4) + (m(6) * (t177 * t68 + t178 * t67) + t197 * mrSges(6,3)) * pkin(4) - t266 * mrSges(7,3) + t206 + t223; -m(6) * t197 * pkin(4) + m(7) * t266 - t130 * mrSges(6,1) + t129 * mrSges(6,2) - t154 * qJD(4) + t210; 0.2e1 * m(7) * (t112 * t126 - t113 * t125) + 0.2e1 * t204; m(6) * t100 + m(7) * t56 + t42 + t6; m(6) * t162 + m(7) * t102 + t16 + t84; 0; 0; 0; t215 + t258; t206; t210; t204; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
