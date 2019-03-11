% Calculate time derivative of joint inertia matrix for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:48
% EndTime: 2019-03-09 06:10:57
% DurationCPUTime: 3.63s
% Computational Cost: add. (6776->322), mult. (14759->460), div. (0->0), fcn. (15205->8), ass. (0->159)
t272 = Ifges(6,1) + Ifges(7,1);
t151 = sin(pkin(10));
t155 = sin(qJ(3));
t152 = cos(pkin(10));
t158 = cos(qJ(3));
t210 = t152 * t158;
t119 = -t155 * t151 + t210;
t120 = t151 * t158 + t155 * t152;
t154 = sin(qJ(4));
t157 = cos(qJ(4));
t169 = t157 * t119 - t120 * t154;
t233 = Ifges(7,4) + Ifges(6,5);
t271 = t169 * t233;
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t270 = t153 ^ 2 + t156 ^ 2;
t222 = Ifges(7,5) * t153;
t224 = Ifges(6,4) * t153;
t269 = t272 * t156 + t222 - t224;
t221 = Ifges(7,5) * t156;
t223 = Ifges(6,4) * t156;
t268 = t153 * t272 - t221 + t223;
t267 = Ifges(6,6) * t156 + t233 * t153;
t266 = Ifges(7,2) + Ifges(6,3);
t203 = qJD(5) * t153;
t104 = t119 * qJD(3);
t105 = t120 * qJD(3);
t70 = t169 * qJD(4) + t104 * t157 - t105 * t154;
t213 = t156 * t70;
t95 = t119 * t154 + t120 * t157;
t166 = t95 * t203 - t213;
t202 = qJD(5) * t156;
t193 = t95 * t202;
t217 = t153 * t70;
t167 = t193 + t217;
t71 = t95 * qJD(4) + t104 * t154 + t157 * t105;
t265 = t233 * t71 + (-Ifges(6,4) + Ifges(7,5)) * t167 - t272 * t166;
t228 = t269 * t95 - t271;
t131 = -t156 * mrSges(6,1) + t153 * mrSges(6,2);
t264 = -mrSges(5,1) + t131;
t263 = t269 * qJD(5);
t262 = Ifges(7,6) * t203 + t202 * t233;
t218 = pkin(3) * qJD(4);
t197 = t157 * t218;
t260 = t270 * t197;
t232 = pkin(7) + qJ(2);
t128 = t232 * t151;
t129 = t232 * t152;
t97 = -t155 * t128 + t158 * t129;
t89 = -t120 * qJD(2) - qJD(3) * t97;
t161 = -t104 * pkin(8) + t89;
t117 = t158 * t128;
t96 = -t129 * t155 - t117;
t91 = -pkin(8) * t120 + t96;
t92 = pkin(8) * t119 + t97;
t255 = -t154 * t92 + t157 * t91;
t88 = -qJD(3) * t117 + qJD(2) * t210 + (-qJD(2) * t151 - qJD(3) * t129) * t155;
t77 = -pkin(8) * t105 + t88;
t24 = qJD(4) * t255 + t154 * t161 + t157 * t77;
t240 = pkin(3) * t105;
t36 = pkin(4) * t71 - pkin(9) * t70 + t240;
t51 = t154 * t91 + t157 * t92;
t190 = -pkin(2) * t152 - pkin(1);
t99 = -pkin(3) * t119 + t190;
t52 = -pkin(4) * t169 - pkin(9) * t95 + t99;
t6 = t153 * t36 + t156 * t24 + t52 * t202 - t51 * t203;
t227 = t153 * t52 + t156 * t51;
t7 = -qJD(5) * t227 - t153 * t24 + t156 * t36;
t259 = -t7 * t153 + t156 * t6;
t2 = qJ(6) * t71 - qJD(6) * t169 + t6;
t26 = -qJ(6) * t169 + t227;
t32 = -t153 * t51 + t156 * t52;
t27 = pkin(5) * t169 - t32;
t4 = -pkin(5) * t71 - t7;
t258 = t153 * t4 + t156 * t2 + t27 * t202 - t26 * t203;
t19 = t167 * mrSges(6,1) - t166 * mrSges(6,2);
t25 = t51 * qJD(4) + t154 * t77 - t157 * t161;
t256 = m(6) * t25 + t19;
t254 = (-mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t270) * t197;
t253 = 0.2e1 * m(6);
t252 = 2 * m(7);
t251 = -2 * mrSges(5,3);
t250 = -2 * Ifges(5,4);
t249 = 0.2e1 * t25;
t248 = -0.2e1 * t255;
t179 = t153 * mrSges(7,1) - t156 * mrSges(7,3);
t121 = t179 * qJD(5);
t247 = 0.2e1 * t121;
t246 = t71 / 0.2e1;
t133 = Ifges(6,2) * t156 + t224;
t243 = -t133 / 0.2e1;
t241 = Ifges(6,6) * t169;
t239 = pkin(3) * t157;
t235 = t25 * t255;
t28 = mrSges(6,1) * t71 + t166 * mrSges(6,3);
t29 = -t71 * mrSges(7,1) - mrSges(7,2) * t166;
t231 = t28 - t29;
t30 = -mrSges(6,2) * t71 - t167 * mrSges(6,3);
t31 = -t167 * mrSges(7,2) + mrSges(7,3) * t71;
t230 = t30 + t31;
t175 = Ifges(7,3) * t153 + t221;
t42 = -Ifges(7,6) * t169 + t175 * t95;
t176 = -Ifges(6,2) * t153 + t223;
t43 = t176 * t95 - t241;
t229 = t42 - t43;
t216 = t153 * t95;
t72 = mrSges(6,2) * t169 - mrSges(6,3) * t216;
t75 = -mrSges(7,2) * t216 - mrSges(7,3) * t169;
t226 = t72 + t75;
t212 = t156 * t95;
t73 = -mrSges(6,1) * t169 - mrSges(6,3) * t212;
t74 = mrSges(7,1) * t169 + mrSges(7,2) * t212;
t225 = -t73 + t74;
t220 = Ifges(6,6) * t153;
t215 = t154 * t255;
t211 = qJD(5) * t95;
t209 = t153 * t157;
t208 = t156 * t157;
t139 = pkin(3) * t154 + pkin(9);
t207 = t260 * t139;
t206 = t260 * pkin(9);
t201 = qJD(6) * t156;
t200 = 0.2e1 * t240;
t199 = m(7) * t201;
t198 = t154 * t218;
t189 = t71 * mrSges(5,1) + t70 * mrSges(5,2);
t188 = -t203 / 0.2e1;
t183 = mrSges(7,2) * t201 + t262;
t180 = t153 * mrSges(6,1) + t156 * mrSges(6,2);
t130 = -t156 * mrSges(7,1) - t153 * mrSges(7,3);
t174 = -pkin(5) * t156 - qJ(6) * t153;
t173 = pkin(5) * t153 - qJ(6) * t156;
t170 = t167 * Ifges(7,6) + t213 * t233 + t266 * t71;
t168 = t264 * t198;
t127 = -pkin(4) + t174;
t102 = -pkin(5) * t203 + qJ(6) * t202 + qJD(6) * t153;
t165 = t225 * t153 + t226 * t156;
t164 = t174 * mrSges(7,2) - t220;
t123 = t175 * qJD(5);
t124 = t176 * qJD(5);
t132 = -Ifges(7,3) * t156 + t222;
t163 = (t132 - t133) * t203 + t268 * t202 + (-t123 + t124) * t156 + t263 * t153;
t162 = m(7) * t174 + t130 + t131;
t160 = t230 * t156 - t231 * t153 + (-t226 * t153 + t225 * t156) * qJD(5) + m(7) * t258 + m(6) * (-t32 * t202 - t203 * t227 + t259);
t122 = t180 * qJD(5);
t14 = -t166 * Ifges(7,5) + Ifges(7,6) * t71 + t167 * Ifges(7,3);
t15 = -t166 * Ifges(6,4) - t167 * Ifges(6,2) + Ifges(6,6) * t71;
t37 = t173 * t95 - t255;
t9 = (pkin(5) * t70 + qJ(6) * t211) * t153 + (-qJ(6) * t70 + (pkin(5) * qJD(5) - qJD(6)) * t95) * t156 + t25;
t159 = t42 * t203 / 0.2e1 + t193 * t243 + t43 * t188 - t24 * mrSges(5,2) + Ifges(5,5) * t70 - Ifges(5,6) * t71 + t37 * t121 - t255 * t122 + t9 * t130 + t264 * t25 + t267 * t246 - (-Ifges(6,6) * t203 + t262) * t169 / 0.2e1 + t265 * t153 / 0.2e1 + (t132 / 0.2e1 + t243) * t217 + (-t124 / 0.2e1 + t123 / 0.2e1) * t216 + t263 * t212 / 0.2e1 + (-Ifges(7,6) * t246 + t15 / 0.2e1 - t14 / 0.2e1) * t156 + ((-t153 * t227 - t156 * t32) * qJD(5) + t259) * mrSges(6,3) + (t95 * t132 + t228) * t202 / 0.2e1 + t258 * mrSges(7,2) + t268 * (t95 * t188 + t213 / 0.2e1);
t142 = mrSges(7,2) * t202;
t140 = -pkin(4) - t239;
t110 = t127 - t239;
t101 = t104 * mrSges(4,2);
t98 = -t102 + t198;
t64 = t180 * t95;
t63 = t179 * t95;
t18 = t167 * mrSges(7,1) + t166 * mrSges(7,3);
t1 = [0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t151 ^ 2 + t152 ^ 2) * qJD(2) + 0.2e1 * t99 * t189 + (mrSges(5,3) * t248 + t229 * t153 + t228 * t156) * t70 + t19 * t248 + t64 * t249 + (t2 * t26 + t27 * t4 + t37 * t9) * t252 + 0.2e1 * m(5) * (t24 * t51 + t99 * t240 - t235) + (t227 * t6 + t32 * t7 - t235) * t253 + 0.2e1 * t190 * (t105 * mrSges(4,1) + t101) + 0.2e1 * (t104 * t119 - t105 * t120) * Ifges(4,4) + 0.2e1 * (-t104 * t96 - t105 * t97 + t119 * t88 - t120 * t89) * mrSges(4,3) + 0.2e1 * t227 * t30 - 0.2e1 * t119 * Ifges(4,2) * t105 + 0.2e1 * t120 * t104 * Ifges(4,1) + 0.2e1 * m(4) * (t88 * t97 + t89 * t96) - (mrSges(5,1) * t200 + t24 * t251 + (t250 - t220) * t70 + ((2 * Ifges(5,2)) + t266) * t71 + t170) * t169 + t51 * t71 * t251 + (mrSges(5,2) * t200 + mrSges(5,3) * t249 + 0.2e1 * Ifges(5,1) * t70 + t265 * t156 + (t14 - t15) * t153 + (t250 + t233 * t156 + (-Ifges(6,6) + Ifges(7,6)) * t153) * t71 + ((t229 + t241) * t156 + (-t228 + t271) * t153) * qJD(5)) * t95 + 0.2e1 * t27 * t29 + 0.2e1 * t26 * t31 + 0.2e1 * t32 * t28 + 0.2e1 * t37 * t18 + 0.2e1 * t9 * t63 + 0.2e1 * t6 * t72 + 0.2e1 * t7 * t73 + 0.2e1 * t4 * t74 + 0.2e1 * t2 * t75; t101 + t231 * t156 + t230 * t153 - (-m(5) * pkin(3) - mrSges(4,1)) * t105 + t165 * qJD(5) + m(7) * (t153 * t2 - t156 * t4 + (t153 * t27 + t156 * t26) * qJD(5)) + m(6) * (t153 * t6 + t156 * t7 + (-t153 * t32 + t156 * t227) * qJD(5)) + t189; 0; t160 * t139 + (m(5) * (t154 * t24 - t157 * t25) + (-t154 * t71 - t157 * t70) * mrSges(5,3) + ((t95 * mrSges(5,3) + t64) * t154 + (mrSges(5,3) * t169 + t165) * t157 + m(7) * (t26 * t208 + t27 * t209) + m(6) * (t208 * t227 - t32 * t209 - t215) + m(5) * (t157 * t51 - t215)) * qJD(4)) * pkin(3) + t159 + m(7) * (t110 * t9 + t37 * t98) - t88 * mrSges(4,2) + t89 * mrSges(4,1) + t98 * t63 + Ifges(4,5) * t104 - Ifges(4,6) * t105 + t110 * t18 + t256 * t140; 0; t110 * t247 + 0.2e1 * t140 * t122 + 0.2e1 * t98 * t130 + 0.2e1 * t168 + (t110 * t98 + t207) * t252 + (t140 * t198 + t207) * t253 + 0.2e1 * t254 + t163; m(7) * (-t102 * t37 + t127 * t9) + t160 * pkin(9) + t159 - t102 * t63 + t127 * t18 - t256 * pkin(4); 0; (t98 - t102) * t130 + (t140 - pkin(4)) * t122 + (t110 + t127) * t121 + t168 + m(7) * (-t102 * t110 + t127 * t98 + t206) + m(6) * (-pkin(4) * t198 + t206) + t254 + t163; -0.2e1 * pkin(4) * t122 + t127 * t247 + 0.2e1 * (-m(7) * t127 - t130) * t102 + t163; -Ifges(6,6) * t217 - pkin(5) * t29 + m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t26) + qJD(6) * t75 + qJ(6) * t31 + t2 * mrSges(7,3) + t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) - t267 * t211 + t170; m(7) * t102 + ((-mrSges(6,2) + mrSges(7,3)) * t156 + (-mrSges(6,1) - mrSges(7,1)) * t153) * qJD(5); t139 * t199 + (-m(7) * t173 - t179 - t180) * t197 + (t162 * t139 + t164) * qJD(5) + t183; t164 * qJD(5) + (t162 * qJD(5) + t199) * pkin(9) + t183; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t4 + t29; m(7) * t203; t142 + (t139 * t202 + t153 * t197) * m(7); m(7) * pkin(9) * t202 + t142; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
