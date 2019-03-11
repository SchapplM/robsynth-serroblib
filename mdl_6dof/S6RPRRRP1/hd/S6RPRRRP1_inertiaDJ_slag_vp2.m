% Calculate time derivative of joint inertia matrix for
% S6RPRRRP1
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:33
% EndTime: 2019-03-09 05:55:40
% DurationCPUTime: 3.38s
% Computational Cost: add. (3919->322), mult. (8472->451), div. (0->0), fcn. (7326->8), ass. (0->162)
t283 = Ifges(6,1) + Ifges(7,1);
t150 = sin(qJ(4));
t151 = sin(qJ(3));
t153 = cos(qJ(4));
t154 = cos(qJ(3));
t118 = t150 * t151 - t153 * t154;
t240 = Ifges(7,4) + Ifges(6,5);
t282 = t240 * t118;
t149 = sin(qJ(5));
t152 = cos(qJ(5));
t209 = t149 ^ 2 + t152 ^ 2;
t230 = Ifges(7,5) * t149;
t232 = Ifges(6,4) * t149;
t281 = t283 * t152 + t230 - t232;
t280 = mrSges(7,2) + mrSges(6,3);
t229 = Ifges(7,5) * t152;
t231 = Ifges(6,4) * t152;
t279 = t283 * t149 - t229 + t231;
t278 = Ifges(6,6) * t152 + t240 * t149;
t277 = Ifges(7,2) + Ifges(6,3);
t119 = t150 * t154 + t153 * t151;
t206 = qJD(5) * t149;
t264 = qJD(3) + qJD(4);
t85 = t264 * t118;
t219 = t152 * t85;
t164 = t119 * t206 + t219;
t205 = qJD(5) * t152;
t192 = t119 * t205;
t221 = t149 * t85;
t165 = t192 - t221;
t86 = t264 * t119;
t276 = t240 * t86 + (-Ifges(6,4) + Ifges(7,5)) * t165 - t283 * t164;
t237 = t281 * t119 + t282;
t129 = -t152 * mrSges(6,1) + t149 * mrSges(6,2);
t275 = -mrSges(5,1) + t129;
t274 = t281 * qJD(5);
t128 = -t152 * mrSges(7,1) - t149 * mrSges(7,3);
t273 = t128 + t129;
t272 = Ifges(7,6) * t206 + t240 * t205;
t207 = qJD(4) * t153;
t199 = pkin(3) * t207;
t270 = t209 * t199;
t269 = t209 * t85;
t136 = sin(pkin(10)) * pkin(1) + pkin(7);
t239 = pkin(8) + t136;
t188 = qJD(3) * t239;
t105 = t151 * t188;
t180 = t154 * t188;
t111 = t239 * t151;
t112 = t239 * t154;
t265 = -t153 * t111 - t112 * t150;
t34 = qJD(4) * t265 - t153 * t105 - t150 * t180;
t201 = pkin(3) * qJD(3) * t151;
t250 = pkin(4) * t86;
t38 = pkin(9) * t85 + t201 + t250;
t195 = -cos(pkin(10)) * pkin(1) - pkin(2);
t126 = -pkin(3) * t154 + t195;
t61 = t118 * pkin(4) - t119 * pkin(9) + t126;
t63 = -t111 * t150 + t112 * t153;
t6 = t149 * t38 + t152 * t34 + t61 * t205 - t63 * t206;
t236 = t149 * t61 + t152 * t63;
t7 = -qJD(5) * t236 - t149 * t34 + t152 * t38;
t268 = -t7 * t149 + t152 * t6;
t2 = qJ(6) * t86 + qJD(6) * t118 + t6;
t22 = qJ(6) * t118 + t236;
t31 = -t149 * t63 + t152 * t61;
t23 = -pkin(5) * t118 - t31;
t4 = -pkin(5) * t86 - t7;
t267 = t149 * t4 + t152 * t2 + t23 * t205 - t22 * t206;
t20 = t165 * mrSges(6,1) - t164 * mrSges(6,2);
t35 = t63 * qJD(4) - t105 * t150 + t153 * t180;
t266 = m(6) * t35 + t20;
t263 = (t280 * t209 - mrSges(5,2)) * t199;
t262 = 2 * m(5);
t261 = 0.2e1 * m(6);
t260 = 2 * m(7);
t259 = -2 * mrSges(5,3);
t258 = -2 * Ifges(5,4);
t257 = 0.2e1 * t35;
t256 = -0.2e1 * t265;
t176 = t149 * mrSges(7,1) - t152 * mrSges(7,3);
t120 = t176 * qJD(5);
t255 = 0.2e1 * t120;
t254 = 0.2e1 * t126;
t253 = m(6) / 0.2e1;
t252 = t86 / 0.2e1;
t131 = Ifges(6,2) * t152 + t232;
t248 = -t131 / 0.2e1;
t246 = pkin(3) * t153;
t242 = t35 * t265;
t172 = Ifges(7,3) * t149 + t229;
t43 = Ifges(7,6) * t118 + t172 * t119;
t173 = -Ifges(6,2) * t149 + t231;
t228 = Ifges(6,6) * t118;
t44 = t173 * t119 + t228;
t238 = t43 - t44;
t235 = t269 * pkin(9);
t217 = t119 * t149;
t76 = -mrSges(6,2) * t118 - mrSges(6,3) * t217;
t79 = -mrSges(7,2) * t217 + mrSges(7,3) * t118;
t234 = t76 + t79;
t216 = t119 * t152;
t77 = mrSges(6,1) * t118 - mrSges(6,3) * t216;
t78 = -mrSges(7,1) * t118 + mrSges(7,2) * t216;
t233 = -t77 + t78;
t227 = Ifges(6,6) * t149;
t225 = t118 * t86;
t224 = t119 * t85;
t220 = t150 * t265;
t215 = t149 * t153;
t213 = t152 * t153;
t138 = pkin(3) * t150 + pkin(9);
t212 = t270 * t138;
t211 = t270 * pkin(9);
t208 = qJD(4) * t150;
t204 = qJD(6) * t152;
t202 = m(7) * t204;
t200 = pkin(3) * t208;
t194 = t118 * t208;
t190 = -t206 / 0.2e1;
t187 = 0.2e1 * t201;
t186 = t270 * t119 - t269 * t138;
t183 = mrSges(7,2) * t204 + t272;
t179 = t118 * t35 - t265 * t86;
t178 = mrSges(4,1) * t151 + mrSges(4,2) * t154;
t177 = t149 * mrSges(6,1) + t152 * mrSges(6,2);
t171 = pkin(5) * t152 + qJ(6) * t149;
t170 = pkin(5) * t149 - qJ(6) * t152;
t167 = t275 * t200;
t166 = t165 * Ifges(7,6) - t240 * t219 + t277 * t86;
t127 = -pkin(4) - t171;
t106 = pkin(5) * t206 - qJ(6) * t205 - qJD(6) * t149;
t163 = t233 * t149 + t234 * t152;
t26 = -t86 * mrSges(7,1) - t164 * mrSges(7,2);
t162 = -t171 * mrSges(7,2) - t227;
t121 = t177 * qJD(5);
t80 = t86 * mrSges(5,1);
t161 = t85 * mrSges(5,2) - t80 + t273 * t86 + (t120 + t121) * t118 - t280 * t269;
t122 = t172 * qJD(5);
t123 = t173 * qJD(5);
t130 = -Ifges(7,3) * t152 + t230;
t160 = (t130 - t131) * t206 + t279 * t205 + (-t122 + t123) * t152 + t274 * t149;
t159 = -m(7) * t171 + t273;
t158 = -m(7) * t170 - t176 - t177;
t157 = t159 * qJD(5) + t202;
t25 = mrSges(6,1) * t86 + t164 * mrSges(6,3);
t27 = -mrSges(6,2) * t86 - t165 * mrSges(6,3);
t28 = -t165 * mrSges(7,2) + mrSges(7,3) * t86;
t156 = (t27 + t28) * t152 + (-t25 + t26) * t149 + (-t234 * t149 + t233 * t152) * qJD(5) + m(6) * (-t31 * t205 - t206 * t236 + t268) + m(7) * t267;
t14 = -t164 * Ifges(7,5) + Ifges(7,6) * t86 + t165 * Ifges(7,3);
t15 = -t164 * Ifges(6,4) - t165 * Ifges(6,2) + Ifges(6,6) * t86;
t37 = t170 * t119 - t265;
t9 = -t170 * t85 + (t171 * qJD(5) - t204) * t119 + t35;
t155 = t43 * t206 / 0.2e1 + t44 * t190 + t192 * t248 - t34 * mrSges(5,2) - Ifges(5,5) * t85 - Ifges(5,6) * t86 + t37 * t120 - t265 * t121 + t9 * t128 + t275 * t35 + t278 * t252 + (-Ifges(6,6) * t206 + t272) * t118 / 0.2e1 + t276 * t149 / 0.2e1 - (t130 / 0.2e1 + t248) * t221 + (-t123 / 0.2e1 + t122 / 0.2e1) * t217 + t274 * t216 / 0.2e1 + (-Ifges(7,6) * t252 + t15 / 0.2e1 - t14 / 0.2e1) * t152 + ((-t149 * t236 - t152 * t31) * qJD(5) + t268) * mrSges(6,3) + (t119 * t130 + t237) * t205 / 0.2e1 + t267 * mrSges(7,2) + t279 * (t119 * t190 - t219 / 0.2e1);
t141 = mrSges(7,2) * t205;
t139 = -pkin(4) - t246;
t110 = t127 - t246;
t99 = t106 + t200;
t65 = t177 * t119;
t64 = t176 * t119;
t19 = t165 * mrSges(7,1) + t164 * mrSges(7,3);
t1 = [t63 * t86 * t259 + t80 * t254 + 0.2e1 * t37 * t19 + 0.2e1 * t2 * t79 + t20 * t256 + 0.2e1 * t22 * t28 + 0.2e1 * t23 * t26 + 0.2e1 * t31 * t25 + 0.2e1 * t236 * t27 + t65 * t257 + 0.2e1 * t4 * t78 + 0.2e1 * t6 * t76 + 0.2e1 * t9 * t64 + 0.2e1 * t7 * t77 + (t2 * t22 + t23 * t4 + t37 * t9) * t260 + (t236 * t6 + t31 * t7 - t242) * t261 + (t126 * t201 + t34 * t63 - t242) * t262 - (mrSges(5,2) * t254 + mrSges(5,3) * t256 + t238 * t149 + t237 * t152) * t85 + (mrSges(5,1) * t187 + t34 * t259 - (t258 - t227) * t85 + ((2 * Ifges(5,2)) + t277) * t86 + t166) * t118 + (mrSges(5,2) * t187 + mrSges(5,3) * t257 - 0.2e1 * Ifges(5,1) * t85 + t276 * t152 + (t14 - t15) * t149 + (t258 + t240 * t152 + (-Ifges(6,6) + Ifges(7,6)) * t149) * t86 + ((-t228 + t238) * t152 + (-t237 - t282) * t149) * qJD(5)) * t119 + 0.2e1 * (t195 * t178 + (-t151 ^ 2 + t154 ^ 2) * Ifges(4,4) + (Ifges(4,1) - Ifges(4,2)) * t151 * t154) * qJD(3); (t64 + t65) * t86 + (t19 + t20) * t118 - t163 * t85 + m(6) * (-t219 * t236 + t31 * t221 + t179) + m(7) * (t118 * t9 - t22 * t219 - t23 * t221 + t37 * t86) + m(5) * (-t63 * t85 + t179) + (m(5) * t34 + t156) * t119; (-t224 + t225) * t262 + 0.4e1 * (t253 + m(7) / 0.2e1) * (-t209 * t224 + t225); t155 + t156 * t138 + (m(5) * (t150 * t34 - t153 * t35) + (-t150 * t86 + t153 * t85) * mrSges(5,3) + ((t119 * mrSges(5,3) + t65) * t150 + (-t118 * mrSges(5,3) + t163) * t153 + m(5) * (t153 * t63 - t220) + m(7) * (t22 * t213 + t23 * t215) + m(6) * (t213 * t236 - t31 * t215 - t220)) * qJD(4)) * pkin(3) + m(7) * (t110 * t9 + t37 * t99) + (Ifges(4,5) * t154 - Ifges(4,6) * t151 + (-mrSges(4,1) * t154 + mrSges(4,2) * t151) * t136) * qJD(3) + t99 * t64 + t110 * t19 + t266 * t139; t161 + m(6) * (t139 * t86 + t186) + 0.2e1 * (t194 * t253 + m(5) * (t119 * t207 - t150 * t85 - t153 * t86 + t194) / 0.2e1) * pkin(3) + m(7) * (t110 * t86 + t118 * t99 + t186) - t178 * qJD(3); t110 * t255 + 0.2e1 * t139 * t121 + 0.2e1 * t99 * t128 + 0.2e1 * t167 + (t110 * t99 + t212) * t260 + (t139 * t200 + t212) * t261 + 0.2e1 * t263 + t160; t155 + m(7) * (t106 * t37 + t127 * t9) + t156 * pkin(9) + t106 * t64 + t127 * t19 - t266 * pkin(4); m(6) * (-t235 - t250) + m(7) * (t106 * t118 + t127 * t86 - t235) + t161; (t106 + t99) * t128 + (t139 - pkin(4)) * t121 + (t110 + t127) * t120 + t167 + m(7) * (t106 * t110 + t127 * t99 + t211) + m(6) * (-pkin(4) * t200 + t211) + t263 + t160; -0.2e1 * pkin(4) * t121 + t127 * t255 + 0.2e1 * (m(7) * t127 + t128) * t106 + t160; Ifges(6,6) * t221 + m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t22) - t4 * mrSges(7,1) - t6 * mrSges(6,2) - pkin(5) * t26 + t2 * mrSges(7,3) + qJD(6) * t79 + qJ(6) * t28 + t7 * mrSges(6,1) - t278 * t119 * qJD(5) + t166; t119 * t157 - t158 * t85; t138 * t202 + t158 * t199 + (t138 * t159 + t162) * qJD(5) + t183; pkin(9) * t157 + qJD(5) * t162 + t183; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t4 + t26; t165 * m(7); t141 + (t138 * t205 + t149 * t199) * m(7); m(7) * pkin(9) * t205 + t141; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
