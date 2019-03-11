% Calculate time derivative of joint inertia matrix for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:14
% EndTime: 2019-03-09 00:01:23
% DurationCPUTime: 4.12s
% Computational Cost: add. (4159->383), mult. (10303->542), div. (0->0), fcn. (9496->10), ass. (0->187)
t312 = Ifges(6,1) + Ifges(7,1);
t168 = sin(qJ(4));
t169 = sin(qJ(3));
t172 = cos(qJ(4));
t173 = cos(qJ(3));
t128 = t168 * t169 - t172 * t173;
t267 = Ifges(7,4) + Ifges(6,5);
t311 = t128 * t267;
t167 = sin(qJ(5));
t171 = cos(qJ(5));
t310 = t167 ^ 2 + t171 ^ 2;
t257 = Ifges(7,5) * t167;
t259 = Ifges(6,4) * t167;
t309 = t171 * t312 + t257 - t259;
t256 = Ifges(7,5) * t171;
t258 = Ifges(6,4) * t171;
t308 = t167 * t312 - t256 + t258;
t307 = Ifges(6,6) * t171 + t167 * t267;
t306 = Ifges(7,2) + Ifges(6,3);
t129 = t168 * t173 + t169 * t172;
t231 = qJD(5) * t167;
t291 = qJD(3) + qJD(4);
t94 = t291 * t128;
t246 = t171 * t94;
t184 = t129 * t231 + t246;
t230 = qJD(5) * t171;
t212 = t129 * t230;
t248 = t167 * t94;
t185 = t212 - t248;
t95 = t291 * t129;
t305 = t267 * t95 + (-Ifges(6,4) + Ifges(7,5)) * t185 - t312 * t184;
t262 = t129 * t309 + t311;
t143 = -mrSges(6,1) * t171 + mrSges(6,2) * t167;
t304 = -mrSges(5,1) + t143;
t303 = t309 * qJD(5);
t142 = -mrSges(7,1) * t171 - mrSges(7,3) * t167;
t302 = t142 + t143;
t301 = Ifges(7,6) * t231 + t230 * t267;
t253 = pkin(3) * qJD(4);
t224 = t172 * t253;
t299 = t310 * t224;
t277 = -pkin(9) - pkin(8);
t148 = t277 * t169;
t149 = t277 * t173;
t111 = t148 * t168 - t149 * t172;
t226 = pkin(3) * qJD(3) * t169;
t43 = pkin(4) * t95 + pkin(10) * t94 + t226;
t216 = qJD(3) * t277;
t138 = t169 * t216;
t203 = t173 * t216;
t292 = t148 * t172 + t149 * t168;
t52 = qJD(4) * t292 + t172 * t138 + t168 * t203;
t155 = -pkin(3) * t173 - pkin(2);
t84 = pkin(4) * t128 - pkin(10) * t129 + t155;
t7 = -t111 * t231 + t167 * t43 + t171 * t52 + t230 * t84;
t294 = t111 * t171 + t167 * t84;
t8 = -qJD(5) * t294 - t167 * t52 + t171 * t43;
t298 = -t8 * t167 + t171 * t7;
t165 = sin(pkin(6));
t170 = sin(qJ(2));
t233 = qJD(2) * t170;
t215 = t165 * t233;
t174 = cos(qJ(2));
t240 = t165 * t174;
t218 = t167 * t240;
t166 = cos(pkin(6));
t241 = t165 * t170;
t117 = t166 * t169 + t173 * t241;
t232 = qJD(2) * t174;
t214 = t165 * t232;
t101 = -qJD(3) * t117 - t169 * t214;
t116 = t166 * t173 - t169 * t241;
t102 = qJD(3) * t116 + t173 * t214;
t188 = t116 * t172 - t117 * t168;
t30 = qJD(4) * t188 + t101 * t168 + t102 * t172;
t74 = t116 * t168 + t117 * t172;
t17 = -qJD(5) * t218 + t167 * t30 - t171 * t215 + t230 * t74;
t249 = t167 * t17;
t59 = t167 * t74 + t171 * t240;
t60 = t171 * t74 - t218;
t297 = t230 * t59 - t231 * t60 + t249;
t2 = qJ(6) * t95 + qJD(6) * t128 + t7;
t39 = qJ(6) * t128 + t294;
t4 = -pkin(5) * t95 - t8;
t45 = -t111 * t167 + t171 * t84;
t40 = -pkin(5) * t128 - t45;
t296 = t167 * t4 + t171 * t2 + t230 * t40 - t231 * t39;
t279 = m(7) / 0.2e1;
t280 = m(6) / 0.2e1;
t295 = 0.2e1 * (t280 + t279) * t297;
t33 = mrSges(6,1) * t185 - mrSges(6,2) * t184;
t53 = qJD(4) * t111 + t138 * t168 - t172 * t203;
t293 = m(6) * t53 + t33;
t290 = m(7) * qJ(6) + mrSges(7,3);
t289 = (-mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t310) * t224;
t288 = 0.2e1 * m(6);
t287 = 0.2e1 * m(7);
t286 = -2 * mrSges(5,3);
t285 = -2 * Ifges(5,4);
t284 = 0.2e1 * t53;
t283 = -0.2e1 * t292;
t197 = t167 * mrSges(7,1) - t171 * mrSges(7,3);
t131 = t197 * qJD(5);
t282 = 0.2e1 * t131;
t281 = m(5) / 0.2e1;
t278 = t95 / 0.2e1;
t145 = Ifges(6,2) * t171 + t259;
t274 = -t145 / 0.2e1;
t272 = pkin(3) * t172;
t31 = qJD(4) * t74 - t101 * t172 + t102 * t168;
t18 = t188 * t31;
t153 = pkin(3) * t168 + pkin(10);
t16 = -qJD(5) * t59 + t167 * t215 + t171 * t30;
t250 = t16 * t171;
t266 = t171 * t224 * t60 + t153 * t250;
t35 = mrSges(6,1) * t95 + mrSges(6,3) * t184;
t36 = -t95 * mrSges(7,1) - mrSges(7,2) * t184;
t265 = t36 - t35;
t37 = -mrSges(6,2) * t95 - mrSges(6,3) * t185;
t38 = -mrSges(7,2) * t185 + mrSges(7,3) * t95;
t264 = t37 + t38;
t193 = Ifges(7,3) * t167 + t256;
t63 = Ifges(7,6) * t128 + t129 * t193;
t194 = -Ifges(6,2) * t167 + t258;
t251 = t128 * Ifges(6,6);
t64 = t129 * t194 + t251;
t263 = t63 - t64;
t244 = t129 * t167;
t86 = -mrSges(6,2) * t128 - mrSges(6,3) * t244;
t89 = -mrSges(7,2) * t244 + mrSges(7,3) * t128;
t261 = t86 + t89;
t243 = t129 * t171;
t87 = mrSges(6,1) * t128 - mrSges(6,3) * t243;
t88 = -mrSges(7,1) * t128 + mrSges(7,2) * t243;
t260 = -t87 + t88;
t255 = Ifges(6,6) * t167;
t252 = t292 * t53;
t247 = t168 * t188;
t245 = t292 * t168;
t239 = t167 * t172;
t238 = t171 * t172;
t237 = t299 * t153;
t236 = t299 * pkin(10);
t229 = qJD(6) * t171;
t228 = 0.2e1 * t169;
t227 = m(7) * t229;
t225 = t168 * t253;
t223 = t59 * t239;
t210 = -t231 / 0.2e1;
t205 = mrSges(7,2) * t229 + t301;
t202 = t165 ^ 2 * t170 * t232;
t200 = -t188 * t53 - t292 * t31;
t199 = -mrSges(4,1) * t173 + mrSges(4,2) * t169;
t198 = t167 * mrSges(6,1) + t171 * mrSges(6,2);
t192 = pkin(5) * t171 + qJ(6) * t167;
t191 = pkin(5) * t167 - qJ(6) * t171;
t187 = t304 * t225;
t186 = Ifges(7,6) * t185 - t246 * t267 + t306 * t95;
t139 = -pkin(4) - t192;
t115 = pkin(5) * t231 - qJ(6) * t230 - qJD(6) * t167;
t183 = -mrSges(7,2) * t192 - t255;
t134 = t193 * qJD(5);
t135 = t194 * qJD(5);
t144 = -Ifges(7,3) * t171 + t257;
t181 = (t144 - t145) * t231 + t308 * t230 + (-t134 + t135) * t171 + t303 * t167;
t180 = -t101 * t169 + t102 * t173 + (-t116 * t173 - t117 * t169) * qJD(3);
t179 = -m(7) * t192 + t302;
t132 = t198 * qJD(5);
t177 = -t30 * mrSges(5,2) + (-t131 - t132) * t188 + (t250 + t249 + (-t167 * t60 + t171 * t59) * qJD(5)) * mrSges(6,3) + (-mrSges(5,1) + t302) * t31 + (t250 + t297) * mrSges(7,2);
t176 = t264 * t171 + t265 * t167 + (-t167 * t261 + t171 * t260) * qJD(5) + m(7) * t296 + m(6) * (-t230 * t45 - t231 * t294 + t298);
t11 = -t191 * t94 + (qJD(5) * t192 - t229) * t129 + t53;
t23 = -Ifges(7,5) * t184 + t95 * Ifges(7,6) + Ifges(7,3) * t185;
t24 = -Ifges(6,4) * t184 - Ifges(6,2) * t185 + t95 * Ifges(6,6);
t51 = t129 * t191 - t292;
t175 = t11 * t142 + t51 * t131 - t292 * t132 - Ifges(5,6) * t95 - Ifges(5,5) * t94 - t52 * mrSges(5,2) + t63 * t231 / 0.2e1 + t64 * t210 + t212 * t274 + t304 * t53 + t307 * t278 + (-Ifges(6,6) * t231 + t301) * t128 / 0.2e1 + t305 * t167 / 0.2e1 - (t144 / 0.2e1 + t274) * t248 + (-t135 / 0.2e1 + t134 / 0.2e1) * t244 + t303 * t243 / 0.2e1 + (t24 / 0.2e1 - t23 / 0.2e1 - Ifges(7,6) * t278) * t171 + ((-t167 * t294 - t171 * t45) * qJD(5) + t298) * mrSges(6,3) + (t129 * t144 + t262) * t230 / 0.2e1 + t296 * mrSges(7,2) + t308 * (t129 * t210 - t246 / 0.2e1);
t157 = mrSges(7,2) * t230;
t154 = -pkin(4) - t272;
t133 = (mrSges(4,1) * t169 + mrSges(4,2) * t173) * qJD(3);
t122 = t139 - t272;
t112 = t115 + t225;
t100 = mrSges(5,1) * t128 + mrSges(5,2) * t129;
t79 = t198 * t129;
t78 = t197 * t129;
t47 = mrSges(5,1) * t95 - mrSges(5,2) * t94;
t32 = mrSges(7,1) * t185 + mrSges(7,3) * t184;
t15 = pkin(10) * t250;
t1 = [0.2e1 * m(5) * (t30 * t74 - t18 - t202) + 0.2e1 * m(4) * (t101 * t116 + t102 * t117 - t202) + 0.2e1 * (m(6) + m(7)) * (t16 * t60 + t17 * t59 - t18); -(t32 + t33) * t188 + t264 * t60 + t265 * t59 + (t78 + t79) * t31 + t260 * t17 + t261 * t16 + (-t128 * t30 + t129 * t31 + t188 * t94 - t74 * t95) * mrSges(5,3) + t180 * mrSges(4,3) + ((-t133 - t47) * t174 + (-t174 * mrSges(3,2) + (-mrSges(3,1) + t100 + t199) * t170) * qJD(2)) * t165 + m(5) * (t111 * t30 + t52 * t74 + (t155 * t233 - t174 * t226) * t165 + t200) + m(6) * (t16 * t294 - t17 * t45 - t59 * t8 + t60 * t7 + t200) + m(7) * (-t11 * t188 + t16 * t39 + t17 * t40 + t2 * t60 + t31 * t51 + t4 * t59) + (-pkin(2) * t215 + pkin(8) * t180) * m(4); t111 * t95 * t286 - 0.2e1 * pkin(2) * t133 + 0.2e1 * t11 * t78 + t33 * t283 + 0.2e1 * t155 * t47 + 0.2e1 * t2 * t89 + 0.2e1 * t51 * t32 + 0.2e1 * t45 * t35 + 0.2e1 * t40 * t36 + 0.2e1 * t294 * t37 + 0.2e1 * t39 * t38 + 0.2e1 * t4 * t88 + t79 * t284 + 0.2e1 * t7 * t86 + 0.2e1 * t8 * t87 + 0.2e1 * m(5) * (t111 * t52 + t155 * t226 - t252) + (t294 * t7 + t45 * t8 - t252) * t288 + (t11 * t51 + t2 * t39 + t4 * t40) * t287 - (mrSges(5,3) * t283 + t263 * t167 + t262 * t171) * t94 + (t52 * t286 - (t285 - t255) * t94 + ((2 * Ifges(5,2)) + t306) * t95 + t186) * t128 + (mrSges(5,3) * t284 - 0.2e1 * Ifges(5,1) * t94 + t305 * t171 + (t23 - t24) * t167 + (t285 + t267 * t171 + (-Ifges(6,6) + Ifges(7,6)) * t167) * t95 + ((-t251 + t263) * t171 + (-t262 - t311) * t167) * qJD(5)) * t129 + ((-Ifges(4,4) * t169 + pkin(3) * t100) * t228 + (0.2e1 * Ifges(4,4) * t173 + (Ifges(4,1) - Ifges(4,2)) * t228) * t173) * qJD(3); 0.2e1 * ((t168 * t30 - t172 * t31) * t281 + ((t223 - t247) * t280 + t223 * t279 + (t172 * t74 - t247) * t281) * qJD(4)) * pkin(3) + m(6) * (t154 * t31 + t266) + m(7) * (-t112 * t188 + t122 * t31 + t266) + t153 * t295 + t101 * mrSges(4,1) - t102 * mrSges(4,2) + t177; (Ifges(4,5) * t173 - Ifges(4,6) * t169 + pkin(8) * t199) * qJD(3) + m(7) * (t11 * t122 + t112 * t51) + t176 * t153 + (m(5) * (t168 * t52 - t172 * t53) + (-t168 * t95 + t172 * t94) * mrSges(5,3) + ((mrSges(5,3) * t129 + t79) * t168 + (-mrSges(5,3) * t128 + t167 * t260 + t171 * t261) * t172 + m(7) * (t238 * t39 + t239 * t40) + m(6) * (t238 * t294 - t239 * t45 - t245) + m(5) * (t111 * t172 - t245)) * qJD(4)) * pkin(3) + t122 * t32 + t112 * t78 + t175 + t293 * t154; 0.2e1 * t112 * t142 + t122 * t282 + 0.2e1 * t154 * t132 + 0.2e1 * t187 + (t112 * t122 + t237) * t287 + (t154 * t225 + t237) * t288 + 0.2e1 * t289 + t181; m(6) * (-pkin(4) * t31 + t15) + m(7) * (-t115 * t188 + t139 * t31 + t15) + pkin(10) * t295 + t177; t176 * pkin(10) + t139 * t32 + t115 * t78 + m(7) * (t11 * t139 + t115 * t51) + t175 - t293 * pkin(4); (t112 + t115) * t142 + (t154 - pkin(4)) * t132 + (t122 + t139) * t131 + t187 + m(7) * (t112 * t139 + t115 * t122 + t236) + m(6) * (-pkin(4) * t225 + t236) + t289 + t181; -0.2e1 * pkin(4) * t132 + t139 * t282 + 0.2e1 * (m(7) * t139 + t142) * t115 + t181; m(7) * qJD(6) * t60 + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t17 + (-mrSges(6,2) + t290) * t16; m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t39) + t2 * mrSges(7,3) - t4 * mrSges(7,1) - pkin(5) * t36 + qJD(6) * t89 + qJ(6) * t38 + Ifges(6,6) * t248 - t7 * mrSges(6,2) + t8 * mrSges(6,1) - t307 * t129 * qJD(5) + t186; t153 * t227 + (-m(7) * t191 - t197 - t198) * t224 + (t153 * t179 + t183) * qJD(5) + t205; t183 * qJD(5) + (qJD(5) * t179 + t227) * pkin(10) + t205; 0.2e1 * t290 * qJD(6); m(7) * t17; m(7) * t4 + t36; t157 + m(7) * (t153 * t230 + t167 * t224); m(7) * pkin(10) * t230 + t157; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
