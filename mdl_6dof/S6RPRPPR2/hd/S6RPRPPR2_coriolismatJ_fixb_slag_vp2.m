% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:20
% EndTime: 2019-03-09 02:41:24
% DurationCPUTime: 2.57s
% Computational Cost: add. (6737->304), mult. (12745->411), div. (0->0), fcn. (13243->8), ass. (0->174)
t170 = sin(pkin(10));
t173 = sin(qJ(3));
t239 = cos(pkin(10));
t265 = cos(qJ(3));
t149 = t170 * t173 - t239 * t265;
t151 = t170 * t265 + t239 * t173;
t215 = -cos(pkin(9)) * pkin(1) - pkin(2);
t152 = -t265 * pkin(3) + t215;
t176 = -t151 * qJ(5) + t152;
t82 = t149 * pkin(4) + t176;
t283 = m(6) * t82 - mrSges(6,2) * t149 - mrSges(6,3) * t151;
t261 = pkin(3) * t170;
t213 = t239 * pkin(3);
t259 = -mrSges(5,1) + mrSges(6,2);
t174 = cos(qJ(6));
t172 = sin(qJ(6));
t163 = sin(pkin(9)) * pkin(1) + pkin(7);
t147 = (-qJ(4) - t163) * t173;
t214 = t265 * t163;
t148 = t265 * qJ(4) + t214;
t280 = t170 * t147 + t239 * t148;
t281 = -t149 * pkin(5) + t280;
t167 = t173 * pkin(3);
t207 = qJ(5) * t149 + t167;
t273 = pkin(4) + pkin(8);
t71 = t273 * t151 + t207;
t37 = -t172 * t71 + t174 * t281;
t242 = t174 * t37;
t38 = t172 * t281 + t174 * t71;
t246 = t172 * t38;
t196 = t242 + t246;
t258 = mrSges(5,3) + mrSges(6,1);
t168 = t172 ^ 2;
t169 = t174 ^ 2;
t219 = t168 + t169;
t249 = Ifges(7,6) * t174;
t251 = Ifges(7,5) * t172;
t279 = Ifges(6,6) + Ifges(5,4) - t251 / 0.2e1 - t249 / 0.2e1;
t277 = 0.2e1 * t151;
t276 = m(5) / 0.2e1;
t275 = m(6) / 0.2e1;
t274 = m(7) / 0.2e1;
t234 = t149 * t174;
t141 = Ifges(7,4) * t234;
t235 = t149 * t172;
t102 = -Ifges(7,2) * t235 + t141;
t272 = -t102 / 0.4e1;
t254 = Ifges(7,4) * t172;
t158 = Ifges(7,1) * t174 - t254;
t103 = t149 * t158;
t271 = t103 / 0.4e1;
t270 = -t149 / 0.2e1;
t269 = t151 / 0.2e1;
t244 = t174 * mrSges(7,1);
t247 = t172 * mrSges(7,2);
t153 = -t244 + t247;
t268 = -t153 / 0.2e1;
t267 = -t172 / 0.2e1;
t266 = t174 / 0.2e1;
t264 = m(6) * t151;
t105 = t219 * t151;
t263 = m(7) * t105;
t262 = m(7) * t151;
t255 = mrSges(7,3) * t149;
t253 = Ifges(7,4) * t174;
t252 = Ifges(7,5) * t151;
t250 = Ifges(7,6) * t151;
t248 = t172 * mrSges(7,1);
t69 = Ifges(7,1) * t235 + t141 + t252;
t245 = t172 * t69;
t243 = t174 * mrSges(7,2);
t201 = Ifges(7,2) * t174 + t254;
t67 = t201 * t149 + t250;
t241 = t174 * t67;
t104 = pkin(4) * t151 + t207;
t232 = t151 * t172;
t107 = -t149 * mrSges(7,1) - mrSges(7,3) * t232;
t231 = t151 * t174;
t109 = t149 * mrSges(7,2) + mrSges(7,3) * t231;
t175 = (-t172 * t37 + t174 * t38) * t274 + t107 * t267 + t109 * t266 + t167 * t276;
t154 = t243 + t248;
t101 = t149 * t154;
t164 = -t213 - pkin(4);
t159 = -pkin(8) + t164;
t160 = qJ(5) + t261;
t236 = t149 * t160;
t192 = (t159 * t105 - t236) * t274 + (-t149 * t170 - t239 * t151) * pkin(3) * t276;
t177 = -t101 / 0.2e1 + t192;
t208 = t168 / 0.2e1 + t169 / 0.2e1;
t205 = t208 * mrSges(7,3);
t143 = t149 * mrSges(6,3);
t144 = t149 * mrSges(5,2);
t220 = t144 - t143;
t78 = t151 * t164 - t236;
t9 = (-t205 + t259) * t151 + 0.2e1 * (t78 / 0.4e1 - t104 / 0.4e1) * m(6) - t175 + t177 + t220;
t240 = t9 * qJD(1);
t108 = -t151 * mrSges(7,2) + mrSges(7,3) * t234;
t223 = t174 * t108;
t106 = t151 * mrSges(7,1) - mrSges(7,3) * t235;
t228 = t172 * t106;
t53 = t273 * t149 + t176;
t96 = t239 * t147 - t170 * t148;
t56 = -t151 * pkin(5) + t96;
t31 = -t172 * t53 - t174 * t56;
t32 = -t172 * t56 + t174 * t53;
t13 = (m(7) * (t172 * t31 - t174 * t32) - t223 + t228 - t283) * t151;
t238 = qJD(1) * t13;
t184 = t223 / 0.2e1 - t228 / 0.2e1;
t14 = t101 * t269 + (-t149 * t205 + t184) * t149;
t237 = qJD(1) * t14;
t182 = (t248 / 0.2e1 + t243 / 0.2e1) * t151;
t15 = t208 * t255 + t182 - t184;
t233 = t15 * qJD(1);
t188 = t247 / 0.2e1 - t244 / 0.2e1;
t181 = t188 * t151;
t225 = t174 * t106;
t227 = t172 * t108;
t186 = t227 / 0.2e1 + t225 / 0.2e1;
t17 = -t181 + t186;
t230 = t17 * qJD(1);
t226 = t172 * t109;
t224 = t174 * t107;
t183 = m(7) * (t105 - t151) * t149;
t26 = t183 / 0.2e1;
t222 = t26 * qJD(1);
t206 = -m(7) * t219 / 0.4e1;
t36 = -t263 / 0.2e1 + (t206 - m(6) / 0.2e1) * t277;
t221 = t36 * qJD(1);
t218 = t78 * t275;
t217 = t265 * mrSges(4,2);
t212 = t232 / 0.2e1;
t211 = t231 / 0.2e1;
t156 = -Ifges(7,2) * t172 + t253;
t202 = Ifges(7,1) * t172 + t253;
t210 = -t202 / 0.4e1 - t156 / 0.4e1;
t209 = t158 / 0.4e1 - t201 / 0.4e1;
t200 = -t249 - t251;
t100 = t151 * t153;
t68 = -Ifges(7,6) * t149 + t201 * t151;
t70 = -Ifges(7,5) * t149 + t202 * t151;
t99 = t153 * t149;
t1 = t82 * t143 - t152 * t144 + t265 ^ 2 * Ifges(4,4) + t215 * t217 + t56 * t99 + t281 * t100 + t37 * t106 + t31 * t107 + t38 * t108 + t32 * t109 + m(7) * (t281 * t56 + t31 * t37 + t32 * t38) + (t68 * t266 + t172 * t70 / 0.2e1 + t279 * t149) * t149 + (t215 * mrSges(4,1) - Ifges(4,4) * t173 + (m(5) * t152 + mrSges(5,1) * t149) * pkin(3) + (Ifges(4,1) - Ifges(4,2)) * t265) * t173 + (t241 / 0.2e1 + t245 / 0.2e1 + mrSges(5,2) * t167 - t82 * mrSges(6,2) + t152 * mrSges(5,1) - t279 * t151 + (Ifges(6,3) - Ifges(5,1) - Ifges(6,2) + Ifges(5,2) - Ifges(7,3)) * t149) * t151 + t283 * t104;
t185 = t226 / 0.2e1 + t224 / 0.2e1;
t197 = t172 * t32 + t174 * t31;
t4 = (t100 / 0.2e1 + t186) * t151 + (-t99 / 0.2e1 + t185) * t149 + ((t197 + t56) * t151 + (t196 - t281) * t149) * t274;
t199 = t1 * qJD(1) + t4 * qJD(2);
t140 = Ifges(7,5) * t234;
t3 = t31 * t108 - t32 * t106 + t140 * t269 + t281 * t101 + ((t69 / 0.2e1 + t102 / 0.2e1 - t31 * mrSges(7,3)) * t174 + (t103 / 0.2e1 - t67 / 0.2e1 - t32 * mrSges(7,3) - t250 / 0.2e1) * t172) * t149;
t195 = t3 * qJD(1) + t14 * qJD(2);
t194 = t4 * qJD(1) + qJD(2) * t183;
t7 = (t258 * t149 - t99) * t149 + (t258 * t151 + t225 + t227) * t151 + m(7) * (-t149 * t281 + t197 * t151) + (m(6) + m(5)) * (-t149 * t280 - t151 * t96);
t193 = qJD(1) * t7 + qJD(2) * t26;
t191 = t159 * t205;
t190 = -t37 * mrSges(7,1) / 0.2e1 + t38 * mrSges(7,2) / 0.2e1;
t187 = t160 * t101 / 0.2e1 + t281 * t268;
t42 = -t160 * t153 + (-t202 / 0.2e1 - t156 / 0.2e1) * t174 + (-t158 / 0.2e1 + t201 / 0.2e1) * t172;
t45 = (t153 / 0.2e1 - t188) * t151;
t6 = (Ifges(7,3) / 0.2e1 - t191) * t149 + (t159 * t108 / 0.2e1 - 0.3e1 / 0.4e1 * t250 - t67 / 0.4e1 + t271 + t209 * t149) * t174 + (-t159 * t106 / 0.2e1 - 0.3e1 / 0.4e1 * t252 + t272 - t69 / 0.4e1 + t210 * t149) * t172 + t187 + t190;
t180 = t6 * qJD(1) - t45 * qJD(2) + t42 * qJD(3);
t110 = mrSges(6,3) + 0.4e1 * (m(7) / 0.4e1 + m(6) / 0.4e1) * t160 + t154;
t12 = t188 * t149 + 0.2e1 * (t281 / 0.4e1 - t246 / 0.4e1 - t242 / 0.4e1) * m(7) - t185;
t84 = t263 / 0.2e1;
t50 = t84 - t262 / 0.2e1;
t178 = qJD(1) * t12 - qJD(2) * t50 + qJD(3) * t110;
t46 = t151 * t268 - t181;
t40 = t264 + t262 / 0.2e1 + t84;
t35 = t264 / 0.2e1 + t84 + (t206 - m(6) / 0.4e1) * t277;
t18 = -t181 - t186;
t16 = t182 + t184 - t219 * t255 / 0.2e1;
t11 = (-mrSges(6,1) + t188) * t149 + t185 + 0.2e1 * t280 * t275 + (t196 + t281) * t274;
t10 = t104 * t275 - t151 * t205 + t175 + t177 + t218;
t5 = t151 * t200 / 0.4e1 - t241 / 0.4e1 + t172 * t272 - t245 / 0.4e1 + t174 * t271 + Ifges(7,5) * t212 + Ifges(7,6) * t211 + Ifges(7,3) * t270 + t184 * t159 + (t210 * t172 + t209 * t174 - t191) * t149 + t187 - t190;
t2 = qJD(3) * t4 + qJD(4) * t26 + qJD(6) * t14;
t8 = [qJD(3) * t1 + qJD(4) * t7 + qJD(5) * t13 + qJD(6) * t3, t2, t10 * qJD(4) + t11 * qJD(5) + t5 * qJD(6) + t199 + (t156 * t211 + t158 * t212 + t70 * t266 + t68 * t267 + (Ifges(7,5) * t174 - Ifges(7,6) * t172) * t270 - mrSges(4,1) * t214 + t56 * t154 + Ifges(4,5) * t265 + (-mrSges(5,3) * t261 + Ifges(6,5) - Ifges(5,6)) * t151 + (-t164 * mrSges(6,1) + mrSges(5,3) * t213 + Ifges(6,4) - Ifges(5,5)) * t149 + (m(5) * t261 - mrSges(5,2) + mrSges(6,3)) * t96 + (-m(5) * t213 + m(6) * t164 + t259) * t280 + (mrSges(4,2) * t163 - Ifges(4,6)) * t173 + (m(6) * t96 + m(7) * t56 - t151 * mrSges(6,1) + t100) * t160 + (m(7) * t196 + t224 + t226) * t159 - t196 * mrSges(7,3)) * qJD(3), qJD(3) * t10 + qJD(5) * t35 + qJD(6) * t18 + t193, qJD(3) * t11 + qJD(4) * t35 + qJD(6) * t16 + t238, t5 * qJD(3) + t18 * qJD(4) + t16 * qJD(5) + (-mrSges(7,1) * t32 - mrSges(7,2) * t31 - Ifges(7,6) * t235 + t140) * qJD(6) + t195; t2, t183 * qJD(3), t40 * qJD(5) + t46 * qJD(6) + t194 + (-t173 * mrSges(4,1) - t101 - t217 + t220 + 0.2e1 * t218 + 0.2e1 * t192 + (-mrSges(7,3) * t219 + t259) * t151) * qJD(3), t222, t40 * qJD(3), qJD(3) * t46 - qJD(6) * t101 + t237; qJD(4) * t9 + qJD(5) * t12 + qJD(6) * t6 - t199, -qJD(5) * t50 - qJD(6) * t45 - t194, qJD(5) * t110 + qJD(6) * t42, t240, t178 (-t154 * t159 + t200) * qJD(6) + t180; -qJD(3) * t9 + qJD(5) * t36 - qJD(6) * t17 - t193, -t222, -t240, 0, t221, qJD(6) * t153 - t230; -qJD(3) * t12 - qJD(4) * t36 - qJD(6) * t15 - t238, t50 * qJD(3), -t178, -t221, 0, -qJD(6) * t154 - t233; -qJD(3) * t6 + qJD(4) * t17 + qJD(5) * t15 - t195, qJD(3) * t45 - t237, -t180, t230, t233, 0;];
Cq  = t8;
