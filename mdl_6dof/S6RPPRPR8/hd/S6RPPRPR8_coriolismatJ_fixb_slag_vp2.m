% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:21
% EndTime: 2019-03-09 01:55:27
% DurationCPUTime: 2.63s
% Computational Cost: add. (6512->317), mult. (11875->425), div. (0->0), fcn. (12438->6), ass. (0->187)
t169 = sin(qJ(6));
t164 = t169 ^ 2;
t170 = cos(qJ(6));
t165 = t170 ^ 2;
t227 = t165 + t164;
t303 = mrSges(7,3) * t227;
t166 = sin(pkin(9));
t167 = cos(pkin(9));
t171 = cos(qJ(4));
t276 = sin(qJ(4));
t148 = t171 * t166 + t276 * t167;
t302 = t148 * t303;
t150 = -t276 * t166 + t167 * t171;
t236 = t148 * t170;
t94 = -t150 * mrSges(7,2) + mrSges(7,3) * t236;
t247 = t170 * t94;
t237 = t148 * t169;
t93 = t150 * mrSges(7,1) - mrSges(7,3) * t237;
t254 = t169 * t93;
t189 = t254 / 0.2e1 - t247 / 0.2e1;
t218 = t165 / 0.2e1 + t164 / 0.2e1;
t298 = mrSges(7,3) * t218;
t301 = -t148 * t298 - t189;
t299 = -mrSges(5,1) + mrSges(6,2);
t271 = mrSges(5,2) - mrSges(6,3);
t297 = -t302 / 0.2e1;
t168 = -pkin(1) - qJ(3);
t270 = -pkin(7) + t168;
t153 = t270 * t166;
t154 = t270 * t167;
t98 = t171 * t153 + t276 * t154;
t72 = -t148 * pkin(5) + t98;
t242 = qJ(5) * t148;
t284 = pkin(4) + pkin(8);
t73 = t284 * t150 + t242;
t40 = -t169 * t73 + t170 * t72;
t250 = t170 * t40;
t41 = t169 * t72 + t170 * t73;
t256 = t169 * t41;
t197 = t250 + t256;
t261 = Ifges(7,6) * t170;
t263 = Ifges(7,5) * t169;
t191 = t263 / 0.2e1 + t261 / 0.2e1;
t296 = -Ifges(6,6) - Ifges(5,4) + t191;
t157 = t166 * pkin(3) + qJ(2);
t241 = qJ(5) * t150;
t55 = t284 * t148 + t157 - t241;
t97 = t276 * t153 - t171 * t154;
t70 = t150 * pkin(5) + t97;
t36 = -t169 * t55 + t170 * t70;
t37 = t169 * t70 + t170 * t55;
t198 = t36 * t169 - t37 * t170;
t202 = -pkin(4) * t148 + t241;
t90 = -t202 + t157;
t295 = m(6) * t90 - m(7) * t198 + t247 - t254;
t109 = t148 ^ 2;
t294 = t150 ^ 2;
t293 = 0.2e1 * t150;
t292 = 2 * qJD(4);
t291 = m(6) / 0.2e1;
t290 = m(6) / 0.4e1;
t289 = -m(7) / 0.2e1;
t288 = m(7) / 0.2e1;
t287 = -Ifges(7,2) / 0.2e1;
t286 = -Ifges(7,3) / 0.2e1;
t285 = t72 / 0.2e1;
t283 = -qJ(5) / 0.2e1;
t282 = -t150 / 0.2e1;
t281 = t150 / 0.2e1;
t251 = t170 * mrSges(7,2);
t258 = t169 * mrSges(7,1);
t156 = t251 + t258;
t280 = -t156 / 0.2e1;
t279 = -t169 / 0.2e1;
t278 = t169 / 0.2e1;
t277 = t170 / 0.2e1;
t100 = pkin(4) * t150 + t242;
t275 = m(6) * t100;
t274 = m(7) * t148;
t273 = t72 * mrSges(7,1);
t272 = t72 * mrSges(7,2);
t269 = mrSges(7,1) * t284;
t268 = mrSges(7,2) * t284;
t267 = mrSges(7,3) * t150;
t266 = Ifges(7,4) * t169;
t265 = Ifges(7,4) * t170;
t264 = Ifges(7,5) * t150;
t262 = Ifges(7,6) * t150;
t252 = t170 * mrSges(7,1);
t257 = t169 * mrSges(7,2);
t207 = t252 - t257;
t194 = -mrSges(6,1) - t207;
t199 = t169 * t37 + t170 * t36;
t228 = t166 ^ 2 + t167 ^ 2;
t217 = m(4) * t228;
t248 = t170 * t93;
t253 = t169 * t94;
t7 = t228 * mrSges(4,3) + (t253 + t248 + (mrSges(6,1) + mrSges(5,3)) * t150) * t150 + (mrSges(5,3) - t194) * t109 + m(7) * (-t148 * t72 + t199 * t150) - t168 * t217 + (m(6) + m(5)) * (-t148 * t98 + t150 * t97);
t260 = qJD(1) * t7;
t134 = Ifges(7,4) * t236;
t67 = Ifges(7,1) * t237 + t134 + t264;
t255 = t169 * t67;
t203 = Ifges(7,2) * t170 + t266;
t65 = t203 * t148 + t262;
t249 = t170 * t65;
t246 = t284 * t93;
t245 = t284 * t94;
t133 = Ifges(7,5) * t236;
t4 = t133 * t281 + t36 * t94 - t37 * t93 + ((t67 / 0.2e1 + t134 / 0.2e1 + t272 - t36 * mrSges(7,3)) * t170 + (-t65 / 0.2e1 + t273 - t262 / 0.2e1 - t37 * mrSges(7,3) + (-t266 / 0.2e1 + (Ifges(7,1) / 0.2e1 + t287) * t170) * t148) * t169) * t148;
t244 = t4 * qJD(1);
t138 = t148 * mrSges(6,3);
t140 = t148 * mrSges(5,2);
t214 = t227 * t284;
t224 = -t275 / 0.2e1;
t176 = (-t150 * t214 - t242) * t288 + t224;
t178 = t224 + (-t169 * t40 + t170 * t41) * t289;
t187 = t251 / 0.2e1 + t258 / 0.2e1;
t9 = t140 - t138 + (t280 - t187) * t148 + (-0.2e1 * t298 + t299) * t150 + t176 + t178;
t243 = t9 * qJD(1);
t139 = t148 * mrSges(6,2);
t99 = -t150 * mrSges(6,3) - t139;
t13 = (-t295 - t99) * t150;
t240 = qJD(1) * t13;
t229 = t148 * mrSges(5,1) - t139;
t14 = t166 * mrSges(4,1) + t167 * mrSges(4,2) + mrSges(3,3) + t271 * t150 + (m(4) + m(3)) * qJ(2) + m(5) * t157 + t229 + t295;
t239 = qJD(1) * t14;
t185 = t109 * t156;
t11 = -t185 / 0.2e1 + t301 * t150 - t187;
t238 = t11 * qJD(1);
t110 = t148 * t150;
t182 = t187 * t150;
t17 = t182 - t301;
t235 = t17 * qJD(1);
t222 = -t257 / 0.2e1;
t190 = t222 + t252 / 0.2e1;
t183 = t190 * t150;
t188 = t253 / 0.2e1 + t248 / 0.2e1;
t19 = t183 + t188;
t234 = t19 * qJD(1);
t177 = (-t227 * t294 - t109) * t288 - t217 / 0.2e1 + 0.2e1 * (t290 + m(5) / 0.4e1) * (-t109 - t294);
t216 = m(7) * t227;
t179 = -m(4) / 0.2e1 - m(5) / 0.2e1 - m(6) / 0.2e1 - t216 / 0.2e1;
t23 = t177 + t179;
t233 = t23 * qJD(1);
t32 = (t218 * m(7) + t291) * t293;
t232 = t32 * qJD(1);
t231 = mrSges(7,2) * t237 / 0.2e1 - mrSges(7,1) * t236 / 0.2e1;
t219 = t236 / 0.2e1;
t230 = mrSges(7,1) * t219 + t148 * t222;
t226 = qJD(4) * t150;
t225 = qJD(6) * t156;
t215 = t227 * t150;
t210 = t216 / 0.4e1;
t206 = Ifges(7,1) * t170 - t266;
t205 = Ifges(7,1) * t169 + t265;
t204 = -Ifges(7,2) * t169 + t265;
t66 = -Ifges(7,6) * t148 + t203 * t150;
t68 = -Ifges(7,5) * t148 + t205 * t150;
t1 = t100 * t99 - t157 * t140 + t40 * t93 + t41 * t94 + m(7) * (t36 * t40 + t37 * t41 - t70 * t72) + (t249 / 0.2e1 + t255 / 0.2e1 + t157 * mrSges(5,1) - t72 * t207 + t296 * t150 - t198 * mrSges(7,3)) * t150 + (-t150 * mrSges(6,2) + t138 + t275) * t90 + (-t36 * mrSges(7,1) + t37 * mrSges(7,2) + t70 * t207 + t66 * t277 + t68 * t278 + (Ifges(6,3) - Ifges(5,1) - Ifges(6,2) + Ifges(5,2) - Ifges(7,3)) * t150 - t296 * t148) * t148;
t3 = m(7) * (-t197 + t72) * t282 + (t207 * t281 + (t199 - t70) * t289 - t188) * t148;
t201 = t1 * qJD(1) - t3 * qJD(2);
t27 = m(7) * (-t148 * t215 + t110);
t196 = -t3 * qJD(1) + t27 * qJD(2);
t195 = -t41 * mrSges(7,2) / 0.2e1 + t40 * mrSges(7,1) / 0.2e1;
t193 = -t148 * mrSges(7,1) - t169 * t267;
t192 = t148 * mrSges(7,2) + t170 * t267;
t186 = t148 * t207;
t44 = t165 * Ifges(7,4) - qJ(5) * t207 + (-t266 + (Ifges(7,1) - Ifges(7,2)) * t170) * t169;
t45 = -t186 / 0.2e1 + t230;
t5 = (-t284 * t298 + t286) * t148 + (0.3e1 / 0.4e1 * t262 - t273 / 0.2e1 + t65 / 0.4e1 + t245 / 0.2e1 + (mrSges(7,2) * t283 + (-Ifges(7,1) / 0.2e1 + Ifges(7,2) / 0.4e1) * t170) * t148) * t170 + (0.3e1 / 0.4e1 * t264 + t272 / 0.2e1 + t67 / 0.4e1 + t134 / 0.4e1 - t246 / 0.2e1 + (0.5e1 / 0.4e1 * t265 + mrSges(7,1) * t283 + (t287 + Ifges(7,1) / 0.4e1) * t169) * t148) * t169 + t195;
t181 = t5 * qJD(1) + t45 * qJD(2) + t44 * qJD(4);
t136 = mrSges(6,3) + (m(6) + m(7)) * qJ(5) + t156;
t16 = 0.2e1 * (t72 / 0.4e1 - t256 / 0.4e1 - t250 / 0.4e1) * m(7) + t190 * t148 + t231;
t51 = 0.2e1 * (t165 / 0.4e1 + t164 / 0.4e1 - 0.1e1 / 0.4e1) * t274;
t180 = qJD(1) * t16 - qJD(2) * t51 + qJD(4) * t136;
t46 = t186 / 0.2e1 + t230;
t43 = t274 / 0.2e1 + 0.2e1 * (t291 + t210) * t148;
t33 = m(6) * t282 - t215 * t288 + (t210 + t290) * t293;
t22 = t177 - t179;
t20 = t183 - t188;
t18 = t182 - t189 + t297;
t15 = m(6) * t98 + m(7) * t285 - t148 * mrSges(6,1) + t192 * t278 + t193 * t277 + t197 * t288 + t231;
t12 = t185 / 0.2e1 - t187 + t110 * t303 / 0.2e1 + t189 * t150;
t10 = t148 * t280 - t150 * t298 + t192 * t277 + t193 * t279 + t176 - t178;
t6 = t207 * t285 - t255 / 0.4e1 + t206 * t219 - t249 / 0.4e1 - t169 * (-Ifges(7,2) * t237 + t134) / 0.4e1 + t156 * t242 / 0.2e1 + t150 * (-t261 - t263) / 0.4e1 - t203 * t236 / 0.4e1 - t245 * t277 - t246 * t279 + t148 * t286 + t191 * t150 + t195 - t297 * t284 - (t205 + t204) * t237 / 0.4e1;
t2 = t3 * qJD(4);
t8 = [qJD(2) * t14 + qJD(3) * t7 + qJD(4) * t1 + qJD(5) * t13 + qJD(6) * t4, qJD(3) * t22 + qJD(6) * t12 - t2 + t239, qJD(2) * t22 + qJD(4) * t10 + qJD(5) * t33 + qJD(6) * t20 + t260, t10 * qJD(3) + t15 * qJD(5) + t6 * qJD(6) + ((-qJ(5) * t70 - t197 * t284) * t288 + (-pkin(4) * t98 - qJ(5) * t97) * t291) * t292 + (t194 * qJ(5) + t204 * t277 + t206 * t278 + Ifges(6,5) - Ifges(5,6)) * t226 + t201 + (-t70 * t156 + t68 * t277 + t66 * t279 + (pkin(4) * mrSges(6,1) + Ifges(6,4) - Ifges(5,5) + (-Ifges(7,5) / 0.2e1 + t269) * t170 + (Ifges(7,6) / 0.2e1 - t268) * t169) * t148 + t299 * t98 + t271 * t97 - t197 * mrSges(7,3)) * qJD(4), qJD(3) * t33 + qJD(4) * t15 + qJD(6) * t18 + t240, t244 + t12 * qJD(2) + t20 * qJD(3) + t6 * qJD(4) + t18 * qJD(5) + (-mrSges(7,1) * t37 - mrSges(7,2) * t36 - Ifges(7,6) * t237 + t133) * qJD(6); qJD(3) * t23 - qJD(6) * t11 - t2 - t239, t27 * qJD(4), t233 (-t229 - t302) * qJD(4) + t43 * qJD(5) + t46 * qJD(6) + (t156 - t271) * t226 + ((-t148 * t214 + t241) * t288 + t202 * t291) * t292 + t196, t43 * qJD(4), t46 * qJD(4) + t150 * t225 - t238; -qJD(2) * t23 - qJD(4) * t9 - qJD(5) * t32 - qJD(6) * t19 - t260, -t233, 0, -t243, -t232, -qJD(6) * t207 - t234; qJD(3) * t9 + qJD(5) * t16 - qJD(6) * t5 - t201, -qJD(5) * t51 - qJD(6) * t45 - t196, t243, qJD(5) * t136 - qJD(6) * t44, t180 ((-Ifges(7,6) + t268) * t170 + (-Ifges(7,5) + t269) * t169) * qJD(6) - t181; qJD(3) * t32 - qJD(4) * t16 - qJD(6) * t17 - t240, t51 * qJD(4), t232, -t180, 0, -t225 - t235; qJD(2) * t11 + qJD(3) * t19 + qJD(4) * t5 + qJD(5) * t17 - t244, qJD(4) * t45 + t238, t234, t181, t235, 0;];
Cq  = t8;
