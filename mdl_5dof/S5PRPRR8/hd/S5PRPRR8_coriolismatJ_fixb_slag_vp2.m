% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:31
% EndTime: 2019-12-05 16:02:36
% DurationCPUTime: 1.76s
% Computational Cost: add. (2513->251), mult. (6396->374), div. (0->0), fcn. (5748->8), ass. (0->161)
t287 = m(6) * pkin(8) + mrSges(6,3);
t169 = sin(qJ(4));
t171 = cos(qJ(5));
t249 = t171 * mrSges(6,1);
t168 = sin(qJ(5));
t253 = t168 * mrSges(6,2);
t209 = t249 - t253;
t286 = t169 * t209;
t172 = cos(qJ(4));
t122 = t209 * t172;
t162 = t168 ^ 2;
t164 = t171 ^ 2;
t227 = t162 + t164;
t284 = t169 * (-0.1e1 + t227);
t160 = Ifges(6,5) * t171;
t256 = Ifges(6,6) * t168;
t283 = t160 - t256;
t264 = pkin(4) * t172;
t149 = pkin(8) * t169 + t264;
t174 = -pkin(2) - pkin(7);
t230 = t172 * t174;
t93 = t171 * t149 - t168 * t230;
t94 = t168 * t149 + t171 * t230;
t282 = -t168 * t93 + t171 * t94;
t161 = Ifges(6,4) * t171;
t257 = Ifges(6,2) * t168;
t281 = t257 - t161;
t259 = Ifges(6,4) * t168;
t145 = Ifges(6,2) * t171 + t259;
t146 = Ifges(6,1) * t168 + t161;
t266 = -t171 / 0.2e1;
t280 = -(Ifges(6,1) * t171 - t259) * t168 / 0.2e1 + t168 * t145 / 0.2e1 + t146 * t266;
t279 = -m(6) * pkin(4) - mrSges(5,1) - t209;
t166 = sin(pkin(5));
t173 = cos(qJ(2));
t170 = sin(qJ(2));
t234 = t170 * t171;
t96 = (t168 * t173 + t169 * t234) * t166;
t245 = t171 * t96;
t238 = t168 * t170;
t95 = (-t169 * t238 + t171 * t173) * t166;
t250 = t168 * t95;
t204 = t245 - t250;
t243 = t166 * t170;
t218 = -t243 / 0.2e1;
t274 = m(6) / 0.2e1;
t278 = (pkin(8) * t204 + t243 * t264) * t274 - t218 * t122 + (t245 / 0.2e1 - t250 / 0.2e1) * mrSges(6,3);
t163 = t169 ^ 2;
t165 = t172 ^ 2;
t277 = m(5) / 0.2e1;
t275 = -m(6) / 0.2e1;
t167 = cos(pkin(5));
t242 = t166 * t173;
t121 = t167 * t172 - t169 * t242;
t63 = t121 * t171 + t166 * t238;
t273 = -t63 / 0.2e1;
t272 = m(6) * t172 * t284;
t120 = t167 * t169 + t172 * t242;
t271 = t120 / 0.2e1;
t270 = -t122 / 0.2e1;
t237 = t168 * t172;
t223 = mrSges(6,3) * t237;
t134 = -mrSges(6,2) * t169 - t223;
t269 = -t134 / 0.2e1;
t248 = t171 * mrSges(6,2);
t254 = t168 * mrSges(6,1);
t144 = t248 + t254;
t268 = -t144 / 0.2e1;
t265 = -t172 / 0.2e1;
t263 = t169 * pkin(4);
t247 = t171 * t63;
t62 = -t121 * t168 + t166 * t234;
t252 = t168 * t62;
t196 = t121 - t247 + t252;
t17 = (-t120 * t284 - t172 * t196) * t274;
t262 = t17 * qJD(4);
t261 = mrSges(5,1) * t172;
t260 = mrSges(5,2) * t169;
t258 = Ifges(6,5) * t169;
t255 = Ifges(6,6) * t169;
t244 = t120 * t172;
t241 = t168 * t134;
t231 = t171 * t172;
t224 = mrSges(6,3) * t231;
t135 = mrSges(6,1) * t169 - t224;
t240 = t168 * t135;
t239 = t168 * t169;
t236 = t169 * t171;
t235 = t169 * t174;
t233 = t171 * t134;
t232 = t171 * t135;
t18 = (t121 * t169 + t242 - t244) * t243 * m(5) + (-t243 * t244 + t62 * t95 + t63 * t96) * m(6);
t229 = t18 * qJD(1);
t193 = -t241 / 0.2e1 - t232 / 0.2e1;
t212 = mrSges(6,3) * (-t164 / 0.2e1 - t162 / 0.2e1);
t178 = (t172 * t212 + t193) * t169 + t122 * t265;
t195 = t253 / 0.2e1 - t249 / 0.2e1;
t23 = t178 + t195;
t228 = t23 * qJD(2);
t226 = t163 + t165;
t222 = pkin(4) * t270;
t221 = t174 * t243;
t217 = t243 / 0.2e1;
t216 = Ifges(5,4) - t160;
t214 = t227 * t172;
t213 = m(5) * t217;
t211 = -t260 + t261;
t210 = -t169 * mrSges(5,1) - t172 * mrSges(5,2);
t207 = -Ifges(6,5) * t168 - Ifges(6,6) * t171;
t180 = t120 * t270 + t62 * t269 + t63 * t135 / 0.2e1;
t200 = t95 * mrSges(6,1) / 0.2e1 - t96 * mrSges(6,2) / 0.2e1;
t4 = (t247 / 0.2e1 - t252 / 0.2e1) * t172 * mrSges(6,3) + t180 + t200;
t142 = -pkin(8) * t172 + qJ(3) + t263;
t85 = t171 * t142 - t168 * t235;
t86 = t168 * t142 + t171 * t235;
t8 = t86 * t135 + (t174 * t122 + (-Ifges(6,4) * t237 + t258) * t168 + (Ifges(6,4) * t231 + t86 * mrSges(6,3) + t255 + (Ifges(6,1) - Ifges(6,2)) * t237) * t171) * t172 + (-t223 - t134) * t85;
t206 = -t4 * qJD(1) - t8 * qJD(2);
t205 = -t168 * t85 + t171 * t86;
t176 = (t165 * t243 + t169 * t204) * t274 + t226 * t213;
t189 = (t168 * t63 + t171 * t62) * t274;
t20 = t213 + t189 - t176;
t199 = mrSges(4,3) - t210;
t26 = t232 + t241 + m(6) * (t168 * t86 + t171 * t85) + (m(5) + m(4)) * qJ(3) + t199;
t203 = qJD(1) * t20 + qJD(2) * t26;
t15 = m(6) * t196 * t120;
t202 = t15 * qJD(1) + t17 * qJD(3);
t201 = -t93 * mrSges(6,1) / 0.2e1 + t94 * mrSges(6,2) / 0.2e1;
t198 = t172 * mrSges(6,1) + mrSges(6,3) * t236;
t197 = -t172 * mrSges(6,2) + mrSges(6,3) * t239;
t194 = t248 / 0.2e1 + t254 / 0.2e1;
t191 = t172 * t144;
t190 = -t145 / 0.4e1 + (Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.4e1) * t171;
t175 = -t62 * t198 / 0.2e1 + t197 * t273 - t121 * t191 / 0.2e1;
t183 = -t121 * t230 + t93 * t62 + t94 * t63;
t123 = t144 * t169;
t185 = t123 / 0.2e1 - t240 / 0.2e1 + t233 / 0.2e1;
t187 = -t205 + t235;
t1 = t183 * t275 + (t187 * t275 + t185) * t120 + t175 + t278;
t10 = (t205 * t274 + t185) * t172 + (-0.2e1 * t230 + t282) * t169 * t274;
t3 = -t191 * t235 - t93 * t135 - t85 * t198 - t94 * t134 - t86 * t197 - m(6) * (t85 * t93 + t86 * t94) - qJ(3) * t211 + (-t216 - t256) * t163 + (-t174 * t123 + t216 * t172 + (m(6) * t174 ^ 2 + t164 * Ifges(6,1) + Ifges(5,1) - Ifges(5,2) - Ifges(6,3)) * t169 + (Ifges(6,6) * t172 + (t257 - 0.2e1 * t161) * t169) * t168) * t172;
t188 = -t1 * qJD(1) - t3 * qJD(2) + t10 * qJD(3);
t186 = t268 + t194;
t184 = -t146 / 0.4e1 - t161 / 0.4e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.4e1) * t168;
t182 = t17 * qJD(1) + t10 * qJD(2) + qJD(3) * t272;
t181 = pkin(8) * t212 + t174 * t268;
t13 = t186 * t120;
t35 = pkin(4) * t144 - t266 * t281 + t280;
t6 = t222 + t193 * pkin(8) + t283 * t169 + (-Ifges(6,3) / 0.2e1 + t190 * t171 + (-0.5e1 / 0.4e1 * t161 + t184) * t168 + t181) * t172 + t201;
t60 = t186 * t172;
t179 = t13 * qJD(1) - t6 * qJD(2) - t60 * qJD(3) + t35 * qJD(4);
t156 = qJ(3) * t242;
t138 = t165 * t221;
t61 = t144 * t265 - t172 * t194;
t22 = t178 - t195;
t19 = t189 + (m(4) + t277) * t243 + t176;
t14 = t120 * t194 + t144 * t271;
t9 = t10 * qJD(4);
t7 = t169 * t160 / 0.4e1 + t222 - Ifges(6,5) * t236 / 0.2e1 + Ifges(6,6) * t239 / 0.2e1 + (-t255 / 0.2e1 + pkin(8) * t269) * t168 + (t258 / 0.4e1 - pkin(8) * t135 / 0.2e1) * t171 - t201 + (Ifges(6,3) / 0.2e1 + t181 + t184 * t168 + (-0.5e1 / 0.4e1 * t259 + t190) * t171) * t172;
t5 = t224 * t273 + t62 * t223 / 0.2e1 - t180 + t200;
t2 = (t120 * t187 + t183) * t274 - t120 * t233 / 0.2e1 + t218 * t260 - t175 + (-t123 + t240) * t271 + (t211 + t261) * t217 + t278;
t11 = [t18 * qJD(2) + t15 * qJD(4), t19 * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t229 + (t96 * t134 + t95 * t135 + ((-mrSges(3,2) + t199) * t173 + (-mrSges(5,3) * t226 - t144 * t165 - mrSges(3,1) + mrSges(4,2)) * t170) * t166 + m(4) * (-pkin(2) * t243 + t156) + 0.2e1 * (t163 * t221 + t138 + t156) * t277 + 0.2e1 * (t85 * t95 + t86 * t96 + t138) * t274) * qJD(2), qJD(2) * t19 + t262, t2 * qJD(2) + (t279 * t121 + (-t287 * t227 + mrSges(5,2)) * t120) * qJD(4) + t14 * qJD(5) + t202, t5 * qJD(2) + t14 * qJD(4) + (-mrSges(6,1) * t63 - mrSges(6,2) * t62) * qJD(5); qJD(3) * t20 - qJD(4) * t1 - qJD(5) * t4 - t229, qJD(3) * t26 - qJD(4) * t3 - qJD(5) * t8, qJD(5) * t22 + t203 + t9, t7 * qJD(5) + t188 + (pkin(4) * t123 + (t171 * t281 / 0.2e1 - Ifges(5,5) + t279 * t174 + t280) * t169 + (-t174 * mrSges(5,2) - pkin(8) * t144 - Ifges(5,6) - t207) * t172 + t287 * t282) * qJD(4), t22 * qJD(3) + t7 * qJD(4) + (-mrSges(6,1) * t86 - mrSges(6,2) * t85 + t172 * t207) * qJD(5) + t206; -qJD(2) * t20 + t262, qJD(5) * t23 - t203 + t9, qJD(4) * t272, (-t286 + m(6) * (pkin(8) * t214 - t263) + mrSges(6,3) * t214 + t210) * qJD(4) + t61 * qJD(5) + t182, t61 * qJD(4) - qJD(5) * t286 + t228; qJD(2) * t1 - qJD(5) * t13 - t202, qJD(5) * t6 - t188, t60 * qJD(5) - t182, -t35 * qJD(5), (-pkin(8) * t209 + t283) * qJD(5) - t179; t4 * qJD(2) + t13 * qJD(4), -qJD(3) * t23 - qJD(4) * t6 - t206, -t60 * qJD(4) - t228, t179, 0;];
Cq = t11;
