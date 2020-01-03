% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:25
% EndTime: 2019-12-31 20:12:31
% DurationCPUTime: 2.72s
% Computational Cost: add. (3547->345), mult. (7036->443), div. (0->0), fcn. (5415->4), ass. (0->165)
t202 = sin(qJ(2));
t203 = cos(qJ(4));
t204 = cos(qJ(2));
t250 = t203 * t204;
t246 = mrSges(5,3) * t250;
t149 = -t202 * mrSges(5,2) - t246;
t191 = t202 * mrSges(6,3);
t150 = -mrSges(6,2) * t250 + t191;
t304 = t149 + t150;
t289 = pkin(3) + pkin(6);
t309 = mrSges(6,2) + mrSges(5,3);
t307 = Ifges(6,4) + Ifges(5,5);
t314 = Ifges(5,6) - Ifges(6,6);
t201 = sin(qJ(4));
t199 = t201 ^ 2;
t200 = t203 ^ 2;
t313 = (t200 / 0.2e1 + t199 / 0.2e1) * t309;
t194 = Ifges(6,5) * t203;
t227 = Ifges(6,1) * t201 - t194;
t101 = Ifges(6,4) * t202 - t204 * t227;
t271 = Ifges(5,4) * t203;
t228 = Ifges(5,1) * t201 + t271;
t103 = Ifges(5,5) * t202 - t204 * t228;
t312 = t101 + t103;
t192 = t203 * mrSges(6,3);
t265 = t201 * mrSges(6,1);
t301 = t192 - t265;
t280 = -t202 / 0.2e1;
t311 = -t204 / 0.2e1;
t276 = t204 / 0.2e1;
t310 = mrSges(5,1) + mrSges(6,1);
t308 = -Ifges(3,4) - Ifges(4,6);
t306 = Ifges(6,2) + Ifges(5,3);
t252 = t202 * t203;
t148 = -t204 * mrSges(5,2) + mrSges(5,3) * t252;
t151 = mrSges(6,2) * t252 + t204 * mrSges(6,3);
t305 = t148 + t151;
t270 = Ifges(6,5) * t201;
t166 = Ifges(6,1) * t203 + t270;
t272 = Ifges(5,4) * t201;
t167 = Ifges(5,1) * t203 - t272;
t303 = t166 + t167;
t259 = qJ(5) * t203;
t224 = pkin(4) * t201 - t259;
t264 = t203 * mrSges(5,2);
t266 = t201 * mrSges(5,1);
t230 = -t264 - t266;
t302 = -m(6) * t224 + t230 + t301;
t193 = Ifges(6,6) * t203;
t300 = -Ifges(5,6) * t203 + t193;
t164 = Ifges(6,3) * t201 + t194;
t158 = qJ(3) + t224;
t275 = m(6) * t158;
t298 = -t275 + t301;
t241 = -Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t190 = m(6) * qJ(5) + mrSges(6,3);
t297 = -m(6) / 0.2e1;
t296 = m(6) / 0.2e1;
t295 = mrSges(6,1) / 0.2e1;
t294 = mrSges(6,3) / 0.2e1;
t205 = -pkin(2) - pkin(7);
t234 = -qJ(3) * t202 - pkin(1);
t130 = t204 * t205 + t234;
t251 = t203 * t130;
t168 = t289 * t202;
t256 = t201 * t168;
t56 = t251 + t256;
t291 = t56 / 0.2e1;
t163 = pkin(4) * t203 + qJ(5) * t201;
t169 = t289 * t204;
t75 = t163 * t204 + t169;
t290 = t75 / 0.2e1;
t233 = t202 * pkin(2) - qJ(3) * t204;
t142 = pkin(7) * t202 + t233;
t57 = -t142 * t201 + t169 * t203;
t52 = -pkin(4) * t204 - t57;
t288 = m(6) * t52;
t287 = -qJ(3) / 0.2e1;
t255 = t201 * t202;
t186 = mrSges(6,2) * t255;
t263 = t204 * mrSges(6,1);
t145 = t186 - t263;
t286 = t145 / 0.2e1;
t285 = -t151 / 0.2e1;
t284 = t169 / 0.2e1;
t283 = -t201 / 0.2e1;
t281 = t201 / 0.2e1;
t278 = -t203 / 0.2e1;
t277 = t203 / 0.2e1;
t274 = m(6) * t201;
t254 = t201 * t204;
t247 = mrSges(5,3) * t254;
t146 = t202 * mrSges(5,1) + t247;
t159 = -pkin(2) * t204 + t234;
t235 = m(4) * t159 + t204 * mrSges(4,2) - t202 * mrSges(4,3);
t245 = mrSges(6,2) * t254;
t147 = -t202 * mrSges(6,1) - t245;
t257 = t201 * t147;
t253 = t202 * qJ(5);
t48 = t56 + t253;
t55 = -t130 * t201 + t168 * t203;
t49 = -pkin(4) * t202 - t55;
t9 = -t146 * t255 + (t257 + t304 * t203 - m(6) * (-t49 * t201 - t48 * t203) - m(5) * (t55 * t201 - t56 * t203) + t235) * t202;
t267 = qJD(1) * t9;
t229 = t203 * mrSges(6,1) + t201 * mrSges(6,3);
t126 = t229 * t202;
t231 = t203 * mrSges(5,1) - t201 * mrSges(5,2);
t127 = t231 * t202;
t128 = t229 * t204;
t144 = t204 * mrSges(5,1) - mrSges(5,3) * t255;
t215 = t168 * mrSges(5,2) + t307 * t311 + (t227 + t228) * t280;
t225 = -Ifges(6,3) * t203 + t270;
t226 = Ifges(5,2) * t203 + t272;
t218 = Ifges(6,6) * t276 + t202 * t225 / 0.2e1 + Ifges(5,6) * t311 + t226 * t280 - t168 * mrSges(5,1);
t242 = Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t58 = t203 * t142 + t201 * t169;
t50 = qJ(5) * t204 + t58;
t74 = (-t163 - t289) * t202;
t97 = Ifges(6,6) * t202 - t204 * t225;
t99 = Ifges(5,6) * t202 - t204 * t226;
t3 = -t75 * t126 - t169 * t127 + t74 * t128 + t55 * t144 + t49 * t145 + t57 * t146 + t52 * t147 + t56 * t148 + t58 * t149 + t50 * t150 + t48 * t151 + t235 * t233 + m(6) * (t48 * t50 + t49 * t52 + t74 * t75) + m(5) * (-t168 * t169 + t55 * t57 + t56 * t58) + (-pkin(1) * mrSges(3,1) - t159 * mrSges(4,2) + t308 * t202 + (-t97 / 0.2e1 + t99 / 0.2e1 + t241 * t202) * t203 + (t101 / 0.2e1 + t103 / 0.2e1 + t242 * t202) * t201) * t202 + (-pkin(1) * mrSges(3,2) - t159 * mrSges(4,3) + t218 * t203 + t215 * t201 + (Ifges(3,1) - Ifges(3,2) + Ifges(4,2) - Ifges(4,3) + t306) * t202 + (-t242 * t201 - t241 * t203 - t308) * t204) * t204;
t262 = t3 * qJD(1);
t129 = t224 * t204;
t165 = -Ifges(5,2) * t201 + t271;
t4 = -t48 * t245 + (m(6) * t75 + t128) * t129 + (-m(6) * t49 + t146 - t147 - t247) * t56 + (-m(6) * t48 - t246 - t304) * t55 + (-t75 * t301 - t169 * t230 + t97 * t281 + t99 * t283 + t49 * t203 * mrSges(6,2) + (t314 * t201 - t307 * t203) * t280 + t312 * t277 + (t303 * t283 + (t165 - t164) * t278) * t204) * t204;
t261 = t4 * qJD(1);
t209 = t257 / 0.2e1 + t146 * t283 + ((t48 - t56) * t203 + (t49 + t55) * t201) * t296 + t304 * t277;
t210 = (t224 * t296 - t192 / 0.2e1 + t264 / 0.2e1 + t266 / 0.2e1 + t265 / 0.2e1) * t202;
t5 = -t204 * t313 - t209 + t210;
t260 = t5 * qJD(1);
t17 = t202 * t150 + m(6) * (t202 * t48 + t254 * t75) + t128 * t254;
t258 = qJD(1) * t17;
t249 = m(6) * t255;
t248 = qJD(5) * t274;
t244 = t295 + mrSges(5,1) / 0.2e1;
t243 = t294 - mrSges(5,2) / 0.2e1;
t238 = t301 * t281;
t232 = pkin(4) * mrSges(6,2) - t307;
t223 = t201 * t50 - t203 * t52;
t222 = t201 * t58 + t203 * t57;
t11 = -qJ(3) * t231 - t158 * t229 + t226 * t283 + (t164 - t227) * t278 + (t228 + t165) * t277 + t298 * t163 + (t225 + t303) * t281;
t206 = t209 * t205 + (-t158 * t129 + t163 * t75) * t296 + t129 * t301 / 0.2e1 + t163 * t128 / 0.2e1 + t231 * t284 + t229 * t290 + ((t291 - t48 / 0.2e1) * mrSges(6,2) + t97 / 0.4e1 - t99 / 0.4e1) * t203 + ((-t55 / 0.2e1 - t49 / 0.2e1) * mrSges(6,2) - t312 / 0.4e1) * t201;
t214 = Ifges(5,2) / 0.4e1 + Ifges(6,3) / 0.4e1 - Ifges(5,1) / 0.4e1 - Ifges(6,1) / 0.4e1;
t207 = (-t194 / 0.4e1 + t165 / 0.4e1 - t164 / 0.4e1 + mrSges(5,1) * t287 - t158 * mrSges(6,1) / 0.2e1 - t214 * t201) * t201 + (-t167 / 0.4e1 - t166 / 0.4e1 + mrSges(5,2) * t287 + t158 * t294 + t214 * t203 + (Ifges(5,4) - 0.3e1 / 0.4e1 * Ifges(6,5)) * t201) * t203 + t205 * t313;
t208 = (-pkin(4) * t52 + qJ(5) * t50) * t297 + pkin(4) * t286 + qJ(5) * t285 - t50 * mrSges(6,3) / 0.2e1 + t52 * t295 - t57 * mrSges(5,1) / 0.2e1 + t58 * mrSges(5,2) / 0.2e1;
t2 = t206 + (t193 / 0.4e1 + (-0.3e1 / 0.4e1 * Ifges(5,6) + Ifges(6,6) / 0.2e1) * t203 + (-0.3e1 / 0.4e1 * Ifges(6,4) - 0.3e1 / 0.4e1 * Ifges(5,5)) * t201) * t202 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + t207) * t204 + t208;
t221 = -t2 * qJD(1) + t11 * qJD(2);
t36 = t264 + t275 - t192 + mrSges(4,3) + t310 * t201 + (m(5) + m(4)) * qJ(3);
t212 = -m(5) * t222 / 0.2e1 + t223 * t297;
t216 = m(5) * t284 + m(6) * t290;
t8 = (t286 - t144 / 0.2e1 + t244 * t204) * t203 + (t285 - t148 / 0.2e1 + t243 * t204) * t201 + t212 + t216;
t220 = qJD(1) * t8 + qJD(2) * t36;
t211 = (-t203 * t75 + (t158 * t204 + t202 * t205) * t201) * t296 + t128 * t278;
t14 = -t186 + (-t238 + t295) * t204 - t288 / 0.2e1 + t211;
t60 = t298 * t203;
t219 = -qJD(1) * t14 - qJD(2) * t60;
t21 = -t191 + 0.2e1 * (t56 / 0.4e1 - t253 / 0.2e1 - t251 / 0.4e1 - t256 / 0.4e1) * m(6);
t217 = qJD(1) * t21 - qJD(4) * t190;
t152 = (m(6) * t205 - mrSges(6,2)) * t201;
t19 = (t56 + 0.2e1 * t253) * t296 + m(6) * t291 + t150;
t15 = -t204 * t238 + t288 / 0.2e1 - t263 / 0.2e1 + t211;
t7 = t145 * t278 + t144 * t277 + (m(4) * pkin(6) + t201 * t243 + t203 * t244 + mrSges(4,1)) * t204 - t212 + t216 + t305 * t281;
t6 = t209 + t210 + t309 * (t199 + t200) * t276;
t1 = t207 * t204 + t206 - t208 + (-t201 * t307 + t300) * t202 / 0.4e1 + t306 * t276 + t241 * t252 + t307 * t255 / 0.2e1;
t10 = [qJD(2) * t3 - qJD(3) * t9 - qJD(4) * t4 + qJD(5) * t17, t7 * qJD(3) + t1 * qJD(4) + t15 * qJD(5) + t262 + (-qJ(3) * t127 - t158 * t126 - t74 * t301 + (t52 * mrSges(6,2) - t57 * mrSges(5,3) + (t144 - t145) * t205 - t215) * t203 + (-t50 * mrSges(6,2) - t58 * mrSges(5,3) + t205 * t305 + t218) * t201 + 0.2e1 * (t158 * t74 + t205 * t223) * t296 + m(5) * (-qJ(3) * t168 + t205 * t222) + (-pkin(2) * mrSges(4,1) - Ifges(4,4) + Ifges(3,5) + t242 * t203 - t241 * t201 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(6)) * t204 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + (-t164 / 0.2e1 + t165 / 0.2e1) * t203 + (t166 / 0.2e1 + t167 / 0.2e1) * t201 + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(6)) * t202) * qJD(2), qJD(2) * t7 + qJD(4) * t6 - t267, t1 * qJD(2) + t6 * qJD(3) + t19 * qJD(5) - t261 + ((-m(6) * pkin(4) - t310) * t56 + (-mrSges(5,2) + t190) * t55 + (t232 * t203 + (qJ(5) * mrSges(6,2) + t314) * t201) * t204) * qJD(4), qJD(2) * t15 + qJD(4) * t19 + t258; qJD(3) * t8 + qJD(4) * t2 + qJD(5) * t14 - t262, qJD(3) * t36 - qJD(4) * t11 + qJD(5) * t60, t220, t152 * qJD(5) - t221 + (-mrSges(6,2) * t259 + t232 * t201 + t302 * t205 + t300) * qJD(4), qJD(4) * t152 - t219; -qJD(2) * t8 - qJD(4) * t5 + t202 * t248 + t267, -t220, 0, qJD(4) * t302 + t248 - t260, (qJD(1) * t202 + qJD(4)) * t274; -qJD(2) * t2 + qJD(3) * t5 - qJD(5) * t21 + t261, t221, t260, t190 * qJD(5), -t217; -qJD(2) * t14 - qJD(3) * t249 + qJD(4) * t21 - t258, t219, -qJD(1) * t249, t217, 0;];
Cq = t10;
