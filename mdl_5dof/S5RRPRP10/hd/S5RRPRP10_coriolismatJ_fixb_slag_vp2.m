% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRP10
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRP10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:27
% EndTime: 2019-12-31 20:09:33
% DurationCPUTime: 2.50s
% Computational Cost: add. (3685->360), mult. (7286->481), div. (0->0), fcn. (5672->4), ass. (0->172)
t189 = sin(qJ(4));
t187 = t189 ^ 2;
t191 = cos(qJ(4));
t188 = t191 ^ 2;
t248 = t188 + t187;
t280 = pkin(3) + pkin(6);
t299 = mrSges(6,3) + mrSges(5,3);
t296 = Ifges(5,5) + Ifges(6,5);
t294 = Ifges(5,6) + Ifges(6,6);
t190 = sin(qJ(2));
t298 = -t190 / 0.2e1;
t192 = cos(qJ(2));
t297 = -t192 / 0.2e1;
t295 = -Ifges(4,6) - Ifges(3,4);
t293 = Ifges(5,3) + Ifges(6,3);
t284 = m(6) * pkin(4);
t235 = mrSges(6,1) + t284;
t193 = -pkin(2) - pkin(7);
t292 = -qJ(5) + t193;
t182 = t191 * mrSges(6,1);
t262 = t189 * mrSges(6,2);
t291 = t182 - t262;
t289 = -m(4) * pkin(6) - mrSges(4,1);
t239 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t264 = Ifges(6,4) * t191;
t216 = Ifges(6,1) * t189 + t264;
t72 = Ifges(6,5) * t190 - t192 * t216;
t266 = Ifges(5,4) * t191;
t217 = Ifges(5,1) * t189 + t266;
t74 = Ifges(5,5) * t190 - t192 * t217;
t288 = t239 * t190 + t72 / 0.2e1 + t74 / 0.2e1;
t287 = -m(5) / 0.2e1;
t286 = -m(6) / 0.2e1;
t285 = m(6) / 0.2e1;
t283 = mrSges(5,2) / 0.2e1;
t282 = mrSges(6,2) / 0.2e1;
t223 = t190 * pkin(2) - qJ(3) * t192;
t124 = pkin(7) * t190 + t223;
t152 = t280 * t192;
t128 = t191 * t152;
t31 = pkin(4) * t192 + t128 + (-qJ(5) * t190 - t124) * t189;
t281 = -t31 / 0.2e1;
t249 = t191 * t192;
t98 = pkin(4) * t249 + t152;
t279 = m(6) * t98;
t253 = t189 * t190;
t130 = t192 * mrSges(6,1) - mrSges(6,3) * t253;
t278 = -t130 / 0.2e1;
t252 = t189 * t192;
t261 = t190 * mrSges(6,1);
t132 = mrSges(6,3) * t252 + t261;
t277 = -t132 / 0.2e1;
t140 = t292 * t191;
t276 = t140 / 0.2e1;
t275 = t152 / 0.2e1;
t274 = -t189 / 0.2e1;
t272 = t191 / 0.2e1;
t271 = t192 / 0.2e1;
t177 = pkin(4) * t189 + qJ(3);
t270 = m(6) * t177;
t269 = pkin(4) * t191;
t268 = mrSges(5,2) + mrSges(6,2);
t267 = Ifges(5,4) * t189;
t265 = Ifges(6,4) * t189;
t133 = t190 * mrSges(5,1) + mrSges(5,3) * t252;
t260 = t190 * mrSges(6,2);
t136 = -mrSges(6,3) * t249 - t260;
t137 = -t190 * mrSges(5,2) - mrSges(5,3) * t249;
t225 = -qJ(3) * t190 - pkin(1);
t109 = t193 * t192 + t225;
t151 = t280 * t190;
t42 = -t109 * t189 + t191 * t151;
t33 = qJ(5) * t252 + t42;
t30 = t190 * pkin(4) + t33;
t43 = t109 * t191 + t151 * t189;
t34 = -qJ(5) * t249 + t43;
t204 = m(6) * (t30 * t189 - t191 * t34);
t143 = -pkin(2) * t192 + t225;
t226 = m(4) * t143 + t192 * mrSges(4,2) - t190 * mrSges(4,3);
t9 = ((-t136 - t137) * t191 + (t132 + t133) * t189 + t204 + m(5) * (t189 * t42 - t191 * t43) - t226) * t190;
t263 = qJD(1) * t9;
t259 = t191 * mrSges(6,2);
t100 = t291 * t190;
t219 = t191 * mrSges(5,1) - t189 * mrSges(5,2);
t101 = t219 * t190;
t102 = t291 * t192;
t131 = t192 * mrSges(5,1) - mrSges(5,3) * t253;
t251 = t190 * t191;
t134 = -t192 * mrSges(6,2) + mrSges(6,3) * t251;
t135 = -t192 * mrSges(5,2) + mrSges(5,3) * t251;
t209 = t151 * mrSges(5,2) + (t216 + t217) * t298 + t296 * t297;
t214 = Ifges(6,2) * t191 + t265;
t215 = Ifges(5,2) * t191 + t267;
t210 = -t151 * mrSges(5,1) + (t214 + t215) * t298 + t294 * t297;
t68 = Ifges(6,6) * t190 - t192 * t214;
t70 = Ifges(5,6) * t190 - t192 * t215;
t237 = t68 / 0.2e1 + t70 / 0.2e1;
t238 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t45 = t191 * t124 + t189 * t152;
t35 = qJ(5) * t251 + t45;
t44 = -t124 * t189 + t128;
t97 = (-t269 - t280) * t190;
t3 = -t98 * t100 - t152 * t101 + t97 * t102 + t30 * t130 + t42 * t131 + t31 * t132 + t44 * t133 + t34 * t134 + t43 * t135 + t35 * t136 + t45 * t137 + t226 * t223 + m(6) * (t30 * t31 + t34 * t35 + t97 * t98) + m(5) * (-t151 * t152 + t42 * t44 + t43 * t45) + (-pkin(1) * mrSges(3,1) - t143 * mrSges(4,2) + t295 * t190 + (t190 * t238 + t237) * t191 + t288 * t189) * t190 + (-pkin(1) * mrSges(3,2) - t143 * mrSges(4,3) + t210 * t191 + t209 * t189 + (Ifges(4,2) - Ifges(4,3) + Ifges(3,1) - Ifges(3,2) + t293) * t190 + (-t239 * t189 - t238 * t191 - t295) * t192) * t192;
t258 = t3 * qJD(1);
t145 = t189 * mrSges(6,1) + t259;
t103 = t192 * t145;
t104 = (-t189 * mrSges(5,1) - t191 * mrSges(5,2)) * t192;
t147 = -Ifges(6,2) * t189 + t264;
t105 = t192 * t147;
t148 = -Ifges(5,2) * t189 + t266;
t106 = t192 * t148;
t149 = Ifges(6,1) * t191 - t265;
t107 = t192 * t149;
t150 = Ifges(5,1) * t191 - t267;
t108 = t192 * t150;
t174 = Ifges(6,6) * t252;
t175 = Ifges(5,6) * t252;
t234 = m(6) * (-t30 + t33);
t4 = -t98 * t103 + t152 * t104 - t43 * t133 + t33 * t136 + t42 * t137 + (t174 / 0.2e1 + t175 / 0.2e1) * t190 + (t234 - t132) * t34 + ((t105 / 0.2e1 + t106 / 0.2e1 + t30 * mrSges(6,3) + t42 * mrSges(5,3) - t288) * t191 + (t108 / 0.2e1 + t107 / 0.2e1 + t34 * mrSges(6,3) + t43 * mrSges(5,3) + (-t102 - t279) * pkin(4) + t237) * t189) * t192;
t257 = t4 * qJD(1);
t221 = t234 / 0.2e1;
t205 = t221 + t277;
t242 = mrSges(6,1) / 0.2e1 + mrSges(5,1) / 0.2e1;
t244 = t284 / 0.2e1;
t207 = t244 + t242;
t227 = t188 / 0.2e1 + t187 / 0.2e1;
t241 = t282 + t283;
t5 = (-t136 / 0.2e1 - t137 / 0.2e1 + t241 * t190) * t191 + (t133 / 0.2e1 + t207 * t190 - t205) * t189 - t299 * t192 * t227;
t256 = t5 * qJD(1);
t250 = t191 * t136;
t254 = t189 * t132;
t15 = (t204 - t250 + t254) * t192;
t255 = qJD(1) * t15;
t247 = qJD(4) * t189;
t246 = qJD(4) * t191;
t123 = (-0.1e1 / 0.2e1 - t227) * m(6);
t245 = t123 * qJD(2);
t243 = t97 * t285;
t240 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t229 = t147 / 0.2e1 + t148 / 0.2e1;
t228 = t149 / 0.2e1 + t150 / 0.2e1;
t222 = mrSges(6,3) * pkin(4) - t296;
t220 = t240 * t191;
t11 = -qJ(3) * t219 - t145 * t269 - t177 * t291 + (-t240 * t189 + t228) * t189 + (-pkin(4) * t270 + t220 + (Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t189 + t229) * t191;
t139 = t292 * t189;
t195 = t139 * mrSges(6,3) / 0.2e1 + t148 / 0.4e1 + t147 / 0.4e1 + (Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1) * t189 + t220 + (-t145 / 0.2e1 - t270 / 0.2e1) * pkin(4);
t196 = t205 * t139 + qJ(3) * t104 / 0.2e1 + t136 * t276 - t177 * t103 / 0.2e1 + t98 * t182 / 0.2e1;
t197 = (-t33 / 0.2e1 + t30 / 0.2e1) * mrSges(6,3) + t105 / 0.4e1 + t106 / 0.4e1 - t72 / 0.4e1 - t74 / 0.4e1 - t152 * mrSges(5,2) / 0.2e1 - t193 * t133 / 0.2e1 - t98 * mrSges(6,2) / 0.2e1;
t198 = (t279 / 0.2e1 + t102 / 0.2e1) * pkin(4) - t107 / 0.4e1 - t108 / 0.4e1 - t68 / 0.4e1 - t70 / 0.4e1 + mrSges(5,1) * t275 + t193 * t137 / 0.2e1;
t202 = mrSges(6,1) * t281 + t35 * t282 - t44 * mrSges(5,1) / 0.2e1 + t45 * t283;
t203 = mrSges(6,3) * t276 - t150 / 0.4e1 - t149 / 0.4e1 + (Ifges(5,2) / 0.4e1 + Ifges(6,2) / 0.4e1) * t191;
t206 = t227 * t193 * mrSges(5,3);
t2 = (m(6) * t281 + t278) * pkin(4) + (-Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + t206) * t192 + ((-0.3e1 / 0.4e1 * Ifges(5,6) - 0.3e1 / 0.4e1 * Ifges(6,6)) * t190 + t203 * t192 + t198) * t191 + ((-0.3e1 / 0.4e1 * Ifges(5,5) - 0.3e1 / 0.4e1 * Ifges(6,5)) * t190 + t195 * t192 + t197) * t189 + t196 + t202;
t213 = t2 * qJD(1) - t11 * qJD(2);
t39 = t270 + mrSges(4,3) + t268 * t191 + (mrSges(5,1) + mrSges(6,1)) * t189 + (m(5) + m(4)) * qJ(3);
t199 = m(5) * t275 + ((-t139 * t191 + t140 * t189) * t190 + t98) * t285;
t200 = (t189 * t45 + t191 * t44) * t287 + (t189 * t35 + t191 * t31) * t286;
t8 = (t278 - t131 / 0.2e1 + t242 * t192) * t191 + (-t134 / 0.2e1 - t135 / 0.2e1 - t241 * t192) * t189 + t199 + t200;
t212 = qJD(1) * t8 + qJD(2) * t39;
t201 = m(6) * ((-t139 * t192 - t30) * t191 + (t140 * t192 - t34) * t189);
t12 = (-t261 / 0.2e1 + t132 / 0.2e1) * t191 + (t260 / 0.2e1 + t136 / 0.2e1) * t189 + t243 - t201 / 0.2e1;
t29 = m(6) * (-t139 * t189 - t140 * t191) + t248 * mrSges(6,3);
t211 = -qJD(1) * t12 + qJD(2) * t29;
t110 = -m(6) * t269 - t291;
t55 = (t235 * t189 + t259) * t192;
t208 = qJD(1) * t55 + qJD(2) * t110;
t122 = t248 * t286 + t285;
t13 = t201 / 0.2e1 + t136 * t274 + t191 * t277 + t243 + (t262 / 0.2e1 - t182 / 0.2e1) * t190;
t7 = (-t241 * t189 + t191 * t242 - t289) * t192 + t199 - t200 + (t134 + t135) * t189 / 0.2e1 + (t130 + t131) * t272;
t6 = t189 * t221 + t250 / 0.2e1 - t254 / 0.2e1 + t137 * t272 + t133 * t274 + (t207 * t189 + t191 * t241) * t190 + t299 * t248 * t271;
t1 = ((-Ifges(5,5) / 0.4e1 - Ifges(6,5) / 0.4e1) * t190 + t197) * t189 + (t195 * t189 + t191 * t203 + t206) * t192 + t196 + ((-Ifges(5,6) / 0.4e1 - Ifges(6,6) / 0.4e1) * t190 + t198) * t191 + t31 * t244 + pkin(4) * t130 / 0.2e1 - t202 + t293 * t271 + t296 * t253 / 0.2e1 + t294 * t251 / 0.2e1;
t10 = [qJD(2) * t3 + qJD(3) * t9 + qJD(4) * t4 + qJD(5) * t15, t7 * qJD(3) + t1 * qJD(4) + t13 * qJD(5) + t258 + (-t177 * t100 + t140 * t130 + t139 * t134 + t97 * t145 + 0.2e1 * (t139 * t35 + t140 * t31 + t177 * t97) * t285 + (-t31 * mrSges(6,3) - t44 * mrSges(5,3) + (m(5) * t44 + t131) * t193 - t209) * t191 + (-t35 * mrSges(6,3) - t45 * mrSges(5,3) + (m(5) * t45 + t135) * t193 + t210) * t189 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) + Ifges(3,5) + t239 * t191 - t238 * t189 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(6)) * t192 + (Ifges(4,5) - Ifges(3,6) + t229 * t191 + t228 * t189 + (mrSges(3,2) - mrSges(4,3)) * pkin(6)) * t190 + (0.2e1 * t151 * t287 + t289 * t190 - t101) * qJ(3)) * qJD(2), qJD(2) * t7 + qJD(4) * t6 + t263, t257 + t1 * qJD(2) + t6 * qJD(3) + (-mrSges(5,1) * t43 - mrSges(5,2) * t42 - mrSges(6,2) * t33 - t235 * t34 + t174 + t175) * qJD(4) + t222 * t192 * t246, qJD(2) * t13 + t255; qJD(3) * t8 + qJD(4) * t2 - qJD(5) * t12 - t258, qJD(3) * t39 - qJD(4) * t11 + qJD(5) * t29, qJD(5) * t122 + t212, (-mrSges(6,2) * t140 - t235 * t139) * qJD(4) + (-mrSges(5,2) * t193 - t294) * t246 + (-mrSges(5,1) * t193 + t222) * t247 + t213, qJD(3) * t122 + t211; -qJD(2) * t8 - qJD(4) * t5 - t263, qJD(5) * t123 - t212, 0, -t256 - t268 * t246 + (-mrSges(5,1) - t235) * t247, t245; -qJD(2) * t2 + qJD(3) * t5 + qJD(5) * t55 - t257, t110 * qJD(5) - t213, t256, 0, t208; qJD(2) * t12 - qJD(4) * t55 - t255, -qJD(3) * t123 - qJD(4) * t110 - t211, -t245, -t208, 0;];
Cq = t10;
