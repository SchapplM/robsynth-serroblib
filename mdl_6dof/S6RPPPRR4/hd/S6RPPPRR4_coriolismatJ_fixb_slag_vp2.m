% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:19
% EndTime: 2019-03-09 01:35:23
% DurationCPUTime: 2.22s
% Computational Cost: add. (4032->318), mult. (7598->469), div. (0->0), fcn. (6239->6), ass. (0->208)
t189 = cos(qJ(6));
t285 = t189 * mrSges(7,2);
t187 = sin(qJ(6));
t289 = t187 * mrSges(7,1);
t152 = t285 + t289;
t190 = cos(qJ(5));
t251 = t190 * t152;
t328 = -t251 / 0.2e1;
t188 = sin(qJ(5));
t306 = t188 / 0.2e1;
t181 = t187 ^ 2;
t183 = t189 ^ 2;
t246 = t181 + t183;
t327 = 0.1e1 - t246;
t224 = t188 * mrSges(6,1) + t190 * mrSges(6,2);
t326 = mrSges(5,3) + t224;
t178 = Ifges(7,5) * t189;
t291 = Ifges(7,6) * t187;
t325 = Ifges(6,4) - t178 / 0.2e1 + t291 / 0.2e1;
t179 = Ifges(7,4) * t189;
t154 = Ifges(7,1) * t187 + t179;
t304 = t189 / 0.2e1;
t324 = t154 * t304;
t186 = cos(pkin(9));
t323 = t186 * qJ(2);
t184 = t190 ^ 2;
t185 = sin(pkin(9));
t170 = t184 * t185;
t182 = t188 ^ 2;
t267 = t182 * t185;
t322 = t170 + t267;
t222 = -Ifges(7,2) * t187 + t179;
t176 = t185 * qJ(2);
t191 = -pkin(1) - pkin(2);
t226 = -t186 * t191 + pkin(3) + t176;
t140 = pkin(7) + t226;
t303 = pkin(5) * t190;
t157 = -pkin(8) * t188 - t303;
t258 = t187 * t190;
t56 = -t140 * t258 + t157 * t189;
t253 = t189 * t190;
t57 = t140 * t253 + t157 * t187;
t219 = -t56 * t187 + t57 * t189;
t256 = t188 * t189;
t128 = -t187 * t185 + t186 * t256;
t272 = t128 * t189;
t259 = t187 * t188;
t126 = t185 * t189 + t186 * t259;
t274 = t126 * t187;
t216 = t272 + t274;
t264 = t186 * t188;
t32 = ((0.2e1 - t246) * t264 - t216) * t190;
t33 = t216 * t188 + (t327 * t184 - t182) * t186;
t263 = t186 * t190;
t36 = (t216 - t264) * t263;
t286 = t189 * mrSges(7,1);
t288 = t187 * mrSges(7,2);
t150 = -t286 + t288;
t177 = t188 * mrSges(6,2);
t142 = mrSges(7,2) * t190 - mrSges(7,3) * t259;
t144 = -mrSges(7,1) * t190 - mrSges(7,3) * t256;
t130 = t188 * t152;
t143 = mrSges(7,2) * t188 + mrSges(7,3) * t258;
t241 = mrSges(7,3) * t253;
t145 = -mrSges(7,1) * t188 + t241;
t305 = -t189 / 0.2e1;
t308 = t187 / 0.2e1;
t201 = t130 / 0.2e1 + t145 * t308 + t143 * t305;
t141 = t185 * t191 - qJ(4) + t323;
t301 = t190 * pkin(8);
t302 = t188 * pkin(5);
t86 = t141 + t301 - t302;
t48 = -t140 * t259 + t189 * t86;
t287 = t187 * t48;
t49 = t140 * t256 + t187 * t86;
t220 = -t189 * t49 + t287;
t314 = t126 / 0.2e1;
t318 = m(7) / 0.2e1;
t192 = (t251 * t306 + ((0.2e1 * t140 * t188 + t220) * t318 + t201) * t190) * t186 + (t126 * t56 - t128 * t57) * t318 + t144 * t314 - t128 * t142 / 0.2e1;
t127 = t185 * t256 + t186 * t187;
t273 = t127 * t189;
t125 = -t185 * t259 + t186 * t189;
t275 = t125 * t187;
t217 = t273 - t275;
t265 = t185 * t190;
t197 = m(7) * (pkin(5) * t265 + t217 * pkin(8));
t151 = -t190 * mrSges(6,1) + t177;
t313 = t151 / 0.2e1;
t6 = (-t273 / 0.2e1 + t275 / 0.2e1) * mrSges(7,3) + (t313 + t177 / 0.2e1 + (t150 / 0.2e1 - mrSges(6,1) / 0.2e1) * t190) * t185 - t197 / 0.2e1 + t192;
t321 = (t36 * qJD(2) + t33 * qJD(3) / 0.2e1 + t32 * qJD(4) / 0.2e1) * m(7) + t6 * qJD(1);
t320 = m(6) / 0.2e1;
t319 = -m(7) / 0.2e1;
t317 = -pkin(8) / 0.2e1;
t316 = mrSges(7,3) / 0.2e1;
t266 = t185 * t188;
t40 = (t217 - t266) * t190;
t315 = m(7) * t40;
t312 = -t152 / 0.2e1;
t311 = -t185 / 0.2e1;
t310 = t186 / 0.2e1;
t309 = -t187 / 0.2e1;
t307 = -t188 / 0.4e1;
t300 = -qJD(2) / 0.2e1;
t297 = m(7) * qJD(5);
t237 = t297 / 0.2e1;
t299 = t32 * t237;
t298 = t33 * t237;
t295 = Ifges(7,4) * t187;
t294 = Ifges(7,5) * t188;
t292 = Ifges(7,2) * t189;
t290 = Ifges(7,6) * t188;
t284 = t190 * mrSges(7,3);
t129 = t190 * t150;
t194 = (t129 * t310 + (-t272 / 0.2e1 - t274 / 0.2e1) * mrSges(7,3)) * t190 + t143 * t314 + t128 * t145 / 0.2e1;
t213 = t125 * mrSges(7,1) / 0.2e1 - t127 * mrSges(7,2) / 0.2e1;
t8 = t194 - t213;
t280 = t8 * qJD(1);
t279 = -mrSges(6,1) + t150;
t254 = t189 * t145;
t261 = t187 * t143;
t22 = -t261 - t254 + m(7) * (-t187 * t49 - t189 * t48) + 0.4e1 * (-m(6) / 0.4e1 - m(5) / 0.4e1) * t141 + t326;
t278 = qJD(1) * t22;
t195 = (t217 * t188 + t170) * t318 + t322 * t320;
t206 = m(7) * (-t126 * t189 + t128 * t187);
t23 = (m(5) + t320) * t185 - t206 / 0.2e1 + t195;
t277 = qJD(1) * t23;
t255 = t189 * t142;
t260 = t187 * t144;
t200 = t328 - t260 / 0.2e1 + t255 / 0.2e1;
t209 = m(7) * t219;
t245 = t182 - t184;
t11 = t245 * t140 * t318 + (t209 / 0.2e1 + t200) * t190 + (t220 * t318 + t201) * t188;
t276 = t11 * qJD(1);
t270 = t140 * t190;
t13 = (t220 * t319 - t201) * t190 + ((t219 - 0.2e1 * t270) * t318 + t200) * t188;
t271 = t13 * qJD(1);
t269 = t141 * t186;
t232 = -t254 / 0.2e1;
t235 = -t261 / 0.2e1;
t208 = t235 + t232;
t212 = t288 / 0.2e1 - t286 / 0.2e1;
t231 = t183 / 0.2e1 + t181 / 0.2e1;
t252 = t190 * t129;
t18 = t252 / 0.2e1 + (-t231 * t284 - t208) * t188 + t212;
t268 = t18 * qJD(1);
t110 = -t190 * t222 - t290;
t262 = t187 * t110;
t257 = t188 * t150;
t250 = t190 * t188;
t225 = mrSges(7,3) * t231;
t21 = t129 * t306 + t184 * t225 + t208 * t190;
t249 = t21 * qJD(1);
t58 = t327 * t245;
t248 = t58 * qJD(4);
t228 = t246 * t190;
t72 = t188 * t228 - t250;
t247 = t72 * qJD(3);
t244 = qJD(5) * t188;
t243 = t72 * t297;
t242 = -t315 / 0.2e1;
t240 = -pkin(5) * t129 / 0.2e1;
t239 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t238 = qJD(2) * t318;
t236 = t289 / 0.2e1;
t234 = -t259 / 0.2e1;
t233 = t256 / 0.2e1;
t230 = t178 - t291;
t229 = t246 * t188;
t223 = Ifges(7,1) * t189 - t295;
t221 = Ifges(7,5) * t187 + Ifges(7,6) * t189;
t168 = Ifges(7,4) * t258;
t112 = -Ifges(7,1) * t253 + t168 - t294;
t5 = -t48 * t143 + (t221 * t306 + mrSges(7,3) * t287 + t110 * t305 + t140 * t129 + t190 * t324 + (t292 * t190 + t112 + t168) * t309) * t190 + (-t241 + t145) * t49;
t218 = -t5 * qJD(1) + t21 * qJD(3);
t215 = -t56 * mrSges(7,1) / 0.2e1 + t57 * mrSges(7,2) / 0.2e1;
t211 = t285 / 0.2e1 + t236;
t210 = t168 / 0.4e1 + t112 / 0.4e1 + t145 * t317;
t153 = t292 + t295;
t207 = t153 * t309 + t324;
t205 = t190 * t221;
t87 = t140 * t170;
t10 = t127 * t143 + t125 * t145 + m(7) * (t125 * t48 + t127 * t49 + t87) + m(6) * (t140 * t267 + t269 + t87) + m(5) * t269 + m(3) * qJ(2) + mrSges(3,3) + t251 * t265 + t322 * mrSges(6,3) + (m(4) * t176 + m(5) * t226 + mrSges(4,1) - mrSges(5,2)) * t185 + (m(4) * t323 + mrSges(4,2) - t326) * t186;
t204 = -t10 * qJD(1) + qJD(3) * t242;
t203 = t312 + t211;
t109 = -Ifges(7,6) * t190 + t188 * t222;
t111 = -Ifges(7,5) * t190 + t188 * t223;
t1 = t141 * t151 + t57 * t143 + t49 * t142 + t56 * t145 + t48 * t144 + m(7) * (t48 * t56 + t49 * t57) + (t112 * t304 - t262 / 0.2e1 - t140 * t251 + t325 * t188) * t188 + (t111 * t305 + t109 * t308 - t140 * t130 - t325 * t190 + (-m(7) * t140 ^ 2 - Ifges(6,1) + Ifges(6,2) + Ifges(7,3)) * t188) * t190;
t202 = t1 * qJD(1) + t11 * qJD(3) + t13 * qJD(4);
t198 = -t72 * qJD(4) + t32 * t300 - t58 * qJD(3) / 0.2e1;
t24 = t203 * t263;
t35 = pkin(5) * t152 + t222 * t305 + t223 * t309 - t207;
t193 = t140 * t312 + pkin(8) * t225 + (t154 / 0.4e1 + t179 / 0.4e1 + t239 * t187) * t187 + (t295 / 0.2e1 + t153 / 0.4e1 - t239 * t189) * t189;
t4 = t240 + t178 * t307 + (-t294 / 0.2e1 + t210) * t189 + (-t110 / 0.4e1 + 0.3e1 / 0.4e1 * t290 + t143 * t317) * t187 + (Ifges(7,3) / 0.2e1 + t193) * t190 + t215;
t63 = t203 * t190;
t65 = t203 * t188;
t196 = t4 * qJD(1) - t24 * qJD(2) - t65 * qJD(3) + t63 * qJD(4) - t35 * qJD(5);
t66 = -t211 * t190 + t328;
t64 = mrSges(7,2) * t233 + t152 * t306 + t188 * t236;
t53 = t58 * t237;
t26 = t206 / 0.2e1 + m(6) * t311 + t195;
t25 = t211 * t263 + t251 * t310;
t19 = -t252 / 0.2e1 + t143 * t234 + t188 * t232 + t212 + t246 * t250 * t316;
t12 = t13 * qJD(5);
t9 = t194 + t213;
t7 = t185 * t313 + t129 * t311 - mrSges(6,2) * t266 / 0.2e1 + mrSges(6,1) * t265 / 0.2e1 + t273 * t316 - mrSges(7,3) * t275 / 0.2e1 + t197 / 0.2e1 + t192;
t3 = -t262 / 0.4e1 + t240 + t230 * t307 + pkin(8) * t235 + Ifges(7,5) * t233 + Ifges(7,6) * t234 + t210 * t189 - t215 + (-Ifges(7,3) / 0.2e1 + t193) * t190;
t2 = t11 * qJD(5) + t21 * qJD(6) + t40 * t238;
t14 = [qJD(2) * t10 + qJD(4) * t22 + qJD(5) * t1 - qJD(6) * t5, 0.2e1 * ((t126 * t125 - t128 * t127) * t318 + (t184 * t319 + (-t182 - t184 + 0.1e1) * t320) * t186 * t185) * qJD(2) + t26 * qJD(4) + t7 * qJD(5) + t9 * qJD(6) - t204, t2, qJD(2) * t26 + qJD(6) * t19 + t12 + t278, t7 * qJD(2) + t3 * qJD(6) + (Ifges(6,5) + (-m(7) * pkin(5) + t279) * t140 + t207) * t244 + t202 + (-mrSges(6,2) * t270 + t109 * t304 + t111 * t308 - pkin(5) * t130 - t205 / 0.2e1 + Ifges(6,6) * t190 + (t209 + t255 - t260) * pkin(8) + t219 * mrSges(7,3)) * qJD(5), t9 * qJD(2) + t19 * qJD(4) + t3 * qJD(5) + (-mrSges(7,1) * t49 - mrSges(7,2) * t48 + t205) * qJD(6) + t218; -t23 * qJD(4) + t6 * qJD(5) + t8 * qJD(6) + t204, t36 * t297, qJD(1) * t242 + t298, -t277 + t299, t25 * qJD(6) + (-t257 + m(7) * (-t246 * t301 + t302) - t246 * t284 + t224) * qJD(5) * t186 + t321, t280 + t25 * qJD(5) + (mrSges(7,1) * t128 - mrSges(7,2) * t126) * qJD(6); t2, qJD(1) * t315 / 0.2e1 + t298, -t243, t53, t276 + t33 * t238 - m(7) * t247 + t248 * t318 + (t129 + m(7) * (-pkin(8) * t229 - t303) - mrSges(7,3) * t229 + t151) * qJD(5) + t64 * qJD(6), t64 * qJD(5) + t129 * qJD(6) + t249; qJD(2) * t23 - qJD(6) * t18 + t12 - t278, t277 + t299, t53, t243, t271 + t66 * qJD(6) + t279 * t244 + (t246 * mrSges(7,3) - mrSges(6,2)) * qJD(5) * t190 + ((pkin(8) * t228 - t302) * qJD(5) - t198) * m(7), t66 * qJD(5) + qJD(6) * t257 - t268; -qJD(2) * t6 + qJD(6) * t4 - t202, -t24 * qJD(6) - t321, -t276 - t65 * qJD(6) + (t247 + t33 * t300 - t248 / 0.2e1) * m(7), m(7) * t198 + t63 * qJD(6) - t271, -t35 * qJD(6) (pkin(8) * t150 + t230) * qJD(6) + t196; -qJD(2) * t8 + qJD(4) * t18 - qJD(5) * t4 - t218, t24 * qJD(5) - t280, t65 * qJD(5) - t249, -t63 * qJD(5) + t268, -t196, 0;];
Cq  = t14;
