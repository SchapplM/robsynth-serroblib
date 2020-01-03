% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:17
% EndTime: 2019-12-31 21:16:33
% DurationCPUTime: 9.65s
% Computational Cost: add. (164879->312), mult. (343419->397), div. (0->0), fcn. (238693->10), ass. (0->124)
t265 = sin(qJ(2));
t268 = cos(qJ(2));
t286 = qJD(1) * qJD(2);
t245 = t265 * qJDD(1) + t268 * t286;
t266 = sin(qJ(1));
t269 = cos(qJ(1));
t252 = -t269 * g(1) - t266 * g(2);
t270 = qJD(1) ^ 2;
t239 = -t270 * pkin(1) + qJDD(1) * pkin(6) + t252;
t289 = t265 * t239;
t290 = pkin(2) * t270;
t202 = qJDD(2) * pkin(2) - t245 * pkin(7) - t289 + (pkin(7) * t286 + t265 * t290 - g(3)) * t268;
t228 = -t265 * g(3) + t268 * t239;
t246 = t268 * qJDD(1) - t265 * t286;
t288 = qJD(1) * t265;
t250 = qJD(2) * pkin(2) - pkin(7) * t288;
t260 = t268 ^ 2;
t203 = t246 * pkin(7) - qJD(2) * t250 - t260 * t290 + t228;
t264 = sin(qJ(3));
t291 = cos(qJ(3));
t185 = t264 * t202 + t291 * t203;
t237 = (t264 * t268 + t291 * t265) * qJD(1);
t210 = t237 * qJD(3) + t264 * t245 - t291 * t246;
t287 = qJD(1) * t268;
t236 = t264 * t288 - t291 * t287;
t221 = t236 * mrSges(4,1) + t237 * mrSges(4,2);
t258 = qJD(2) + qJD(3);
t230 = t258 * mrSges(4,1) - t237 * mrSges(4,3);
t257 = qJDD(2) + qJDD(3);
t211 = -t236 * qJD(3) + t291 * t245 + t264 * t246;
t251 = t266 * g(1) - t269 * g(2);
t280 = -qJDD(1) * pkin(1) - t251;
t212 = -t246 * pkin(2) + t250 * t288 + (-pkin(7) * t260 - pkin(6)) * t270 + t280;
t174 = (t236 * t258 - t211) * qJ(4) + (t237 * t258 + t210) * pkin(3) + t212;
t220 = t236 * pkin(3) - t237 * qJ(4);
t256 = t258 ^ 2;
t177 = -t256 * pkin(3) + t257 * qJ(4) - t236 * t220 + t185;
t261 = sin(pkin(9));
t262 = cos(pkin(9));
t226 = t262 * t237 + t261 * t258;
t163 = -0.2e1 * qJD(4) * t226 + t262 * t174 - t261 * t177;
t196 = t262 * t211 + t261 * t257;
t225 = -t261 * t237 + t262 * t258;
t161 = (t225 * t236 - t196) * pkin(8) + (t225 * t226 + t210) * pkin(4) + t163;
t164 = 0.2e1 * qJD(4) * t225 + t261 * t174 + t262 * t177;
t195 = -t261 * t211 + t262 * t257;
t215 = t236 * pkin(4) - t226 * pkin(8);
t224 = t225 ^ 2;
t162 = -t224 * pkin(4) + t195 * pkin(8) - t236 * t215 + t164;
t263 = sin(qJ(5));
t267 = cos(qJ(5));
t159 = t267 * t161 - t263 * t162;
t193 = t267 * t225 - t263 * t226;
t173 = t193 * qJD(5) + t263 * t195 + t267 * t196;
t194 = t263 * t225 + t267 * t226;
t182 = -t193 * mrSges(6,1) + t194 * mrSges(6,2);
t232 = qJD(5) + t236;
t186 = -t232 * mrSges(6,2) + t193 * mrSges(6,3);
t209 = qJDD(5) + t210;
t154 = m(6) * t159 + t209 * mrSges(6,1) - t173 * mrSges(6,3) - t194 * t182 + t232 * t186;
t160 = t263 * t161 + t267 * t162;
t172 = -t194 * qJD(5) + t267 * t195 - t263 * t196;
t187 = t232 * mrSges(6,1) - t194 * mrSges(6,3);
t155 = m(6) * t160 - t209 * mrSges(6,2) + t172 * mrSges(6,3) + t193 * t182 - t232 * t187;
t146 = t267 * t154 + t263 * t155;
t198 = -t225 * mrSges(5,1) + t226 * mrSges(5,2);
t213 = -t236 * mrSges(5,2) + t225 * mrSges(5,3);
t144 = m(5) * t163 + t210 * mrSges(5,1) - t196 * mrSges(5,3) - t226 * t198 + t236 * t213 + t146;
t214 = t236 * mrSges(5,1) - t226 * mrSges(5,3);
t282 = -t263 * t154 + t267 * t155;
t145 = m(5) * t164 - t210 * mrSges(5,2) + t195 * mrSges(5,3) + t225 * t198 - t236 * t214 + t282;
t283 = -t261 * t144 + t262 * t145;
t137 = m(4) * t185 - t257 * mrSges(4,2) - t210 * mrSges(4,3) - t236 * t221 - t258 * t230 + t283;
t184 = t291 * t202 - t264 * t203;
t229 = -t258 * mrSges(4,2) - t236 * mrSges(4,3);
t176 = -t257 * pkin(3) - t256 * qJ(4) + t237 * t220 + qJDD(4) - t184;
t165 = -t195 * pkin(4) - t224 * pkin(8) + t226 * t215 + t176;
t278 = m(6) * t165 - t172 * mrSges(6,1) + t173 * mrSges(6,2) - t193 * t186 + t194 * t187;
t274 = -m(5) * t176 + t195 * mrSges(5,1) - t196 * mrSges(5,2) + t225 * t213 - t226 * t214 - t278;
t150 = m(4) * t184 + t257 * mrSges(4,1) - t211 * mrSges(4,3) - t237 * t221 + t258 * t229 + t274;
t128 = t264 * t137 + t291 * t150;
t227 = -t268 * g(3) - t289;
t234 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t265 + Ifges(3,2) * t268) * qJD(1);
t235 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t265 + Ifges(3,4) * t268) * qJD(1);
t178 = Ifges(6,5) * t194 + Ifges(6,6) * t193 + Ifges(6,3) * t232;
t180 = Ifges(6,1) * t194 + Ifges(6,4) * t193 + Ifges(6,5) * t232;
t147 = -mrSges(6,1) * t165 + mrSges(6,3) * t160 + Ifges(6,4) * t173 + Ifges(6,2) * t172 + Ifges(6,6) * t209 - t194 * t178 + t232 * t180;
t179 = Ifges(6,4) * t194 + Ifges(6,2) * t193 + Ifges(6,6) * t232;
t148 = mrSges(6,2) * t165 - mrSges(6,3) * t159 + Ifges(6,1) * t173 + Ifges(6,4) * t172 + Ifges(6,5) * t209 + t193 * t178 - t232 * t179;
t188 = Ifges(5,5) * t226 + Ifges(5,6) * t225 + Ifges(5,3) * t236;
t190 = Ifges(5,1) * t226 + Ifges(5,4) * t225 + Ifges(5,5) * t236;
t130 = -mrSges(5,1) * t176 + mrSges(5,3) * t164 + Ifges(5,4) * t196 + Ifges(5,2) * t195 + Ifges(5,6) * t210 - pkin(4) * t278 + pkin(8) * t282 + t267 * t147 + t263 * t148 - t226 * t188 + t236 * t190;
t189 = Ifges(5,4) * t226 + Ifges(5,2) * t225 + Ifges(5,6) * t236;
t132 = mrSges(5,2) * t176 - mrSges(5,3) * t163 + Ifges(5,1) * t196 + Ifges(5,4) * t195 + Ifges(5,5) * t210 - pkin(8) * t146 - t263 * t147 + t267 * t148 + t225 * t188 - t236 * t189;
t217 = Ifges(4,4) * t237 - Ifges(4,2) * t236 + Ifges(4,6) * t258;
t218 = Ifges(4,1) * t237 - Ifges(4,4) * t236 + Ifges(4,5) * t258;
t275 = -mrSges(4,1) * t184 + mrSges(4,2) * t185 - Ifges(4,5) * t211 + Ifges(4,6) * t210 - Ifges(4,3) * t257 - pkin(3) * t274 - qJ(4) * t283 - t262 * t130 - t261 * t132 - t237 * t217 - t236 * t218;
t292 = mrSges(3,1) * t227 - mrSges(3,2) * t228 + Ifges(3,5) * t245 + Ifges(3,6) * t246 + Ifges(3,3) * qJDD(2) + pkin(2) * t128 + (t265 * t234 - t268 * t235) * qJD(1) - t275;
t139 = t262 * t144 + t261 * t145;
t244 = (-mrSges(3,1) * t268 + mrSges(3,2) * t265) * qJD(1);
t249 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t287;
t126 = m(3) * t227 + qJDD(2) * mrSges(3,1) - t245 * mrSges(3,3) + qJD(2) * t249 - t244 * t288 + t128;
t248 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t288;
t284 = t291 * t137 - t264 * t150;
t127 = m(3) * t228 - qJDD(2) * mrSges(3,2) + t246 * mrSges(3,3) - qJD(2) * t248 + t244 * t287 + t284;
t285 = -t265 * t126 + t268 * t127;
t216 = Ifges(4,5) * t237 - Ifges(4,6) * t236 + Ifges(4,3) * t258;
t120 = mrSges(4,2) * t212 - mrSges(4,3) * t184 + Ifges(4,1) * t211 - Ifges(4,4) * t210 + Ifges(4,5) * t257 - qJ(4) * t139 - t261 * t130 + t262 * t132 - t236 * t216 - t258 * t217;
t277 = -mrSges(6,1) * t159 + mrSges(6,2) * t160 - Ifges(6,5) * t173 - Ifges(6,6) * t172 - Ifges(6,3) * t209 - t194 * t179 + t193 * t180;
t272 = -mrSges(5,1) * t163 + mrSges(5,2) * t164 - Ifges(5,5) * t196 - Ifges(5,6) * t195 - pkin(4) * t146 - t226 * t189 + t225 * t190 + t277;
t124 = (-Ifges(4,2) - Ifges(5,3)) * t210 + Ifges(4,6) * t257 + t258 * t218 - t237 * t216 + Ifges(4,4) * t211 - mrSges(4,1) * t212 + mrSges(4,3) * t185 + t272 - pkin(3) * t139;
t233 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t265 + Ifges(3,6) * t268) * qJD(1);
t238 = -t270 * pkin(6) + t280;
t276 = m(4) * t212 + t210 * mrSges(4,1) + t211 * mrSges(4,2) + t236 * t229 + t237 * t230 + t139;
t116 = -mrSges(3,1) * t238 + mrSges(3,3) * t228 + Ifges(3,4) * t245 + Ifges(3,2) * t246 + Ifges(3,6) * qJDD(2) - pkin(2) * t276 + pkin(7) * t284 + qJD(2) * t235 + t264 * t120 + t291 * t124 - t233 * t288;
t119 = mrSges(3,2) * t238 - mrSges(3,3) * t227 + Ifges(3,1) * t245 + Ifges(3,4) * t246 + Ifges(3,5) * qJDD(2) - pkin(7) * t128 - qJD(2) * t234 + t291 * t120 - t264 * t124 + t233 * t287;
t273 = -m(3) * t238 + t246 * mrSges(3,1) - t245 * mrSges(3,2) - t248 * t288 + t249 * t287 - t276;
t279 = mrSges(2,1) * t251 - mrSges(2,2) * t252 + Ifges(2,3) * qJDD(1) + pkin(1) * t273 + pkin(6) * t285 + t268 * t116 + t265 * t119;
t133 = m(2) * t251 + qJDD(1) * mrSges(2,1) - t270 * mrSges(2,2) + t273;
t123 = t268 * t126 + t265 * t127;
t121 = m(2) * t252 - t270 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t285;
t117 = mrSges(2,1) * g(3) + mrSges(2,3) * t252 + t270 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t123 - t292;
t114 = -mrSges(2,2) * g(3) - mrSges(2,3) * t251 + Ifges(2,5) * qJDD(1) - t270 * Ifges(2,6) - pkin(6) * t123 - t265 * t116 + t268 * t119;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t114 - t266 * t117 - pkin(5) * (t266 * t121 + t269 * t133), t114, t119, t120, t132, t148; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t266 * t114 + t269 * t117 + pkin(5) * (t269 * t121 - t266 * t133), t117, t116, t124, t130, t147; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t279, t279, t292, -t275, Ifges(5,3) * t210 - t272, -t277;];
m_new = t1;
