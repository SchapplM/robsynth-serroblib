% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR12
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:07
% EndTime: 2019-12-31 20:29:17
% DurationCPUTime: 4.31s
% Computational Cost: add. (49076->307), mult. (100984->381), div. (0->0), fcn. (59471->8), ass. (0->116)
t287 = sin(qJ(1));
t291 = cos(qJ(1));
t263 = -t291 * g(1) - t287 * g(2);
t293 = qJD(1) ^ 2;
t236 = -t293 * pkin(1) + qJDD(1) * pkin(6) + t263;
t286 = sin(qJ(2));
t290 = cos(qJ(2));
t215 = -t290 * g(3) - t286 * t236;
t216 = -t286 * g(3) + t290 * t236;
t227 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t286 - Ifges(4,3) * t290) * qJD(1);
t230 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t286 + Ifges(3,2) * t290) * qJD(1);
t250 = (-mrSges(4,1) * t290 - mrSges(4,3) * t286) * qJD(1);
t311 = qJD(1) * qJD(2);
t309 = t290 * t311;
t252 = t286 * qJDD(1) + t309;
t310 = t286 * t311;
t253 = t290 * qJDD(1) - t310;
t249 = (-pkin(2) * t290 - qJ(3) * t286) * qJD(1);
t292 = qJD(2) ^ 2;
t312 = qJD(1) * t290;
t318 = 2 * qJD(3);
t194 = -t292 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t318 + t249 * t312 + t216;
t313 = qJD(1) * t286;
t261 = -qJD(2) * pkin(3) - pkin(7) * t313;
t316 = t290 ^ 2 * t293;
t182 = -pkin(3) * t316 - t253 * pkin(7) + qJD(2) * t261 + t194;
t201 = -qJDD(2) * pkin(2) - t292 * qJ(3) + t249 * t313 + qJDD(3) - t215;
t183 = (-t252 + t309) * pkin(7) + (-t286 * t290 * t293 - qJDD(2)) * pkin(3) + t201;
t285 = sin(qJ(4));
t289 = cos(qJ(4));
t174 = t289 * t182 + t285 * t183;
t234 = (-t285 * t290 + t286 * t289) * qJD(1);
t202 = -t234 * qJD(4) - t285 * t252 - t289 * t253;
t233 = (t285 * t286 + t289 * t290) * qJD(1);
t211 = t233 * mrSges(5,1) + t234 * mrSges(5,2);
t274 = -qJD(2) + qJD(4);
t218 = t274 * mrSges(5,1) - t234 * mrSges(5,3);
t273 = -qJDD(2) + qJDD(4);
t262 = t287 * g(1) - t291 * g(2);
t235 = -qJDD(1) * pkin(1) - t293 * pkin(6) - t262;
t305 = -t253 * pkin(2) + t235 + (-t252 - t309) * qJ(3);
t176 = -pkin(2) * t310 + t253 * pkin(3) - pkin(7) * t316 - t305 + (t261 + t318) * t313;
t203 = -t233 * qJD(4) + t289 * t252 - t285 * t253;
t168 = t176 + (t233 * t274 - t203) * pkin(8) + (t234 * t274 - t202) * pkin(4);
t212 = t233 * pkin(4) - t234 * pkin(8);
t272 = t274 ^ 2;
t170 = -t272 * pkin(4) + t273 * pkin(8) - t233 * t212 + t174;
t284 = sin(qJ(5));
t288 = cos(qJ(5));
t166 = t288 * t168 - t284 * t170;
t213 = -t284 * t234 + t288 * t274;
t179 = t213 * qJD(5) + t288 * t203 + t284 * t273;
t214 = t288 * t234 + t284 * t274;
t189 = -t213 * mrSges(6,1) + t214 * mrSges(6,2);
t200 = qJDD(5) - t202;
t226 = qJD(5) + t233;
t204 = -t226 * mrSges(6,2) + t213 * mrSges(6,3);
t162 = m(6) * t166 + t200 * mrSges(6,1) - t179 * mrSges(6,3) - t214 * t189 + t226 * t204;
t167 = t284 * t168 + t288 * t170;
t178 = -t214 * qJD(5) - t284 * t203 + t288 * t273;
t205 = t226 * mrSges(6,1) - t214 * mrSges(6,3);
t163 = m(6) * t167 - t200 * mrSges(6,2) + t178 * mrSges(6,3) + t213 * t189 - t226 * t205;
t307 = -t284 * t162 + t288 * t163;
t149 = m(5) * t174 - t273 * mrSges(5,2) + t202 * mrSges(5,3) - t233 * t211 - t274 * t218 + t307;
t173 = -t285 * t182 + t289 * t183;
t217 = -t274 * mrSges(5,2) - t233 * mrSges(5,3);
t169 = -t273 * pkin(4) - t272 * pkin(8) + t234 * t212 - t173;
t302 = -m(6) * t169 + t178 * mrSges(6,1) - t179 * mrSges(6,2) + t213 * t204 - t214 * t205;
t158 = m(5) * t173 + t273 * mrSges(5,1) - t203 * mrSges(5,3) - t234 * t211 + t274 * t217 + t302;
t143 = t285 * t149 + t289 * t158;
t184 = Ifges(6,5) * t214 + Ifges(6,6) * t213 + Ifges(6,3) * t226;
t186 = Ifges(6,1) * t214 + Ifges(6,4) * t213 + Ifges(6,5) * t226;
t156 = -mrSges(6,1) * t169 + mrSges(6,3) * t167 + Ifges(6,4) * t179 + Ifges(6,2) * t178 + Ifges(6,6) * t200 - t214 * t184 + t226 * t186;
t185 = Ifges(6,4) * t214 + Ifges(6,2) * t213 + Ifges(6,6) * t226;
t157 = mrSges(6,2) * t169 - mrSges(6,3) * t166 + Ifges(6,1) * t179 + Ifges(6,4) * t178 + Ifges(6,5) * t200 + t213 * t184 - t226 * t185;
t207 = Ifges(5,4) * t234 - Ifges(5,2) * t233 + Ifges(5,6) * t274;
t208 = Ifges(5,1) * t234 - Ifges(5,4) * t233 + Ifges(5,5) * t274;
t300 = mrSges(5,1) * t173 - mrSges(5,2) * t174 + Ifges(5,5) * t203 + Ifges(5,6) * t202 + Ifges(5,3) * t273 + pkin(4) * t302 + pkin(8) * t307 + t288 * t156 + t284 * t157 + t234 * t207 + t233 * t208;
t296 = -mrSges(4,1) * t201 + mrSges(4,3) * t194 + Ifges(4,4) * t252 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t253 - pkin(3) * t143 - t300;
t260 = mrSges(4,2) * t312 + qJD(2) * mrSges(4,3);
t301 = -m(4) * t201 + qJDD(2) * mrSges(4,1) + qJD(2) * t260 - t143;
t144 = t289 * t149 - t285 * t158;
t258 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t313;
t304 = m(4) * t194 + qJDD(2) * mrSges(4,3) + qJD(2) * t258 + t250 * t312 + t144;
t231 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t286 - Ifges(4,5) * t290) * qJD(1);
t314 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t286 + Ifges(3,4) * t290) * qJD(1) + t231;
t320 = -(t314 * t290 + (t227 - t230) * t286) * qJD(1) + mrSges(3,1) * t215 - mrSges(3,2) * t216 + Ifges(3,5) * t252 + Ifges(3,6) * t253 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t252 * mrSges(4,2) - t250 * t313 + t301) + qJ(3) * (t253 * mrSges(4,2) + t304) + t296;
t317 = mrSges(3,3) + mrSges(4,2);
t152 = t288 * t162 + t284 * t163;
t251 = (-mrSges(3,1) * t290 + mrSges(3,2) * t286) * qJD(1);
t257 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t313;
t139 = m(3) * t216 - qJDD(2) * mrSges(3,2) - qJD(2) * t257 + t251 * t312 + t317 * t253 + t304;
t259 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t312;
t140 = m(3) * t215 + qJDD(2) * mrSges(3,1) + qJD(2) * t259 - t317 * t252 + (-t250 - t251) * t313 + t301;
t308 = t290 * t139 - t286 * t140;
t150 = -m(5) * t176 + t202 * mrSges(5,1) - t203 * mrSges(5,2) - t233 * t217 - t234 * t218 - t152;
t188 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t313 + t305;
t147 = m(4) * t188 - t253 * mrSges(4,1) - t252 * mrSges(4,3) - t258 * t313 - t260 * t312 + t150;
t228 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t286 + Ifges(3,6) * t290) * qJD(1);
t229 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t286 - Ifges(4,6) * t290) * qJD(1);
t206 = Ifges(5,5) * t234 - Ifges(5,6) * t233 + Ifges(5,3) * t274;
t136 = mrSges(5,2) * t176 - mrSges(5,3) * t173 + Ifges(5,1) * t203 + Ifges(5,4) * t202 + Ifges(5,5) * t273 - pkin(8) * t152 - t284 * t156 + t288 * t157 - t233 * t206 - t274 * t207;
t297 = mrSges(6,1) * t166 - mrSges(6,2) * t167 + Ifges(6,5) * t179 + Ifges(6,6) * t178 + Ifges(6,3) * t200 + t214 * t185 - t213 * t186;
t137 = -mrSges(5,1) * t176 + mrSges(5,3) * t174 + Ifges(5,4) * t203 + Ifges(5,2) * t202 + Ifges(5,6) * t273 - pkin(4) * t152 - t234 * t206 + t274 * t208 - t297;
t298 = -mrSges(4,1) * t188 + mrSges(4,2) * t194 - pkin(3) * t150 - pkin(7) * t144 - t285 * t136 - t289 * t137;
t129 = -mrSges(3,1) * t235 + mrSges(3,3) * t216 - pkin(2) * t147 + (Ifges(3,2) + Ifges(4,3)) * t253 + (Ifges(3,4) - Ifges(4,5)) * t252 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t314 * qJD(2) + (-t228 - t229) * t313 + t298;
t299 = mrSges(4,2) * t201 - mrSges(4,3) * t188 + Ifges(4,1) * t252 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t253 - pkin(7) * t143 + qJD(2) * t227 + t289 * t136 - t285 * t137 + t229 * t312;
t131 = mrSges(3,2) * t235 - mrSges(3,3) * t215 + Ifges(3,1) * t252 + Ifges(3,4) * t253 + Ifges(3,5) * qJDD(2) - qJ(3) * t147 - qJD(2) * t230 + t228 * t312 + t299;
t295 = -m(3) * t235 + t253 * mrSges(3,1) - t252 * mrSges(3,2) - t257 * t313 + t259 * t312 - t147;
t303 = mrSges(2,1) * t262 - mrSges(2,2) * t263 + Ifges(2,3) * qJDD(1) + pkin(1) * t295 + pkin(6) * t308 + t290 * t129 + t286 * t131;
t145 = m(2) * t262 + qJDD(1) * mrSges(2,1) - t293 * mrSges(2,2) + t295;
t134 = t286 * t139 + t290 * t140;
t132 = m(2) * t263 - t293 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t308;
t127 = mrSges(2,1) * g(3) + mrSges(2,3) * t263 + t293 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t134 - t320;
t126 = -mrSges(2,2) * g(3) - mrSges(2,3) * t262 + Ifges(2,5) * qJDD(1) - t293 * Ifges(2,6) - pkin(6) * t134 - t286 * t129 + t290 * t131;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t291 * t126 - t287 * t127 - pkin(5) * (t287 * t132 + t291 * t145), t126, t131, t299, t136, t157; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t287 * t126 + t291 * t127 + pkin(5) * (t291 * t132 - t287 * t145), t127, t129, t296 + (-t286 * t227 - t290 * t231) * qJD(1), t137, t156; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t303, t303, t320, Ifges(4,5) * t252 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t253 - qJD(2) * t231 + t229 * t313 - t298, t300, t297;];
m_new = t1;
