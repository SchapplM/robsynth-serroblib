% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR2
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
% m [6x1]
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:42
% EndTime: 2022-01-20 11:30:42
% DurationCPUTime: 0.47s
% Computational Cost: add. (5275->110), mult. (6167->147), div. (0->0), fcn. (3316->10), ass. (0->59)
t257 = qJD(1) + qJD(2);
t253 = qJD(3) + t257;
t261 = sin(qJ(5));
t279 = t253 * t261;
t265 = cos(qJ(5));
t278 = t253 * t265;
t264 = sin(qJ(1));
t268 = cos(qJ(1));
t274 = t264 * g(1) - t268 * g(2);
t246 = qJDD(1) * pkin(1) + t274;
t272 = -t268 * g(1) - t264 * g(2);
t247 = -qJD(1) ^ 2 * pkin(1) + t272;
t263 = sin(qJ(2));
t267 = cos(qJ(2));
t231 = t267 * t246 - t263 * t247;
t256 = qJDD(1) + qJDD(2);
t228 = t256 * pkin(2) + t231;
t232 = t263 * t246 + t267 * t247;
t255 = t257 ^ 2;
t229 = -t255 * pkin(2) + t232;
t262 = sin(qJ(3));
t266 = cos(qJ(3));
t223 = t266 * t228 - t262 * t229;
t251 = t253 ^ 2;
t252 = qJDD(3) + t256;
t220 = t252 * pkin(3) + t223;
t224 = t262 * t228 + t266 * t229;
t221 = -t251 * pkin(3) + t224;
t259 = sin(pkin(9));
t260 = cos(pkin(9));
t217 = t259 * t220 + t260 * t221;
t214 = -t251 * pkin(4) + t252 * pkin(8) + t217;
t258 = -g(3) + qJDD(4);
t211 = -t261 * t214 + t265 * t258;
t238 = (-mrSges(6,1) * t265 + mrSges(6,2) * t261) * t253;
t275 = qJD(5) * t253;
t239 = t261 * t252 + t265 * t275;
t245 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t278;
t209 = m(6) * t211 + qJDD(5) * mrSges(6,1) - t239 * mrSges(6,3) + qJD(5) * t245 - t238 * t279;
t212 = t265 * t214 + t261 * t258;
t240 = t265 * t252 - t261 * t275;
t244 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t279;
t210 = m(6) * t212 - qJDD(5) * mrSges(6,2) + t240 * mrSges(6,3) - qJD(5) * t244 + t238 * t278;
t273 = -t261 * t209 + t265 * t210;
t201 = m(5) * t217 - t251 * mrSges(5,1) - t252 * mrSges(5,2) + t273;
t216 = t260 * t220 - t259 * t221;
t213 = -t252 * pkin(4) - t251 * pkin(8) - t216;
t271 = -m(6) * t213 + t240 * mrSges(6,1) - t239 * mrSges(6,2) - t244 * t279 + t245 * t278;
t206 = m(5) * t216 + t252 * mrSges(5,1) - t251 * mrSges(5,2) + t271;
t276 = t259 * t201 + t260 * t206;
t197 = m(4) * t223 + t252 * mrSges(4,1) - t251 * mrSges(4,2) + t276;
t198 = m(4) * t224 - t251 * mrSges(4,1) - t252 * mrSges(4,2) + t260 * t201 - t259 * t206;
t277 = t266 * t197 + t262 * t198;
t233 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t261 + Ifges(6,6) * t265) * t253;
t234 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t261 + Ifges(6,2) * t265) * t253;
t235 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t261 + Ifges(6,4) * t265) * t253;
t270 = -mrSges(4,2) * t224 - mrSges(5,2) * t217 + pkin(3) * t276 + pkin(8) * t273 + t261 * (mrSges(6,2) * t213 - mrSges(6,3) * t211 + Ifges(6,1) * t239 + Ifges(6,4) * t240 + Ifges(6,5) * qJDD(5) - qJD(5) * t234 + t233 * t278) + t265 * (-mrSges(6,1) * t213 + mrSges(6,3) * t212 + Ifges(6,4) * t239 + Ifges(6,2) * t240 + Ifges(6,6) * qJDD(5) + qJD(5) * t235 - t233 * t279) + pkin(4) * t271 + mrSges(5,1) * t216 + mrSges(4,1) * t223 + (Ifges(5,3) + Ifges(4,3)) * t252;
t269 = mrSges(3,1) * t231 - mrSges(3,2) * t232 + Ifges(3,3) * t256 + pkin(2) * t277 + t270;
t1 = [t269 + Ifges(2,3) * qJDD(1) + pkin(1) * (t263 * (m(3) * t232 - t255 * mrSges(3,1) - t256 * mrSges(3,2) - t262 * t197 + t266 * t198) + t267 * (m(3) * t231 + t256 * mrSges(3,1) - t255 * mrSges(3,2) + t277)) + mrSges(2,1) * t274 - mrSges(2,2) * t272; t269; t270; m(5) * t258 + t265 * t209 + t261 * t210; mrSges(6,1) * t211 - mrSges(6,2) * t212 + Ifges(6,5) * t239 + Ifges(6,6) * t240 + Ifges(6,3) * qJDD(5) + (t234 * t261 - t235 * t265) * t253;];
tauJ = t1;
