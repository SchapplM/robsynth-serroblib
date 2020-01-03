% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:10
% EndTime: 2019-12-31 16:36:11
% DurationCPUTime: 0.48s
% Computational Cost: add. (1652->143), mult. (3103->191), div. (0->0), fcn. (2011->10), ass. (0->67)
t245 = sin(pkin(8));
t247 = cos(pkin(8));
t238 = t245 * g(1) - t247 * g(2);
t244 = -g(3) + qJDD(1);
t246 = sin(pkin(4));
t248 = cos(pkin(4));
t271 = t238 * t248 + t244 * t246;
t239 = -t247 * g(1) - t245 * g(2);
t251 = sin(qJ(2));
t254 = cos(qJ(2));
t209 = -t251 * t239 + t271 * t254;
t222 = -t246 * t238 + t248 * t244;
t253 = cos(qJ(3));
t268 = t253 * t222;
t210 = t254 * t239 + t271 * t251;
t256 = qJD(2) ^ 2;
t208 = -t256 * pkin(2) + qJDD(2) * pkin(6) + t210;
t250 = sin(qJ(3));
t205 = t253 * t208 + t250 * t222;
t267 = qJD(2) * t250;
t266 = t253 * qJD(2);
t265 = qJD(2) * qJD(3);
t264 = t250 * t265;
t263 = t253 * t265;
t234 = (-t253 * mrSges(4,1) + t250 * mrSges(4,2)) * qJD(2);
t237 = t253 * qJDD(2) - t264;
t240 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t267;
t235 = (-t253 * pkin(3) - t250 * pkin(7)) * qJD(2);
t255 = qJD(3) ^ 2;
t202 = -t255 * pkin(3) + qJDD(3) * pkin(7) + t235 * t266 + t205;
t207 = -qJDD(2) * pkin(2) - t256 * pkin(6) - t209;
t236 = t250 * qJDD(2) + t263;
t203 = (-t236 - t263) * pkin(7) + (-t237 + t264) * pkin(3) + t207;
t249 = sin(qJ(4));
t252 = cos(qJ(4));
t199 = -t249 * t202 + t252 * t203;
t232 = t252 * qJD(3) - t249 * t267;
t217 = t232 * qJD(4) + t249 * qJDD(3) + t252 * t236;
t233 = t249 * qJD(3) + t252 * t267;
t218 = -t232 * mrSges(5,1) + t233 * mrSges(5,2);
t243 = qJD(4) - t266;
t220 = -t243 * mrSges(5,2) + t232 * mrSges(5,3);
t229 = qJDD(4) - t237;
t197 = m(5) * t199 + t229 * mrSges(5,1) - t217 * mrSges(5,3) - t233 * t218 + t243 * t220;
t200 = t252 * t202 + t249 * t203;
t216 = -t233 * qJD(4) + t252 * qJDD(3) - t249 * t236;
t221 = t243 * mrSges(5,1) - t233 * mrSges(5,3);
t198 = m(5) * t200 - t229 * mrSges(5,2) + t216 * mrSges(5,3) + t232 * t218 - t243 * t221;
t261 = -t249 * t197 + t252 * t198;
t191 = m(4) * t205 - qJDD(3) * mrSges(4,2) + t237 * mrSges(4,3) - qJD(3) * t240 + t234 * t266 + t261;
t204 = -t250 * t208 + t268;
t241 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t266;
t201 = -qJDD(3) * pkin(3) - t255 * pkin(7) - t268 + (qJD(2) * t235 + t208) * t250;
t259 = -m(5) * t201 + t216 * mrSges(5,1) - t217 * mrSges(5,2) + t232 * t220 - t233 * t221;
t195 = m(4) * t204 + qJDD(3) * mrSges(4,1) - t236 * mrSges(4,3) + qJD(3) * t241 - t234 * t267 + t259;
t262 = t253 * t191 - t250 * t195;
t192 = t252 * t197 + t249 * t198;
t258 = -m(4) * t207 + t237 * mrSges(4,1) - t236 * mrSges(4,2) - t240 * t267 + t241 * t266 - t192;
t212 = Ifges(5,4) * t233 + Ifges(5,2) * t232 + Ifges(5,6) * t243;
t213 = Ifges(5,1) * t233 + Ifges(5,4) * t232 + Ifges(5,5) * t243;
t257 = mrSges(5,1) * t199 - mrSges(5,2) * t200 + Ifges(5,5) * t217 + Ifges(5,6) * t216 + Ifges(5,3) * t229 + t233 * t212 - t232 * t213;
t226 = Ifges(4,5) * qJD(3) + (t250 * Ifges(4,1) + t253 * Ifges(4,4)) * qJD(2);
t225 = Ifges(4,6) * qJD(3) + (t250 * Ifges(4,4) + t253 * Ifges(4,2)) * qJD(2);
t211 = Ifges(5,5) * t233 + Ifges(5,6) * t232 + Ifges(5,3) * t243;
t194 = mrSges(5,2) * t201 - mrSges(5,3) * t199 + Ifges(5,1) * t217 + Ifges(5,4) * t216 + Ifges(5,5) * t229 + t232 * t211 - t243 * t212;
t193 = -mrSges(5,1) * t201 + mrSges(5,3) * t200 + Ifges(5,4) * t217 + Ifges(5,2) * t216 + Ifges(5,6) * t229 - t233 * t211 + t243 * t213;
t1 = [m(2) * t244 + t248 * (m(3) * t222 + t250 * t191 + t253 * t195) + (t251 * (m(3) * t210 - t256 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t262) + t254 * (m(3) * t209 + qJDD(2) * mrSges(3,1) - t256 * mrSges(3,2) + t258)) * t246; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t209 - mrSges(3,2) * t210 + t250 * (mrSges(4,2) * t207 - mrSges(4,3) * t204 + Ifges(4,1) * t236 + Ifges(4,4) * t237 + Ifges(4,5) * qJDD(3) - pkin(7) * t192 - qJD(3) * t225 - t249 * t193 + t252 * t194) + t253 * (-mrSges(4,1) * t207 + mrSges(4,3) * t205 + Ifges(4,4) * t236 + Ifges(4,2) * t237 + Ifges(4,6) * qJDD(3) - pkin(3) * t192 + qJD(3) * t226 - t257) + pkin(2) * t258 + pkin(6) * t262; Ifges(4,5) * t236 + Ifges(4,6) * t237 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t204 - mrSges(4,2) * t205 + t249 * t194 + t252 * t193 + pkin(3) * t259 + pkin(7) * t261 + (t250 * t225 - t253 * t226) * qJD(2); t257;];
tauJ = t1;
