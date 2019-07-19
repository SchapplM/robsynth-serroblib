% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:33
% EndTime: 2019-07-18 18:16:33
% DurationCPUTime: 0.32s
% Computational Cost: add. (3558->118), mult. (4057->132), div. (0->0), fcn. (1656->6), ass. (0->54)
t251 = -m(4) - m(5);
t250 = -pkin(2) - pkin(3);
t249 = -mrSges(3,1) - mrSges(4,1);
t248 = Ifges(4,4) + Ifges(3,5);
t247 = Ifges(3,6) - Ifges(4,6);
t233 = sin(qJ(1));
t236 = cos(qJ(1));
t216 = t233 * g(1) - t236 * g(2);
t214 = qJDD(1) * pkin(1) + t216;
t217 = -t236 * g(1) - t233 * g(2);
t237 = qJD(1) ^ 2;
t215 = -t237 * pkin(1) + t217;
t232 = sin(qJ(2));
t235 = cos(qJ(2));
t210 = t232 * t214 + t235 * t215;
t230 = (qJD(1) + qJD(2));
t228 = t230 ^ 2;
t229 = qJDD(1) + qJDD(2);
t242 = t229 * qJ(3) + (2 * qJD(3) * t230) + t210;
t204 = t250 * t228 + t242;
t209 = t235 * t214 - t232 * t215;
t239 = -t228 * qJ(3) + qJDD(3) - t209;
t206 = t250 * t229 + t239;
t231 = sin(qJ(4));
t234 = cos(qJ(4));
t202 = -t231 * t204 + t234 * t206;
t226 = qJD(4) - t230;
t223 = t226 ^ 2;
t224 = qJDD(4) - t229;
t200 = m(5) * t202 + t224 * mrSges(5,1) - (t223 * mrSges(5,2));
t203 = t234 * t204 + t231 * t206;
t201 = m(5) * t203 - t223 * mrSges(5,1) - t224 * mrSges(5,2);
t207 = -t228 * pkin(2) + t242;
t240 = m(4) * t207 + t229 * mrSges(4,3) - t231 * t200 + t234 * t201;
t194 = m(3) * t210 - t229 * mrSges(3,2) + t249 * t228 + t240;
t208 = -t229 * pkin(2) + t239;
t241 = t234 * t200 + t231 * t201;
t238 = -m(4) * t208 + t229 * mrSges(4,1) + t228 * mrSges(4,3) - t241;
t196 = m(3) * t209 + t229 * mrSges(3,1) - t228 * mrSges(3,2) + t238;
t189 = t232 * t194 + t235 * t196;
t187 = m(2) * t216 + qJDD(1) * mrSges(2,1) - t237 * mrSges(2,2) + t189;
t243 = t235 * t194 - t232 * t196;
t188 = m(2) * t217 - t237 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t243;
t246 = t236 * t187 + t233 * t188;
t245 = -m(3) + t251;
t244 = -t233 * t187 + t236 * t188;
t220 = t251 * g(3);
t199 = mrSges(5,2) * g(3) - mrSges(5,3) * t202 + Ifges(5,5) * t224 - (t223 * Ifges(5,6));
t198 = -mrSges(5,1) * g(3) + mrSges(5,3) * t203 + t223 * Ifges(5,5) + Ifges(5,6) * t224;
t191 = mrSges(4,2) * t208 - mrSges(3,3) * t209 - qJ(3) * t220 - t231 * t198 + t234 * t199 + t248 * t229 - t247 * t228 + (-mrSges(3,2) + mrSges(4,3)) * g(3);
t190 = mrSges(4,2) * t207 + mrSges(3,3) * t210 - pkin(2) * t220 - t234 * t198 - t231 * t199 + t247 * t229 + t248 * t228 + (m(5) * pkin(3) - t249) * g(3);
t183 = -mrSges(2,2) * g(3) - mrSges(2,3) * t216 + Ifges(2,5) * qJDD(1) - t237 * Ifges(2,6) - pkin(5) * t189 - t232 * t190 + t235 * t191;
t182 = Ifges(2,6) * qJDD(1) + t237 * Ifges(2,5) + mrSges(2,3) * t217 + t232 * t191 + t235 * t190 + pkin(5) * t243 + (-pkin(1) * t245 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t244; -m(1) * g(2) + t246; (-m(1) - m(2) + t245) * g(3); -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t246 - t233 * t182 + t236 * t183; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t244 + t236 * t182 + t233 * t183; pkin(1) * t189 + mrSges(2,1) * t216 - mrSges(2,2) * t217 + mrSges(3,1) * t209 - mrSges(3,2) * t210 + pkin(2) * t238 + qJ(3) * (-t228 * mrSges(4,1) + t240) - mrSges(4,1) * t208 + mrSges(4,3) * t207 - pkin(3) * t241 - mrSges(5,1) * t202 + mrSges(5,2) * t203 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(5,3) * t224 + Ifges(2,3) * qJDD(1) + (Ifges(3,3) + Ifges(4,2)) * t229;];
tauB  = t1;
