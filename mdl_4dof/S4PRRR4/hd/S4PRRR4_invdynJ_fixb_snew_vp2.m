% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:34
% EndTime: 2019-12-31 16:32:34
% DurationCPUTime: 0.33s
% Computational Cost: add. (1262->133), mult. (2502->181), div. (0->0), fcn. (1522->8), ass. (0->60)
t245 = qJD(2) ^ 2;
t237 = sin(pkin(7));
t238 = cos(pkin(7));
t225 = t237 * g(1) - t238 * g(2);
t226 = -t238 * g(1) - t237 * g(2);
t241 = sin(qJ(2));
t244 = cos(qJ(2));
t255 = t241 * t225 + t244 * t226;
t211 = -t245 * pkin(2) + qJDD(2) * pkin(5) + t255;
t236 = -g(3) + qJDD(1);
t240 = sin(qJ(3));
t243 = cos(qJ(3));
t203 = -t240 * t211 + t243 * t236;
t252 = qJD(2) * qJD(3);
t251 = t243 * t252;
t223 = t240 * qJDD(2) + t251;
t195 = (-t223 + t251) * pkin(6) + (t240 * t243 * t245 + qJDD(3)) * pkin(3) + t203;
t204 = t243 * t211 + t240 * t236;
t224 = t243 * qJDD(2) - t240 * t252;
t253 = t240 * qJD(2);
t229 = qJD(3) * pkin(3) - pkin(6) * t253;
t235 = t243 ^ 2;
t196 = -t235 * t245 * pkin(3) + t224 * pkin(6) - qJD(3) * t229 + t204;
t239 = sin(qJ(4));
t242 = cos(qJ(4));
t193 = t242 * t195 - t239 * t196;
t217 = (-t240 * t239 + t243 * t242) * qJD(2);
t202 = t217 * qJD(4) + t242 * t223 + t239 * t224;
t218 = (t243 * t239 + t240 * t242) * qJD(2);
t209 = -t217 * mrSges(5,1) + t218 * mrSges(5,2);
t234 = qJD(3) + qJD(4);
t212 = -t234 * mrSges(5,2) + t217 * mrSges(5,3);
t233 = qJDD(3) + qJDD(4);
t190 = m(5) * t193 + t233 * mrSges(5,1) - t202 * mrSges(5,3) - t218 * t209 + t234 * t212;
t194 = t239 * t195 + t242 * t196;
t201 = -t218 * qJD(4) - t239 * t223 + t242 * t224;
t213 = t234 * mrSges(5,1) - t218 * mrSges(5,3);
t191 = m(5) * t194 - t233 * mrSges(5,2) + t201 * mrSges(5,3) + t217 * t209 - t234 * t213;
t184 = t242 * t190 + t239 * t191;
t254 = qJD(2) * t243;
t250 = -t239 * t190 + t242 * t191;
t249 = t244 * t225 - t241 * t226;
t248 = -qJDD(2) * pkin(2) - t249;
t206 = Ifges(5,4) * t218 + Ifges(5,2) * t217 + Ifges(5,6) * t234;
t207 = Ifges(5,1) * t218 + Ifges(5,4) * t217 + Ifges(5,5) * t234;
t247 = mrSges(5,1) * t193 - mrSges(5,2) * t194 + Ifges(5,5) * t202 + Ifges(5,6) * t201 + Ifges(5,3) * t233 + t218 * t206 - t217 * t207;
t197 = t229 * t253 - t224 * pkin(3) + (-pkin(6) * t235 - pkin(5)) * t245 + t248;
t246 = m(5) * t197 - t201 * mrSges(5,1) + t202 * mrSges(5,2) - t217 * t212 + t218 * t213;
t228 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t254;
t227 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t253;
t222 = (-t243 * mrSges(4,1) + t240 * mrSges(4,2)) * qJD(2);
t216 = Ifges(4,5) * qJD(3) + (t240 * Ifges(4,1) + t243 * Ifges(4,4)) * qJD(2);
t215 = Ifges(4,6) * qJD(3) + (t240 * Ifges(4,4) + t243 * Ifges(4,2)) * qJD(2);
t210 = -t245 * pkin(5) + t248;
t205 = Ifges(5,5) * t218 + Ifges(5,6) * t217 + Ifges(5,3) * t234;
t186 = mrSges(5,2) * t197 - mrSges(5,3) * t193 + Ifges(5,1) * t202 + Ifges(5,4) * t201 + Ifges(5,5) * t233 + t217 * t205 - t234 * t206;
t185 = -mrSges(5,1) * t197 + mrSges(5,3) * t194 + Ifges(5,4) * t202 + Ifges(5,2) * t201 + Ifges(5,6) * t233 - t218 * t205 + t234 * t207;
t183 = m(4) * t204 - qJDD(3) * mrSges(4,2) + t224 * mrSges(4,3) - qJD(3) * t227 + t222 * t254 + t250;
t182 = m(4) * t203 + qJDD(3) * mrSges(4,1) - t223 * mrSges(4,3) + qJD(3) * t228 - t222 * t253 + t184;
t1 = [t243 * t182 + t240 * t183 + (m(2) + m(3)) * t236; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t249 - mrSges(3,2) * t255 + t240 * (mrSges(4,2) * t210 - mrSges(4,3) * t203 + Ifges(4,1) * t223 + Ifges(4,4) * t224 + Ifges(4,5) * qJDD(3) - pkin(6) * t184 - qJD(3) * t215 - t239 * t185 + t242 * t186) + t243 * (-mrSges(4,1) * t210 + mrSges(4,3) * t204 + Ifges(4,4) * t223 + Ifges(4,2) * t224 + Ifges(4,6) * qJDD(3) - pkin(3) * t246 + pkin(6) * t250 + qJD(3) * t216 + t242 * t185 + t239 * t186) + pkin(5) * (-t240 * t182 + t243 * t183) + (-m(4) * t210 + t224 * mrSges(4,1) - t223 * mrSges(4,2) - t246 + (-t227 * t240 + t228 * t243) * qJD(2)) * pkin(2); mrSges(4,1) * t203 - mrSges(4,2) * t204 + Ifges(4,5) * t223 + Ifges(4,6) * t224 + Ifges(4,3) * qJDD(3) + pkin(3) * t184 + (t240 * t215 - t243 * t216) * qJD(2) + t247; t247;];
tauJ = t1;
