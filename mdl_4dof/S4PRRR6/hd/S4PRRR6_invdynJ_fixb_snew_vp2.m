% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRR6
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:44
% EndTime: 2019-12-31 16:34:45
% DurationCPUTime: 0.33s
% Computational Cost: add. (1319->138), mult. (2594->186), div. (0->0), fcn. (1567->8), ass. (0->60)
t237 = sin(pkin(7));
t238 = cos(pkin(7));
t227 = -t238 * g(1) - t237 * g(2);
t236 = -g(3) + qJDD(1);
t241 = sin(qJ(2));
t244 = cos(qJ(2));
t211 = t244 * t227 + t241 * t236;
t245 = qJD(2) ^ 2;
t207 = -t245 * pkin(2) + qJDD(2) * pkin(5) + t211;
t226 = -t237 * g(1) + t238 * g(2);
t240 = sin(qJ(3));
t243 = cos(qJ(3));
t198 = -t240 * t207 + t243 * t226;
t253 = qJD(2) * qJD(3);
t252 = t243 * t253;
t224 = t240 * qJDD(2) + t252;
t189 = (-t224 + t252) * pkin(6) + (t240 * t243 * t245 + qJDD(3)) * pkin(3) + t198;
t199 = t243 * t207 + t240 * t226;
t225 = t243 * qJDD(2) - t240 * t253;
t254 = t240 * qJD(2);
t230 = qJD(3) * pkin(3) - pkin(6) * t254;
t235 = t243 ^ 2;
t190 = -t235 * t245 * pkin(3) + t225 * pkin(6) - qJD(3) * t230 + t199;
t239 = sin(qJ(4));
t242 = cos(qJ(4));
t187 = t242 * t189 - t239 * t190;
t215 = (-t240 * t239 + t243 * t242) * qJD(2);
t197 = t215 * qJD(4) + t242 * t224 + t239 * t225;
t216 = (t243 * t239 + t240 * t242) * qJD(2);
t204 = -t215 * mrSges(5,1) + t216 * mrSges(5,2);
t234 = qJD(3) + qJD(4);
t208 = -t234 * mrSges(5,2) + t215 * mrSges(5,3);
t233 = qJDD(3) + qJDD(4);
t184 = m(5) * t187 + t233 * mrSges(5,1) - t197 * mrSges(5,3) - t216 * t204 + t234 * t208;
t188 = t239 * t189 + t242 * t190;
t196 = -t216 * qJD(4) - t239 * t224 + t242 * t225;
t209 = t234 * mrSges(5,1) - t216 * mrSges(5,3);
t185 = m(5) * t188 - t233 * mrSges(5,2) + t196 * mrSges(5,3) + t215 * t204 - t234 * t209;
t178 = t242 * t184 + t239 * t185;
t255 = qJD(2) * t243;
t223 = (-t243 * mrSges(4,1) + t240 * mrSges(4,2)) * qJD(2);
t228 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t254;
t229 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t255;
t250 = -t239 * t184 + t242 * t185;
t251 = -t240 * (m(4) * t198 + qJDD(3) * mrSges(4,1) - t224 * mrSges(4,3) + qJD(3) * t229 - t223 * t254 + t178) + t243 * (m(4) * t199 - qJDD(3) * mrSges(4,2) + t225 * mrSges(4,3) - qJD(3) * t228 + t223 * t255 + t250);
t210 = -t241 * t227 + t244 * t236;
t249 = -qJDD(2) * pkin(2) - t210;
t191 = t230 * t254 - t225 * pkin(3) + (-pkin(6) * t235 - pkin(5)) * t245 + t249;
t248 = m(5) * t191 - t196 * mrSges(5,1) + t197 * mrSges(5,2) - t215 * t208 + t216 * t209;
t201 = Ifges(5,4) * t216 + Ifges(5,2) * t215 + Ifges(5,6) * t234;
t202 = Ifges(5,1) * t216 + Ifges(5,4) * t215 + Ifges(5,5) * t234;
t247 = mrSges(5,1) * t187 - mrSges(5,2) * t188 + Ifges(5,5) * t197 + Ifges(5,6) * t196 + Ifges(5,3) * t233 + t216 * t201 - t215 * t202;
t206 = -t245 * pkin(5) + t249;
t246 = -m(4) * t206 + t225 * mrSges(4,1) - t224 * mrSges(4,2) - t228 * t254 + t229 * t255 - t248;
t214 = Ifges(4,5) * qJD(3) + (t240 * Ifges(4,1) + t243 * Ifges(4,4)) * qJD(2);
t213 = Ifges(4,6) * qJD(3) + (t240 * Ifges(4,4) + t243 * Ifges(4,2)) * qJD(2);
t200 = Ifges(5,5) * t216 + Ifges(5,6) * t215 + Ifges(5,3) * t234;
t180 = mrSges(5,2) * t191 - mrSges(5,3) * t187 + Ifges(5,1) * t197 + Ifges(5,4) * t196 + Ifges(5,5) * t233 + t215 * t200 - t234 * t201;
t179 = -mrSges(5,1) * t191 + mrSges(5,3) * t188 + Ifges(5,4) * t197 + Ifges(5,2) * t196 + Ifges(5,6) * t233 - t216 * t200 + t234 * t202;
t1 = [m(2) * t236 + t241 * (m(3) * t211 - t245 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t251) + t244 * (m(3) * t210 + qJDD(2) * mrSges(3,1) - t245 * mrSges(3,2) + t246); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t210 - mrSges(3,2) * t211 + t240 * (mrSges(4,2) * t206 - mrSges(4,3) * t198 + Ifges(4,1) * t224 + Ifges(4,4) * t225 + Ifges(4,5) * qJDD(3) - pkin(6) * t178 - qJD(3) * t213 - t239 * t179 + t242 * t180) + t243 * (-mrSges(4,1) * t206 + mrSges(4,3) * t199 + Ifges(4,4) * t224 + Ifges(4,2) * t225 + Ifges(4,6) * qJDD(3) - pkin(3) * t248 + pkin(6) * t250 + qJD(3) * t214 + t242 * t179 + t239 * t180) + pkin(2) * t246 + pkin(5) * t251; mrSges(4,1) * t198 - mrSges(4,2) * t199 + Ifges(4,5) * t224 + Ifges(4,6) * t225 + Ifges(4,3) * qJDD(3) + pkin(3) * t178 + (t240 * t213 - t243 * t214) * qJD(2) + t247; t247;];
tauJ = t1;
