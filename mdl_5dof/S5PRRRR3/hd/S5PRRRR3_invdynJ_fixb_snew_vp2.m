% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:13
% EndTime: 2019-12-05 17:06:13
% DurationCPUTime: 0.28s
% Computational Cost: add. (2312->99), mult. (2923->132), div. (0->0), fcn. (1922->10), ass. (0->54)
t246 = qJD(2) + qJD(3);
t243 = qJD(4) + t246;
t250 = sin(qJ(5));
t267 = t243 * t250;
t254 = cos(qJ(5));
t266 = t243 * t254;
t248 = sin(pkin(9));
t249 = cos(pkin(9));
t237 = t248 * g(1) - t249 * g(2);
t238 = -t249 * g(1) - t248 * g(2);
t253 = sin(qJ(2));
t257 = cos(qJ(2));
t261 = t257 * t237 - t253 * t238;
t222 = qJDD(2) * pkin(2) + t261;
t264 = t253 * t237 + t257 * t238;
t223 = -qJD(2) ^ 2 * pkin(2) + t264;
t252 = sin(qJ(3));
t256 = cos(qJ(3));
t217 = t256 * t222 - t252 * t223;
t245 = qJDD(2) + qJDD(3);
t214 = t245 * pkin(3) + t217;
t218 = t252 * t222 + t256 * t223;
t244 = t246 ^ 2;
t215 = -t244 * pkin(3) + t218;
t251 = sin(qJ(4));
t255 = cos(qJ(4));
t211 = t251 * t214 + t255 * t215;
t241 = t243 ^ 2;
t242 = qJDD(4) + t245;
t208 = -t241 * pkin(4) + t242 * pkin(8) + t211;
t247 = -g(3) + qJDD(1);
t205 = -t250 * t208 + t254 * t247;
t229 = (-mrSges(6,1) * t254 + mrSges(6,2) * t250) * t243;
t263 = qJD(5) * t243;
t230 = t250 * t242 + t254 * t263;
t233 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t266;
t203 = m(6) * t205 + qJDD(5) * mrSges(6,1) - t230 * mrSges(6,3) + qJD(5) * t233 - t229 * t267;
t206 = t254 * t208 + t250 * t247;
t231 = t254 * t242 - t250 * t263;
t232 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t267;
t204 = m(6) * t206 - qJDD(5) * mrSges(6,2) + t231 * mrSges(6,3) - qJD(5) * t232 + t229 * t266;
t262 = -t250 * t203 + t254 * t204;
t195 = m(5) * t211 - t241 * mrSges(5,1) - t242 * mrSges(5,2) + t262;
t210 = t255 * t214 - t251 * t215;
t207 = -t242 * pkin(4) - t241 * pkin(8) - t210;
t259 = -m(6) * t207 + t231 * mrSges(6,1) - t230 * mrSges(6,2) - t232 * t267 + t233 * t266;
t200 = m(5) * t210 + t242 * mrSges(5,1) - t241 * mrSges(5,2) + t259;
t265 = t251 * t195 + t255 * t200;
t224 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t250 + Ifges(6,6) * t254) * t243;
t225 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t250 + Ifges(6,2) * t254) * t243;
t226 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t250 + Ifges(6,4) * t254) * t243;
t260 = -mrSges(5,2) * t211 + pkin(8) * t262 + t250 * (mrSges(6,2) * t207 - mrSges(6,3) * t205 + Ifges(6,1) * t230 + Ifges(6,4) * t231 + Ifges(6,5) * qJDD(5) - qJD(5) * t225 + t224 * t266) + t254 * (-mrSges(6,1) * t207 + mrSges(6,3) * t206 + Ifges(6,4) * t230 + Ifges(6,2) * t231 + Ifges(6,6) * qJDD(5) + qJD(5) * t226 - t224 * t267) + pkin(4) * t259 + mrSges(5,1) * t210 + Ifges(5,3) * t242;
t258 = mrSges(4,1) * t217 - mrSges(4,2) * t218 + Ifges(4,3) * t245 + pkin(3) * t265 + t260;
t1 = [t254 * t203 + t250 * t204 + (m(2) + m(3) + m(4) + m(5)) * t247; t258 + Ifges(3,3) * qJDD(2) - mrSges(3,2) * t264 + mrSges(3,1) * t261 + pkin(2) * (t252 * (m(4) * t218 - t244 * mrSges(4,1) - t245 * mrSges(4,2) + t255 * t195 - t251 * t200) + t256 * (m(4) * t217 + t245 * mrSges(4,1) - t244 * mrSges(4,2) + t265)); t258; t260; mrSges(6,1) * t205 - mrSges(6,2) * t206 + Ifges(6,5) * t230 + Ifges(6,6) * t231 + Ifges(6,3) * qJDD(5) + (t225 * t250 - t226 * t254) * t243;];
tauJ = t1;
