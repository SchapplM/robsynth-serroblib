% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPPR5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:45
% EndTime: 2019-12-31 16:39:45
% DurationCPUTime: 0.42s
% Computational Cost: add. (3126->149), mult. (5170->180), div. (0->0), fcn. (1848->6), ass. (0->61)
t266 = -pkin(1) - pkin(2);
t265 = mrSges(2,1) + mrSges(3,1);
t264 = Ifges(3,4) + Ifges(2,5);
t263 = Ifges(2,6) - Ifges(3,6);
t246 = sin(qJ(1));
t248 = cos(qJ(1));
t235 = -t248 * g(1) - t246 * g(2);
t249 = qJD(1) ^ 2;
t254 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t235;
t222 = -t249 * pkin(1) + t254;
t219 = t266 * t249 + t254;
t234 = t246 * g(1) - t248 * g(2);
t253 = -t249 * qJ(2) + qJDD(2) - t234;
t221 = t266 * qJDD(1) + t253;
t243 = sin(pkin(6));
t244 = cos(pkin(6));
t216 = t244 * t219 + t243 * t221;
t214 = -t249 * pkin(3) - qJDD(1) * pkin(5) + t216;
t241 = g(3) + qJDD(3);
t245 = sin(qJ(4));
t247 = cos(qJ(4));
t211 = -t245 * t214 + t247 * t241;
t229 = (mrSges(5,1) * t247 - mrSges(5,2) * t245) * qJD(1);
t259 = qJD(1) * qJD(4);
t230 = -t245 * qJDD(1) - t247 * t259;
t260 = qJD(1) * t247;
t233 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t260;
t261 = qJD(1) * t245;
t209 = m(5) * t211 + qJDD(4) * mrSges(5,1) - t230 * mrSges(5,3) + qJD(4) * t233 + t229 * t261;
t212 = t247 * t214 + t245 * t241;
t231 = -t247 * qJDD(1) + t245 * t259;
t232 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t261;
t210 = m(5) * t212 - qJDD(4) * mrSges(5,2) + t231 * mrSges(5,3) - qJD(4) * t232 - t229 * t260;
t256 = -t245 * t209 + t247 * t210;
t202 = m(4) * t216 - t249 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t256;
t215 = -t243 * t219 + t244 * t221;
t213 = qJDD(1) * pkin(3) - t249 * pkin(5) - t215;
t250 = -m(5) * t213 + t231 * mrSges(5,1) - t230 * mrSges(5,2) + t232 * t261 - t233 * t260;
t207 = m(4) * t215 - qJDD(1) * mrSges(4,1) - t249 * mrSges(4,2) + t250;
t257 = t244 * t202 - t243 * t207;
t255 = m(3) * t222 + qJDD(1) * mrSges(3,3) + t257;
t197 = m(2) * t235 - qJDD(1) * mrSges(2,2) - t265 * t249 + t255;
t200 = t243 * t202 + t244 * t207;
t223 = -qJDD(1) * pkin(1) + t253;
t251 = -m(3) * t223 + qJDD(1) * mrSges(3,1) + t249 * mrSges(3,3) - t200;
t198 = m(2) * t234 + qJDD(1) * mrSges(2,1) - t249 * mrSges(2,2) + t251;
t262 = t246 * t197 + t248 * t198;
t258 = t248 * t197 - t246 * t198;
t204 = t247 * t209 + t245 * t210;
t252 = -m(4) * t241 - t204;
t226 = (Ifges(5,5) * qJD(4)) + (-Ifges(5,1) * t245 - Ifges(5,4) * t247) * qJD(1);
t225 = (Ifges(5,6) * qJD(4)) + (-Ifges(5,4) * t245 - Ifges(5,2) * t247) * qJD(1);
t224 = (Ifges(5,3) * qJD(4)) + (-Ifges(5,5) * t245 - Ifges(5,6) * t247) * qJD(1);
t206 = mrSges(5,2) * t213 - mrSges(5,3) * t211 + Ifges(5,1) * t230 + Ifges(5,4) * t231 + Ifges(5,5) * qJDD(4) - qJD(4) * t225 - t224 * t260;
t205 = -mrSges(5,1) * t213 + mrSges(5,3) * t212 + Ifges(5,4) * t230 + Ifges(5,2) * t231 + Ifges(5,6) * qJDD(4) + qJD(4) * t226 + t224 * t261;
t203 = -m(3) * g(3) + t252;
t199 = -mrSges(4,1) * t241 - mrSges(5,1) * t211 + mrSges(5,2) * t212 + mrSges(4,3) * t216 + t249 * Ifges(4,5) - Ifges(5,5) * t230 - Ifges(4,6) * qJDD(1) - Ifges(5,6) * t231 - Ifges(5,3) * qJDD(4) - pkin(3) * t204 + (t225 * t245 - t226 * t247) * qJD(1);
t193 = mrSges(4,2) * t241 - mrSges(4,3) * t215 - Ifges(4,5) * qJDD(1) - t249 * Ifges(4,6) - pkin(5) * t204 - t245 * t205 + t247 * t206;
t192 = mrSges(3,2) * t223 - mrSges(2,3) * t234 - qJ(2) * t203 - qJ(3) * t200 + t244 * t193 - t243 * t199 - t263 * t249 + t264 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t191 = mrSges(3,2) * t222 + mrSges(2,3) * t235 - pkin(1) * t203 - pkin(2) * t252 + t265 * g(3) - qJ(3) * t257 + t263 * qJDD(1) - t243 * t193 - t244 * t199 + t264 * t249;
t1 = [-m(1) * g(1) + t258; -m(1) * g(2) + t262; (-m(1) - m(2) - m(3)) * g(3) + t252; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t262 - t246 * t191 + t248 * t192; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t258 + t248 * t191 + t246 * t192; pkin(1) * t251 + qJ(2) * (-t249 * mrSges(3,1) + t255) + mrSges(2,1) * t234 - mrSges(2,2) * t235 - pkin(2) * t200 - mrSges(3,1) * t223 + mrSges(3,3) * t222 - mrSges(4,1) * t215 + mrSges(4,2) * t216 - t245 * t206 - t247 * t205 - pkin(3) * t250 - pkin(5) * t256 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB = t1;
