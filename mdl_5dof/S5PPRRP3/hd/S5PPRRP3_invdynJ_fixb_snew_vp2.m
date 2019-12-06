% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:39
% EndTime: 2019-12-05 15:10:40
% DurationCPUTime: 0.50s
% Computational Cost: add. (909->124), mult. (1569->158), div. (0->0), fcn. (880->8), ass. (0->57)
t253 = sin(qJ(4));
t255 = cos(qJ(4));
t273 = Ifges(5,4) - Ifges(6,5);
t279 = t253 * (Ifges(5,1) + Ifges(6,1)) + t255 * t273;
t278 = -t253 * t273 - t255 * (Ifges(5,2) + Ifges(6,3));
t272 = Ifges(6,4) + Ifges(5,5);
t271 = Ifges(5,6) - Ifges(6,6);
t275 = t253 * (t278 * qJD(3) - t271 * qJD(4)) + t255 * (t279 * qJD(3) + t272 * qJD(4));
t274 = mrSges(5,3) + mrSges(6,2);
t264 = qJD(3) * qJD(4);
t236 = qJDD(3) * t255 - t253 * t264;
t270 = mrSges(6,1) * t236;
t250 = sin(pkin(7));
t252 = cos(pkin(7));
t240 = -g(1) * t252 - g(2) * t250;
t248 = -g(3) + qJDD(1);
t249 = sin(pkin(8));
t251 = cos(pkin(8));
t217 = t240 * t249 - t248 * t251;
t269 = t217 * t255;
t218 = t240 * t251 + t248 * t249;
t239 = -g(1) * t250 + g(2) * t252 + qJDD(2);
t254 = sin(qJ(3));
t256 = cos(qJ(3));
t213 = t256 * t218 + t254 * t239;
t258 = qJD(3) ^ 2;
t211 = -pkin(3) * t258 + qJDD(3) * pkin(6) + t213;
t208 = t255 * t211 + t253 * t217;
t266 = qJD(3) * t253;
t265 = qJD(3) * t255;
t234 = (-t255 * mrSges(5,1) + t253 * mrSges(5,2)) * qJD(3);
t241 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t266;
t232 = (-t255 * pkin(4) - t253 * qJ(5)) * qJD(3);
t257 = qJD(4) ^ 2;
t204 = -pkin(4) * t257 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t232 * t265 + t208;
t233 = (-t255 * mrSges(6,1) - t253 * mrSges(6,3)) * qJD(3);
t242 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t266;
t262 = m(6) * t204 + qJDD(4) * mrSges(6,3) + qJD(4) * t242 + t233 * t265;
t199 = m(5) * t208 - qJDD(4) * mrSges(5,2) - qJD(4) * t241 + t234 * t265 + t274 * t236 + t262;
t207 = -t211 * t253 + t269;
t235 = qJDD(3) * t253 + t255 * t264;
t243 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t265;
t205 = -qJDD(4) * pkin(4) - qJ(5) * t257 - t269 + qJDD(5) + (qJD(3) * t232 + t211) * t253;
t244 = mrSges(6,2) * t265 + qJD(4) * mrSges(6,3);
t261 = -m(6) * t205 + qJDD(4) * mrSges(6,1) + qJD(4) * t244;
t200 = m(5) * t207 + qJDD(4) * mrSges(5,1) + qJD(4) * t243 - t274 * t235 + (-t233 - t234) * t266 + t261;
t263 = t255 * t199 - t200 * t253;
t212 = -t254 * t218 + t239 * t256;
t210 = -qJDD(3) * pkin(3) - pkin(6) * t258 - t212;
t206 = -pkin(4) * t236 - qJ(5) * t235 + (-0.2e1 * qJD(5) * t253 + (pkin(4) * t253 - qJ(5) * t255) * qJD(4)) * qJD(3) + t210;
t260 = m(6) * t206 - t235 * mrSges(6,3) - t242 * t266 - t244 * t265;
t259 = -m(5) * t210 + t236 * mrSges(5,1) - t241 * t266 + t243 * t265 - t260;
t202 = mrSges(6,2) * t235 + t233 * t266 - t261;
t201 = t260 - t270;
t197 = m(4) * t212 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t258 - mrSges(5,2) * t235 + t259 + t270;
t196 = m(4) * t213 - mrSges(4,1) * t258 - qJDD(3) * mrSges(4,2) + t263;
t1 = [m(2) * t248 + t249 * (m(3) * t218 + t196 * t256 - t197 * t254) + t251 * (-t199 * t253 - t200 * t255 + (-m(3) - m(4)) * t217); m(3) * t239 + t196 * t254 + t197 * t256; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t212 - mrSges(4,2) * t213 + t253 * (mrSges(5,2) * t210 + mrSges(6,2) * t205 - mrSges(5,3) * t207 - mrSges(6,3) * t206 - qJ(5) * t201) + t255 * (-mrSges(5,1) * t210 - mrSges(6,1) * t206 + mrSges(6,2) * t204 + mrSges(5,3) * t208 - pkin(4) * t201) + pkin(3) * t259 + pkin(6) * t263 + (pkin(3) * mrSges(6,1) - t278) * t236 + (-pkin(3) * mrSges(5,2) + t279) * t235 + (t253 * t272 + t255 * t271) * qJDD(4) + t275 * qJD(4); mrSges(5,1) * t207 - mrSges(5,2) * t208 - mrSges(6,1) * t205 + mrSges(6,3) * t204 - pkin(4) * t202 + qJ(5) * t262 + (qJ(5) * mrSges(6,2) + t271) * t236 + t272 * t235 + (Ifges(5,3) + Ifges(6,2)) * qJDD(4) - t275 * qJD(3); t202;];
tauJ = t1;
