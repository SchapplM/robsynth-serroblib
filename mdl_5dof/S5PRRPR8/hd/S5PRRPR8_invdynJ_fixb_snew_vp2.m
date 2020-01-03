% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPR8
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:27
% EndTime: 2019-12-31 17:42:27
% DurationCPUTime: 0.34s
% Computational Cost: add. (2504->104), mult. (3141->141), div. (0->0), fcn. (1902->10), ass. (0->54)
t246 = qJD(2) + qJD(3);
t252 = sin(qJ(5));
t266 = t246 * t252;
t255 = cos(qJ(5));
t265 = t246 * t255;
t249 = sin(pkin(8));
t251 = cos(pkin(8));
t239 = -t251 * g(1) - t249 * g(2);
t247 = -g(3) + qJDD(1);
t254 = sin(qJ(2));
t257 = cos(qJ(2));
t225 = -t254 * t239 + t257 * t247;
t223 = qJDD(2) * pkin(2) + t225;
t226 = t257 * t239 + t254 * t247;
t258 = qJD(2) ^ 2;
t224 = -t258 * pkin(2) + t226;
t253 = sin(qJ(3));
t256 = cos(qJ(3));
t218 = t256 * t223 - t253 * t224;
t244 = t246 ^ 2;
t245 = qJDD(2) + qJDD(3);
t215 = t245 * pkin(3) + t218;
t219 = t253 * t223 + t256 * t224;
t216 = -t244 * pkin(3) + t219;
t248 = sin(pkin(9));
t250 = cos(pkin(9));
t212 = t248 * t215 + t250 * t216;
t209 = -t244 * pkin(4) + t245 * pkin(7) + t212;
t238 = -t249 * g(1) + t251 * g(2) + qJDD(4);
t206 = -t252 * t209 + t255 * t238;
t232 = (-mrSges(6,1) * t255 + mrSges(6,2) * t252) * t246;
t262 = qJD(5) * t246;
t233 = t252 * t245 + t255 * t262;
t237 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t265;
t203 = m(6) * t206 + qJDD(5) * mrSges(6,1) - t233 * mrSges(6,3) + qJD(5) * t237 - t232 * t266;
t207 = t255 * t209 + t252 * t238;
t234 = t255 * t245 - t252 * t262;
t236 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t266;
t204 = m(6) * t207 - qJDD(5) * mrSges(6,2) + t234 * mrSges(6,3) - qJD(5) * t236 + t232 * t265;
t261 = -t252 * t203 + t255 * t204;
t196 = m(5) * t212 - t244 * mrSges(5,1) - t245 * mrSges(5,2) + t261;
t211 = t250 * t215 - t248 * t216;
t208 = -t245 * pkin(4) - t244 * pkin(7) - t211;
t260 = -m(6) * t208 + t234 * mrSges(6,1) - t233 * mrSges(6,2) - t236 * t266 + t237 * t265;
t201 = m(5) * t211 + t245 * mrSges(5,1) - t244 * mrSges(5,2) + t260;
t263 = t248 * t196 + t250 * t201;
t192 = m(4) * t218 + t245 * mrSges(4,1) - t244 * mrSges(4,2) + t263;
t193 = m(4) * t219 - t244 * mrSges(4,1) - t245 * mrSges(4,2) + t250 * t196 - t248 * t201;
t264 = t256 * t192 + t253 * t193;
t227 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t252 + Ifges(6,6) * t255) * t246;
t228 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t252 + Ifges(6,2) * t255) * t246;
t229 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t252 + Ifges(6,4) * t255) * t246;
t259 = -mrSges(4,2) * t219 - mrSges(5,2) * t212 + pkin(3) * t263 + pkin(7) * t261 + t252 * (mrSges(6,2) * t208 - mrSges(6,3) * t206 + Ifges(6,1) * t233 + Ifges(6,4) * t234 + Ifges(6,5) * qJDD(5) - qJD(5) * t228 + t227 * t265) + t255 * (-mrSges(6,1) * t208 + mrSges(6,3) * t207 + Ifges(6,4) * t233 + Ifges(6,2) * t234 + Ifges(6,6) * qJDD(5) + qJD(5) * t229 - t227 * t266) + pkin(4) * t260 + mrSges(5,1) * t211 + mrSges(4,1) * t218 + (Ifges(5,3) + Ifges(4,3)) * t245;
t1 = [m(2) * t247 + t254 * (m(3) * t226 - t258 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t253 * t192 + t256 * t193) + t257 * (m(3) * t225 + qJDD(2) * mrSges(3,1) - t258 * mrSges(3,2) + t264); mrSges(3,1) * t225 - mrSges(3,2) * t226 + Ifges(3,3) * qJDD(2) + pkin(2) * t264 + t259; t259; m(5) * t238 + t255 * t203 + t252 * t204; mrSges(6,1) * t206 - mrSges(6,2) * t207 + Ifges(6,5) * t233 + Ifges(6,6) * t234 + Ifges(6,3) * qJDD(5) + (t228 * t252 - t229 * t255) * t246;];
tauJ = t1;
