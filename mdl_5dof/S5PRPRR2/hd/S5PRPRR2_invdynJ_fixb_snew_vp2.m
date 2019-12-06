% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:56
% EndTime: 2019-12-05 15:44:56
% DurationCPUTime: 0.31s
% Computational Cost: add. (2170->105), mult. (2879->141), div. (0->0), fcn. (1742->10), ass. (0->54)
t243 = qJD(2) + qJD(4);
t249 = sin(qJ(5));
t263 = t243 * t249;
t252 = cos(qJ(5));
t262 = t243 * t252;
t246 = sin(pkin(8));
t248 = cos(pkin(8));
t237 = -t248 * g(1) - t246 * g(2);
t244 = -g(3) + qJDD(1);
t251 = sin(qJ(2));
t254 = cos(qJ(2));
t223 = -t251 * t237 + t254 * t244;
t221 = qJDD(2) * pkin(2) + t223;
t224 = t254 * t237 + t251 * t244;
t255 = qJD(2) ^ 2;
t222 = -t255 * pkin(2) + t224;
t245 = sin(pkin(9));
t247 = cos(pkin(9));
t216 = t247 * t221 - t245 * t222;
t214 = qJDD(2) * pkin(3) + t216;
t217 = t245 * t221 + t247 * t222;
t215 = -t255 * pkin(3) + t217;
t250 = sin(qJ(4));
t253 = cos(qJ(4));
t211 = t250 * t214 + t253 * t215;
t241 = t243 ^ 2;
t242 = qJDD(2) + qJDD(4);
t208 = -t241 * pkin(4) + t242 * pkin(7) + t211;
t236 = -t246 * g(1) + t248 * g(2) + qJDD(3);
t205 = -t249 * t208 + t252 * t236;
t230 = (-mrSges(6,1) * t252 + mrSges(6,2) * t249) * t243;
t259 = qJD(5) * t243;
t231 = t249 * t242 + t252 * t259;
t235 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t262;
t202 = m(6) * t205 + qJDD(5) * mrSges(6,1) - t231 * mrSges(6,3) + qJD(5) * t235 - t230 * t263;
t206 = t252 * t208 + t249 * t236;
t232 = t252 * t242 - t249 * t259;
t234 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t263;
t203 = m(6) * t206 - qJDD(5) * mrSges(6,2) + t232 * mrSges(6,3) - qJD(5) * t234 + t230 * t262;
t258 = -t249 * t202 + t252 * t203;
t195 = m(5) * t211 - t241 * mrSges(5,1) - t242 * mrSges(5,2) + t258;
t210 = t253 * t214 - t250 * t215;
t207 = -t242 * pkin(4) - t241 * pkin(7) - t210;
t256 = -m(6) * t207 + t232 * mrSges(6,1) - t231 * mrSges(6,2) - t234 * t263 + t235 * t262;
t200 = m(5) * t210 + t242 * mrSges(5,1) - t241 * mrSges(5,2) + t256;
t260 = t250 * t195 + t253 * t200;
t192 = m(4) * t216 + qJDD(2) * mrSges(4,1) - t255 * mrSges(4,2) + t260;
t193 = m(4) * t217 - t255 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t253 * t195 - t250 * t200;
t261 = t247 * t192 + t245 * t193;
t225 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t249 + Ifges(6,6) * t252) * t243;
t226 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t249 + Ifges(6,2) * t252) * t243;
t227 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t249 + Ifges(6,4) * t252) * t243;
t257 = -mrSges(5,2) * t211 + pkin(7) * t258 + t249 * (mrSges(6,2) * t207 - mrSges(6,3) * t205 + Ifges(6,1) * t231 + Ifges(6,4) * t232 + Ifges(6,5) * qJDD(5) - qJD(5) * t226 + t225 * t262) + t252 * (-mrSges(6,1) * t207 + mrSges(6,3) * t206 + Ifges(6,4) * t231 + Ifges(6,2) * t232 + Ifges(6,6) * qJDD(5) + qJD(5) * t227 - t225 * t263) + pkin(4) * t256 + mrSges(5,1) * t210 + Ifges(5,3) * t242;
t1 = [m(2) * t244 + t251 * (m(3) * t224 - t255 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t245 * t192 + t247 * t193) + t254 * (m(3) * t223 + qJDD(2) * mrSges(3,1) - t255 * mrSges(3,2) + t261); pkin(2) * t261 - mrSges(3,2) * t224 + mrSges(3,1) * t223 + pkin(3) * t260 + mrSges(4,1) * t216 - mrSges(4,2) * t217 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + t257; t252 * t202 + t249 * t203 + (m(4) + m(5)) * t236; t257; mrSges(6,1) * t205 - mrSges(6,2) * t206 + Ifges(6,5) * t231 + Ifges(6,6) * t232 + Ifges(6,3) * qJDD(5) + (t226 * t249 - t227 * t252) * t243;];
tauJ = t1;
