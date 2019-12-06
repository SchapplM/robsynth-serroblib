% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRRR1
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:44
% EndTime: 2019-12-05 15:12:44
% DurationCPUTime: 0.22s
% Computational Cost: add. (1580->96), mult. (2058->132), div. (0->0), fcn. (1430->10), ass. (0->51)
t225 = qJD(3) + qJD(4);
t231 = sin(qJ(5));
t244 = t225 * t231;
t234 = cos(qJ(5));
t243 = t225 * t234;
t228 = sin(pkin(8));
t230 = cos(pkin(8));
t221 = -t230 * g(1) - t228 * g(2);
t226 = -g(3) + qJDD(1);
t227 = sin(pkin(9));
t229 = cos(pkin(9));
t208 = -t227 * t221 + t229 * t226;
t209 = t229 * t221 + t227 * t226;
t233 = sin(qJ(3));
t236 = cos(qJ(3));
t203 = t236 * t208 - t233 * t209;
t201 = qJDD(3) * pkin(3) + t203;
t204 = t233 * t208 + t236 * t209;
t237 = qJD(3) ^ 2;
t202 = -t237 * pkin(3) + t204;
t232 = sin(qJ(4));
t235 = cos(qJ(4));
t198 = t232 * t201 + t235 * t202;
t223 = t225 ^ 2;
t224 = qJDD(3) + qJDD(4);
t195 = -t223 * pkin(4) + t224 * pkin(7) + t198;
t220 = -t228 * g(1) + t230 * g(2) + qJDD(2);
t192 = -t231 * t195 + t234 * t220;
t215 = (-mrSges(6,1) * t234 + mrSges(6,2) * t231) * t225;
t241 = qJD(5) * t225;
t216 = t231 * t224 + t234 * t241;
t219 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t243;
t189 = m(6) * t192 + qJDD(5) * mrSges(6,1) - t216 * mrSges(6,3) + qJD(5) * t219 - t215 * t244;
t193 = t234 * t195 + t231 * t220;
t217 = t234 * t224 - t231 * t241;
t218 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t244;
t190 = m(6) * t193 - qJDD(5) * mrSges(6,2) + t217 * mrSges(6,3) - qJD(5) * t218 + t215 * t243;
t240 = -t231 * t189 + t234 * t190;
t182 = m(5) * t198 - t223 * mrSges(5,1) - t224 * mrSges(5,2) + t240;
t197 = t235 * t201 - t232 * t202;
t194 = -t224 * pkin(4) - t223 * pkin(7) - t197;
t238 = -m(6) * t194 + t217 * mrSges(6,1) - t216 * mrSges(6,2) - t218 * t244 + t219 * t243;
t187 = m(5) * t197 + t224 * mrSges(5,1) - t223 * mrSges(5,2) + t238;
t242 = t232 * t182 + t235 * t187;
t210 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t231 + Ifges(6,6) * t234) * t225;
t211 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t231 + Ifges(6,2) * t234) * t225;
t212 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t231 + Ifges(6,4) * t234) * t225;
t239 = -mrSges(5,2) * t198 + pkin(7) * t240 + t231 * (mrSges(6,2) * t194 - mrSges(6,3) * t192 + Ifges(6,1) * t216 + Ifges(6,4) * t217 + Ifges(6,5) * qJDD(5) - qJD(5) * t211 + t210 * t243) + t234 * (-mrSges(6,1) * t194 + mrSges(6,3) * t193 + Ifges(6,4) * t216 + Ifges(6,2) * t217 + Ifges(6,6) * qJDD(5) + qJD(5) * t212 - t210 * t244) + pkin(4) * t238 + mrSges(5,1) * t197 + Ifges(5,3) * t224;
t180 = m(4) * t204 - t237 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t235 * t182 - t232 * t187;
t179 = m(4) * t203 + qJDD(3) * mrSges(4,1) - t237 * mrSges(4,2) + t242;
t1 = [m(2) * t226 + t227 * (m(3) * t209 - t233 * t179 + t236 * t180) + t229 * (m(3) * t208 + t236 * t179 + t233 * t180); t234 * t189 + t231 * t190 + (m(3) + m(4) + m(5)) * t220; mrSges(4,1) * t203 - mrSges(4,2) * t204 + Ifges(4,3) * qJDD(3) + pkin(3) * t242 + t239; t239; mrSges(6,1) * t192 - mrSges(6,2) * t193 + Ifges(6,5) * t216 + Ifges(6,6) * t217 + Ifges(6,3) * qJDD(5) + (t211 * t231 - t212 * t234) * t225;];
tauJ = t1;
