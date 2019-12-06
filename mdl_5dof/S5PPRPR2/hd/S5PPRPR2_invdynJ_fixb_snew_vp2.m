% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PPRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:04
% EndTime: 2019-12-05 15:03:05
% DurationCPUTime: 0.18s
% Computational Cost: add. (626->90), mult. (973->118), div. (0->0), fcn. (582->8), ass. (0->44)
t227 = -pkin(3) - pkin(6);
t211 = sin(pkin(7));
t213 = cos(pkin(7));
t203 = -t213 * g(1) - t211 * g(2);
t207 = -g(3) + qJDD(1);
t210 = sin(pkin(8));
t212 = cos(pkin(8));
t191 = -t210 * t203 + t212 * t207;
t192 = t212 * t203 + t210 * t207;
t215 = sin(qJ(3));
t217 = cos(qJ(3));
t187 = t215 * t191 + t217 * t192;
t214 = sin(qJ(5));
t226 = qJD(3) * t214;
t216 = cos(qJ(5));
t225 = qJD(3) * t216;
t224 = qJD(3) * qJD(5);
t186 = t217 * t191 - t215 * t192;
t218 = qJD(3) ^ 2;
t222 = -t218 * qJ(4) + qJDD(4) - t186;
t183 = t227 * qJDD(3) + t222;
t202 = -t211 * g(1) + t213 * g(2) + qJDD(2);
t179 = t216 * t183 - t214 * t202;
t199 = (t214 * mrSges(6,1) + t216 * mrSges(6,2)) * qJD(3);
t201 = t216 * qJDD(3) - t214 * t224;
t204 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t226;
t177 = m(6) * t179 + qJDD(5) * mrSges(6,1) - t201 * mrSges(6,3) + qJD(5) * t204 - t199 * t225;
t180 = t214 * t183 + t216 * t202;
t200 = -t214 * qJDD(3) - t216 * t224;
t205 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t225;
t178 = m(6) * t180 - qJDD(5) * mrSges(6,2) + t200 * mrSges(6,3) - qJD(5) * t205 - t199 * t226;
t223 = t216 * t177 + t214 * t178;
t221 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t187;
t185 = -qJDD(3) * pkin(3) + t222;
t220 = -m(5) * t185 + t218 * mrSges(5,3) - t223;
t182 = t227 * t218 + t221;
t184 = t218 * pkin(3) - t221;
t219 = -m(5) * t184 + m(6) * t182 - t200 * mrSges(6,1) + t218 * mrSges(5,2) + t201 * mrSges(6,2) + qJDD(3) * mrSges(5,3) + t204 * t226 + t205 * t225;
t195 = Ifges(6,5) * qJD(5) + (t216 * Ifges(6,1) - t214 * Ifges(6,4)) * qJD(3);
t194 = Ifges(6,6) * qJD(5) + (t216 * Ifges(6,4) - t214 * Ifges(6,2)) * qJD(3);
t176 = m(4) * t187 - t218 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t219;
t175 = qJDD(3) * mrSges(5,2) - t220;
t174 = m(4) * t186 - t218 * mrSges(4,2) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + t220;
t1 = [m(2) * t207 + t210 * (m(3) * t192 - t215 * t174 + t217 * t176) + t212 * (m(3) * t191 + t217 * t174 + t215 * t176); -t214 * t177 + t216 * t178 + (m(3) + m(4) + m(5)) * t202; mrSges(4,1) * t186 - mrSges(4,2) * t187 + mrSges(5,2) * t185 - mrSges(5,3) * t184 + t216 * (mrSges(6,2) * t182 - mrSges(6,3) * t179 + Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * qJDD(5) - qJD(5) * t194) - t214 * (-mrSges(6,1) * t182 + mrSges(6,3) * t180 + Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * qJDD(5) + qJD(5) * t195) - pkin(6) * t223 - pkin(3) * t175 + qJ(4) * t219 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3); t175; mrSges(6,1) * t179 - mrSges(6,2) * t180 + Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * qJDD(5) + (t216 * t194 + t214 * t195) * qJD(3);];
tauJ = t1;
