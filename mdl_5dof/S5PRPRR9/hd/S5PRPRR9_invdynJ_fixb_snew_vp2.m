% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRR9
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:41
% EndTime: 2019-12-31 17:39:41
% DurationCPUTime: 0.20s
% Computational Cost: add. (1301->100), mult. (1793->129), div. (0->0), fcn. (874->8), ass. (0->51)
t238 = -pkin(2) - pkin(3);
t217 = -qJD(2) + qJD(4);
t222 = sin(qJ(5));
t237 = t217 * t222;
t225 = cos(qJ(5));
t236 = t217 * t225;
t228 = qJD(2) ^ 2;
t220 = sin(pkin(8));
t221 = cos(pkin(8));
t212 = t220 * g(1) - t221 * g(2);
t213 = -t221 * g(1) - t220 * g(2);
t224 = sin(qJ(2));
t227 = cos(qJ(2));
t235 = t224 * t212 + t227 * t213;
t232 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t235;
t195 = t238 * t228 + t232;
t233 = t227 * t212 - t224 * t213;
t230 = -t228 * qJ(3) + qJDD(3) - t233;
t196 = t238 * qJDD(2) + t230;
t223 = sin(qJ(4));
t226 = cos(qJ(4));
t192 = t226 * t195 + t223 * t196;
t234 = qJD(5) * t217;
t215 = t217 ^ 2;
t216 = -qJDD(2) + qJDD(4);
t190 = -t215 * pkin(4) + t216 * pkin(7) + t192;
t219 = g(3) - qJDD(1);
t187 = -t222 * t190 + t225 * t219;
t204 = (-mrSges(6,1) * t225 + mrSges(6,2) * t222) * t217;
t205 = t222 * t216 + t225 * t234;
t211 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t236;
t185 = m(6) * t187 + qJDD(5) * mrSges(6,1) - t205 * mrSges(6,3) + qJD(5) * t211 - t204 * t237;
t188 = t225 * t190 + t222 * t219;
t206 = t225 * t216 - t222 * t234;
t210 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t237;
t186 = m(6) * t188 - qJDD(5) * mrSges(6,2) + t206 * mrSges(6,3) - qJD(5) * t210 + t204 * t236;
t179 = -t222 * t185 + t225 * t186;
t178 = m(5) * t192 - t215 * mrSges(5,1) - t216 * mrSges(5,2) + t179;
t191 = -t223 * t195 + t226 * t196;
t189 = -t216 * pkin(4) - t215 * pkin(7) - t191;
t183 = -m(6) * t189 + t206 * mrSges(6,1) - t205 * mrSges(6,2) - t210 * t237 + t211 * t236;
t182 = m(5) * t191 + t216 * mrSges(5,1) - t215 * mrSges(5,2) + t183;
t231 = t223 * t178 + t226 * t182;
t199 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t222 + Ifges(6,6) * t225) * t217;
t200 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t222 + Ifges(6,2) * t225) * t217;
t201 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t222 + Ifges(6,4) * t225) * t217;
t229 = mrSges(5,1) * t191 - mrSges(5,2) * t192 + Ifges(5,3) * t216 + pkin(4) * t183 + pkin(7) * t179 + t225 * (-mrSges(6,1) * t189 + mrSges(6,3) * t188 + Ifges(6,4) * t205 + Ifges(6,2) * t206 + Ifges(6,6) * qJDD(5) + qJD(5) * t201 - t199 * t237) + t222 * (mrSges(6,2) * t189 - mrSges(6,3) * t187 + Ifges(6,1) * t205 + Ifges(6,4) * t206 + Ifges(6,5) * qJDD(5) - qJD(5) * t200 + t199 * t236);
t198 = -qJDD(2) * pkin(2) + t230;
t197 = -t228 * pkin(2) + t232;
t177 = m(4) * t198 - qJDD(2) * mrSges(4,1) - t228 * mrSges(4,3) + t231;
t1 = [-t225 * t185 - t222 * t186 + (-m(2) - m(3) - m(4) - m(5)) * t219; -pkin(2) * t177 + qJ(3) * (m(4) * t197 - t228 * mrSges(4,1) + t226 * t178 - t223 * t182) + mrSges(3,1) * t233 - mrSges(3,2) * t235 - mrSges(4,1) * t198 + mrSges(4,3) * t197 - pkin(3) * t231 + (qJ(3) * mrSges(4,3) + Ifges(4,2) + Ifges(3,3)) * qJDD(2) - t229; t177; t229; mrSges(6,1) * t187 - mrSges(6,2) * t188 + Ifges(6,5) * t205 + Ifges(6,6) * t206 + Ifges(6,3) * qJDD(5) + (t200 * t222 - t201 * t225) * t217;];
tauJ = t1;
