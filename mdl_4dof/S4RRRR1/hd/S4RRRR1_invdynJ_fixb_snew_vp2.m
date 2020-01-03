% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:11
% EndTime: 2019-12-31 17:22:11
% DurationCPUTime: 0.23s
% Computational Cost: add. (1761->91), mult. (1957->125), div. (0->0), fcn. (1004->8), ass. (0->47)
t204 = qJD(1) + qJD(2);
t200 = qJD(3) + t204;
t205 = sin(qJ(4));
t222 = t200 * t205;
t209 = cos(qJ(4));
t221 = t200 * t209;
t208 = sin(qJ(1));
t212 = cos(qJ(1));
t218 = t208 * g(1) - t212 * g(2);
t194 = qJDD(1) * pkin(1) + t218;
t216 = -t212 * g(1) - t208 * g(2);
t195 = -qJD(1) ^ 2 * pkin(1) + t216;
t207 = sin(qJ(2));
t211 = cos(qJ(2));
t179 = t211 * t194 - t207 * t195;
t203 = qJDD(1) + qJDD(2);
t176 = t203 * pkin(2) + t179;
t180 = t207 * t194 + t211 * t195;
t202 = t204 ^ 2;
t177 = -t202 * pkin(2) + t180;
t206 = sin(qJ(3));
t210 = cos(qJ(3));
t173 = t206 * t176 + t210 * t177;
t198 = t200 ^ 2;
t199 = qJDD(3) + t203;
t170 = -t198 * pkin(3) + t199 * pkin(7) + t173;
t167 = -t209 * g(3) - t205 * t170;
t168 = -t205 * g(3) + t209 * t170;
t186 = (-mrSges(5,1) * t209 + mrSges(5,2) * t205) * t200;
t219 = qJD(4) * t200;
t187 = t205 * t199 + t209 * t219;
t188 = t209 * t199 - t205 * t219;
t192 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t222;
t193 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t221;
t217 = -t205 * (m(5) * t167 + qJDD(4) * mrSges(5,1) - t187 * mrSges(5,3) + qJD(4) * t193 - t186 * t222) + t209 * (m(5) * t168 - qJDD(4) * mrSges(5,2) + t188 * mrSges(5,3) - qJD(4) * t192 + t186 * t221);
t158 = m(4) * t173 - t198 * mrSges(4,1) - t199 * mrSges(4,2) + t217;
t172 = t210 * t176 - t206 * t177;
t169 = -t199 * pkin(3) - t198 * pkin(7) - t172;
t214 = -m(5) * t169 + t188 * mrSges(5,1) - t187 * mrSges(5,2) - t192 * t222 + t193 * t221;
t163 = m(4) * t172 + t199 * mrSges(4,1) - t198 * mrSges(4,2) + t214;
t220 = t206 * t158 + t210 * t163;
t181 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t205 + Ifges(5,6) * t209) * t200;
t182 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t205 + Ifges(5,2) * t209) * t200;
t183 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t205 + Ifges(5,4) * t209) * t200;
t215 = -mrSges(4,2) * t173 + pkin(7) * t217 + t205 * (mrSges(5,2) * t169 - mrSges(5,3) * t167 + Ifges(5,1) * t187 + Ifges(5,4) * t188 + Ifges(5,5) * qJDD(4) - qJD(4) * t182 + t181 * t221) + t209 * (-mrSges(5,1) * t169 + mrSges(5,3) * t168 + Ifges(5,4) * t187 + Ifges(5,2) * t188 + Ifges(5,6) * qJDD(4) + qJD(4) * t183 - t181 * t222) + pkin(3) * t214 + mrSges(4,1) * t172 + Ifges(4,3) * t199;
t213 = mrSges(3,1) * t179 - mrSges(3,2) * t180 + Ifges(3,3) * t203 + pkin(2) * t220 + t215;
t1 = [t213 + pkin(1) * (t207 * (m(3) * t180 - t202 * mrSges(3,1) - t203 * mrSges(3,2) + t210 * t158 - t206 * t163) + t211 * (m(3) * t179 + t203 * mrSges(3,1) - t202 * mrSges(3,2) + t220)) - mrSges(2,2) * t216 + mrSges(2,1) * t218 + Ifges(2,3) * qJDD(1); t213; t215; mrSges(5,1) * t167 - mrSges(5,2) * t168 + Ifges(5,5) * t187 + Ifges(5,6) * t188 + Ifges(5,3) * qJDD(4) + (t182 * t205 - t183 * t209) * t200;];
tauJ = t1;
