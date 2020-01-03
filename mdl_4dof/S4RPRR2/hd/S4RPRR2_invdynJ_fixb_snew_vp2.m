% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RPRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:10
% EndTime: 2019-12-31 16:48:10
% DurationCPUTime: 0.19s
% Computational Cost: add. (1061->93), mult. (1492->127), div. (0->0), fcn. (764->8), ass. (0->47)
t198 = qJD(1) + qJD(3);
t202 = sin(qJ(4));
t217 = t198 * t202;
t205 = cos(qJ(4));
t216 = t198 * t205;
t204 = sin(qJ(1));
t207 = cos(qJ(1));
t213 = t204 * g(1) - t207 * g(2);
t190 = qJDD(1) * pkin(1) + t213;
t208 = qJD(1) ^ 2;
t211 = -t207 * g(1) - t204 * g(2);
t191 = -t208 * pkin(1) + t211;
t200 = sin(pkin(7));
t201 = cos(pkin(7));
t177 = t201 * t190 - t200 * t191;
t175 = qJDD(1) * pkin(2) + t177;
t178 = t200 * t190 + t201 * t191;
t176 = -t208 * pkin(2) + t178;
t203 = sin(qJ(3));
t206 = cos(qJ(3));
t172 = t203 * t175 + t206 * t176;
t196 = t198 ^ 2;
t197 = qJDD(1) + qJDD(3);
t169 = -t196 * pkin(3) + t197 * pkin(6) + t172;
t199 = -g(3) + qJDD(2);
t166 = -t202 * t169 + t205 * t199;
t184 = (-mrSges(5,1) * t205 + mrSges(5,2) * t202) * t198;
t214 = qJD(4) * t198;
t185 = t202 * t197 + t205 * t214;
t193 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t216;
t164 = m(5) * t166 + qJDD(4) * mrSges(5,1) - t185 * mrSges(5,3) + qJD(4) * t193 - t184 * t217;
t167 = t205 * t169 + t202 * t199;
t186 = t205 * t197 - t202 * t214;
t192 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t217;
t165 = m(5) * t167 - qJDD(4) * mrSges(5,2) + t186 * mrSges(5,3) - qJD(4) * t192 + t184 * t216;
t212 = -t202 * t164 + t205 * t165;
t156 = m(4) * t172 - t196 * mrSges(4,1) - t197 * mrSges(4,2) + t212;
t171 = t206 * t175 - t203 * t176;
t168 = -t197 * pkin(3) - t196 * pkin(6) - t171;
t209 = -m(5) * t168 + t186 * mrSges(5,1) - t185 * mrSges(5,2) - t192 * t217 + t193 * t216;
t161 = m(4) * t171 + t197 * mrSges(4,1) - t196 * mrSges(4,2) + t209;
t215 = t203 * t156 + t206 * t161;
t179 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t202 + Ifges(5,6) * t205) * t198;
t180 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t202 + Ifges(5,2) * t205) * t198;
t181 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t202 + Ifges(5,4) * t205) * t198;
t210 = -mrSges(4,2) * t172 + pkin(6) * t212 + t202 * (mrSges(5,2) * t168 - mrSges(5,3) * t166 + Ifges(5,1) * t185 + Ifges(5,4) * t186 + Ifges(5,5) * qJDD(4) - qJD(4) * t180 + t179 * t216) + t205 * (-mrSges(5,1) * t168 + mrSges(5,3) * t167 + Ifges(5,4) * t185 + Ifges(5,2) * t186 + Ifges(5,6) * qJDD(4) + qJD(4) * t181 - t179 * t217) + pkin(3) * t209 + mrSges(4,1) * t171 + Ifges(4,3) * t197;
t1 = [pkin(1) * (t200 * (m(3) * t178 - t208 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t206 * t156 - t203 * t161) + t201 * (m(3) * t177 + qJDD(1) * mrSges(3,1) - t208 * mrSges(3,2) + t215)) - mrSges(2,2) * t211 + mrSges(2,1) * t213 + mrSges(3,1) * t177 - mrSges(3,2) * t178 + pkin(2) * t215 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * qJDD(1) + t210; t205 * t164 + t202 * t165 + (m(3) + m(4)) * t199; t210; mrSges(5,1) * t166 - mrSges(5,2) * t167 + Ifges(5,5) * t185 + Ifges(5,6) * t186 + Ifges(5,3) * qJDD(4) + (t180 * t202 - t181 * t205) * t198;];
tauJ = t1;
