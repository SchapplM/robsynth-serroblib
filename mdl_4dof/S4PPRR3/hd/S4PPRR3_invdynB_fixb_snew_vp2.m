% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PPRR3
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:24
% EndTime: 2019-12-31 16:17:24
% DurationCPUTime: 0.32s
% Computational Cost: add. (2024->125), mult. (3420->160), div. (0->0), fcn. (1848->6), ass. (0->53)
t213 = -mrSges(2,2) + mrSges(3,3);
t195 = sin(pkin(6));
t196 = cos(pkin(6));
t189 = g(1) * t195 - g(2) * t196;
t187 = qJDD(2) - t189;
t190 = -g(1) * t196 - g(2) * t195;
t198 = sin(qJ(3));
t200 = cos(qJ(3));
t176 = t198 * t187 + t200 * t190;
t201 = qJD(3) ^ 2;
t174 = -pkin(3) * t201 + qJDD(3) * pkin(5) + t176;
t194 = g(3) - qJDD(1);
t197 = sin(qJ(4));
t199 = cos(qJ(4));
t171 = -t174 * t197 + t194 * t199;
t184 = (-mrSges(5,1) * t199 + mrSges(5,2) * t197) * qJD(3);
t209 = qJD(3) * qJD(4);
t185 = qJDD(3) * t197 + t199 * t209;
t210 = qJD(3) * t199;
t192 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t210;
t211 = qJD(3) * t197;
t169 = m(5) * t171 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t185 + qJD(4) * t192 - t184 * t211;
t172 = t174 * t199 + t194 * t197;
t186 = qJDD(3) * t199 - t197 * t209;
t191 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t211;
t170 = m(5) * t172 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t186 - qJD(4) * t191 + t184 * t210;
t206 = -t169 * t197 + t199 * t170;
t162 = m(4) * t176 - mrSges(4,1) * t201 - qJDD(3) * mrSges(4,2) + t206;
t175 = t187 * t200 - t190 * t198;
t173 = -qJDD(3) * pkin(3) - pkin(5) * t201 - t175;
t202 = -m(5) * t173 + t186 * mrSges(5,1) - mrSges(5,2) * t185 - t191 * t211 + t192 * t210;
t167 = m(4) * t175 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t201 + t202;
t160 = t162 * t198 + t167 * t200;
t203 = -m(3) * t187 - t160;
t158 = m(2) * t189 + t203;
t207 = t200 * t162 - t167 * t198;
t205 = m(3) * t190 + t207;
t159 = m(2) * t190 + t205;
t212 = t196 * t158 + t195 * t159;
t208 = -t158 * t195 + t196 * t159;
t164 = t199 * t169 + t197 * t170;
t204 = -m(3) * t194 - t164;
t179 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t197 + Ifges(5,4) * t199) * qJD(3);
t178 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t197 + Ifges(5,2) * t199) * qJD(3);
t177 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t197 + Ifges(5,6) * t199) * qJD(3);
t166 = mrSges(5,2) * t173 - mrSges(5,3) * t171 + Ifges(5,1) * t185 + Ifges(5,4) * t186 + Ifges(5,5) * qJDD(4) - qJD(4) * t178 + t177 * t210;
t165 = -mrSges(5,1) * t173 + mrSges(5,3) * t172 + Ifges(5,4) * t185 + Ifges(5,2) * t186 + Ifges(5,6) * qJDD(4) + qJD(4) * t179 - t177 * t211;
t163 = -m(4) * t194 + t204;
t154 = -mrSges(4,1) * t194 - mrSges(5,1) * t171 + mrSges(5,2) * t172 + mrSges(4,3) * t176 + t201 * Ifges(4,5) - Ifges(5,5) * t185 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t186 - Ifges(5,3) * qJDD(4) - pkin(3) * t164 + (-t178 * t197 + t179 * t199) * qJD(3);
t153 = mrSges(4,2) * t194 - mrSges(4,3) * t175 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t201 - pkin(5) * t164 - t165 * t197 + t166 * t199;
t152 = mrSges(3,2) * t187 - mrSges(2,3) * t189 - pkin(4) * t160 - qJ(2) * t163 + t200 * t153 - t198 * t154 + t194 * t213;
t151 = -t198 * t153 - t200 * t154 + pkin(2) * t164 - pkin(4) * t207 - pkin(1) * t163 + (m(4) * pkin(2) + mrSges(2,1) + mrSges(3,1)) * t194 + (mrSges(3,2) + mrSges(2,3)) * t190;
t1 = [-m(1) * g(1) + t208; -m(1) * g(2) + t212; -m(1) * g(3) + (-m(2) - m(4)) * t194 + t204; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t212 - t195 * t151 + t196 * t152; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t208 + t196 * t151 + t195 * t152; -mrSges(1,1) * g(2) + mrSges(2,1) * t189 - mrSges(3,1) * t187 - mrSges(4,1) * t175 + mrSges(1,2) * g(1) + mrSges(4,2) * t176 - Ifges(4,3) * qJDD(3) + pkin(1) * t203 - pkin(2) * t160 - pkin(3) * t202 - pkin(5) * t206 + qJ(2) * t205 - t199 * t165 - t197 * t166 + t190 * t213;];
tauB = t1;
