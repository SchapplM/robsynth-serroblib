% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4PPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:43
% EndTime: 2019-12-31 16:19:44
% DurationCPUTime: 0.30s
% Computational Cost: add. (1971->124), mult. (3228->159), div. (0->0), fcn. (1696->6), ass. (0->52)
t224 = mrSges(2,2) - mrSges(3,3);
t207 = sin(pkin(6));
t208 = cos(pkin(6));
t200 = t207 * g(1) - t208 * g(2);
t198 = qJDD(2) - t200;
t206 = -g(3) + qJDD(1);
t210 = sin(qJ(3));
t212 = cos(qJ(3));
t188 = t210 * t198 + t212 * t206;
t213 = qJD(3) ^ 2;
t186 = -t213 * pkin(3) + qJDD(3) * pkin(5) + t188;
t201 = t208 * g(1) + t207 * g(2);
t209 = sin(qJ(4));
t211 = cos(qJ(4));
t183 = -t209 * t186 - t211 * t201;
t195 = (-mrSges(5,1) * t211 + mrSges(5,2) * t209) * qJD(3);
t220 = qJD(3) * qJD(4);
t196 = t209 * qJDD(3) + t211 * t220;
t221 = qJD(3) * t211;
t203 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t221;
t222 = qJD(3) * t209;
t181 = m(5) * t183 + qJDD(4) * mrSges(5,1) - t196 * mrSges(5,3) + qJD(4) * t203 - t195 * t222;
t184 = t211 * t186 - t209 * t201;
t197 = t211 * qJDD(3) - t209 * t220;
t202 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t222;
t182 = m(5) * t184 - qJDD(4) * mrSges(5,2) + t197 * mrSges(5,3) - qJD(4) * t202 + t195 * t221;
t216 = -t209 * t181 + t211 * t182;
t173 = m(4) * t188 - t213 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t216;
t187 = t212 * t198 - t210 * t206;
t185 = -qJDD(3) * pkin(3) - t213 * pkin(5) - t187;
t214 = -m(5) * t185 + t197 * mrSges(5,1) - t196 * mrSges(5,2) - t202 * t222 + t203 * t221;
t177 = m(4) * t187 + qJDD(3) * mrSges(4,1) - t213 * mrSges(4,2) + t214;
t168 = t210 * t173 + t212 * t177;
t215 = -m(3) * t198 - t168;
t166 = m(2) * t200 + t215;
t174 = t211 * t181 + t209 * t182;
t219 = m(4) * t201 - t174;
t171 = (-m(2) - m(3)) * t201 - t219;
t223 = t208 * t166 + t207 * t171;
t218 = -t207 * t166 + t208 * t171;
t217 = t212 * t173 - t210 * t177;
t167 = m(3) * t206 + t217;
t191 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t209 + Ifges(5,4) * t211) * qJD(3);
t190 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t209 + Ifges(5,2) * t211) * qJD(3);
t189 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t209 + Ifges(5,6) * t211) * qJD(3);
t176 = mrSges(5,2) * t185 - mrSges(5,3) * t183 + Ifges(5,1) * t196 + Ifges(5,4) * t197 + Ifges(5,5) * qJDD(4) - qJD(4) * t190 + t189 * t221;
t175 = -mrSges(5,1) * t185 + mrSges(5,3) * t184 + Ifges(5,4) * t196 + Ifges(5,2) * t197 + Ifges(5,6) * qJDD(4) + qJD(4) * t191 - t189 * t222;
t164 = mrSges(4,1) * t201 - mrSges(5,1) * t183 + mrSges(5,2) * t184 + mrSges(4,3) * t188 + t213 * Ifges(4,5) - Ifges(5,5) * t196 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t197 - Ifges(5,3) * qJDD(4) - pkin(3) * t174 + (-t190 * t209 + t191 * t211) * qJD(3);
t163 = -mrSges(4,2) * t201 - mrSges(4,3) * t187 + Ifges(4,5) * qJDD(3) - t213 * Ifges(4,6) - pkin(5) * t174 - t209 * t175 + t211 * t176;
t162 = mrSges(3,1) * t198 + mrSges(4,1) * t187 - mrSges(4,2) * t188 - mrSges(2,3) * t200 + Ifges(4,3) * qJDD(3) + pkin(2) * t168 + pkin(3) * t214 + pkin(5) * t216 - qJ(2) * t167 + t211 * t175 + t209 * t176 + t224 * t206;
t161 = -t210 * t163 - t212 * t164 - pkin(2) * t219 - pkin(4) * t217 - pkin(1) * t167 + (-mrSges(2,1) + mrSges(3,2)) * t206 + (-mrSges(3,1) - mrSges(2,3)) * t201;
t1 = [-m(1) * g(1) + t218; -m(1) * g(2) + t223; -m(1) * g(3) + m(2) * t206 + t167; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t223 - t207 * t161 + t208 * t162; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t218 + t208 * t161 + t207 * t162; pkin(1) * t215 - qJ(2) * t219 + t212 * t163 - t210 * t164 - pkin(4) * t168 + mrSges(2,1) * t200 + mrSges(3,2) * t198 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-qJ(2) * m(3) + t224) * t201;];
tauB = t1;
