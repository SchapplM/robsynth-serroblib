% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_invdynJ_fixb_snew_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:08:43
% EndTime: 2019-12-05 18:08:45
% DurationCPUTime: 0.47s
% Computational Cost: add. (1821->178), mult. (3217->230), div. (0->0), fcn. (2261->8), ass. (0->64)
t216 = cos(qJ(3));
t225 = qJD(1) * t216;
t212 = sin(qJ(3));
t224 = t212 * qJD(1);
t223 = qJD(1) * qJD(3);
t213 = sin(qJ(1));
t217 = cos(qJ(1));
t222 = -g(1) * t217 - g(2) * t213;
t221 = g(1) * t213 - g(2) * t217;
t211 = sin(qJ(4));
t215 = cos(qJ(4));
t202 = qJD(3) * t215 - t211 * t224;
t206 = qJDD(1) * t216 - t212 * t223;
t203 = qJD(3) * t211 + t215 * t224;
t205 = qJDD(1) * t212 + t216 * t223;
t185 = -qJD(4) * t203 + qJDD(3) * t215 - t205 * t211;
t186 = qJD(4) * t202 + qJDD(3) * t211 + t205 * t215;
t209 = qJD(4) - t225;
t210 = sin(qJ(5));
t214 = cos(qJ(5));
t189 = t203 * t214 + t209 * t210;
t201 = qJDD(4) - t206;
t169 = -qJD(5) * t189 - t186 * t210 + t201 * t214;
t188 = -t203 * t210 + t209 * t214;
t170 = qJD(5) * t188 + t186 * t214 + t201 * t210;
t195 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t222;
t191 = -g(3) * t212 + t195 * t216;
t218 = qJD(1) ^ 2;
t200 = -qJ(2) * t218 + qJDD(2) - t221;
t178 = t191 * t215 + t200 * t211;
t190 = g(3) * t216 + t195 * t212;
t171 = -t178 * t210 + t190 * t214;
t172 = t178 * t214 + t190 * t210;
t199 = qJD(5) - t202;
t174 = Ifges(6,4) * t189 + Ifges(6,2) * t188 + Ifges(6,6) * t199;
t175 = Ifges(6,1) * t189 + Ifges(6,4) * t188 + Ifges(6,5) * t199;
t184 = qJDD(5) - t185;
t220 = mrSges(6,1) * t171 - mrSges(6,2) * t172 + Ifges(6,5) * t170 + Ifges(6,6) * t169 + Ifges(6,3) * t184 + t174 * t189 - t175 * t188;
t173 = Ifges(6,5) * t189 + Ifges(6,6) * t188 + Ifges(6,3) * t199;
t177 = t191 * t211 - t215 * t200;
t165 = -mrSges(6,1) * t177 + mrSges(6,3) * t172 + Ifges(6,4) * t170 + Ifges(6,2) * t169 + Ifges(6,6) * t184 - t173 * t189 + t175 * t199;
t166 = mrSges(6,2) * t177 - mrSges(6,3) * t171 + Ifges(6,1) * t170 + Ifges(6,4) * t169 + Ifges(6,5) * t184 + t173 * t188 - t174 * t199;
t182 = Ifges(5,4) * t203 + Ifges(5,2) * t202 + Ifges(5,6) * t209;
t183 = Ifges(5,1) * t203 + Ifges(5,4) * t202 + Ifges(5,5) * t209;
t219 = -mrSges(5,1) * t177 - mrSges(5,2) * t178 + Ifges(5,5) * t186 + Ifges(5,6) * t185 + Ifges(5,3) * t201 + t214 * t165 + t210 * t166 + t203 * t182 - t202 * t183;
t208 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t225;
t207 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t224;
t204 = (-mrSges(4,1) * t216 + mrSges(4,2) * t212) * qJD(1);
t198 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t212 + Ifges(4,4) * t216) * qJD(1);
t197 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t212 + Ifges(4,2) * t216) * qJD(1);
t193 = mrSges(5,1) * t209 - mrSges(5,3) * t203;
t192 = -mrSges(5,2) * t209 + mrSges(5,3) * t202;
t187 = -mrSges(5,1) * t202 + mrSges(5,2) * t203;
t181 = Ifges(5,5) * t203 + Ifges(5,6) * t202 + Ifges(5,3) * t209;
t180 = mrSges(6,1) * t199 - mrSges(6,3) * t189;
t179 = -mrSges(6,2) * t199 + mrSges(6,3) * t188;
t176 = -mrSges(6,1) * t188 + mrSges(6,2) * t189;
t168 = m(6) * t172 - mrSges(6,2) * t184 + mrSges(6,3) * t169 + t176 * t188 - t180 * t199;
t167 = m(6) * t171 + mrSges(6,1) * t184 - mrSges(6,3) * t170 - t176 * t189 + t179 * t199;
t164 = mrSges(5,1) * t201 + mrSges(6,1) * t169 - mrSges(6,2) * t170 - mrSges(5,3) * t186 + t179 * t188 - t180 * t189 - t187 * t203 + t192 * t209 + (-m(5) - m(6)) * t177;
t163 = -mrSges(5,1) * t190 + mrSges(5,3) * t178 + Ifges(5,4) * t186 + Ifges(5,2) * t185 + Ifges(5,6) * t201 - t181 * t203 + t183 * t209 - t220;
t162 = m(5) * t178 - mrSges(5,2) * t201 + mrSges(5,3) * t185 - t167 * t210 + t168 * t214 + t187 * t202 - t193 * t209;
t161 = mrSges(5,2) * t190 + mrSges(5,3) * t177 + Ifges(5,1) * t186 + Ifges(5,4) * t185 + Ifges(5,5) * t201 - t165 * t210 + t166 * t214 + t181 * t202 - t182 * t209;
t1 = [mrSges(2,1) * t221 - mrSges(2,2) * t222 - mrSges(3,1) * t200 + mrSges(3,3) * t195 + t212 * (mrSges(4,2) * t200 + mrSges(4,3) * t190 + Ifges(4,1) * t205 + Ifges(4,4) * t206 + Ifges(4,5) * qJDD(3) - qJD(3) * t197 + t161 * t215 - t163 * t211) + t216 * (-mrSges(4,1) * t200 + mrSges(4,3) * t191 + Ifges(4,4) * t205 + Ifges(4,2) * t206 + Ifges(4,6) * qJDD(3) + qJD(3) * t198 - t219) + qJ(2) * (m(3) * t195 - t218 * mrSges(3,1) + t216 * (m(4) * t191 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t206 - qJD(3) * t207 + t162 * t215 - t164 * t211 + t204 * t225) - t212 * (-t204 * t224 + qJDD(3) * mrSges(4,1) + mrSges(5,1) * t185 - mrSges(5,2) * t186 - mrSges(4,3) * t205 + qJD(3) * t208 - t167 * t214 - t168 * t210 + t192 * t202 - t193 * t203 + (-m(4) - m(5)) * t190)) + (mrSges(3,3) * qJ(2) + Ifges(3,2) + Ifges(2,3)) * qJDD(1); -qJDD(1) * mrSges(3,1) - mrSges(4,1) * t206 + mrSges(4,2) * t205 - mrSges(3,3) * t218 + t162 * t211 + t164 * t215 + (m(3) + m(4)) * t200 + (t207 * t212 - t208 * t216) * qJD(1); -mrSges(4,1) * t190 - mrSges(4,2) * t191 + Ifges(4,5) * t205 + Ifges(4,6) * t206 + Ifges(4,3) * qJDD(3) + t161 * t211 + t163 * t215 + (t212 * t197 - t216 * t198) * qJD(1); t219; t220;];
tauJ = t1;
