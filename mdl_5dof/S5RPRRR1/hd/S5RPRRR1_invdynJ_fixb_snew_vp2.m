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
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
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
% StartTime: 2019-07-18 13:24:59
% EndTime: 2019-07-18 13:25:01
% DurationCPUTime: 0.46s
% Computational Cost: add. (1821->178), mult. (3217->230), div. (0->0), fcn. (2261->8), ass. (0->64)
t212 = sin(qJ(3));
t225 = qJD(1) * t212;
t216 = cos(qJ(3));
t224 = t216 * qJD(1);
t223 = qJD(1) * qJD(3);
t213 = sin(qJ(1));
t217 = cos(qJ(1));
t222 = -t217 * g(1) - t213 * g(2);
t221 = t213 * g(1) - t217 * g(2);
t211 = sin(qJ(4));
t215 = cos(qJ(4));
t202 = t215 * qJD(3) - t211 * t225;
t206 = t216 * qJDD(1) - t212 * t223;
t203 = t211 * qJD(3) + t215 * t225;
t205 = t212 * qJDD(1) + t216 * t223;
t185 = -t203 * qJD(4) + t215 * qJDD(3) - t211 * t205;
t186 = t202 * qJD(4) + t211 * qJDD(3) + t215 * t205;
t209 = qJD(4) - t224;
t210 = sin(qJ(5));
t214 = cos(qJ(5));
t189 = t214 * t203 + t210 * t209;
t201 = qJDD(4) - t206;
t169 = -t189 * qJD(5) - t210 * t186 + t214 * t201;
t188 = -t210 * t203 + t214 * t209;
t170 = t188 * qJD(5) + t214 * t186 + t210 * t201;
t195 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t222;
t191 = -t212 * g(3) + t216 * t195;
t218 = qJD(1) ^ 2;
t200 = -t218 * qJ(2) + qJDD(2) - t221;
t178 = t215 * t191 + t211 * t200;
t190 = t216 * g(3) + t212 * t195;
t171 = -t210 * t178 + t214 * t190;
t172 = t214 * t178 + t210 * t190;
t199 = qJD(5) - t202;
t174 = Ifges(6,4) * t189 + Ifges(6,2) * t188 + Ifges(6,6) * t199;
t175 = Ifges(6,1) * t189 + Ifges(6,4) * t188 + Ifges(6,5) * t199;
t184 = qJDD(5) - t185;
t220 = mrSges(6,1) * t171 - mrSges(6,2) * t172 + Ifges(6,5) * t170 + Ifges(6,6) * t169 + Ifges(6,3) * t184 + t189 * t174 - t188 * t175;
t173 = Ifges(6,5) * t189 + Ifges(6,6) * t188 + Ifges(6,3) * t199;
t177 = t211 * t191 - t215 * t200;
t165 = -mrSges(6,1) * t177 + mrSges(6,3) * t172 + Ifges(6,4) * t170 + Ifges(6,2) * t169 + Ifges(6,6) * t184 - t189 * t173 + t199 * t175;
t166 = mrSges(6,2) * t177 - mrSges(6,3) * t171 + Ifges(6,1) * t170 + Ifges(6,4) * t169 + Ifges(6,5) * t184 + t188 * t173 - t199 * t174;
t182 = Ifges(5,4) * t203 + Ifges(5,2) * t202 + Ifges(5,6) * t209;
t183 = Ifges(5,1) * t203 + Ifges(5,4) * t202 + Ifges(5,5) * t209;
t219 = -mrSges(5,1) * t177 - mrSges(5,2) * t178 + Ifges(5,5) * t186 + Ifges(5,6) * t185 + Ifges(5,3) * t201 + t214 * t165 + t210 * t166 + t203 * t182 - t202 * t183;
t208 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t224;
t207 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t225;
t204 = (-t216 * mrSges(4,1) + t212 * mrSges(4,2)) * qJD(1);
t198 = Ifges(4,5) * qJD(3) + (t212 * Ifges(4,1) + t216 * Ifges(4,4)) * qJD(1);
t197 = Ifges(4,6) * qJD(3) + (t212 * Ifges(4,4) + t216 * Ifges(4,2)) * qJD(1);
t193 = t209 * mrSges(5,1) - t203 * mrSges(5,3);
t192 = -t209 * mrSges(5,2) + t202 * mrSges(5,3);
t187 = -t202 * mrSges(5,1) + t203 * mrSges(5,2);
t181 = Ifges(5,5) * t203 + Ifges(5,6) * t202 + Ifges(5,3) * t209;
t180 = t199 * mrSges(6,1) - t189 * mrSges(6,3);
t179 = -t199 * mrSges(6,2) + t188 * mrSges(6,3);
t176 = -t188 * mrSges(6,1) + t189 * mrSges(6,2);
t168 = m(6) * t172 - t184 * mrSges(6,2) + t169 * mrSges(6,3) + t188 * t176 - t199 * t180;
t167 = m(6) * t171 + t184 * mrSges(6,1) - t170 * mrSges(6,3) - t189 * t176 + t199 * t179;
t164 = t201 * mrSges(5,1) + t169 * mrSges(6,1) - t170 * mrSges(6,2) - t186 * mrSges(5,3) + t188 * t179 - t189 * t180 - t203 * t187 + t209 * t192 + (-m(5) - m(6)) * t177;
t163 = -mrSges(5,1) * t190 + mrSges(5,3) * t178 + Ifges(5,4) * t186 + Ifges(5,2) * t185 + Ifges(5,6) * t201 - t203 * t181 + t209 * t183 - t220;
t162 = m(5) * t178 - t201 * mrSges(5,2) + t185 * mrSges(5,3) - t210 * t167 + t214 * t168 + t202 * t187 - t209 * t193;
t161 = mrSges(5,2) * t190 + mrSges(5,3) * t177 + Ifges(5,1) * t186 + Ifges(5,4) * t185 + Ifges(5,5) * t201 - t210 * t165 + t214 * t166 + t202 * t181 - t209 * t182;
t1 = [mrSges(2,1) * t221 - mrSges(2,2) * t222 - mrSges(3,1) * t200 + mrSges(3,3) * t195 + t212 * (mrSges(4,2) * t200 + mrSges(4,3) * t190 + Ifges(4,1) * t205 + Ifges(4,4) * t206 + Ifges(4,5) * qJDD(3) - qJD(3) * t197 + t215 * t161 - t211 * t163) + t216 * (-mrSges(4,1) * t200 + mrSges(4,3) * t191 + Ifges(4,4) * t205 + Ifges(4,2) * t206 + Ifges(4,6) * qJDD(3) + qJD(3) * t198 - t219) + qJ(2) * (m(3) * t195 - t218 * mrSges(3,1) + t216 * (m(4) * t191 - qJDD(3) * mrSges(4,2) + t206 * mrSges(4,3) - qJD(3) * t207 + t215 * t162 - t211 * t164 + t204 * t224) - t212 * (-t204 * t225 + qJDD(3) * mrSges(4,1) + t185 * mrSges(5,1) - t186 * mrSges(5,2) - t205 * mrSges(4,3) + qJD(3) * t208 - t214 * t167 - t210 * t168 + t202 * t192 - t203 * t193 + (-m(4) - m(5)) * t190)) + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1); -qJDD(1) * mrSges(3,1) - t206 * mrSges(4,1) + t205 * mrSges(4,2) - t218 * mrSges(3,3) + t211 * t162 + t215 * t164 + (m(3) + m(4)) * t200 + (t207 * t212 - t208 * t216) * qJD(1); -mrSges(4,1) * t190 - mrSges(4,2) * t191 + Ifges(4,5) * t205 + Ifges(4,6) * t206 + Ifges(4,3) * qJDD(3) + t211 * t161 + t215 * t163 + (t212 * t197 - t216 * t198) * qJD(1); t219; t220;];
tauJ  = t1;
