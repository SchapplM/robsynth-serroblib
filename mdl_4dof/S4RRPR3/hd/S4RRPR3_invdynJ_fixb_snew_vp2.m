% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR3
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:31
% EndTime: 2019-12-31 17:01:31
% DurationCPUTime: 0.25s
% Computational Cost: add. (1298->92), mult. (1704->126), div. (0->0), fcn. (874->8), ass. (0->46)
t200 = qJD(1) + qJD(2);
t204 = sin(qJ(4));
t218 = t200 * t204;
t207 = cos(qJ(4));
t217 = t200 * t207;
t206 = sin(qJ(1));
t209 = cos(qJ(1));
t214 = t206 * g(1) - t209 * g(2);
t191 = qJDD(1) * pkin(1) + t214;
t212 = -t209 * g(1) - t206 * g(2);
t192 = -qJD(1) ^ 2 * pkin(1) + t212;
t205 = sin(qJ(2));
t208 = cos(qJ(2));
t178 = t208 * t191 - t205 * t192;
t199 = qJDD(1) + qJDD(2);
t175 = t199 * pkin(2) + t178;
t179 = t205 * t191 + t208 * t192;
t198 = t200 ^ 2;
t176 = -t198 * pkin(2) + t179;
t202 = sin(pkin(7));
t203 = cos(pkin(7));
t172 = t202 * t175 + t203 * t176;
t169 = -t198 * pkin(3) + t199 * pkin(6) + t172;
t201 = -g(3) + qJDD(3);
t166 = -t204 * t169 + t207 * t201;
t185 = (-mrSges(5,1) * t207 + mrSges(5,2) * t204) * t200;
t215 = qJD(4) * t200;
t186 = t204 * t199 + t207 * t215;
t194 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t217;
t164 = m(5) * t166 + qJDD(4) * mrSges(5,1) - t186 * mrSges(5,3) + qJD(4) * t194 - t185 * t218;
t167 = t207 * t169 + t204 * t201;
t187 = t207 * t199 - t204 * t215;
t193 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t218;
t165 = m(5) * t167 - qJDD(4) * mrSges(5,2) + t187 * mrSges(5,3) - qJD(4) * t193 + t185 * t217;
t213 = -t204 * t164 + t207 * t165;
t156 = m(4) * t172 - t198 * mrSges(4,1) - t199 * mrSges(4,2) + t213;
t171 = t203 * t175 - t202 * t176;
t168 = -t199 * pkin(3) - t198 * pkin(6) - t171;
t211 = -m(5) * t168 + t187 * mrSges(5,1) - t186 * mrSges(5,2) - t193 * t218 + t194 * t217;
t161 = m(4) * t171 + t199 * mrSges(4,1) - t198 * mrSges(4,2) + t211;
t216 = t202 * t156 + t203 * t161;
t180 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t204 + Ifges(5,6) * t207) * t200;
t181 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t204 + Ifges(5,2) * t207) * t200;
t182 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t204 + Ifges(5,4) * t207) * t200;
t210 = -mrSges(3,2) * t179 - mrSges(4,2) * t172 + pkin(2) * t216 + pkin(6) * t213 + t204 * (mrSges(5,2) * t168 - mrSges(5,3) * t166 + Ifges(5,1) * t186 + Ifges(5,4) * t187 + Ifges(5,5) * qJDD(4) - qJD(4) * t181 + t180 * t217) + t207 * (-mrSges(5,1) * t168 + mrSges(5,3) * t167 + Ifges(5,4) * t186 + Ifges(5,2) * t187 + Ifges(5,6) * qJDD(4) + qJD(4) * t182 - t180 * t218) + pkin(3) * t211 + mrSges(4,1) * t171 + mrSges(3,1) * t178 + (Ifges(4,3) + Ifges(3,3)) * t199;
t1 = [t210 + Ifges(2,3) * qJDD(1) + pkin(1) * (t205 * (m(3) * t179 - t198 * mrSges(3,1) - t199 * mrSges(3,2) + t203 * t156 - t202 * t161) + t208 * (m(3) * t178 + t199 * mrSges(3,1) - t198 * mrSges(3,2) + t216)) - mrSges(2,2) * t212 + mrSges(2,1) * t214; t210; m(4) * t201 + t207 * t164 + t204 * t165; mrSges(5,1) * t166 - mrSges(5,2) * t167 + Ifges(5,5) * t186 + Ifges(5,6) * t187 + Ifges(5,3) * qJDD(4) + (t181 * t204 - t182 * t207) * t200;];
tauJ = t1;
