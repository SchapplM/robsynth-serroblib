% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR5
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynJ_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:30
% EndTime: 2019-12-31 17:03:30
% DurationCPUTime: 0.19s
% Computational Cost: add. (862->89), mult. (1028->117), div. (0->0), fcn. (452->6), ass. (0->41)
t211 = -pkin(2) - pkin(6);
t194 = qJD(1) + qJD(2);
t195 = sin(qJ(4));
t210 = t194 * t195;
t198 = cos(qJ(4));
t209 = t194 * t198;
t197 = sin(qJ(1));
t200 = cos(qJ(1));
t207 = t197 * g(1) - t200 * g(2);
t182 = qJDD(1) * pkin(1) + t207;
t206 = -t200 * g(1) - t197 * g(2);
t183 = -qJD(1) ^ 2 * pkin(1) + t206;
t196 = sin(qJ(2));
t199 = cos(qJ(2));
t169 = t196 * t182 + t199 * t183;
t208 = qJD(4) * t194;
t168 = t199 * t182 - t196 * t183;
t193 = qJDD(1) + qJDD(2);
t192 = t194 ^ 2;
t204 = -t192 * qJ(3) + qJDD(3) - t168;
t163 = t211 * t193 + t204;
t159 = t195 * g(3) + t198 * t163;
t160 = -t198 * g(3) + t195 * t163;
t176 = (mrSges(5,1) * t195 + mrSges(5,2) * t198) * t194;
t177 = -t195 * t193 - t198 * t208;
t178 = t198 * t193 - t195 * t208;
t184 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t210;
t185 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t209;
t153 = t198 * (m(5) * t159 + qJDD(4) * mrSges(5,1) - t178 * mrSges(5,3) + qJD(4) * t184 - t176 * t209) + t195 * (m(5) * t160 - qJDD(4) * mrSges(5,2) + t177 * mrSges(5,3) - qJD(4) * t185 - t176 * t210);
t205 = t193 * qJ(3) + 0.2e1 * qJD(3) * t194 + t169;
t166 = -t193 * pkin(2) + t204;
t203 = -m(4) * t166 + t192 * mrSges(4,3) - t153;
t162 = t211 * t192 + t205;
t164 = t192 * pkin(2) - t205;
t202 = -m(4) * t164 + m(5) * t162 - t177 * mrSges(5,1) + t192 * mrSges(4,2) + t178 * mrSges(5,2) + t193 * mrSges(4,3) + t184 * t210 + t185 * t209;
t152 = t193 * mrSges(4,2) - t203;
t170 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t198 - Ifges(5,6) * t195) * t194;
t171 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t198 - Ifges(5,2) * t195) * t194;
t172 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t198 - Ifges(5,4) * t195) * t194;
t201 = -mrSges(3,2) * t169 - mrSges(4,3) * t164 - pkin(2) * t152 - pkin(6) * t153 - t195 * (-mrSges(5,1) * t162 + mrSges(5,3) * t160 + Ifges(5,4) * t178 + Ifges(5,2) * t177 + Ifges(5,6) * qJDD(4) + qJD(4) * t172 - t170 * t209) + t198 * (mrSges(5,2) * t162 - mrSges(5,3) * t159 + Ifges(5,1) * t178 + Ifges(5,4) * t177 + Ifges(5,5) * qJDD(4) - qJD(4) * t171 - t170 * t210) + qJ(3) * t202 + mrSges(4,2) * t166 + mrSges(3,1) * t168 + (Ifges(3,3) + Ifges(4,1)) * t193;
t1 = [mrSges(2,1) * t207 - mrSges(2,2) * t206 + t201 + Ifges(2,3) * qJDD(1) + pkin(1) * (t196 * (m(3) * t169 - t192 * mrSges(3,1) - t193 * mrSges(3,2) + t202) + t199 * (m(3) * t168 - t192 * mrSges(3,2) + (mrSges(3,1) - mrSges(4,2)) * t193 + t203)); t201; t152; mrSges(5,1) * t159 - mrSges(5,2) * t160 + Ifges(5,5) * t178 + Ifges(5,6) * t177 + Ifges(5,3) * qJDD(4) + (t171 * t198 + t172 * t195) * t194;];
tauJ = t1;
