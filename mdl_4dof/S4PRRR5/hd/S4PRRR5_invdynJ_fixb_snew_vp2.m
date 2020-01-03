% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PRRR5
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PRRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR5_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:38
% EndTime: 2019-12-31 16:33:38
% DurationCPUTime: 0.19s
% Computational Cost: add. (840->85), mult. (1055->119), div. (0->0), fcn. (598->8), ass. (0->43)
t192 = qJD(2) + qJD(3);
t196 = sin(qJ(4));
t209 = t192 * t196;
t199 = cos(qJ(4));
t208 = t192 * t199;
t194 = sin(pkin(7));
t195 = cos(pkin(7));
t186 = -t195 * g(1) - t194 * g(2);
t193 = -g(3) + qJDD(1);
t198 = sin(qJ(2));
t201 = cos(qJ(2));
t172 = -t198 * t186 + t201 * t193;
t170 = qJDD(2) * pkin(2) + t172;
t173 = t201 * t186 + t198 * t193;
t202 = qJD(2) ^ 2;
t171 = -t202 * pkin(2) + t173;
t197 = sin(qJ(3));
t200 = cos(qJ(3));
t167 = t197 * t170 + t200 * t171;
t190 = t192 ^ 2;
t191 = qJDD(2) + qJDD(3);
t164 = -t190 * pkin(3) + t191 * pkin(6) + t167;
t185 = -t194 * g(1) + t195 * g(2);
t161 = -t196 * t164 + t199 * t185;
t162 = t199 * t164 + t196 * t185;
t179 = (-mrSges(5,1) * t199 + mrSges(5,2) * t196) * t192;
t206 = qJD(4) * t192;
t180 = t196 * t191 + t199 * t206;
t181 = t199 * t191 - t196 * t206;
t183 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t209;
t184 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t208;
t205 = -t196 * (m(5) * t161 + qJDD(4) * mrSges(5,1) - t180 * mrSges(5,3) + qJD(4) * t184 - t179 * t209) + t199 * (m(5) * t162 - qJDD(4) * mrSges(5,2) + t181 * mrSges(5,3) - qJD(4) * t183 + t179 * t208);
t152 = m(4) * t167 - t190 * mrSges(4,1) - t191 * mrSges(4,2) + t205;
t166 = t200 * t170 - t197 * t171;
t163 = -t191 * pkin(3) - t190 * pkin(6) - t166;
t203 = -m(5) * t163 + t181 * mrSges(5,1) - t180 * mrSges(5,2) - t183 * t209 + t184 * t208;
t157 = m(4) * t166 + t191 * mrSges(4,1) - t190 * mrSges(4,2) + t203;
t207 = t197 * t152 + t200 * t157;
t174 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t196 + Ifges(5,6) * t199) * t192;
t175 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t196 + Ifges(5,2) * t199) * t192;
t176 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t196 + Ifges(5,4) * t199) * t192;
t204 = -mrSges(4,2) * t167 + pkin(6) * t205 + t196 * (mrSges(5,2) * t163 - mrSges(5,3) * t161 + Ifges(5,1) * t180 + Ifges(5,4) * t181 + Ifges(5,5) * qJDD(4) - qJD(4) * t175 + t174 * t208) + t199 * (-mrSges(5,1) * t163 + mrSges(5,3) * t162 + Ifges(5,4) * t180 + Ifges(5,2) * t181 + Ifges(5,6) * qJDD(4) + qJD(4) * t176 - t174 * t209) + pkin(3) * t203 + mrSges(4,1) * t166 + Ifges(4,3) * t191;
t1 = [m(2) * t193 + t198 * (m(3) * t173 - t202 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t200 * t152 - t197 * t157) + t201 * (m(3) * t172 + qJDD(2) * mrSges(3,1) - t202 * mrSges(3,2) + t207); mrSges(3,1) * t172 - mrSges(3,2) * t173 + Ifges(3,3) * qJDD(2) + pkin(2) * t207 + t204; t204; mrSges(5,1) * t161 - mrSges(5,2) * t162 + Ifges(5,5) * t180 + Ifges(5,6) * t181 + Ifges(5,3) * qJDD(4) + (t175 * t196 - t176 * t199) * t192;];
tauJ = t1;
