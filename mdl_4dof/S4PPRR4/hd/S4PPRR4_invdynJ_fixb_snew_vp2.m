% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4PPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:30
% EndTime: 2019-12-31 16:18:31
% DurationCPUTime: 0.16s
% Computational Cost: add. (388->73), mult. (641->105), div. (0->0), fcn. (406->8), ass. (0->37)
t165 = sin(pkin(6));
t167 = cos(pkin(6));
t160 = -t167 * g(1) - t165 * g(2);
t163 = -g(3) + qJDD(1);
t164 = sin(pkin(7));
t166 = cos(pkin(7));
t149 = -t164 * t160 + t166 * t163;
t150 = t166 * t160 + t164 * t163;
t169 = sin(qJ(3));
t171 = cos(qJ(3));
t146 = t169 * t149 + t171 * t150;
t168 = sin(qJ(4));
t177 = qJD(3) * t168;
t170 = cos(qJ(4));
t176 = qJD(3) * t170;
t175 = qJD(3) * qJD(4);
t172 = qJD(3) ^ 2;
t144 = -t172 * pkin(3) + qJDD(3) * pkin(5) + t146;
t159 = -t165 * g(1) + t167 * g(2) + qJDD(2);
t141 = -t168 * t144 + t170 * t159;
t156 = (-t170 * mrSges(5,1) + t168 * mrSges(5,2)) * qJD(3);
t157 = t168 * qJDD(3) + t170 * t175;
t162 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t176;
t139 = m(5) * t141 + qJDD(4) * mrSges(5,1) - t157 * mrSges(5,3) + qJD(4) * t162 - t156 * t177;
t142 = t170 * t144 + t168 * t159;
t158 = t170 * qJDD(3) - t168 * t175;
t161 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t177;
t140 = m(5) * t142 - qJDD(4) * mrSges(5,2) + t158 * mrSges(5,3) - qJD(4) * t161 + t156 * t176;
t174 = -t168 * t139 + t170 * t140;
t145 = t171 * t149 - t169 * t150;
t143 = -qJDD(3) * pkin(3) - t172 * pkin(5) - t145;
t173 = -m(5) * t143 + t158 * mrSges(5,1) - t157 * mrSges(5,2) - t161 * t177 + t162 * t176;
t153 = Ifges(5,5) * qJD(4) + (t168 * Ifges(5,1) + t170 * Ifges(5,4)) * qJD(3);
t152 = Ifges(5,6) * qJD(4) + (t168 * Ifges(5,4) + t170 * Ifges(5,2)) * qJD(3);
t137 = m(4) * t145 + qJDD(3) * mrSges(4,1) - t172 * mrSges(4,2) + t173;
t136 = m(4) * t146 - t172 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t174;
t1 = [m(2) * t163 + t164 * (m(3) * t150 + t171 * t136 - t169 * t137) + t166 * (m(3) * t149 + t169 * t136 + t171 * t137); t170 * t139 + t168 * t140 + (m(3) + m(4)) * t159; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t145 - mrSges(4,2) * t146 + t168 * (mrSges(5,2) * t143 - mrSges(5,3) * t141 + Ifges(5,1) * t157 + Ifges(5,4) * t158 + Ifges(5,5) * qJDD(4) - qJD(4) * t152) + t170 * (-mrSges(5,1) * t143 + mrSges(5,3) * t142 + Ifges(5,4) * t157 + Ifges(5,2) * t158 + Ifges(5,6) * qJDD(4) + qJD(4) * t153) + pkin(3) * t173 + pkin(5) * t174; mrSges(5,1) * t141 - mrSges(5,2) * t142 + Ifges(5,5) * t157 + Ifges(5,6) * t158 + Ifges(5,3) * qJDD(4) + (t168 * t152 - t170 * t153) * qJD(3);];
tauJ = t1;
