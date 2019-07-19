% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:33
% EndTime: 2019-07-18 18:16:33
% DurationCPUTime: 0.17s
% Computational Cost: add. (972->56), mult. (1103->64), div. (0->0), fcn. (462->6), ass. (0->35)
t181 = -pkin(2) - pkin(3);
t170 = sin(qJ(1));
t173 = cos(qJ(1));
t180 = g(1) * t170 - g(2) * t173;
t152 = qJDD(1) * pkin(1) + t180;
t178 = -g(1) * t173 - g(2) * t170;
t153 = -qJD(1) ^ 2 * pkin(1) + t178;
t169 = sin(qJ(2));
t172 = cos(qJ(2));
t148 = t152 * t169 + t153 * t172;
t167 = qJD(1) + qJD(2);
t166 = qJDD(1) + qJDD(2);
t147 = t172 * t152 - t153 * t169;
t179 = qJ(3) * t166 + 0.2e1 * qJD(3) * t167 + t148;
t165 = t167 ^ 2;
t140 = t165 * t181 + t179;
t175 = -t165 * qJ(3) + qJDD(3) - t147;
t143 = t166 * t181 + t175;
t168 = sin(qJ(4));
t171 = cos(qJ(4));
t138 = -t140 * t168 + t143 * t171;
t163 = qJD(4) - t167;
t160 = t163 ^ 2;
t161 = qJDD(4) - t166;
t135 = m(5) * t138 + mrSges(5,1) * t161 - mrSges(5,2) * t160;
t139 = t140 * t171 + t143 * t168;
t136 = m(5) * t139 - mrSges(5,1) * t160 - mrSges(5,2) * t161;
t133 = t135 * t171 + t136 * t168;
t144 = -pkin(2) * t165 + t179;
t177 = m(4) * t144 + mrSges(4,3) * t166 - t168 * t135 + t136 * t171;
t176 = mrSges(5,1) * t138 - mrSges(5,2) * t139 + Ifges(5,3) * t161;
t145 = -pkin(2) * t166 + t175;
t132 = m(4) * t145 - mrSges(4,1) * t166 - mrSges(4,3) * t165 + t133;
t174 = -mrSges(4,1) * t145 - mrSges(3,2) * t148 - pkin(3) * t133 + qJ(3) * (-mrSges(4,1) * t165 + t177) - pkin(2) * t132 + mrSges(4,3) * t144 + mrSges(3,1) * t147 - t176 + (Ifges(3,3) + Ifges(4,2)) * t166;
t1 = [t174 + pkin(1) * (t169 * (m(3) * t148 - t166 * mrSges(3,2) + (-mrSges(3,1) - mrSges(4,1)) * t165 + t177) + t172 * (m(3) * t147 + mrSges(3,1) * t166 - mrSges(3,2) * t165 - t132)) + mrSges(2,1) * t180 - mrSges(2,2) * t178 + Ifges(2,3) * qJDD(1); t174; t132; t176;];
tauJ  = t1;
