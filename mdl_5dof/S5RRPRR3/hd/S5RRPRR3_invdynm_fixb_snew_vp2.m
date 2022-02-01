% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:06
% EndTime: 2022-01-20 10:34:09
% DurationCPUTime: 3.51s
% Computational Cost: add. (62348->181), mult. (74879->231), div. (0->0), fcn. (40054->10), ass. (0->84)
t174 = sin(qJ(1));
t178 = cos(qJ(1));
t153 = t174 * g(1) - t178 * g(2);
t150 = qJDD(1) * pkin(1) + t153;
t154 = -t178 * g(1) - t174 * g(2);
t179 = qJD(1) ^ 2;
t151 = -t179 * pkin(1) + t154;
t173 = sin(qJ(2));
t177 = cos(qJ(2));
t135 = t177 * t150 - t173 * t151;
t165 = qJDD(1) + qJDD(2);
t132 = t165 * pkin(2) + t135;
t136 = t173 * t150 + t177 * t151;
t166 = qJD(1) + qJD(2);
t164 = t166 ^ 2;
t133 = -t164 * pkin(2) + t136;
t169 = sin(pkin(9));
t170 = cos(pkin(9));
t127 = t170 * t132 - t169 * t133;
t124 = t165 * pkin(3) + t127;
t128 = t169 * t132 + t170 * t133;
t125 = -t164 * pkin(3) + t128;
t172 = sin(qJ(4));
t176 = cos(qJ(4));
t121 = t172 * t124 + t176 * t125;
t160 = qJD(4) + t166;
t158 = t160 ^ 2;
t159 = qJDD(4) + t165;
t118 = -t158 * pkin(4) + t159 * pkin(8) + t121;
t168 = -g(3) + qJDD(3);
t171 = sin(qJ(5));
t175 = cos(qJ(5));
t115 = -t171 * t118 + t175 * t168;
t116 = t175 * t118 + t171 * t168;
t138 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t171 + Ifges(6,2) * t175) * t160;
t139 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t171 + Ifges(6,4) * t175) * t160;
t193 = qJD(5) * t160;
t143 = t171 * t159 + t175 * t193;
t144 = t175 * t159 - t171 * t193;
t196 = mrSges(6,1) * t115 - mrSges(6,2) * t116 + Ifges(6,5) * t143 + Ifges(6,6) * t144 + Ifges(6,3) * qJDD(5) + (t138 * t171 - t139 * t175) * t160;
t142 = (-mrSges(6,1) * t175 + mrSges(6,2) * t171) * t160;
t194 = t160 * t175;
t149 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t194;
t195 = t160 * t171;
t113 = m(6) * t115 + qJDD(5) * mrSges(6,1) - t143 * mrSges(6,3) + qJD(5) * t149 - t142 * t195;
t148 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t195;
t114 = m(6) * t116 - qJDD(5) * mrSges(6,2) + t144 * mrSges(6,3) - qJD(5) * t148 + t142 * t194;
t188 = -t171 * t113 + t175 * t114;
t100 = m(5) * t121 - t158 * mrSges(5,1) - t159 * mrSges(5,2) + t188;
t120 = t176 * t124 - t172 * t125;
t117 = -t159 * pkin(4) - t158 * pkin(8) - t120;
t184 = -m(6) * t117 + t144 * mrSges(6,1) - t143 * mrSges(6,2) - t148 * t195 + t149 * t194;
t108 = m(5) * t120 + t159 * mrSges(5,1) - t158 * mrSges(5,2) + t184;
t97 = t172 * t100 + t176 * t108;
t93 = m(4) * t127 + t165 * mrSges(4,1) - t164 * mrSges(4,2) + t97;
t189 = t176 * t100 - t172 * t108;
t94 = m(4) * t128 - t164 * mrSges(4,1) - t165 * mrSges(4,2) + t189;
t88 = t169 * t94 + t170 * t93;
t85 = m(3) * t135 + t165 * mrSges(3,1) - t164 * mrSges(3,2) + t88;
t191 = -t169 * t93 + t170 * t94;
t86 = m(3) * t136 - t164 * mrSges(3,1) - t165 * mrSges(3,2) + t191;
t79 = t173 * t86 + t177 * t85;
t102 = t175 * t113 + t171 * t114;
t192 = m(5) * t168 + t102;
t190 = -t173 * t85 + t177 * t86;
t187 = m(4) * t168 + t192;
t137 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t171 + Ifges(6,6) * t175) * t160;
t105 = -mrSges(6,1) * t117 + mrSges(6,3) * t116 + Ifges(6,4) * t143 + Ifges(6,2) * t144 + Ifges(6,6) * qJDD(5) + qJD(5) * t139 - t137 * t195;
t106 = mrSges(6,2) * t117 - mrSges(6,3) * t115 + Ifges(6,1) * t143 + Ifges(6,4) * t144 + Ifges(6,5) * qJDD(5) - qJD(5) * t138 + t137 * t194;
t185 = mrSges(5,1) * t120 - mrSges(5,2) * t121 + Ifges(5,3) * t159 + pkin(4) * t184 + pkin(8) * t188 + t175 * t105 + t171 * t106;
t182 = mrSges(4,1) * t127 - mrSges(4,2) * t128 + Ifges(4,3) * t165 + pkin(3) * t97 + t185;
t181 = mrSges(3,1) * t135 - mrSges(3,2) * t136 + Ifges(3,3) * t165 + pkin(2) * t88 + t182;
t180 = mrSges(2,1) * t153 - mrSges(2,2) * t154 + Ifges(2,3) * qJDD(1) + pkin(1) * t79 + t181;
t95 = -mrSges(5,1) * t168 + mrSges(5,3) * t121 + t158 * Ifges(5,5) + Ifges(5,6) * t159 - pkin(4) * t102 - t196;
t89 = mrSges(5,2) * t168 - mrSges(5,3) * t120 + Ifges(5,5) * t159 - t158 * Ifges(5,6) - pkin(8) * t102 - t171 * t105 + t175 * t106;
t81 = mrSges(4,2) * t168 - mrSges(4,3) * t127 + Ifges(4,5) * t165 - t164 * Ifges(4,6) - pkin(7) * t97 - t172 * t95 + t176 * t89;
t80 = -mrSges(4,1) * t168 + mrSges(4,3) * t128 + t164 * Ifges(4,5) + Ifges(4,6) * t165 - pkin(3) * t192 + pkin(7) * t189 + t172 * t89 + t176 * t95;
t77 = m(2) * t154 - t179 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t190;
t76 = m(2) * t153 + qJDD(1) * mrSges(2,1) - t179 * mrSges(2,2) + t79;
t75 = -mrSges(3,2) * g(3) - mrSges(3,3) * t135 + Ifges(3,5) * t165 - t164 * Ifges(3,6) - qJ(3) * t88 - t169 * t80 + t170 * t81;
t74 = mrSges(3,1) * g(3) + mrSges(3,3) * t136 + t164 * Ifges(3,5) + Ifges(3,6) * t165 - pkin(2) * t187 + qJ(3) * t191 + t169 * t81 + t170 * t80;
t73 = -mrSges(2,2) * g(3) - mrSges(2,3) * t153 + Ifges(2,5) * qJDD(1) - t179 * Ifges(2,6) - pkin(6) * t79 - t173 * t74 + t177 * t75;
t72 = Ifges(2,6) * qJDD(1) + t179 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t154 + t173 * t75 + t177 * t74 - pkin(1) * (-m(3) * g(3) + t187) + pkin(6) * t190;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t178 * t73 - t174 * t72 - pkin(5) * (t174 * t77 + t178 * t76), t73, t75, t81, t89, t106; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t174 * t73 + t178 * t72 + pkin(5) * (-t174 * t76 + t178 * t77), t72, t74, t80, t95, t105; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t180, t180, t181, t182, t185, t196;];
m_new = t1;
