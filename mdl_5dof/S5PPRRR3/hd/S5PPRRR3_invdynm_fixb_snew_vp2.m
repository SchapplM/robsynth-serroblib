% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:28
% EndTime: 2019-12-05 15:16:34
% DurationCPUTime: 2.96s
% Computational Cost: add. (37038->204), mult. (67265->265), div. (0->0), fcn. (44032->10), ass. (0->92)
t171 = sin(pkin(8));
t173 = cos(pkin(8));
t161 = -t173 * g(1) - t171 * g(2);
t170 = sin(pkin(9));
t172 = cos(pkin(9));
t196 = g(3) - qJDD(1);
t146 = t172 * t161 - t170 * t196;
t160 = t171 * g(1) - t173 * g(2);
t159 = qJDD(2) - t160;
t176 = sin(qJ(3));
t179 = cos(qJ(3));
t136 = t179 * t146 + t176 * t159;
t180 = qJD(3) ^ 2;
t131 = -t180 * pkin(3) + qJDD(3) * pkin(6) + t136;
t145 = t170 * t161 + t172 * t196;
t175 = sin(qJ(4));
t178 = cos(qJ(4));
t121 = -t175 * t131 + t178 * t145;
t193 = qJD(3) * qJD(4);
t192 = t178 * t193;
t156 = t175 * qJDD(3) + t192;
t118 = (-t156 + t192) * pkin(7) + (t175 * t178 * t180 + qJDD(4)) * pkin(4) + t121;
t122 = t178 * t131 + t175 * t145;
t157 = t178 * qJDD(3) - t175 * t193;
t195 = qJD(3) * t175;
t164 = qJD(4) * pkin(4) - pkin(7) * t195;
t169 = t178 ^ 2;
t119 = -t169 * t180 * pkin(4) + t157 * pkin(7) - qJD(4) * t164 + t122;
t174 = sin(qJ(5));
t177 = cos(qJ(5));
t116 = t177 * t118 - t174 * t119;
t150 = (-t174 * t175 + t177 * t178) * qJD(3);
t129 = t150 * qJD(5) + t177 * t156 + t174 * t157;
t151 = (t174 * t178 + t175 * t177) * qJD(3);
t138 = -t150 * mrSges(6,1) + t151 * mrSges(6,2);
t168 = qJD(4) + qJD(5);
t143 = -t168 * mrSges(6,2) + t150 * mrSges(6,3);
t167 = qJDD(4) + qJDD(5);
t113 = m(6) * t116 + t167 * mrSges(6,1) - t129 * mrSges(6,3) - t151 * t138 + t168 * t143;
t117 = t174 * t118 + t177 * t119;
t128 = -t151 * qJD(5) - t174 * t156 + t177 * t157;
t144 = t168 * mrSges(6,1) - t151 * mrSges(6,3);
t114 = m(6) * t117 - t167 * mrSges(6,2) + t128 * mrSges(6,3) + t150 * t138 - t168 * t144;
t105 = t177 * t113 + t174 * t114;
t148 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t175 + Ifges(5,2) * t178) * qJD(3);
t149 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t175 + Ifges(5,4) * t178) * qJD(3);
t133 = Ifges(6,4) * t151 + Ifges(6,2) * t150 + Ifges(6,6) * t168;
t134 = Ifges(6,1) * t151 + Ifges(6,4) * t150 + Ifges(6,5) * t168;
t184 = -mrSges(6,1) * t116 + mrSges(6,2) * t117 - Ifges(6,5) * t129 - Ifges(6,6) * t128 - Ifges(6,3) * t167 - t151 * t133 + t150 * t134;
t197 = mrSges(5,1) * t121 - mrSges(5,2) * t122 + Ifges(5,5) * t156 + Ifges(5,6) * t157 + Ifges(5,3) * qJDD(4) + pkin(4) * t105 + (t175 * t148 - t178 * t149) * qJD(3) - t184;
t194 = qJD(3) * t178;
t135 = -t176 * t146 + t179 * t159;
t188 = -qJDD(3) * pkin(3) - t135;
t130 = -t180 * pkin(6) + t188;
t162 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t195;
t163 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t194;
t120 = t164 * t195 - t157 * pkin(4) + (-pkin(7) * t169 - pkin(6)) * t180 + t188;
t185 = m(6) * t120 - t128 * mrSges(6,1) + t129 * mrSges(6,2) - t150 * t143 + t151 * t144;
t109 = -m(5) * t130 + t157 * mrSges(5,1) - t156 * mrSges(5,2) - t162 * t195 + t163 * t194 - t185;
t108 = m(4) * t135 + qJDD(3) * mrSges(4,1) - t180 * mrSges(4,2) + t109;
t155 = (-mrSges(5,1) * t178 + mrSges(5,2) * t175) * qJD(3);
t103 = m(5) * t121 + qJDD(4) * mrSges(5,1) - t156 * mrSges(5,3) + qJD(4) * t163 - t155 * t195 + t105;
t190 = -t174 * t113 + t177 * t114;
t104 = m(5) * t122 - qJDD(4) * mrSges(5,2) + t157 * mrSges(5,3) - qJD(4) * t162 + t155 * t194 + t190;
t101 = -t175 * t103 + t178 * t104;
t97 = m(4) * t136 - t180 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t101;
t94 = -t176 * t108 + t179 * t97;
t91 = m(3) * t146 + t94;
t100 = t178 * t103 + t175 * t104;
t98 = (-m(3) - m(4)) * t145 - t100;
t191 = -t170 * t98 + t172 * t91;
t93 = t179 * t108 + t176 * t97;
t186 = -m(3) * t159 - t93;
t132 = Ifges(6,5) * t151 + Ifges(6,6) * t150 + Ifges(6,3) * t168;
t106 = -mrSges(6,1) * t120 + mrSges(6,3) * t117 + Ifges(6,4) * t129 + Ifges(6,2) * t128 + Ifges(6,6) * t167 - t151 * t132 + t168 * t134;
t107 = mrSges(6,2) * t120 - mrSges(6,3) * t116 + Ifges(6,1) * t129 + Ifges(6,4) * t128 + Ifges(6,5) * t167 + t150 * t132 - t168 * t133;
t147 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t175 + Ifges(5,6) * t178) * qJD(3);
t88 = -mrSges(5,1) * t130 + mrSges(5,3) * t122 + Ifges(5,4) * t156 + Ifges(5,2) * t157 + Ifges(5,6) * qJDD(4) - pkin(4) * t185 + pkin(7) * t190 + qJD(4) * t149 + t177 * t106 + t174 * t107 - t147 * t195;
t95 = mrSges(5,2) * t130 - mrSges(5,3) * t121 + Ifges(5,1) * t156 + Ifges(5,4) * t157 + Ifges(5,5) * qJDD(4) - pkin(7) * t105 - qJD(4) * t148 - t174 * t106 + t177 * t107 + t147 * t194;
t83 = mrSges(4,2) * t145 - mrSges(4,3) * t135 + Ifges(4,5) * qJDD(3) - t180 * Ifges(4,6) - pkin(6) * t100 - t175 * t88 + t178 * t95;
t87 = -mrSges(4,1) * t145 + mrSges(4,3) * t136 + t180 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t100 - t197;
t80 = mrSges(3,2) * t159 + mrSges(3,3) * t145 - pkin(5) * t93 - t176 * t87 + t179 * t83;
t182 = mrSges(4,1) * t135 - mrSges(4,2) * t136 + Ifges(4,3) * qJDD(3) + pkin(3) * t109 + pkin(6) * t101 + t175 * t95 + t178 * t88;
t82 = -mrSges(3,1) * t159 + mrSges(3,3) * t146 - pkin(2) * t93 - t182;
t187 = mrSges(2,1) * t160 - mrSges(2,2) * t161 + pkin(1) * t186 + qJ(2) * t191 + t170 * t80 + t172 * t82;
t183 = mrSges(3,1) * t145 + mrSges(3,2) * t146 - pkin(2) * (-m(4) * t145 - t100) - pkin(5) * t94 - t176 * t83 - t179 * t87;
t90 = m(2) * t160 + t186;
t86 = t170 * t91 + t172 * t98;
t84 = m(2) * t161 + t191;
t78 = mrSges(2,1) * t196 + mrSges(2,3) * t161 - pkin(1) * t86 + t183;
t77 = -mrSges(2,2) * t196 - mrSges(2,3) * t160 - qJ(2) * t86 - t170 * t82 + t172 * t80;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t173 * t77 - t171 * t78 - qJ(1) * (t171 * t84 + t173 * t90), t77, t80, t83, t95, t107; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t171 * t77 + t173 * t78 + qJ(1) * (-t171 * t90 + t173 * t84), t78, t82, t87, t88, t106; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t187, t187, -t183, t182, t197, -t184;];
m_new = t1;
