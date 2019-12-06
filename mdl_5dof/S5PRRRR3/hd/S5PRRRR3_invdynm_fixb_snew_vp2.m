% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:10
% EndTime: 2019-12-05 17:06:13
% DurationCPUTime: 2.72s
% Computational Cost: add. (48826->170), mult. (61359->220), div. (0->0), fcn. (40054->10), ass. (0->82)
t165 = sin(pkin(9));
t166 = cos(pkin(9));
t151 = t165 * g(1) - t166 * g(2);
t152 = -t166 * g(1) - t165 * g(2);
t170 = sin(qJ(2));
t174 = cos(qJ(2));
t135 = t174 * t151 - t170 * t152;
t132 = qJDD(2) * pkin(2) + t135;
t136 = t170 * t151 + t174 * t152;
t175 = qJD(2) ^ 2;
t133 = -t175 * pkin(2) + t136;
t169 = sin(qJ(3));
t173 = cos(qJ(3));
t127 = t173 * t132 - t169 * t133;
t161 = qJDD(2) + qJDD(3);
t124 = t161 * pkin(3) + t127;
t128 = t169 * t132 + t173 * t133;
t162 = qJD(2) + qJD(3);
t160 = t162 ^ 2;
t125 = -t160 * pkin(3) + t128;
t168 = sin(qJ(4));
t172 = cos(qJ(4));
t121 = t168 * t124 + t172 * t125;
t157 = qJD(4) + t162;
t155 = t157 ^ 2;
t156 = qJDD(4) + t161;
t118 = -t155 * pkin(4) + t156 * pkin(8) + t121;
t164 = -g(3) + qJDD(1);
t167 = sin(qJ(5));
t171 = cos(qJ(5));
t115 = -t167 * t118 + t171 * t164;
t116 = t171 * t118 + t167 * t164;
t138 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t167 + Ifges(6,2) * t171) * t157;
t139 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t167 + Ifges(6,4) * t171) * t157;
t189 = qJD(5) * t157;
t143 = t167 * t156 + t171 * t189;
t144 = t171 * t156 - t167 * t189;
t192 = mrSges(6,1) * t115 - mrSges(6,2) * t116 + Ifges(6,5) * t143 + Ifges(6,6) * t144 + Ifges(6,3) * qJDD(5) + (t138 * t167 - t139 * t171) * t157;
t142 = (-mrSges(6,1) * t171 + mrSges(6,2) * t167) * t157;
t190 = t157 * t171;
t146 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t190;
t191 = t157 * t167;
t113 = m(6) * t115 + qJDD(5) * mrSges(6,1) - t143 * mrSges(6,3) + qJD(5) * t146 - t142 * t191;
t145 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t191;
t114 = m(6) * t116 - qJDD(5) * mrSges(6,2) + t144 * mrSges(6,3) - qJD(5) * t145 + t142 * t190;
t184 = -t167 * t113 + t171 * t114;
t100 = m(5) * t121 - t155 * mrSges(5,1) - t156 * mrSges(5,2) + t184;
t120 = t172 * t124 - t168 * t125;
t117 = -t156 * pkin(4) - t155 * pkin(8) - t120;
t180 = -m(6) * t117 + t144 * mrSges(6,1) - t143 * mrSges(6,2) - t145 * t191 + t146 * t190;
t108 = m(5) * t120 + t156 * mrSges(5,1) - t155 * mrSges(5,2) + t180;
t97 = t168 * t100 + t172 * t108;
t93 = m(4) * t127 + t161 * mrSges(4,1) - t160 * mrSges(4,2) + t97;
t185 = t172 * t100 - t168 * t108;
t94 = m(4) * t128 - t160 * mrSges(4,1) - t161 * mrSges(4,2) + t185;
t88 = t169 * t94 + t173 * t93;
t85 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t175 * mrSges(3,2) + t88;
t187 = -t169 * t93 + t173 * t94;
t86 = m(3) * t136 - t175 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t187;
t79 = t170 * t86 + t174 * t85;
t102 = t171 * t113 + t167 * t114;
t188 = m(5) * t164 + t102;
t186 = -t170 * t85 + t174 * t86;
t183 = m(4) * t164 + t188;
t137 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t167 + Ifges(6,6) * t171) * t157;
t105 = -mrSges(6,1) * t117 + mrSges(6,3) * t116 + Ifges(6,4) * t143 + Ifges(6,2) * t144 + Ifges(6,6) * qJDD(5) + qJD(5) * t139 - t137 * t191;
t106 = mrSges(6,2) * t117 - mrSges(6,3) * t115 + Ifges(6,1) * t143 + Ifges(6,4) * t144 + Ifges(6,5) * qJDD(5) - qJD(5) * t138 + t137 * t190;
t181 = mrSges(5,1) * t120 - mrSges(5,2) * t121 + Ifges(5,3) * t156 + pkin(4) * t180 + pkin(8) * t184 + t171 * t105 + t167 * t106;
t178 = mrSges(4,1) * t127 - mrSges(4,2) * t128 + Ifges(4,3) * t161 + pkin(3) * t97 + t181;
t177 = mrSges(3,1) * t135 - mrSges(3,2) * t136 + Ifges(3,3) * qJDD(2) + pkin(2) * t88 + t178;
t176 = mrSges(2,1) * t151 - mrSges(2,2) * t152 + pkin(1) * t79 + t177;
t95 = -mrSges(5,1) * t164 + mrSges(5,3) * t121 + t155 * Ifges(5,5) + Ifges(5,6) * t156 - pkin(4) * t102 - t192;
t89 = mrSges(5,2) * t164 - mrSges(5,3) * t120 + Ifges(5,5) * t156 - t155 * Ifges(5,6) - pkin(8) * t102 - t167 * t105 + t171 * t106;
t81 = mrSges(4,2) * t164 - mrSges(4,3) * t127 + Ifges(4,5) * t161 - t160 * Ifges(4,6) - pkin(7) * t97 - t168 * t95 + t172 * t89;
t80 = -mrSges(4,1) * t164 + mrSges(4,3) * t128 + t160 * Ifges(4,5) + Ifges(4,6) * t161 - pkin(3) * t188 + pkin(7) * t185 + t168 * t89 + t172 * t95;
t77 = m(2) * t152 + t186;
t76 = m(2) * t151 + t79;
t75 = mrSges(3,2) * t164 - mrSges(3,3) * t135 + Ifges(3,5) * qJDD(2) - t175 * Ifges(3,6) - pkin(6) * t88 - t169 * t80 + t173 * t81;
t74 = -mrSges(3,1) * t164 + mrSges(3,3) * t136 + t175 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t183 + pkin(6) * t187 + t169 * t81 + t173 * t80;
t73 = mrSges(2,2) * t164 - mrSges(2,3) * t151 - pkin(5) * t79 - t170 * t74 + t174 * t75;
t72 = -mrSges(2,1) * t164 + mrSges(2,3) * t152 + t170 * t75 + t174 * t74 - pkin(1) * (m(3) * t164 + t183) + pkin(5) * t186;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t166 * t73 - t165 * t72 - qJ(1) * (t165 * t77 + t166 * t76), t73, t75, t81, t89, t106; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t165 * t73 + t166 * t72 + qJ(1) * (-t165 * t76 + t166 * t77), t72, t74, t80, t95, t105; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t176, t176, t177, t178, t181, t192;];
m_new = t1;
