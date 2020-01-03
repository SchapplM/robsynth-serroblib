% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:27
% EndTime: 2020-01-03 12:07:30
% DurationCPUTime: 3.07s
% Computational Cost: add. (64642->182), mult. (74879->231), div. (0->0), fcn. (40054->10), ass. (0->83)
t172 = sin(qJ(1));
t176 = cos(qJ(1));
t153 = -t176 * g(2) - t172 * g(3);
t149 = qJDD(1) * pkin(1) + t153;
t152 = -t172 * g(2) + t176 * g(3);
t177 = qJD(1) ^ 2;
t150 = -t177 * pkin(1) + t152;
t171 = sin(qJ(2));
t175 = cos(qJ(2));
t134 = t175 * t149 - t171 * t150;
t163 = qJDD(1) + qJDD(2);
t131 = t163 * pkin(2) + t134;
t135 = t171 * t149 + t175 * t150;
t164 = qJD(1) + qJD(2);
t162 = t164 ^ 2;
t132 = -t162 * pkin(2) + t135;
t170 = sin(qJ(3));
t174 = cos(qJ(3));
t126 = t174 * t131 - t170 * t132;
t158 = qJDD(3) + t163;
t123 = t158 * pkin(3) + t126;
t127 = t170 * t131 + t174 * t132;
t159 = qJD(3) + t164;
t157 = t159 ^ 2;
t124 = -t157 * pkin(3) + t127;
t167 = sin(pkin(9));
t168 = cos(pkin(9));
t120 = t167 * t123 + t168 * t124;
t117 = -t157 * pkin(4) + t158 * pkin(8) + t120;
t166 = -g(1) + qJDD(4);
t169 = sin(qJ(5));
t173 = cos(qJ(5));
t114 = -t169 * t117 + t173 * t166;
t115 = t173 * t117 + t169 * t166;
t137 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t169 + Ifges(6,2) * t173) * t159;
t138 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t169 + Ifges(6,4) * t173) * t159;
t190 = qJD(5) * t159;
t142 = t169 * t158 + t173 * t190;
t143 = t173 * t158 - t169 * t190;
t193 = mrSges(6,1) * t114 - mrSges(6,2) * t115 + Ifges(6,5) * t142 + Ifges(6,6) * t143 + Ifges(6,3) * qJDD(5) + (t137 * t169 - t138 * t173) * t159;
t119 = t168 * t123 - t167 * t124;
t116 = -t158 * pkin(4) - t157 * pkin(8) - t119;
t192 = t159 * t169;
t147 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t192;
t191 = t159 * t173;
t148 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t191;
t182 = -m(6) * t116 + t143 * mrSges(6,1) - t142 * mrSges(6,2) - t147 * t192 + t148 * t191;
t107 = m(5) * t119 + t158 * mrSges(5,1) - t157 * mrSges(5,2) + t182;
t141 = (-mrSges(6,1) * t173 + mrSges(6,2) * t169) * t159;
t112 = m(6) * t114 + qJDD(5) * mrSges(6,1) - t142 * mrSges(6,3) + qJD(5) * t148 - t141 * t192;
t113 = m(6) * t115 - qJDD(5) * mrSges(6,2) + t143 * mrSges(6,3) - qJD(5) * t147 + t141 * t191;
t185 = -t169 * t112 + t173 * t113;
t99 = m(5) * t120 - t157 * mrSges(5,1) - t158 * mrSges(5,2) + t185;
t96 = t168 * t107 + t167 * t99;
t92 = m(4) * t126 + t158 * mrSges(4,1) - t157 * mrSges(4,2) + t96;
t186 = -t167 * t107 + t168 * t99;
t93 = m(4) * t127 - t157 * mrSges(4,1) - t158 * mrSges(4,2) + t186;
t87 = t170 * t93 + t174 * t92;
t84 = m(3) * t134 + t163 * mrSges(3,1) - t162 * mrSges(3,2) + t87;
t188 = -t170 * t92 + t174 * t93;
t85 = m(3) * t135 - t162 * mrSges(3,1) - t163 * mrSges(3,2) + t188;
t78 = t171 * t85 + t175 * t84;
t101 = t173 * t112 + t169 * t113;
t189 = m(5) * t166 + t101;
t187 = -t171 * t84 + t175 * t85;
t136 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t169 + Ifges(6,6) * t173) * t159;
t104 = -mrSges(6,1) * t116 + mrSges(6,3) * t115 + Ifges(6,4) * t142 + Ifges(6,2) * t143 + Ifges(6,6) * qJDD(5) + qJD(5) * t138 - t136 * t192;
t105 = mrSges(6,2) * t116 - mrSges(6,3) * t114 + Ifges(6,1) * t142 + Ifges(6,4) * t143 + Ifges(6,5) * qJDD(5) - qJD(5) * t137 + t136 * t191;
t183 = mrSges(5,1) * t119 - mrSges(5,2) * t120 + Ifges(5,3) * t158 + pkin(4) * t182 + pkin(8) * t185 + t173 * t104 + t169 * t105;
t180 = mrSges(4,1) * t126 - mrSges(4,2) * t127 + Ifges(4,3) * t158 + pkin(3) * t96 + t183;
t179 = mrSges(3,1) * t134 - mrSges(3,2) * t135 + Ifges(3,3) * t163 + pkin(2) * t87 + t180;
t178 = mrSges(2,1) * t153 - mrSges(2,2) * t152 + Ifges(2,3) * qJDD(1) + pkin(1) * t78 + t179;
t94 = -mrSges(5,1) * t166 + mrSges(5,3) * t120 + t157 * Ifges(5,5) + Ifges(5,6) * t158 - pkin(4) * t101 - t193;
t88 = mrSges(5,2) * t166 - mrSges(5,3) * t119 + Ifges(5,5) * t158 - t157 * Ifges(5,6) - pkin(8) * t101 - t169 * t104 + t173 * t105;
t80 = -mrSges(4,2) * g(1) - mrSges(4,3) * t126 + Ifges(4,5) * t158 - t157 * Ifges(4,6) - qJ(4) * t96 - t167 * t94 + t168 * t88;
t79 = mrSges(4,1) * g(1) + mrSges(4,3) * t127 + t157 * Ifges(4,5) + Ifges(4,6) * t158 - pkin(3) * t189 + qJ(4) * t186 + t167 * t88 + t168 * t94;
t76 = m(2) * t153 + qJDD(1) * mrSges(2,1) - t177 * mrSges(2,2) + t78;
t75 = m(2) * t152 - t177 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t187;
t74 = -mrSges(3,2) * g(1) - mrSges(3,3) * t134 + Ifges(3,5) * t163 - t162 * Ifges(3,6) - pkin(7) * t87 - t170 * t79 + t174 * t80;
t73 = Ifges(3,6) * t163 + t162 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t135 + t170 * t80 + t174 * t79 - pkin(2) * (-m(4) * g(1) + t189) + pkin(7) * t188;
t72 = -mrSges(2,2) * g(1) - mrSges(2,3) * t153 + Ifges(2,5) * qJDD(1) - t177 * Ifges(2,6) - pkin(6) * t78 - t171 * t73 + t175 * t74;
t71 = Ifges(2,6) * qJDD(1) + t177 * Ifges(2,5) + mrSges(2,3) * t152 + t171 * t74 + t175 * t73 - pkin(1) * t189 + pkin(6) * t187 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(1);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t178, t72, t74, t80, t88, t105; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t172 * t72 + t176 * t71 - pkin(5) * (t172 * t76 - t176 * t75), t71, t73, t79, t94, t104; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t176 * t72 + t172 * t71 + pkin(5) * (t172 * t75 + t176 * t76), t178, t179, t180, t183, t193;];
m_new = t1;
