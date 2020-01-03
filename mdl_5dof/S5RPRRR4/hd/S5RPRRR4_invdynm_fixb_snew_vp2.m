% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:09
% EndTime: 2020-01-03 11:52:13
% DurationCPUTime: 3.56s
% Computational Cost: add. (57834->181), mult. (74879->231), div. (0->0), fcn. (40054->10), ass. (0->84)
t175 = sin(qJ(1));
t179 = cos(qJ(1));
t155 = -t179 * g(2) - t175 * g(3);
t151 = qJDD(1) * pkin(1) + t155;
t154 = -t175 * g(2) + t179 * g(3);
t180 = qJD(1) ^ 2;
t152 = -t180 * pkin(1) + t154;
t170 = sin(pkin(9));
t171 = cos(pkin(9));
t136 = t171 * t151 - t170 * t152;
t133 = qJDD(1) * pkin(2) + t136;
t137 = t170 * t151 + t171 * t152;
t134 = -t180 * pkin(2) + t137;
t174 = sin(qJ(3));
t178 = cos(qJ(3));
t128 = t178 * t133 - t174 * t134;
t165 = qJDD(1) + qJDD(3);
t125 = t165 * pkin(3) + t128;
t129 = t174 * t133 + t178 * t134;
t166 = qJD(1) + qJD(3);
t164 = t166 ^ 2;
t126 = -t164 * pkin(3) + t129;
t173 = sin(qJ(4));
t177 = cos(qJ(4));
t122 = t173 * t125 + t177 * t126;
t160 = qJD(4) + t166;
t158 = t160 ^ 2;
t159 = qJDD(4) + t165;
t119 = -t158 * pkin(4) + t159 * pkin(8) + t122;
t169 = -g(1) + qJDD(2);
t172 = sin(qJ(5));
t176 = cos(qJ(5));
t116 = -t172 * t119 + t176 * t169;
t117 = t176 * t119 + t172 * t169;
t139 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t172 + Ifges(6,2) * t176) * t160;
t140 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t172 + Ifges(6,4) * t176) * t160;
t194 = qJD(5) * t160;
t144 = t172 * t159 + t176 * t194;
t145 = t176 * t159 - t172 * t194;
t197 = mrSges(6,1) * t116 - mrSges(6,2) * t117 + Ifges(6,5) * t144 + Ifges(6,6) * t145 + Ifges(6,3) * qJDD(5) + (t139 * t172 - t140 * t176) * t160;
t143 = (-mrSges(6,1) * t176 + mrSges(6,2) * t172) * t160;
t195 = t160 * t176;
t150 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t195;
t196 = t160 * t172;
t114 = m(6) * t116 + qJDD(5) * mrSges(6,1) - t144 * mrSges(6,3) + qJD(5) * t150 - t143 * t196;
t149 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t196;
t115 = m(6) * t117 - qJDD(5) * mrSges(6,2) + t145 * mrSges(6,3) - qJD(5) * t149 + t143 * t195;
t189 = -t172 * t114 + t176 * t115;
t101 = m(5) * t122 - t158 * mrSges(5,1) - t159 * mrSges(5,2) + t189;
t121 = t177 * t125 - t173 * t126;
t118 = -t159 * pkin(4) - t158 * pkin(8) - t121;
t185 = -m(6) * t118 + t145 * mrSges(6,1) - t144 * mrSges(6,2) - t149 * t196 + t150 * t195;
t109 = m(5) * t121 + t159 * mrSges(5,1) - t158 * mrSges(5,2) + t185;
t98 = t173 * t101 + t177 * t109;
t94 = m(4) * t128 + t165 * mrSges(4,1) - t164 * mrSges(4,2) + t98;
t190 = t177 * t101 - t173 * t109;
t95 = m(4) * t129 - t164 * mrSges(4,1) - t165 * mrSges(4,2) + t190;
t89 = t174 * t95 + t178 * t94;
t86 = m(3) * t136 + qJDD(1) * mrSges(3,1) - t180 * mrSges(3,2) + t89;
t191 = -t174 * t94 + t178 * t95;
t87 = m(3) * t137 - t180 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t191;
t80 = t170 * t87 + t171 * t86;
t103 = t176 * t114 + t172 * t115;
t193 = m(5) * t169 + t103;
t192 = -t170 * t86 + t171 * t87;
t188 = m(4) * t169 + t193;
t138 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t172 + Ifges(6,6) * t176) * t160;
t106 = -mrSges(6,1) * t118 + mrSges(6,3) * t117 + Ifges(6,4) * t144 + Ifges(6,2) * t145 + Ifges(6,6) * qJDD(5) + qJD(5) * t140 - t138 * t196;
t107 = mrSges(6,2) * t118 - mrSges(6,3) * t116 + Ifges(6,1) * t144 + Ifges(6,4) * t145 + Ifges(6,5) * qJDD(5) - qJD(5) * t139 + t138 * t195;
t186 = mrSges(5,1) * t121 - mrSges(5,2) * t122 + Ifges(5,3) * t159 + pkin(4) * t185 + pkin(8) * t189 + t176 * t106 + t172 * t107;
t183 = mrSges(4,1) * t128 - mrSges(4,2) * t129 + Ifges(4,3) * t165 + pkin(3) * t98 + t186;
t182 = mrSges(3,1) * t136 - mrSges(3,2) * t137 + Ifges(3,3) * qJDD(1) + pkin(2) * t89 + t183;
t181 = mrSges(2,1) * t155 - mrSges(2,2) * t154 + Ifges(2,3) * qJDD(1) + pkin(1) * t80 + t182;
t96 = -mrSges(5,1) * t169 + mrSges(5,3) * t122 + t158 * Ifges(5,5) + Ifges(5,6) * t159 - pkin(4) * t103 - t197;
t90 = mrSges(5,2) * t169 - mrSges(5,3) * t121 + Ifges(5,5) * t159 - t158 * Ifges(5,6) - pkin(8) * t103 - t172 * t106 + t176 * t107;
t82 = mrSges(4,2) * t169 - mrSges(4,3) * t128 + Ifges(4,5) * t165 - t164 * Ifges(4,6) - pkin(7) * t98 - t173 * t96 + t177 * t90;
t81 = -mrSges(4,1) * t169 + mrSges(4,3) * t129 + t164 * Ifges(4,5) + Ifges(4,6) * t165 - pkin(3) * t193 + pkin(7) * t190 + t173 * t90 + t177 * t96;
t78 = m(2) * t155 + qJDD(1) * mrSges(2,1) - t180 * mrSges(2,2) + t80;
t77 = m(2) * t154 - t180 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t192;
t76 = mrSges(3,2) * t169 - mrSges(3,3) * t136 + Ifges(3,5) * qJDD(1) - t180 * Ifges(3,6) - pkin(6) * t89 - t174 * t81 + t178 * t82;
t75 = -mrSges(3,1) * t169 + mrSges(3,3) * t137 + t180 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t188 + pkin(6) * t191 + t174 * t82 + t178 * t81;
t74 = -mrSges(2,2) * g(1) - mrSges(2,3) * t155 + Ifges(2,5) * qJDD(1) - t180 * Ifges(2,6) - qJ(2) * t80 - t170 * t75 + t171 * t76;
t73 = Ifges(2,6) * qJDD(1) + t180 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t154 + t170 * t76 + t171 * t75 - pkin(1) * (m(3) * t169 + t188) + qJ(2) * t192;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t181, t74, t76, t82, t90, t107; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t175 * t74 + t179 * t73 - pkin(5) * (t175 * t78 - t179 * t77), t73, t75, t81, t96, t106; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t179 * t74 + t175 * t73 + pkin(5) * (t175 * t77 + t179 * t78), t181, t182, t183, t186, t197;];
m_new = t1;
