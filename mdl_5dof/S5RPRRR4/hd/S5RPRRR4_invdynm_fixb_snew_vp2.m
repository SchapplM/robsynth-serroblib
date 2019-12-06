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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:14:26
% EndTime: 2019-12-05 18:14:29
% DurationCPUTime: 3.08s
% Computational Cost: add. (57834->181), mult. (74879->231), div. (0->0), fcn. (40054->10), ass. (0->84)
t179 = sin(qJ(1));
t183 = cos(qJ(1));
t157 = t183 * g(2) + t179 * g(3);
t153 = qJDD(1) * pkin(1) + t157;
t156 = t179 * g(2) - t183 * g(3);
t184 = qJD(1) ^ 2;
t154 = -t184 * pkin(1) + t156;
t174 = sin(pkin(9));
t175 = cos(pkin(9));
t138 = t175 * t153 - t174 * t154;
t135 = qJDD(1) * pkin(2) + t138;
t139 = t174 * t153 + t175 * t154;
t136 = -t184 * pkin(2) + t139;
t178 = sin(qJ(3));
t182 = cos(qJ(3));
t130 = t182 * t135 - t178 * t136;
t169 = qJDD(1) + qJDD(3);
t127 = t169 * pkin(3) + t130;
t131 = t178 * t135 + t182 * t136;
t170 = qJD(1) + qJD(3);
t168 = t170 ^ 2;
t128 = -t168 * pkin(3) + t131;
t177 = sin(qJ(4));
t181 = cos(qJ(4));
t124 = t177 * t127 + t181 * t128;
t162 = qJD(4) + t170;
t160 = t162 ^ 2;
t161 = qJDD(4) + t169;
t121 = -t160 * pkin(4) + t161 * pkin(8) + t124;
t173 = -g(1) + qJDD(2);
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t118 = -t176 * t121 + t180 * t173;
t119 = t180 * t121 + t176 * t173;
t141 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t176 + Ifges(6,2) * t180) * t162;
t142 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t176 + Ifges(6,4) * t180) * t162;
t198 = qJD(5) * t162;
t146 = t176 * t161 + t180 * t198;
t147 = t180 * t161 - t176 * t198;
t201 = mrSges(6,1) * t118 - mrSges(6,2) * t119 + Ifges(6,5) * t146 + Ifges(6,6) * t147 + Ifges(6,3) * qJDD(5) + (t141 * t176 - t142 * t180) * t162;
t145 = (-mrSges(6,1) * t180 + mrSges(6,2) * t176) * t162;
t199 = t162 * t180;
t152 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t199;
t200 = t162 * t176;
t116 = m(6) * t118 + qJDD(5) * mrSges(6,1) - t146 * mrSges(6,3) + qJD(5) * t152 - t145 * t200;
t151 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t200;
t117 = m(6) * t119 - qJDD(5) * mrSges(6,2) + t147 * mrSges(6,3) - qJD(5) * t151 + t145 * t199;
t193 = -t176 * t116 + t180 * t117;
t103 = m(5) * t124 - t160 * mrSges(5,1) - t161 * mrSges(5,2) + t193;
t123 = t181 * t127 - t177 * t128;
t120 = -t161 * pkin(4) - t160 * pkin(8) - t123;
t189 = -m(6) * t120 + t147 * mrSges(6,1) - t146 * mrSges(6,2) - t151 * t200 + t152 * t199;
t111 = m(5) * t123 + t161 * mrSges(5,1) - t160 * mrSges(5,2) + t189;
t100 = t177 * t103 + t181 * t111;
t96 = m(4) * t130 + t169 * mrSges(4,1) - t168 * mrSges(4,2) + t100;
t194 = t181 * t103 - t177 * t111;
t97 = m(4) * t131 - t168 * mrSges(4,1) - t169 * mrSges(4,2) + t194;
t91 = t178 * t97 + t182 * t96;
t88 = m(3) * t138 + qJDD(1) * mrSges(3,1) - t184 * mrSges(3,2) + t91;
t195 = -t178 * t96 + t182 * t97;
t89 = m(3) * t139 - t184 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t195;
t82 = t174 * t89 + t175 * t88;
t105 = t180 * t116 + t176 * t117;
t197 = m(5) * t173 + t105;
t196 = -t174 * t88 + t175 * t89;
t192 = m(4) * t173 + t197;
t140 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t176 + Ifges(6,6) * t180) * t162;
t108 = -mrSges(6,1) * t120 + mrSges(6,3) * t119 + Ifges(6,4) * t146 + Ifges(6,2) * t147 + Ifges(6,6) * qJDD(5) + qJD(5) * t142 - t140 * t200;
t109 = mrSges(6,2) * t120 - mrSges(6,3) * t118 + Ifges(6,1) * t146 + Ifges(6,4) * t147 + Ifges(6,5) * qJDD(5) - qJD(5) * t141 + t140 * t199;
t190 = mrSges(5,1) * t123 - mrSges(5,2) * t124 + Ifges(5,3) * t161 + pkin(4) * t189 + pkin(8) * t193 + t180 * t108 + t176 * t109;
t187 = mrSges(4,1) * t130 - mrSges(4,2) * t131 + Ifges(4,3) * t169 + pkin(3) * t100 + t190;
t186 = mrSges(3,1) * t138 - mrSges(3,2) * t139 + Ifges(3,3) * qJDD(1) + pkin(2) * t91 + t187;
t185 = mrSges(2,1) * t157 - mrSges(2,2) * t156 + Ifges(2,3) * qJDD(1) + pkin(1) * t82 + t186;
t98 = -mrSges(5,1) * t173 + mrSges(5,3) * t124 + t160 * Ifges(5,5) + Ifges(5,6) * t161 - pkin(4) * t105 - t201;
t92 = mrSges(5,2) * t173 - mrSges(5,3) * t123 + Ifges(5,5) * t161 - t160 * Ifges(5,6) - pkin(8) * t105 - t176 * t108 + t180 * t109;
t84 = mrSges(4,2) * t173 - mrSges(4,3) * t130 + Ifges(4,5) * t169 - t168 * Ifges(4,6) - pkin(7) * t100 - t177 * t98 + t181 * t92;
t83 = -mrSges(4,1) * t173 + mrSges(4,3) * t131 + t168 * Ifges(4,5) + Ifges(4,6) * t169 - pkin(3) * t197 + pkin(7) * t194 + t177 * t92 + t181 * t98;
t80 = m(2) * t157 + qJDD(1) * mrSges(2,1) - t184 * mrSges(2,2) + t82;
t79 = m(2) * t156 - t184 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t196;
t78 = mrSges(3,2) * t173 - mrSges(3,3) * t138 + Ifges(3,5) * qJDD(1) - t184 * Ifges(3,6) - pkin(6) * t91 - t178 * t83 + t182 * t84;
t77 = -mrSges(3,1) * t173 + mrSges(3,3) * t139 + t184 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t192 + pkin(6) * t195 + t178 * t84 + t182 * t83;
t76 = -mrSges(2,2) * g(1) - mrSges(2,3) * t157 + Ifges(2,5) * qJDD(1) - t184 * Ifges(2,6) - qJ(2) * t82 - t174 * t77 + t175 * t78;
t75 = Ifges(2,6) * qJDD(1) + t184 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t156 + t174 * t78 + t175 * t77 - pkin(1) * (m(3) * t173 + t192) + qJ(2) * t196;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t185, t76, t78, t84, t92, t109; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t179 * t76 - t183 * t75 - pkin(5) * (-t179 * t80 + t183 * t79), t75, t77, t83, t98, t108; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t183 * t76 - t179 * t75 + pkin(5) * (-t179 * t79 - t183 * t80), t185, t186, t187, t190, t201;];
m_new = t1;
