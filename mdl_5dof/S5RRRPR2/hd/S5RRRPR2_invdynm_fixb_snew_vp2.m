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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:40:49
% EndTime: 2019-12-05 18:40:52
% DurationCPUTime: 3.14s
% Computational Cost: add. (64642->182), mult. (74879->231), div. (0->0), fcn. (40054->10), ass. (0->83)
t176 = sin(qJ(1));
t180 = cos(qJ(1));
t155 = t180 * g(2) + t176 * g(3);
t151 = qJDD(1) * pkin(1) + t155;
t154 = t176 * g(2) - t180 * g(3);
t181 = qJD(1) ^ 2;
t152 = -t181 * pkin(1) + t154;
t175 = sin(qJ(2));
t179 = cos(qJ(2));
t136 = t179 * t151 - t175 * t152;
t167 = qJDD(1) + qJDD(2);
t133 = t167 * pkin(2) + t136;
t137 = t175 * t151 + t179 * t152;
t168 = qJD(1) + qJD(2);
t166 = t168 ^ 2;
t134 = -t166 * pkin(2) + t137;
t174 = sin(qJ(3));
t178 = cos(qJ(3));
t128 = t178 * t133 - t174 * t134;
t160 = qJDD(3) + t167;
t125 = t160 * pkin(3) + t128;
t129 = t174 * t133 + t178 * t134;
t161 = qJD(3) + t168;
t159 = t161 ^ 2;
t126 = -t159 * pkin(3) + t129;
t171 = sin(pkin(9));
t172 = cos(pkin(9));
t122 = t171 * t125 + t172 * t126;
t119 = -t159 * pkin(4) + t160 * pkin(8) + t122;
t170 = -g(1) + qJDD(4);
t173 = sin(qJ(5));
t177 = cos(qJ(5));
t116 = -t173 * t119 + t177 * t170;
t117 = t177 * t119 + t173 * t170;
t139 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t173 + Ifges(6,2) * t177) * t161;
t140 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t173 + Ifges(6,4) * t177) * t161;
t194 = qJD(5) * t161;
t144 = t173 * t160 + t177 * t194;
t145 = t177 * t160 - t173 * t194;
t197 = mrSges(6,1) * t116 - mrSges(6,2) * t117 + Ifges(6,5) * t144 + Ifges(6,6) * t145 + Ifges(6,3) * qJDD(5) + (t139 * t173 - t140 * t177) * t161;
t143 = (-mrSges(6,1) * t177 + mrSges(6,2) * t173) * t161;
t195 = t161 * t177;
t150 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t195;
t196 = t161 * t173;
t114 = m(6) * t116 + qJDD(5) * mrSges(6,1) - t144 * mrSges(6,3) + qJD(5) * t150 - t143 * t196;
t149 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t196;
t115 = m(6) * t117 - qJDD(5) * mrSges(6,2) + t145 * mrSges(6,3) - qJD(5) * t149 + t143 * t195;
t189 = -t173 * t114 + t177 * t115;
t101 = m(5) * t122 - t159 * mrSges(5,1) - t160 * mrSges(5,2) + t189;
t121 = t172 * t125 - t171 * t126;
t118 = -t160 * pkin(4) - t159 * pkin(8) - t121;
t186 = -m(6) * t118 + t145 * mrSges(6,1) - t144 * mrSges(6,2) - t149 * t196 + t150 * t195;
t109 = m(5) * t121 + t160 * mrSges(5,1) - t159 * mrSges(5,2) + t186;
t98 = t171 * t101 + t172 * t109;
t94 = m(4) * t128 + t160 * mrSges(4,1) - t159 * mrSges(4,2) + t98;
t190 = t172 * t101 - t171 * t109;
t95 = m(4) * t129 - t159 * mrSges(4,1) - t160 * mrSges(4,2) + t190;
t89 = t174 * t95 + t178 * t94;
t86 = m(3) * t136 + t167 * mrSges(3,1) - t166 * mrSges(3,2) + t89;
t192 = -t174 * t94 + t178 * t95;
t87 = m(3) * t137 - t166 * mrSges(3,1) - t167 * mrSges(3,2) + t192;
t80 = t175 * t87 + t179 * t86;
t103 = t177 * t114 + t173 * t115;
t193 = m(5) * t170 + t103;
t191 = -t175 * t86 + t179 * t87;
t138 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t173 + Ifges(6,6) * t177) * t161;
t106 = -mrSges(6,1) * t118 + mrSges(6,3) * t117 + Ifges(6,4) * t144 + Ifges(6,2) * t145 + Ifges(6,6) * qJDD(5) + qJD(5) * t140 - t138 * t196;
t107 = mrSges(6,2) * t118 - mrSges(6,3) * t116 + Ifges(6,1) * t144 + Ifges(6,4) * t145 + Ifges(6,5) * qJDD(5) - qJD(5) * t139 + t138 * t195;
t187 = mrSges(5,1) * t121 - mrSges(5,2) * t122 + Ifges(5,3) * t160 + pkin(4) * t186 + pkin(8) * t189 + t177 * t106 + t173 * t107;
t184 = mrSges(4,1) * t128 - mrSges(4,2) * t129 + Ifges(4,3) * t160 + pkin(3) * t98 + t187;
t183 = mrSges(3,1) * t136 - mrSges(3,2) * t137 + Ifges(3,3) * t167 + pkin(2) * t89 + t184;
t182 = mrSges(2,1) * t155 - mrSges(2,2) * t154 + Ifges(2,3) * qJDD(1) + pkin(1) * t80 + t183;
t96 = -mrSges(5,1) * t170 + mrSges(5,3) * t122 + t159 * Ifges(5,5) + Ifges(5,6) * t160 - pkin(4) * t103 - t197;
t90 = mrSges(5,2) * t170 - mrSges(5,3) * t121 + Ifges(5,5) * t160 - t159 * Ifges(5,6) - pkin(8) * t103 - t173 * t106 + t177 * t107;
t82 = -mrSges(4,2) * g(1) - mrSges(4,3) * t128 + Ifges(4,5) * t160 - t159 * Ifges(4,6) - qJ(4) * t98 - t171 * t96 + t172 * t90;
t81 = mrSges(4,1) * g(1) + mrSges(4,3) * t129 + t159 * Ifges(4,5) + Ifges(4,6) * t160 - pkin(3) * t193 + qJ(4) * t190 + t171 * t90 + t172 * t96;
t78 = m(2) * t155 + qJDD(1) * mrSges(2,1) - t181 * mrSges(2,2) + t80;
t77 = m(2) * t154 - t181 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t191;
t76 = -mrSges(3,2) * g(1) - mrSges(3,3) * t136 + Ifges(3,5) * t167 - t166 * Ifges(3,6) - pkin(7) * t89 - t174 * t81 + t178 * t82;
t75 = Ifges(3,6) * t167 + t166 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t137 + t174 * t82 + t178 * t81 - pkin(2) * (-m(4) * g(1) + t193) + pkin(7) * t192;
t74 = -mrSges(2,2) * g(1) - mrSges(2,3) * t155 + Ifges(2,5) * qJDD(1) - t181 * Ifges(2,6) - pkin(6) * t80 - t175 * t75 + t179 * t76;
t73 = Ifges(2,6) * qJDD(1) + t181 * Ifges(2,5) + mrSges(2,3) * t154 + t175 * t76 + t179 * t75 - pkin(1) * t193 + pkin(6) * t191 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(1);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t182, t74, t76, t82, t90, t107; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t176 * t74 - t180 * t73 - pkin(5) * (-t176 * t78 + t180 * t77), t73, t75, t81, t96, t106; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t180 * t74 - t176 * t73 + pkin(5) * (-t176 * t77 - t180 * t78), t182, t183, t184, t187, t197;];
m_new = t1;
