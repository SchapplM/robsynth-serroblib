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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:30:23
% EndTime: 2019-12-05 18:30:26
% DurationCPUTime: 3.08s
% Computational Cost: add. (62348->181), mult. (74879->231), div. (0->0), fcn. (40054->10), ass. (0->84)
t178 = sin(qJ(1));
t182 = cos(qJ(1));
t156 = t182 * g(2) + t178 * g(3);
t152 = qJDD(1) * pkin(1) + t156;
t155 = t178 * g(2) - g(3) * t182;
t183 = qJD(1) ^ 2;
t153 = -pkin(1) * t183 + t155;
t177 = sin(qJ(2));
t181 = cos(qJ(2));
t137 = t181 * t152 - t153 * t177;
t169 = qJDD(1) + qJDD(2);
t134 = pkin(2) * t169 + t137;
t138 = t177 * t152 + t181 * t153;
t170 = qJD(1) + qJD(2);
t168 = t170 ^ 2;
t135 = -pkin(2) * t168 + t138;
t173 = sin(pkin(9));
t174 = cos(pkin(9));
t129 = t174 * t134 - t135 * t173;
t126 = pkin(3) * t169 + t129;
t130 = t173 * t134 + t174 * t135;
t127 = -pkin(3) * t168 + t130;
t176 = sin(qJ(4));
t180 = cos(qJ(4));
t123 = t176 * t126 + t180 * t127;
t162 = qJD(4) + t170;
t160 = t162 ^ 2;
t161 = qJDD(4) + t169;
t120 = -pkin(4) * t160 + pkin(8) * t161 + t123;
t172 = -g(1) + qJDD(3);
t175 = sin(qJ(5));
t179 = cos(qJ(5));
t117 = -t120 * t175 + t172 * t179;
t118 = t120 * t179 + t172 * t175;
t140 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t175 + Ifges(6,2) * t179) * t162;
t141 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t175 + Ifges(6,4) * t179) * t162;
t197 = qJD(5) * t162;
t145 = t161 * t175 + t179 * t197;
t146 = t161 * t179 - t175 * t197;
t200 = mrSges(6,1) * t117 - mrSges(6,2) * t118 + Ifges(6,5) * t145 + Ifges(6,6) * t146 + Ifges(6,3) * qJDD(5) + (t140 * t175 - t141 * t179) * t162;
t144 = (-mrSges(6,1) * t179 + mrSges(6,2) * t175) * t162;
t198 = t162 * t179;
t151 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t198;
t199 = t162 * t175;
t115 = m(6) * t117 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t145 + qJD(5) * t151 - t144 * t199;
t150 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t199;
t116 = m(6) * t118 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t146 - qJD(5) * t150 + t144 * t198;
t192 = -t115 * t175 + t179 * t116;
t102 = m(5) * t123 - mrSges(5,1) * t160 - mrSges(5,2) * t161 + t192;
t122 = t126 * t180 - t127 * t176;
t119 = -pkin(4) * t161 - pkin(8) * t160 - t122;
t188 = -m(6) * t119 + t146 * mrSges(6,1) - mrSges(6,2) * t145 - t150 * t199 + t151 * t198;
t110 = m(5) * t122 + mrSges(5,1) * t161 - mrSges(5,2) * t160 + t188;
t99 = t102 * t176 + t110 * t180;
t95 = m(4) * t129 + mrSges(4,1) * t169 - mrSges(4,2) * t168 + t99;
t193 = t102 * t180 - t110 * t176;
t96 = m(4) * t130 - mrSges(4,1) * t168 - mrSges(4,2) * t169 + t193;
t90 = t173 * t96 + t174 * t95;
t87 = m(3) * t137 + mrSges(3,1) * t169 - mrSges(3,2) * t168 + t90;
t195 = -t173 * t95 + t174 * t96;
t88 = m(3) * t138 - mrSges(3,1) * t168 - mrSges(3,2) * t169 + t195;
t81 = t177 * t88 + t181 * t87;
t104 = t179 * t115 + t175 * t116;
t196 = m(5) * t172 + t104;
t194 = -t177 * t87 + t181 * t88;
t191 = m(4) * t172 + t196;
t139 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t175 + Ifges(6,6) * t179) * t162;
t107 = -mrSges(6,1) * t119 + mrSges(6,3) * t118 + Ifges(6,4) * t145 + Ifges(6,2) * t146 + Ifges(6,6) * qJDD(5) + qJD(5) * t141 - t139 * t199;
t108 = mrSges(6,2) * t119 - mrSges(6,3) * t117 + Ifges(6,1) * t145 + Ifges(6,4) * t146 + Ifges(6,5) * qJDD(5) - qJD(5) * t140 + t139 * t198;
t189 = mrSges(5,1) * t122 - mrSges(5,2) * t123 + Ifges(5,3) * t161 + pkin(4) * t188 + pkin(8) * t192 + t107 * t179 + t108 * t175;
t186 = mrSges(4,1) * t129 - mrSges(4,2) * t130 + Ifges(4,3) * t169 + pkin(3) * t99 + t189;
t185 = mrSges(3,1) * t137 - mrSges(3,2) * t138 + Ifges(3,3) * t169 + pkin(2) * t90 + t186;
t184 = mrSges(2,1) * t156 - mrSges(2,2) * t155 + Ifges(2,3) * qJDD(1) + pkin(1) * t81 + t185;
t97 = -mrSges(5,1) * t172 + mrSges(5,3) * t123 + t160 * Ifges(5,5) + Ifges(5,6) * t161 - pkin(4) * t104 - t200;
t91 = mrSges(5,2) * t172 - mrSges(5,3) * t122 + Ifges(5,5) * t161 - Ifges(5,6) * t160 - pkin(8) * t104 - t107 * t175 + t108 * t179;
t83 = mrSges(4,2) * t172 - mrSges(4,3) * t129 + Ifges(4,5) * t169 - Ifges(4,6) * t168 - pkin(7) * t99 - t176 * t97 + t180 * t91;
t82 = -mrSges(4,1) * t172 + mrSges(4,3) * t130 + t168 * Ifges(4,5) + Ifges(4,6) * t169 - pkin(3) * t196 + pkin(7) * t193 + t176 * t91 + t180 * t97;
t79 = m(2) * t156 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t183 + t81;
t78 = m(2) * t155 - mrSges(2,1) * t183 - qJDD(1) * mrSges(2,2) + t194;
t77 = -mrSges(3,2) * g(1) - mrSges(3,3) * t137 + Ifges(3,5) * t169 - Ifges(3,6) * t168 - qJ(3) * t90 - t173 * t82 + t174 * t83;
t76 = mrSges(3,1) * g(1) + mrSges(3,3) * t138 + t168 * Ifges(3,5) + Ifges(3,6) * t169 - pkin(2) * t191 + qJ(3) * t195 + t173 * t83 + t174 * t82;
t75 = -mrSges(2,2) * g(1) - mrSges(2,3) * t156 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t183 - pkin(6) * t81 - t177 * t76 + t181 * t77;
t74 = Ifges(2,6) * qJDD(1) + t183 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t155 + t177 * t77 + t181 * t76 - pkin(1) * (-m(3) * g(1) + t191) + pkin(6) * t194;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t184, t75, t77, t83, t91, t108; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t178 * t75 - t182 * t74 - pkin(5) * (-t178 * t79 + t182 * t78), t74, t76, t82, t97, t107; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t182 * t75 - t178 * t74 + pkin(5) * (-t178 * t78 - t182 * t79), t184, t185, t186, t189, t200;];
m_new = t1;
