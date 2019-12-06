% Calculate vector of cutting torques with Newton-Euler for
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:57
% EndTime: 2019-12-05 15:01:01
% DurationCPUTime: 2.75s
% Computational Cost: add. (32983->182), mult. (63549->236), div. (0->0), fcn. (43606->10), ass. (0->91)
t187 = qJD(3) ^ 2;
t179 = sin(pkin(7));
t182 = cos(pkin(7));
t166 = -t182 * g(1) - t179 * g(2);
t176 = -g(3) + qJDD(1);
t178 = sin(pkin(8));
t181 = cos(pkin(8));
t153 = -t178 * t166 + t181 * t176;
t154 = t181 * t166 + t178 * t176;
t184 = sin(qJ(3));
t186 = cos(qJ(3));
t137 = t184 * t153 + t186 * t154;
t134 = -t187 * pkin(3) + qJDD(3) * qJ(4) + t137;
t177 = sin(pkin(9));
t165 = t179 * g(1) - t182 * g(2);
t164 = qJDD(2) - t165;
t180 = cos(pkin(9));
t209 = qJD(3) * qJD(4);
t211 = t180 * t164 - 0.2e1 * t177 * t209;
t214 = pkin(4) * t180;
t126 = (-pkin(6) * qJDD(3) + t187 * t214 - t134) * t177 + t211;
t129 = t177 * t164 + (t134 + 0.2e1 * t209) * t180;
t208 = qJDD(3) * t180;
t173 = t180 ^ 2;
t212 = t173 * t187;
t127 = -pkin(4) * t212 + pkin(6) * t208 + t129;
t183 = sin(qJ(5));
t185 = cos(qJ(5));
t124 = t185 * t126 - t183 * t127;
t197 = -t177 * t183 + t180 * t185;
t155 = t197 * qJD(3);
t198 = t177 * t185 + t180 * t183;
t156 = t198 * qJD(3);
t142 = -t155 * mrSges(6,1) + t156 * mrSges(6,2);
t147 = t155 * qJD(5) + t198 * qJDD(3);
t151 = -qJD(5) * mrSges(6,2) + t155 * mrSges(6,3);
t119 = m(6) * t124 + qJDD(5) * mrSges(6,1) - t147 * mrSges(6,3) + qJD(5) * t151 - t156 * t142;
t125 = t183 * t126 + t185 * t127;
t146 = -t156 * qJD(5) + t197 * qJDD(3);
t152 = qJD(5) * mrSges(6,1) - t156 * mrSges(6,3);
t120 = m(6) * t125 - qJDD(5) * mrSges(6,2) + t146 * mrSges(6,3) - qJD(5) * t152 + t155 * t142;
t111 = t185 * t119 + t183 * t120;
t128 = -t177 * t134 + t211;
t139 = Ifges(6,4) * t156 + Ifges(6,2) * t155 + Ifges(6,6) * qJD(5);
t140 = Ifges(6,1) * t156 + Ifges(6,4) * t155 + Ifges(6,5) * qJD(5);
t191 = -mrSges(6,1) * t124 + mrSges(6,2) * t125 - Ifges(6,5) * t147 - Ifges(6,6) * t146 - Ifges(6,3) * qJDD(5) - t156 * t139 + t155 * t140;
t202 = Ifges(5,4) * t177 + Ifges(5,2) * t180;
t203 = Ifges(5,1) * t177 + Ifges(5,4) * t180;
t215 = -mrSges(5,1) * t128 + mrSges(5,2) * t129 - pkin(4) * t111 - (t177 * t202 - t180 * t203) * t187 + t191;
t213 = mrSges(5,2) * t177;
t196 = mrSges(5,3) * qJDD(3) + t187 * (-mrSges(5,1) * t180 + t213);
t109 = m(5) * t128 - t196 * t177 + t111;
t204 = -t183 * t119 + t185 * t120;
t110 = m(5) * t129 + t196 * t180 + t204;
t205 = -t177 * t109 + t180 * t110;
t100 = m(4) * t137 - t187 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t205;
t136 = t186 * t153 - t184 * t154;
t200 = qJDD(4) - t136;
t132 = -qJDD(3) * pkin(3) - t187 * qJ(4) + t200;
t172 = t177 ^ 2;
t130 = (-pkin(3) - t214) * qJDD(3) + (-qJ(4) + (-t172 - t173) * pkin(6)) * t187 + t200;
t192 = m(6) * t130 - t146 * mrSges(6,1) + t147 * mrSges(6,2) - t155 * t151 + t156 * t152;
t190 = -m(5) * t132 + mrSges(5,1) * t208 - t192 + (t172 * t187 + t212) * mrSges(5,3);
t115 = m(4) * t136 - t187 * mrSges(4,2) + (mrSges(4,1) - t213) * qJDD(3) + t190;
t95 = t184 * t100 + t186 * t115;
t104 = t180 * t109 + t177 * t110;
t201 = Ifges(5,5) * t177 + Ifges(5,6) * t180;
t210 = t187 * t201;
t91 = m(3) * t153 + t95;
t206 = t186 * t100 - t184 * t115;
t92 = m(3) * t154 + t206;
t207 = -t178 * t91 + t181 * t92;
t195 = (-m(3) - m(4)) * t164 - t104;
t138 = Ifges(6,5) * t156 + Ifges(6,6) * t155 + Ifges(6,3) * qJD(5);
t112 = -mrSges(6,1) * t130 + mrSges(6,3) * t125 + Ifges(6,4) * t147 + Ifges(6,2) * t146 + Ifges(6,6) * qJDD(5) + qJD(5) * t140 - t156 * t138;
t113 = mrSges(6,2) * t130 - mrSges(6,3) * t124 + Ifges(6,1) * t147 + Ifges(6,4) * t146 + Ifges(6,5) * qJDD(5) - qJD(5) * t139 + t155 * t138;
t94 = -mrSges(5,1) * t132 + mrSges(5,3) * t129 - pkin(4) * t192 + pkin(6) * t204 + t202 * qJDD(3) + t185 * t112 + t183 * t113 - t177 * t210;
t97 = mrSges(5,2) * t132 - mrSges(5,3) * t128 - pkin(6) * t111 + t203 * qJDD(3) - t183 * t112 + t185 * t113 + t180 * t210;
t85 = mrSges(4,2) * t164 - mrSges(4,3) * t136 + Ifges(4,5) * qJDD(3) - t187 * Ifges(4,6) - qJ(4) * t104 - t177 * t94 + t180 * t97;
t89 = (Ifges(4,6) - t201) * qJDD(3) + t187 * Ifges(4,5) - mrSges(4,1) * t164 + mrSges(4,3) * t137 - pkin(3) * t104 + t215;
t81 = -mrSges(3,1) * t164 + mrSges(3,3) * t154 + t184 * t85 + t186 * t89 - pkin(2) * (m(4) * t164 + t104) + pkin(5) * t206;
t84 = mrSges(3,2) * t164 - mrSges(3,3) * t153 - pkin(5) * t95 - t184 * t89 + t186 * t85;
t194 = mrSges(2,1) * t165 - mrSges(2,2) * t166 + pkin(1) * t195 + qJ(2) * t207 + t178 * t84 + t181 * t81;
t193 = -mrSges(4,2) * t137 + qJ(4) * t205 + t177 * t97 + t180 * t94 + pkin(3) * (-qJDD(3) * t213 + t190) + mrSges(4,1) * t136 + Ifges(4,3) * qJDD(3);
t189 = mrSges(3,1) * t153 - mrSges(3,2) * t154 + pkin(2) * t95 + t193;
t101 = m(2) * t165 + t195;
t88 = t178 * t92 + t181 * t91;
t86 = m(2) * t166 + t207;
t82 = -mrSges(2,1) * t176 + mrSges(2,3) * t166 - pkin(1) * t88 - t189;
t79 = mrSges(2,2) * t176 - mrSges(2,3) * t165 - qJ(2) * t88 - t178 * t81 + t181 * t84;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t182 * t79 - t179 * t82 - qJ(1) * (t182 * t101 + t179 * t86), t79, t84, t85, t97, t113; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t179 * t79 + t182 * t82 + qJ(1) * (-t179 * t101 + t182 * t86), t82, t81, t89, t94, t112; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t194, t194, t189, t193, t201 * qJDD(3) - t215, -t191;];
m_new = t1;
