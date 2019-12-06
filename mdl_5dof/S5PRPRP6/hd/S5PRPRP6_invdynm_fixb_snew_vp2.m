% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRP6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:23
% EndTime: 2019-12-05 15:40:26
% DurationCPUTime: 1.06s
% Computational Cost: add. (8911->210), mult. (15573->253), div. (0->0), fcn. (7013->6), ass. (0->82)
t184 = sin(pkin(7));
t213 = cos(pkin(7));
t165 = -t213 * g(1) - t184 * g(2);
t181 = -g(3) + qJDD(1);
t186 = sin(qJ(2));
t188 = cos(qJ(2));
t131 = t188 * t165 + t186 * t181;
t220 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t131;
t219 = m(3) + m(4);
t218 = -pkin(2) - pkin(6);
t217 = mrSges(3,1) - mrSges(4,2);
t216 = -mrSges(5,3) - mrSges(6,2);
t215 = (Ifges(3,5) - Ifges(4,4));
t214 = -Ifges(3,6) + Ifges(4,5);
t164 = t184 * g(1) - t213 * g(2);
t185 = sin(qJ(4));
t212 = t185 * t164;
t130 = -t186 * t165 + t188 * t181;
t190 = qJD(2) ^ 2;
t199 = -t190 * qJ(3) + qJDD(3) - t130;
t126 = t218 * qJDD(2) + t199;
t187 = cos(qJ(4));
t122 = t185 * t126 - t187 * t164;
t209 = qJD(2) * qJD(4);
t159 = t185 * qJDD(2) + t187 * t209;
t210 = qJD(2) * t187;
t167 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t210;
t157 = (mrSges(6,1) * t185 - mrSges(6,3) * t187) * qJD(2);
t203 = qJD(2) * (-t157 - (mrSges(5,1) * t185 + mrSges(5,2) * t187) * qJD(2));
t156 = (pkin(4) * t185 - qJ(5) * t187) * qJD(2);
t189 = qJD(4) ^ 2;
t211 = qJD(2) * t185;
t117 = -t189 * pkin(4) + qJDD(4) * qJ(5) + (2 * qJD(5) * qJD(4)) - t156 * t211 + t122;
t168 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t210;
t207 = m(6) * t117 + qJDD(4) * mrSges(6,3) + qJD(4) * t168;
t106 = m(5) * t122 - qJDD(4) * mrSges(5,2) - qJD(4) * t167 + t216 * t159 + t185 * t203 + t207;
t121 = t187 * t126 + t212;
t160 = t187 * qJDD(2) - t185 * t209;
t166 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t211;
t119 = -qJDD(4) * pkin(4) - t189 * qJ(5) - t212 + qJDD(5) + (qJD(2) * t156 - t126) * t187;
t169 = -mrSges(6,2) * t211 + qJD(4) * mrSges(6,3);
t202 = -m(6) * t119 + qJDD(4) * mrSges(6,1) + qJD(4) * t169;
t107 = m(5) * t121 + qJDD(4) * mrSges(5,1) + qJD(4) * t166 + t216 * t160 + t187 * t203 + t202;
t99 = t187 * t106 - t185 * t107;
t125 = t218 * t190 - t220;
t114 = t159 * pkin(4) - t160 * qJ(5) + (-0.2e1 * qJD(5) * t187 + (pkin(4) * t187 + qJ(5) * t185) * qJD(4)) * qJD(2) + t125;
t108 = m(6) * t114 + t159 * mrSges(6,1) - t160 * mrSges(6,3) - t168 * t210 + t169 * t211;
t103 = -m(5) * t125 - t159 * mrSges(5,1) - t160 * mrSges(5,2) - t166 * t211 - t167 * t210 - t108;
t127 = (t190 * pkin(2)) + t220;
t102 = -m(4) * t127 + (t190 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t103;
t101 = m(3) * t131 - (t190 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t102;
t129 = -qJDD(2) * pkin(2) + t199;
t98 = t185 * t106 + t187 * t107;
t197 = -m(4) * t129 + t190 * mrSges(4,3) - t98;
t93 = m(3) * t130 - t190 * mrSges(3,2) + t217 * qJDD(2) + t197;
t205 = t188 * t101 - t186 * t93;
t138 = (Ifges(6,2) * qJD(4)) + (Ifges(6,4) * t187 + Ifges(6,6) * t185) * qJD(2);
t204 = qJD(2) * (-(Ifges(5,3) * qJD(4)) - (Ifges(5,5) * t187 - Ifges(5,6) * t185) * qJD(2) - t138);
t201 = -mrSges(6,1) * t114 + mrSges(6,2) * t117;
t140 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t187 + Ifges(6,5) * t185) * qJD(2);
t141 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t187 - Ifges(5,4) * t185) * qJD(2);
t91 = -mrSges(5,1) * t125 + mrSges(5,3) * t122 - pkin(4) * t108 + (Ifges(5,4) - Ifges(6,5)) * t160 + (-Ifges(5,2) - Ifges(6,3)) * t159 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + (t140 + t141) * qJD(4) + t187 * t204 + t201;
t139 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t187 - Ifges(5,2) * t185) * qJD(2);
t136 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t187 + Ifges(6,3) * t185) * qJD(2);
t198 = mrSges(6,2) * t119 - mrSges(6,3) * t114 + Ifges(6,1) * t160 + Ifges(6,4) * qJDD(4) + Ifges(6,5) * t159 + qJD(4) * t136;
t92 = mrSges(5,2) * t125 - mrSges(5,3) * t121 + Ifges(5,1) * t160 - Ifges(5,4) * t159 + Ifges(5,5) * qJDD(4) - qJ(5) * t108 - qJD(4) * t139 + t185 * t204 + t198;
t194 = -mrSges(4,1) * t127 - pkin(3) * t103 - pkin(6) * t99 - t185 * t92 - t187 * t91;
t97 = -m(4) * t164 + t99;
t84 = mrSges(3,3) * t131 - pkin(2) * t97 - t214 * qJDD(2) + t217 * t164 + (t215 * t190) + t194;
t195 = -mrSges(6,1) * t119 + mrSges(6,3) * t117 + Ifges(6,4) * t160 + Ifges(6,2) * qJDD(4) + Ifges(6,6) * t159 - t136 * t210 + t140 * t211;
t193 = -mrSges(5,2) * t122 + qJ(5) * (-t159 * mrSges(6,2) - t157 * t211 + t207) + pkin(4) * (-t160 * mrSges(6,2) - t157 * t210 + t202) + mrSges(5,1) * t121 + t141 * t211 + t139 * t210 - Ifges(5,6) * t159 + Ifges(5,5) * t160 + Ifges(5,3) * qJDD(4) + t195;
t191 = mrSges(4,1) * t129 + pkin(3) * t98 + t193;
t86 = t191 + t214 * t190 + (-mrSges(3,2) + mrSges(4,3)) * t164 + t215 * qJDD(2) - qJ(3) * t97 - mrSges(3,3) * t130;
t200 = -mrSges(2,2) * t165 + pkin(5) * t205 + t186 * t86 + t188 * t84 + mrSges(2,1) * t164 + pkin(1) * (t219 * t164 - t99);
t196 = mrSges(4,2) * t129 - mrSges(4,3) * t127 + Ifges(4,1) * qJDD(2) - pkin(6) * t98 - t185 * t91 + t187 * t92;
t192 = mrSges(3,1) * t130 - mrSges(3,2) * t131 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) + t197) + qJ(3) * t102 + t196;
t95 = (m(2) + t219) * t164 - t99;
t89 = t186 * t101 + t188 * t93;
t87 = m(2) * t165 + t205;
t82 = -mrSges(2,1) * t181 + mrSges(2,3) * t165 - pkin(1) * t89 - t192;
t81 = mrSges(2,2) * t181 - mrSges(2,3) * t164 - pkin(5) * t89 - t186 * t84 + t188 * t86;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t213 * t81 - t184 * t82 - qJ(1) * (t184 * t87 + t213 * t95), t81, t86, t196, t92, -t138 * t211 + t198; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t184 * t81 + t213 * t82 + qJ(1) * (-t184 * t95 + t213 * t87), t82, t84, -mrSges(4,3) * t164 + Ifges(4,4) * qJDD(2) - (t190 * Ifges(4,5)) - t191, t91, t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t200, t200, t192, mrSges(4,2) * t164 + t190 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t194, t193, Ifges(6,5) * t160 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t159 - qJD(4) * t140 + t138 * t210 - t201;];
m_new = t1;
