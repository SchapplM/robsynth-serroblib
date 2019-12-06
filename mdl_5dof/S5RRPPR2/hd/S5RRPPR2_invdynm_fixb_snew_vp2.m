% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:10
% EndTime: 2019-12-05 18:20:15
% DurationCPUTime: 3.35s
% Computational Cost: add. (55820->203), mult. (73380->265), div. (0->0), fcn. (40823->10), ass. (0->101)
t231 = 2 * qJD(4);
t195 = sin(qJ(1));
t198 = cos(qJ(1));
t173 = t198 * g(2) + t195 * g(3);
t166 = qJDD(1) * pkin(1) + t173;
t172 = t195 * g(2) - t198 * g(3);
t199 = qJD(1) ^ 2;
t167 = -t199 * pkin(1) + t172;
t194 = sin(qJ(2));
t197 = cos(qJ(2));
t150 = t197 * t166 - t194 * t167;
t185 = qJDD(1) + qJDD(2);
t144 = t185 * pkin(2) + t150;
t151 = t194 * t166 + t197 * t167;
t186 = (qJD(1) + qJD(2));
t184 = t186 ^ 2;
t145 = -t184 * pkin(2) + t151;
t190 = sin(pkin(8));
t192 = cos(pkin(8));
t140 = t190 * t144 + t192 * t145;
t137 = -t184 * pkin(3) + t185 * qJ(4) + t140;
t230 = (t186 * t231) + t137;
t189 = sin(pkin(9));
t225 = t186 * t189;
t191 = cos(pkin(9));
t223 = t191 * t186;
t188 = -g(1) + qJDD(3);
t133 = t189 * t188 + t230 * t191;
t213 = -pkin(4) * t191 - pkin(7) * t189;
t162 = t213 * t186;
t131 = t162 * t223 + t133;
t139 = t192 * t144 - t190 * t145;
t206 = -t184 * qJ(4) + qJDD(4) - t139;
t134 = (-pkin(3) + t213) * t185 + t206;
t193 = sin(qJ(5));
t196 = cos(qJ(5));
t127 = -t193 * t131 + t196 * t134;
t170 = qJD(5) - t223;
t219 = t193 * t225;
t153 = -t170 * mrSges(6,2) - mrSges(6,3) * t219;
t155 = (mrSges(6,1) * t193 + mrSges(6,2) * t196) * t225;
t220 = qJD(5) * t186;
t157 = (t185 * t196 - t193 * t220) * t189;
t224 = t191 * t185;
t169 = qJDD(5) - t224;
t218 = t196 * t225;
t125 = m(6) * t127 + t169 * mrSges(6,1) - t157 * mrSges(6,3) + t170 * t153 - t155 * t218;
t128 = t196 * t131 + t193 * t134;
t154 = t170 * mrSges(6,1) - mrSges(6,3) * t218;
t156 = (-t185 * t193 - t196 * t220) * t189;
t126 = m(6) * t128 - t169 * mrSges(6,2) + t156 * mrSges(6,3) - t170 * t154 - t155 * t219;
t119 = -t193 * t125 + t196 * t126;
t222 = t191 * t188;
t130 = -t222 + (t137 + (t231 + t162) * t186) * t189;
t146 = Ifges(6,3) * t170 + (Ifges(6,5) * t196 - Ifges(6,6) * t193) * t225;
t148 = Ifges(6,5) * t170 + (Ifges(6,1) * t196 - Ifges(6,4) * t193) * t225;
t120 = -mrSges(6,1) * t130 + mrSges(6,3) * t128 + Ifges(6,4) * t157 + Ifges(6,2) * t156 + Ifges(6,6) * t169 - t146 * t218 + t170 * t148;
t147 = Ifges(6,6) * t170 + (Ifges(6,4) * t196 - Ifges(6,2) * t193) * t225;
t121 = mrSges(6,2) * t130 - mrSges(6,3) * t127 + Ifges(6,1) * t157 + Ifges(6,4) * t156 + Ifges(6,5) * t169 - t146 * t219 - t170 * t147;
t132 = -t230 * t189 + t222;
t207 = -m(6) * t130 + t156 * mrSges(6,1) - t157 * mrSges(6,2);
t209 = -t153 * t193 - t154 * t196;
t212 = Ifges(5,1) * t189 + Ifges(5,4) * t191;
t229 = -((Ifges(5,4) * t189 + Ifges(5,2) * t191) * t225 - t212 * t223) * t186 - mrSges(5,1) * t132 + mrSges(5,2) * t133 - pkin(4) * (t209 * t225 + t207) - pkin(7) * t119 - t196 * t120 - t193 * t121;
t228 = mrSges(5,2) * t189;
t158 = (-mrSges(5,1) * t191 + t228) * t186;
t226 = mrSges(5,3) * t185;
t116 = m(5) * t133 + (t158 * t186 + t226) * t191 + t119;
t123 = m(5) * t132 + (-t226 + (-t158 + t209) * t186) * t189 + t207;
t214 = t191 * t116 - t189 * t123;
t109 = m(4) * t140 - t184 * mrSges(4,1) - t185 * mrSges(4,2) + t214;
t118 = t196 * t125 + t193 * t126;
t136 = -t185 * pkin(3) + t206;
t204 = -m(5) * t136 + mrSges(5,1) * t224 - t118 + (t189 ^ 2 + t191 ^ 2) * mrSges(5,3) * t184;
t113 = m(4) * t139 - t184 * mrSges(4,2) + (mrSges(4,1) - t228) * t185 + t204;
t215 = t192 * t109 - t190 * t113;
t100 = m(3) * t151 - t184 * mrSges(3,1) - t185 * mrSges(3,2) + t215;
t102 = t190 * t109 + t192 * t113;
t99 = m(3) * t150 + t185 * mrSges(3,1) - t184 * mrSges(3,2) + t102;
t94 = t194 * t100 + t197 * t99;
t111 = t189 * t116 + t191 * t123;
t217 = m(4) * t188 + t111;
t216 = t197 * t100 - t194 * t99;
t211 = Ifges(5,5) * t189 + Ifges(5,6) * t191;
t210 = t147 * t196 + t148 * t193;
t159 = t211 * t186;
t104 = mrSges(5,2) * t136 - mrSges(5,3) * t132 - pkin(7) * t118 - t193 * t120 + t196 * t121 + t159 * t223 + t212 * t185;
t203 = mrSges(6,1) * t127 - mrSges(6,2) * t128 + Ifges(6,5) * t157 + Ifges(6,6) * t156 + Ifges(6,3) * t169;
t106 = Ifges(5,2) * t224 - mrSges(5,1) * t136 + mrSges(5,3) * t133 - pkin(4) * t118 + (Ifges(5,4) * t185 + (-t159 - t210) * t186) * t189 - t203;
t205 = -mrSges(4,2) * t140 + qJ(4) * t214 + t189 * t104 + t191 * t106 + pkin(3) * (-t185 * t228 + t204) + mrSges(4,1) * t139 + Ifges(4,3) * t185;
t202 = mrSges(3,1) * t150 - mrSges(3,2) * t151 + Ifges(3,3) * t185 + pkin(2) * t102 + t205;
t200 = mrSges(2,1) * t173 - mrSges(2,2) * t172 + Ifges(2,3) * qJDD(1) + pkin(1) * t94 + t202;
t95 = -mrSges(4,1) * t188 + mrSges(4,3) * t140 + t184 * Ifges(4,5) - pkin(3) * t111 + (Ifges(4,6) - t211) * t185 + t229;
t92 = m(2) * t173 + qJDD(1) * mrSges(2,1) - t199 * mrSges(2,2) + t94;
t91 = m(2) * t172 - t199 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t216;
t90 = mrSges(4,2) * t188 - mrSges(4,3) * t139 + Ifges(4,5) * t185 - t184 * Ifges(4,6) - qJ(4) * t111 + t191 * t104 - t189 * t106;
t89 = -mrSges(3,2) * g(1) - mrSges(3,3) * t150 + Ifges(3,5) * t185 - t184 * Ifges(3,6) - qJ(3) * t102 - t190 * t95 + t192 * t90;
t88 = mrSges(3,1) * g(1) + mrSges(3,3) * t151 + t184 * Ifges(3,5) + Ifges(3,6) * t185 - pkin(2) * t217 + qJ(3) * t215 + t190 * t90 + t192 * t95;
t87 = -mrSges(2,2) * g(1) - mrSges(2,3) * t173 + Ifges(2,5) * qJDD(1) - t199 * Ifges(2,6) - pkin(6) * t94 - t194 * t88 + t197 * t89;
t86 = Ifges(2,6) * qJDD(1) + t199 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t172 + t194 * t89 + t197 * t88 - pkin(1) * (-m(3) * g(1) + t217) + pkin(6) * t216;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t200, t87, t89, t90, t104, t121; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t195 * t87 - t198 * t86 - pkin(5) * (-t195 * t92 + t198 * t91), t86, t88, t95, t106, t120; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t198 * t87 - t195 * t86 + pkin(5) * (-t195 * t91 - t198 * t92), t200, t202, t205, t211 * t185 - t229, t210 * t225 + t203;];
m_new = t1;
