% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:44
% EndTime: 2019-12-31 16:53:48
% DurationCPUTime: 2.32s
% Computational Cost: add. (23543->215), mult. (55284->274), div. (0->0), fcn. (37003->8), ass. (0->97)
t193 = qJD(1) ^ 2;
t184 = sin(pkin(7));
t185 = cos(pkin(7));
t187 = sin(qJ(3));
t190 = cos(qJ(3));
t202 = t184 * t187 - t185 * t190;
t164 = t202 * qJD(1);
t188 = sin(qJ(1));
t191 = cos(qJ(1));
t173 = -t191 * g(1) - t188 * g(2);
t166 = -t193 * pkin(1) + qJDD(1) * qJ(2) + t173;
t214 = qJD(1) * qJD(2);
t212 = -t185 * g(3) - 0.2e1 * t184 * t214;
t155 = -t184 * t166 + t212;
t156 = -t184 * g(3) + (t166 + 0.2e1 * t214) * t185;
t220 = pkin(2) * t185;
t140 = (-pkin(5) * qJDD(1) + t193 * t220 - t166) * t184 + t212;
t213 = qJDD(1) * t185;
t181 = t185 ^ 2;
t218 = t181 * t193;
t141 = -pkin(2) * t218 + pkin(5) * t213 + t156;
t126 = t187 * t140 + t190 * t141;
t203 = t184 * t190 + t185 * t187;
t165 = t203 * qJD(1);
t151 = t164 * pkin(3) - t165 * pkin(6);
t192 = qJD(3) ^ 2;
t122 = -t192 * pkin(3) + qJDD(3) * pkin(6) - t164 * t151 + t126;
t180 = t184 ^ 2;
t172 = t188 * g(1) - t191 * g(2);
t208 = qJDD(2) - t172;
t152 = (-pkin(1) - t220) * qJDD(1) + (-qJ(2) + (-t180 - t181) * pkin(5)) * t193 + t208;
t215 = t165 * qJD(3);
t153 = -t202 * qJDD(1) - t215;
t216 = t164 * qJD(3);
t154 = t203 * qJDD(1) - t216;
t123 = (-t154 + t216) * pkin(6) + (-t153 + t215) * pkin(3) + t152;
t186 = sin(qJ(4));
t189 = cos(qJ(4));
t120 = t189 * t122 + t186 * t123;
t125 = t190 * t140 - t187 * t141;
t121 = -qJDD(3) * pkin(3) - t192 * pkin(6) + t165 * t151 - t125;
t157 = t189 * qJD(3) - t186 * t165;
t158 = t186 * qJD(3) + t189 * t165;
t162 = qJD(4) + t164;
t127 = Ifges(5,5) * t158 + Ifges(5,6) * t157 + Ifges(5,3) * t162;
t129 = Ifges(5,1) * t158 + Ifges(5,4) * t157 + Ifges(5,5) * t162;
t132 = -t158 * qJD(4) + t189 * qJDD(3) - t186 * t154;
t133 = t157 * qJD(4) + t186 * qJDD(3) + t189 * t154;
t150 = qJDD(4) - t153;
t108 = -mrSges(5,1) * t121 + mrSges(5,3) * t120 + Ifges(5,4) * t133 + Ifges(5,2) * t132 + Ifges(5,6) * t150 - t158 * t127 + t162 * t129;
t119 = -t186 * t122 + t189 * t123;
t128 = Ifges(5,4) * t158 + Ifges(5,2) * t157 + Ifges(5,6) * t162;
t109 = mrSges(5,2) * t121 - mrSges(5,3) * t119 + Ifges(5,1) * t133 + Ifges(5,4) * t132 + Ifges(5,5) * t150 + t157 * t127 - t162 * t128;
t143 = Ifges(4,4) * t165 - Ifges(4,2) * t164 + Ifges(4,6) * qJD(3);
t144 = Ifges(4,1) * t165 - Ifges(4,4) * t164 + Ifges(4,5) * qJD(3);
t137 = -t162 * mrSges(5,2) + t157 * mrSges(5,3);
t138 = t162 * mrSges(5,1) - t158 * mrSges(5,3);
t199 = -m(5) * t121 + t132 * mrSges(5,1) - t133 * mrSges(5,2) + t157 * t137 - t158 * t138;
t135 = -t157 * mrSges(5,1) + t158 * mrSges(5,2);
t115 = m(5) * t119 + t150 * mrSges(5,1) - t133 * mrSges(5,3) - t158 * t135 + t162 * t137;
t116 = m(5) * t120 - t150 * mrSges(5,2) + t132 * mrSges(5,3) + t157 * t135 - t162 * t138;
t209 = -t186 * t115 + t189 * t116;
t197 = -mrSges(4,1) * t125 + mrSges(4,2) * t126 - Ifges(4,5) * t154 - Ifges(4,6) * t153 - Ifges(4,3) * qJDD(3) - pkin(3) * t199 - pkin(6) * t209 - t189 * t108 - t186 * t109 - t165 * t143 - t164 * t144;
t206 = Ifges(3,4) * t184 + Ifges(3,2) * t185;
t207 = Ifges(3,1) * t184 + Ifges(3,4) * t185;
t146 = t164 * mrSges(4,1) + t165 * mrSges(4,2);
t160 = qJD(3) * mrSges(4,1) - t165 * mrSges(4,3);
t102 = m(4) * t126 - qJDD(3) * mrSges(4,2) + t153 * mrSges(4,3) - qJD(3) * t160 - t164 * t146 + t209;
t159 = -qJD(3) * mrSges(4,2) - t164 * mrSges(4,3);
t111 = m(4) * t125 + qJDD(3) * mrSges(4,1) - t154 * mrSges(4,3) + qJD(3) * t159 - t165 * t146 + t199;
t97 = t187 * t102 + t190 * t111;
t221 = -mrSges(3,1) * t155 + mrSges(3,2) * t156 - pkin(2) * t97 - (t184 * t206 - t185 * t207) * t193 + t197;
t219 = mrSges(3,2) * t184;
t104 = t189 * t115 + t186 * t116;
t205 = Ifges(3,5) * t184 + Ifges(3,6) * t185;
t217 = t193 * t205;
t201 = mrSges(3,3) * qJDD(1) + t193 * (-mrSges(3,1) * t185 + t219);
t95 = m(3) * t155 - t201 * t184 + t97;
t210 = t190 * t102 - t187 * t111;
t96 = m(3) * t156 + t201 * t185 + t210;
t211 = -t184 * t95 + t185 * t96;
t163 = -qJDD(1) * pkin(1) - t193 * qJ(2) + t208;
t198 = m(4) * t152 - t153 * mrSges(4,1) + t154 * mrSges(4,2) + t164 * t159 + t165 * t160 + t104;
t196 = -m(3) * t163 + mrSges(3,1) * t213 - t198 + (t180 * t193 + t218) * mrSges(3,3);
t142 = Ifges(4,5) * t165 - Ifges(4,6) * t164 + Ifges(4,3) * qJD(3);
t92 = mrSges(4,2) * t152 - mrSges(4,3) * t125 + Ifges(4,1) * t154 + Ifges(4,4) * t153 + Ifges(4,5) * qJDD(3) - pkin(6) * t104 - qJD(3) * t143 - t186 * t108 + t189 * t109 - t164 * t142;
t195 = mrSges(5,1) * t119 - mrSges(5,2) * t120 + Ifges(5,5) * t133 + Ifges(5,6) * t132 + Ifges(5,3) * t150 + t158 * t128 - t157 * t129;
t93 = -mrSges(4,1) * t152 + mrSges(4,3) * t126 + Ifges(4,4) * t154 + Ifges(4,2) * t153 + Ifges(4,6) * qJDD(3) - pkin(3) * t104 + qJD(3) * t144 - t165 * t142 - t195;
t86 = -mrSges(3,1) * t163 + mrSges(3,3) * t156 - pkin(2) * t198 + pkin(5) * t210 + t206 * qJDD(1) - t184 * t217 + t187 * t92 + t190 * t93;
t88 = mrSges(3,2) * t163 - mrSges(3,3) * t155 - pkin(5) * t97 + t207 * qJDD(1) + t185 * t217 - t187 * t93 + t190 * t92;
t200 = -mrSges(2,2) * t173 + qJ(2) * t211 + t184 * t88 + t185 * t86 + mrSges(2,1) * t172 + Ifges(2,3) * qJDD(1) + pkin(1) * (-qJDD(1) * t219 + t196);
t98 = t196 - t193 * mrSges(2,2) + m(2) * t172 + (mrSges(2,1) - t219) * qJDD(1);
t91 = t184 * t96 + t185 * t95;
t89 = m(2) * t173 - t193 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t211;
t84 = t193 * Ifges(2,5) + mrSges(2,3) * t173 - pkin(1) * t91 + mrSges(2,1) * g(3) + (Ifges(2,6) - t205) * qJDD(1) + t221;
t83 = -mrSges(2,2) * g(3) - mrSges(2,3) * t172 + Ifges(2,5) * qJDD(1) - t193 * Ifges(2,6) - qJ(2) * t91 - t184 * t86 + t185 * t88;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t191 * t83 - t188 * t84 - pkin(4) * (t188 * t89 + t191 * t98), t83, t88, t92, t109; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t188 * t83 + t191 * t84 + pkin(4) * (-t188 * t98 + t191 * t89), t84, t86, t93, t108; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t200, t200, t205 * qJDD(1) - t221, -t197, t195;];
m_new = t1;
