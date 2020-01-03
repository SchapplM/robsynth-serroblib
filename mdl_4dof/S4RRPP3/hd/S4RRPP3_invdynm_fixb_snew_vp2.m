% Calculate vector of cutting torques with Newton-Euler for
% S4RRPP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:24
% EndTime: 2019-12-31 16:57:27
% DurationCPUTime: 1.41s
% Computational Cost: add. (12640->234), mult. (28633->289), div. (0->0), fcn. (16444->6), ass. (0->86)
t192 = cos(qJ(2));
t191 = sin(qJ(1));
t193 = cos(qJ(1));
t179 = -t193 * g(1) - t191 * g(2);
t195 = qJD(1) ^ 2;
t167 = -t195 * pkin(1) + qJDD(1) * pkin(5) + t179;
t190 = sin(qJ(2));
t216 = t190 * t167;
t152 = -t192 * g(3) - t216;
t153 = -t190 * g(3) + t192 * t167;
t163 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t190 + Ifges(3,2) * t192) * qJD(1);
t164 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t190 + Ifges(3,4) * t192) * qJD(1);
t211 = qJD(1) * qJD(2);
t172 = t190 * qJDD(1) + t192 * t211;
t173 = t192 * qJDD(1) - t190 * t211;
t189 = sin(pkin(6));
t217 = cos(pkin(6));
t161 = (t189 * t192 + t217 * t190) * qJD(1);
t219 = pkin(2) * t195;
t123 = qJDD(2) * pkin(2) - t172 * qJ(3) - t216 + (qJ(3) * t211 + t190 * t219 - g(3)) * t192;
t213 = qJD(1) * t190;
t175 = qJD(2) * pkin(2) - qJ(3) * t213;
t188 = t192 ^ 2;
t124 = t173 * qJ(3) - qJD(2) * t175 - t188 * t219 + t153;
t203 = t217 * t123 - t189 * t124;
t220 = -2 * qJD(3);
t117 = t161 * t220 + t203;
t212 = qJD(1) * t192;
t160 = t189 * t213 - t217 * t212;
t118 = t189 * t123 + t217 * t124 + t160 * t220;
t130 = Ifges(4,4) * t161 - Ifges(4,2) * t160 + Ifges(4,6) * qJD(2);
t132 = Ifges(4,1) * t161 - Ifges(4,4) * t160 + Ifges(4,5) * qJD(2);
t137 = t160 * mrSges(5,1) - t161 * mrSges(5,3);
t148 = t189 * t172 - t217 * t173;
t149 = t217 * t172 + t189 * t173;
t136 = t160 * pkin(3) - t161 * qJ(4);
t194 = qJD(2) ^ 2;
t111 = -t194 * pkin(3) + qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) - t160 * t136 + t118;
t113 = -qJDD(2) * pkin(3) - t194 * qJ(4) + qJDD(4) + ((2 * qJD(3)) + t136) * t161 - t203;
t127 = Ifges(5,5) * t161 + Ifges(5,6) * qJD(2) + Ifges(5,3) * t160;
t131 = Ifges(5,1) * t161 + Ifges(5,4) * qJD(2) + Ifges(5,5) * t160;
t200 = mrSges(5,1) * t113 - mrSges(5,3) * t111 - Ifges(5,4) * t149 - Ifges(5,2) * qJDD(2) - Ifges(5,6) * t148 + t161 * t127 - t160 * t131;
t157 = -t160 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t207 = -m(5) * t113 + qJDD(2) * mrSges(5,1) + qJD(2) * t157;
t156 = -qJD(2) * mrSges(5,1) + t161 * mrSges(5,2);
t210 = m(5) * t111 + qJDD(2) * mrSges(5,3) + qJD(2) * t156;
t198 = mrSges(4,2) * t118 - t160 * t132 - qJ(4) * (-t148 * mrSges(5,2) - t160 * t137 + t210) - pkin(3) * (-t149 * mrSges(5,2) - t161 * t137 + t207) - mrSges(4,1) * t117 - t161 * t130 + Ifges(4,6) * t148 - Ifges(4,5) * t149 - Ifges(4,3) * qJDD(2) + t200;
t155 = qJD(2) * mrSges(4,1) - t161 * mrSges(4,3);
t214 = -t160 * mrSges(4,1) - t161 * mrSges(4,2) - t137;
t218 = -mrSges(4,3) - mrSges(5,2);
t102 = m(4) * t118 - qJDD(2) * mrSges(4,2) - qJD(2) * t155 + t218 * t148 + t214 * t160 + t210;
t154 = -qJD(2) * mrSges(4,2) - t160 * mrSges(4,3);
t103 = m(4) * t117 + qJDD(2) * mrSges(4,1) + qJD(2) * t154 + t218 * t149 + t214 * t161 + t207;
t96 = t189 * t102 + t217 * t103;
t221 = mrSges(3,1) * t152 - mrSges(3,2) * t153 + Ifges(3,5) * t172 + Ifges(3,6) * t173 + Ifges(3,3) * qJDD(2) + pkin(2) * t96 + (t190 * t163 - t192 * t164) * qJD(1) - t198;
t129 = Ifges(5,4) * t161 + Ifges(5,2) * qJD(2) + Ifges(5,6) * t160;
t215 = -Ifges(4,5) * t161 + Ifges(4,6) * t160 - Ifges(4,3) * qJD(2) - t129;
t178 = t191 * g(1) - t193 * g(2);
t171 = (-mrSges(3,1) * t192 + mrSges(3,2) * t190) * qJD(1);
t177 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t212;
t94 = m(3) * t152 + qJDD(2) * mrSges(3,1) - t172 * mrSges(3,3) + qJD(2) * t177 - t171 * t213 + t96;
t176 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t213;
t208 = t217 * t102 - t189 * t103;
t95 = m(3) * t153 - qJDD(2) * mrSges(3,2) + t173 * mrSges(3,3) - qJD(2) * t176 + t171 * t212 + t208;
t209 = -t190 * t94 + t192 * t95;
t204 = -qJDD(1) * pkin(1) - t178;
t125 = -t173 * pkin(2) + qJDD(3) + t175 * t213 + (-qJ(3) * t188 - pkin(5)) * t195 + t204;
t115 = -0.2e1 * qJD(4) * t161 + (qJD(2) * t160 - t149) * qJ(4) + (qJD(2) * t161 + t148) * pkin(3) + t125;
t206 = -mrSges(5,1) * t115 + mrSges(5,2) * t111;
t104 = m(5) * t115 + t148 * mrSges(5,1) - t149 * mrSges(5,3) - t161 * t156 + t160 * t157;
t202 = mrSges(5,2) * t113 - mrSges(5,3) * t115 + Ifges(5,1) * t149 + Ifges(5,4) * qJDD(2) + Ifges(5,5) * t148 + qJD(2) * t127;
t166 = -t195 * pkin(5) + t204;
t199 = m(4) * t125 + t148 * mrSges(4,1) + t149 * mrSges(4,2) + t160 * t154 + t161 * t155 + t104;
t197 = -m(3) * t166 + t173 * mrSges(3,1) - t172 * mrSges(3,2) - t176 * t213 + t177 * t212 - t199;
t162 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t190 + Ifges(3,6) * t192) * qJD(1);
t91 = -mrSges(4,1) * t125 + mrSges(4,3) * t118 - pkin(3) * t104 + t215 * t161 + (Ifges(4,4) - Ifges(5,5)) * t149 + (-Ifges(4,2) - Ifges(5,3)) * t148 + (Ifges(4,6) - Ifges(5,6)) * qJDD(2) + (t131 + t132) * qJD(2) + t206;
t92 = mrSges(4,2) * t125 - mrSges(4,3) * t117 + Ifges(4,1) * t149 - Ifges(4,4) * t148 + Ifges(4,5) * qJDD(2) - qJ(4) * t104 - qJD(2) * t130 + t215 * t160 + t202;
t85 = -mrSges(3,1) * t166 + mrSges(3,3) * t153 + Ifges(3,4) * t172 + Ifges(3,2) * t173 + Ifges(3,6) * qJDD(2) - pkin(2) * t199 + qJ(3) * t208 + qJD(2) * t164 - t162 * t213 + t189 * t92 + t217 * t91;
t87 = mrSges(3,2) * t166 - mrSges(3,3) * t152 + Ifges(3,1) * t172 + Ifges(3,4) * t173 + Ifges(3,5) * qJDD(2) - qJ(3) * t96 - qJD(2) * t163 + t162 * t212 - t189 * t91 + t217 * t92;
t201 = mrSges(2,1) * t178 - mrSges(2,2) * t179 + Ifges(2,3) * qJDD(1) + pkin(1) * t197 + pkin(5) * t209 + t190 * t87 + t192 * t85;
t97 = m(2) * t178 + qJDD(1) * mrSges(2,1) - t195 * mrSges(2,2) + t197;
t90 = t190 * t95 + t192 * t94;
t88 = m(2) * t179 - t195 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t209;
t83 = mrSges(2,1) * g(3) + mrSges(2,3) * t179 + t195 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t90 - t221;
t82 = -mrSges(2,2) * g(3) - mrSges(2,3) * t178 + Ifges(2,5) * qJDD(1) - t195 * Ifges(2,6) - pkin(5) * t90 - t190 * t85 + t192 * t87;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t193 * t82 - t191 * t83 - pkin(4) * (t191 * t88 + t193 * t97), t82, t87, t92, -t160 * t129 + t202; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t191 * t82 + t193 * t83 + pkin(4) * (-t191 * t97 + t193 * t88), t83, t85, t91, -t200; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t201, t201, t221, -t198, Ifges(5,5) * t149 + Ifges(5,6) * qJDD(2) + Ifges(5,3) * t148 - qJD(2) * t131 + t161 * t129 - t206;];
m_new = t1;
