% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:45
% EndTime: 2019-12-31 18:40:48
% DurationCPUTime: 2.07s
% Computational Cost: add. (30104->221), mult. (40655->275), div. (0->0), fcn. (19685->8), ass. (0->86)
t199 = sin(qJ(1));
t202 = cos(qJ(1));
t178 = t199 * g(1) - t202 * g(2);
t171 = qJDD(1) * pkin(1) + t178;
t179 = -t202 * g(1) - t199 * g(2);
t204 = qJD(1) ^ 2;
t172 = -t204 * pkin(1) + t179;
t195 = sin(pkin(8));
t196 = cos(pkin(8));
t141 = t196 * t171 - t195 * t172;
t138 = qJDD(1) * pkin(2) + t141;
t142 = t195 * t171 + t196 * t172;
t139 = -t204 * pkin(2) + t142;
t198 = sin(qJ(3));
t201 = cos(qJ(3));
t134 = t198 * t138 + t201 * t139;
t187 = qJD(1) + qJD(3);
t185 = t187 ^ 2;
t186 = qJDD(1) + qJDD(3);
t131 = -t185 * pkin(3) + t186 * pkin(7) + t134;
t197 = sin(qJ(4));
t194 = -g(3) + qJDD(2);
t200 = cos(qJ(4));
t222 = t200 * t194;
t127 = -t197 * t131 + t222;
t128 = t200 * t131 + t197 * t194;
t146 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t197 - Ifges(6,3) * t200) * t187;
t149 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t197 + Ifges(5,2) * t200) * t187;
t162 = (-mrSges(6,1) * t200 - mrSges(6,3) * t197) * t187;
t219 = qJD(4) * t187;
t164 = t197 * t186 + t200 * t219;
t165 = t200 * t186 - t197 * t219;
t161 = (-pkin(4) * t200 - qJ(5) * t197) * t187;
t203 = qJD(4) ^ 2;
t223 = t187 * t200;
t124 = -t203 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t161 * t223 + t128;
t126 = -qJDD(4) * pkin(4) - t203 * qJ(5) - t222 + qJDD(5) + (t161 * t187 + t131) * t197;
t211 = -mrSges(6,1) * t126 + mrSges(6,3) * t124 + Ifges(6,4) * t164 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t165;
t176 = mrSges(6,2) * t223 + qJD(4) * mrSges(6,3);
t213 = -m(6) * t126 + qJDD(4) * mrSges(6,1) + qJD(4) * t176;
t224 = t187 * t197;
t174 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t224;
t214 = m(6) * t124 + qJDD(4) * mrSges(6,3) + qJD(4) * t174 + t162 * t223;
t150 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t197 - Ifges(6,5) * t200) * t187;
t220 = t150 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t197 + Ifges(5,4) * t200) * t187;
t227 = -((t146 - t149) * t197 + t220 * t200) * t187 + mrSges(5,1) * t127 - mrSges(5,2) * t128 + Ifges(5,5) * t164 + Ifges(5,6) * t165 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t164 * mrSges(6,2) - t162 * t224 + t213) + qJ(5) * (t165 * mrSges(6,2) + t214) + t211;
t225 = mrSges(5,3) + mrSges(6,2);
t163 = (-mrSges(5,1) * t200 + mrSges(5,2) * t197) * t187;
t173 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t224;
t116 = m(5) * t128 - qJDD(4) * mrSges(5,2) - qJD(4) * t173 + t163 * t223 + t225 * t165 + t214;
t175 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t223;
t117 = m(5) * t127 + qJDD(4) * mrSges(5,1) + qJD(4) * t175 + (-t162 - t163) * t224 - t225 * t164 + t213;
t215 = t200 * t116 - t197 * t117;
t107 = m(4) * t134 - t185 * mrSges(4,1) - t186 * mrSges(4,2) + t215;
t133 = t201 * t138 - t198 * t139;
t130 = -t186 * pkin(3) - t185 * pkin(7) - t133;
t121 = -t165 * pkin(4) - t164 * qJ(5) + (-0.2e1 * qJD(5) * t197 + (pkin(4) * t197 - qJ(5) * t200) * qJD(4)) * t187 + t130;
t118 = m(6) * t121 - t165 * mrSges(6,1) - t164 * mrSges(6,3) - t174 * t224 - t176 * t223;
t207 = -m(5) * t130 + t165 * mrSges(5,1) - t164 * mrSges(5,2) - t173 * t224 + t175 * t223 - t118;
t111 = m(4) * t133 + t186 * mrSges(4,1) - t185 * mrSges(4,2) + t207;
t100 = t198 * t107 + t201 * t111;
t97 = m(3) * t141 + qJDD(1) * mrSges(3,1) - t204 * mrSges(3,2) + t100;
t216 = t201 * t107 - t198 * t111;
t98 = m(3) * t142 - t204 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t216;
t91 = t195 * t98 + t196 * t97;
t109 = t197 * t116 + t200 * t117;
t218 = m(4) * t194 + t109;
t217 = -t195 * t97 + t196 * t98;
t212 = -mrSges(6,1) * t121 + mrSges(6,2) * t124;
t147 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t197 + Ifges(5,6) * t200) * t187;
t148 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t197 - Ifges(6,6) * t200) * t187;
t103 = -mrSges(5,1) * t130 + mrSges(5,3) * t128 - pkin(4) * t118 + (-t147 - t148) * t224 + (Ifges(5,2) + Ifges(6,3)) * t165 + (Ifges(5,4) - Ifges(6,5)) * t164 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t220 * qJD(4) + t212;
t209 = mrSges(6,2) * t126 - mrSges(6,3) * t121 + Ifges(6,1) * t164 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t165 + qJD(4) * t146 + t148 * t223;
t104 = mrSges(5,2) * t130 - mrSges(5,3) * t127 + Ifges(5,1) * t164 + Ifges(5,4) * t165 + Ifges(5,5) * qJDD(4) - qJ(5) * t118 - qJD(4) * t149 + t147 * t223 + t209;
t210 = mrSges(4,1) * t133 - mrSges(4,2) * t134 + Ifges(4,3) * t186 + pkin(3) * t207 + pkin(7) * t215 + t200 * t103 + t197 * t104;
t208 = mrSges(3,1) * t141 - mrSges(3,2) * t142 + Ifges(3,3) * qJDD(1) + pkin(2) * t100 + t210;
t206 = mrSges(2,1) * t178 - mrSges(2,2) * t179 + Ifges(2,3) * qJDD(1) + pkin(1) * t91 + t208;
t93 = -mrSges(4,1) * t194 + mrSges(4,3) * t134 + t185 * Ifges(4,5) + Ifges(4,6) * t186 - pkin(3) * t109 - t227;
t92 = mrSges(4,2) * t194 - mrSges(4,3) * t133 + Ifges(4,5) * t186 - t185 * Ifges(4,6) - pkin(7) * t109 - t197 * t103 + t200 * t104;
t89 = m(2) * t179 - t204 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t217;
t88 = m(2) * t178 + qJDD(1) * mrSges(2,1) - t204 * mrSges(2,2) + t91;
t87 = mrSges(3,2) * t194 - mrSges(3,3) * t141 + Ifges(3,5) * qJDD(1) - t204 * Ifges(3,6) - pkin(6) * t100 - t198 * t93 + t201 * t92;
t86 = -mrSges(3,1) * t194 + mrSges(3,3) * t142 + t204 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t218 + pkin(6) * t216 + t198 * t92 + t201 * t93;
t85 = -mrSges(2,2) * g(3) - mrSges(2,3) * t178 + Ifges(2,5) * qJDD(1) - t204 * Ifges(2,6) - qJ(2) * t91 - t195 * t86 + t196 * t87;
t84 = Ifges(2,6) * qJDD(1) + t204 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t179 + t195 * t87 + t196 * t86 - pkin(1) * (m(3) * t194 + t218) + qJ(2) * t217;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t202 * t85 - t199 * t84 - pkin(5) * (t199 * t89 + t202 * t88), t85, t87, t92, t104, t209; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t199 * t85 + t202 * t84 + pkin(5) * (-t199 * t88 + t202 * t89), t84, t86, t93, t103, (-t197 * t146 - t200 * t150) * t187 + t211; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t206, t206, t208, t210, t227, Ifges(6,5) * t164 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t165 - qJD(4) * t150 + t148 * t224 - t212;];
m_new = t1;
