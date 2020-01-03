% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:11
% EndTime: 2019-12-31 21:49:14
% DurationCPUTime: 2.18s
% Computational Cost: add. (39065->223), mult. (40655->276), div. (0->0), fcn. (19685->8), ass. (0->87)
t198 = sin(qJ(1));
t202 = cos(qJ(1));
t179 = t198 * g(1) - t202 * g(2);
t175 = qJDD(1) * pkin(1) + t179;
t180 = -t202 * g(1) - t198 * g(2);
t204 = qJD(1) ^ 2;
t176 = -t204 * pkin(1) + t180;
t197 = sin(qJ(2));
t201 = cos(qJ(2));
t141 = t201 * t175 - t197 * t176;
t188 = qJDD(1) + qJDD(2);
t138 = t188 * pkin(2) + t141;
t142 = t197 * t175 + t201 * t176;
t189 = qJD(1) + qJD(2);
t187 = t189 ^ 2;
t139 = -t187 * pkin(2) + t142;
t196 = sin(qJ(3));
t200 = cos(qJ(3));
t134 = t196 * t138 + t200 * t139;
t185 = qJD(3) + t189;
t183 = t185 ^ 2;
t184 = qJDD(3) + t188;
t131 = -t183 * pkin(3) + t184 * pkin(8) + t134;
t195 = sin(qJ(4));
t199 = cos(qJ(4));
t224 = t199 * g(3);
t127 = -t195 * t131 - t224;
t128 = -t195 * g(3) + t199 * t131;
t146 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t195 - Ifges(6,3) * t199) * t185;
t149 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t195 + Ifges(5,2) * t199) * t185;
t162 = (-mrSges(6,1) * t199 - mrSges(6,3) * t195) * t185;
t218 = qJD(4) * t185;
t164 = t195 * t184 + t199 * t218;
t165 = t199 * t184 - t195 * t218;
t161 = (-pkin(4) * t199 - qJ(5) * t195) * t185;
t203 = qJD(4) ^ 2;
t221 = t185 * t199;
t124 = -t203 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t161 * t221 + t128;
t126 = -qJDD(4) * pkin(4) + t224 - t203 * qJ(5) + qJDD(5) + (t161 * t185 + t131) * t195;
t211 = -mrSges(6,1) * t126 + mrSges(6,3) * t124 + Ifges(6,4) * t164 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t165;
t174 = mrSges(6,2) * t221 + qJD(4) * mrSges(6,3);
t213 = -m(6) * t126 + qJDD(4) * mrSges(6,1) + qJD(4) * t174;
t222 = t185 * t195;
t172 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t222;
t214 = m(6) * t124 + qJDD(4) * mrSges(6,3) + qJD(4) * t172 + t162 * t221;
t150 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t195 - Ifges(6,5) * t199) * t185;
t219 = t150 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t195 + Ifges(5,4) * t199) * t185;
t226 = -(t219 * t199 + (t146 - t149) * t195) * t185 + mrSges(5,1) * t127 - mrSges(5,2) * t128 + Ifges(5,5) * t164 + Ifges(5,6) * t165 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t164 * mrSges(6,2) - t162 * t222 + t213) + qJ(5) * (t165 * mrSges(6,2) + t214) + t211;
t223 = mrSges(5,3) + mrSges(6,2);
t163 = (-mrSges(5,1) * t199 + mrSges(5,2) * t195) * t185;
t171 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t222;
t116 = m(5) * t128 - qJDD(4) * mrSges(5,2) - qJD(4) * t171 + t163 * t221 + t223 * t165 + t214;
t173 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t221;
t117 = m(5) * t127 + qJDD(4) * mrSges(5,1) + qJD(4) * t173 + (-t162 - t163) * t222 - t223 * t164 + t213;
t215 = t199 * t116 - t195 * t117;
t107 = m(4) * t134 - t183 * mrSges(4,1) - t184 * mrSges(4,2) + t215;
t133 = t200 * t138 - t196 * t139;
t130 = -t184 * pkin(3) - t183 * pkin(8) - t133;
t121 = -t165 * pkin(4) - t164 * qJ(5) + (-0.2e1 * qJD(5) * t195 + (pkin(4) * t195 - qJ(5) * t199) * qJD(4)) * t185 + t130;
t118 = m(6) * t121 - t165 * mrSges(6,1) - t164 * mrSges(6,3) - t172 * t222 - t174 * t221;
t207 = -m(5) * t130 + t165 * mrSges(5,1) - t164 * mrSges(5,2) - t171 * t222 + t173 * t221 - t118;
t111 = m(4) * t133 + t184 * mrSges(4,1) - t183 * mrSges(4,2) + t207;
t100 = t196 * t107 + t200 * t111;
t97 = m(3) * t141 + t188 * mrSges(3,1) - t187 * mrSges(3,2) + t100;
t216 = t200 * t107 - t196 * t111;
t98 = m(3) * t142 - t187 * mrSges(3,1) - t188 * mrSges(3,2) + t216;
t91 = t197 * t98 + t201 * t97;
t109 = t195 * t116 + t199 * t117;
t217 = -t197 * t97 + t201 * t98;
t212 = -mrSges(6,1) * t121 + mrSges(6,2) * t124;
t147 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t195 + Ifges(5,6) * t199) * t185;
t148 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t195 - Ifges(6,6) * t199) * t185;
t103 = -mrSges(5,1) * t130 + mrSges(5,3) * t128 - pkin(4) * t118 + (-t147 - t148) * t222 + (Ifges(5,2) + Ifges(6,3)) * t165 + (Ifges(5,4) - Ifges(6,5)) * t164 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t219 * qJD(4) + t212;
t209 = mrSges(6,2) * t126 - mrSges(6,3) * t121 + Ifges(6,1) * t164 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t165 + qJD(4) * t146 + t148 * t221;
t104 = mrSges(5,2) * t130 - mrSges(5,3) * t127 + Ifges(5,1) * t164 + Ifges(5,4) * t165 + Ifges(5,5) * qJDD(4) - qJ(5) * t118 - qJD(4) * t149 + t147 * t221 + t209;
t210 = mrSges(4,1) * t133 - mrSges(4,2) * t134 + Ifges(4,3) * t184 + pkin(3) * t207 + pkin(8) * t215 + t199 * t103 + t195 * t104;
t208 = mrSges(3,1) * t141 - mrSges(3,2) * t142 + Ifges(3,3) * t188 + pkin(2) * t100 + t210;
t206 = mrSges(2,1) * t179 - mrSges(2,2) * t180 + Ifges(2,3) * qJDD(1) + pkin(1) * t91 + t208;
t93 = mrSges(4,1) * g(3) + mrSges(4,3) * t134 + t183 * Ifges(4,5) + Ifges(4,6) * t184 - pkin(3) * t109 - t226;
t92 = -mrSges(4,2) * g(3) - mrSges(4,3) * t133 + Ifges(4,5) * t184 - t183 * Ifges(4,6) - pkin(8) * t109 - t195 * t103 + t199 * t104;
t89 = m(2) * t180 - t204 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t217;
t88 = m(2) * t179 + qJDD(1) * mrSges(2,1) - t204 * mrSges(2,2) + t91;
t87 = -mrSges(3,2) * g(3) - mrSges(3,3) * t141 + Ifges(3,5) * t188 - t187 * Ifges(3,6) - pkin(7) * t100 - t196 * t93 + t200 * t92;
t86 = Ifges(3,6) * t188 + t187 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t142 + t196 * t92 + t200 * t93 - pkin(2) * (-m(4) * g(3) + t109) + pkin(7) * t216;
t85 = -mrSges(2,2) * g(3) - mrSges(2,3) * t179 + Ifges(2,5) * qJDD(1) - t204 * Ifges(2,6) - pkin(6) * t91 - t197 * t86 + t201 * t87;
t84 = Ifges(2,6) * qJDD(1) + t204 * Ifges(2,5) + mrSges(2,3) * t180 + t197 * t87 + t201 * t86 - pkin(1) * t109 + pkin(6) * t217 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t202 * t85 - t198 * t84 - pkin(5) * (t198 * t89 + t202 * t88), t85, t87, t92, t104, t209; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t198 * t85 + t202 * t84 + pkin(5) * (-t198 * t88 + t202 * t89), t84, t86, t93, t103, (-t195 * t146 - t199 * t150) * t185 + t211; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t206, t206, t208, t210, t226, Ifges(6,5) * t164 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t165 - qJD(4) * t150 + t148 * t222 - t212;];
m_new = t1;
