% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR8
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR8_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:50
% EndTime: 2019-12-31 17:07:53
% DurationCPUTime: 1.39s
% Computational Cost: add. (10176->237), mult. (21113->298), div. (0->0), fcn. (10664->6), ass. (0->89)
t209 = sin(qJ(1));
t212 = cos(qJ(1));
t191 = -t212 * g(1) - t209 * g(2);
t214 = qJD(1) ^ 2;
t167 = -t214 * pkin(1) + qJDD(1) * pkin(5) + t191;
t208 = sin(qJ(2));
t211 = cos(qJ(2));
t147 = -t211 * g(3) - t208 * t167;
t148 = -t208 * g(3) + t211 * t167;
t158 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t208 - Ifges(4,3) * t211) * qJD(1);
t161 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t208 + Ifges(3,2) * t211) * qJD(1);
t178 = (-mrSges(4,1) * t211 - mrSges(4,3) * t208) * qJD(1);
t228 = qJD(1) * qJD(2);
t227 = t211 * t228;
t180 = t208 * qJDD(1) + t227;
t181 = t211 * qJDD(1) - t208 * t228;
t177 = (-pkin(2) * t211 - qJ(3) * t208) * qJD(1);
t213 = qJD(2) ^ 2;
t229 = qJD(1) * t211;
t236 = 2 * qJD(3);
t132 = -t213 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t236 + t177 * t229 + t148;
t230 = qJD(1) * t208;
t189 = -qJD(2) * pkin(3) - pkin(6) * t230;
t206 = t211 ^ 2;
t125 = -t206 * t214 * pkin(3) - t181 * pkin(6) + qJD(2) * t189 + t132;
t138 = -qJDD(2) * pkin(2) - t213 * qJ(3) + t177 * t230 + qJDD(3) - t147;
t126 = (-t180 + t227) * pkin(6) + (-t208 * t211 * t214 - qJDD(2)) * pkin(3) + t138;
t207 = sin(qJ(4));
t210 = cos(qJ(4));
t121 = -t207 * t125 + t210 * t126;
t164 = (-t207 * t208 - t210 * t211) * qJD(1);
t140 = t164 * qJD(4) + t210 * t180 - t207 * t181;
t165 = (-t207 * t211 + t208 * t210) * qJD(1);
t146 = -t164 * mrSges(5,1) + t165 * mrSges(5,2);
t199 = -qJD(2) + qJD(4);
t149 = -t199 * mrSges(5,2) + t164 * mrSges(5,3);
t198 = -qJDD(2) + qJDD(4);
t116 = m(5) * t121 + t198 * mrSges(5,1) - t140 * mrSges(5,3) - t165 * t146 + t199 * t149;
t122 = t210 * t125 + t207 * t126;
t139 = -t165 * qJD(4) - t207 * t180 - t210 * t181;
t150 = t199 * mrSges(5,1) - t165 * mrSges(5,3);
t117 = m(5) * t122 - t198 * mrSges(5,2) + t139 * mrSges(5,3) + t164 * t146 - t199 * t150;
t107 = t210 * t116 + t207 * t117;
t142 = Ifges(5,4) * t165 + Ifges(5,2) * t164 + Ifges(5,6) * t199;
t143 = Ifges(5,1) * t165 + Ifges(5,4) * t164 + Ifges(5,5) * t199;
t224 = mrSges(5,1) * t121 - mrSges(5,2) * t122 + Ifges(5,5) * t140 + Ifges(5,6) * t139 + Ifges(5,3) * t198 + t165 * t142 - t164 * t143;
t217 = -mrSges(4,1) * t138 + mrSges(4,3) * t132 + Ifges(4,4) * t180 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t181 - pkin(3) * t107 - t224;
t188 = mrSges(4,2) * t229 + qJD(2) * mrSges(4,3);
t220 = -m(4) * t138 + qJDD(2) * mrSges(4,1) + qJD(2) * t188 - t107;
t108 = -t207 * t116 + t210 * t117;
t186 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t230;
t222 = m(4) * t132 + qJDD(2) * mrSges(4,3) + qJD(2) * t186 + t178 * t229 + t108;
t162 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t208 - Ifges(4,5) * t211) * qJD(1);
t231 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t208 + Ifges(3,4) * t211) * qJD(1) + t162;
t238 = -(t231 * t211 + (t158 - t161) * t208) * qJD(1) + mrSges(3,1) * t147 - mrSges(3,2) * t148 + Ifges(3,5) * t180 + Ifges(3,6) * t181 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t180 * mrSges(4,2) - t178 * t230 + t220) + qJ(3) * (t181 * mrSges(4,2) + t222) + t217;
t235 = t214 * pkin(5);
t234 = mrSges(3,3) + mrSges(4,2);
t233 = qJ(3) * t211;
t190 = t209 * g(1) - t212 * g(2);
t179 = (-mrSges(3,1) * t211 + mrSges(3,2) * t208) * qJD(1);
t185 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t230;
t103 = m(3) * t148 - qJDD(2) * mrSges(3,2) - qJD(2) * t185 + t179 * t229 + t234 * t181 + t222;
t187 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t229;
t104 = m(3) * t147 + qJDD(2) * mrSges(3,1) + qJD(2) * t187 - t234 * t180 + (-t178 - t179) * t230 + t220;
t226 = t211 * t103 - t208 * t104;
t225 = qJDD(1) * pkin(1) + t190;
t223 = -t180 * qJ(3) - t225;
t124 = (-pkin(6) * t206 + pkin(5)) * t214 + (pkin(2) + pkin(3)) * t181 + (qJD(2) * t233 + (-pkin(2) * qJD(2) + t189 + t236) * t208) * qJD(1) - t223;
t118 = -m(5) * t124 + t139 * mrSges(5,1) - t140 * mrSges(5,2) + t164 * t149 - t165 * t150;
t127 = -t181 * pkin(2) - t235 + (-0.2e1 * qJD(3) * t208 + (pkin(2) * t208 - t233) * qJD(2)) * qJD(1) + t223;
t114 = m(4) * t127 - t181 * mrSges(4,1) - t180 * mrSges(4,3) - t186 * t230 - t188 * t229 + t118;
t166 = -t225 - t235;
t216 = -m(3) * t166 + t181 * mrSges(3,1) - t180 * mrSges(3,2) - t185 * t230 + t187 * t229 - t114;
t159 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t208 + Ifges(3,6) * t211) * qJD(1);
t160 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t208 - Ifges(4,6) * t211) * qJD(1);
t141 = Ifges(5,5) * t165 + Ifges(5,6) * t164 + Ifges(5,3) * t199;
t111 = -mrSges(5,1) * t124 + mrSges(5,3) * t122 + Ifges(5,4) * t140 + Ifges(5,2) * t139 + Ifges(5,6) * t198 - t165 * t141 + t199 * t143;
t112 = mrSges(5,2) * t124 - mrSges(5,3) * t121 + Ifges(5,1) * t140 + Ifges(5,4) * t139 + Ifges(5,5) * t198 + t164 * t141 - t199 * t142;
t218 = -mrSges(4,1) * t127 + mrSges(4,2) * t132 - pkin(3) * t118 - pkin(6) * t108 - t210 * t111 - t207 * t112;
t96 = -mrSges(3,1) * t166 + mrSges(3,3) * t148 - pkin(2) * t114 + (Ifges(3,2) + Ifges(4,3)) * t181 + (Ifges(3,4) - Ifges(4,5)) * t180 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t231 * qJD(2) + (-t159 - t160) * t230 + t218;
t219 = mrSges(4,2) * t138 - mrSges(4,3) * t127 + Ifges(4,1) * t180 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t181 - pkin(6) * t107 + qJD(2) * t158 - t207 * t111 + t210 * t112 + t160 * t229;
t98 = mrSges(3,2) * t166 - mrSges(3,3) * t147 + Ifges(3,1) * t180 + Ifges(3,4) * t181 + Ifges(3,5) * qJDD(2) - qJ(3) * t114 - qJD(2) * t161 + t159 * t229 + t219;
t221 = mrSges(2,1) * t190 - mrSges(2,2) * t191 + Ifges(2,3) * qJDD(1) + pkin(1) * t216 + pkin(5) * t226 + t208 * t98 + t211 * t96;
t109 = m(2) * t190 + qJDD(1) * mrSges(2,1) - t214 * mrSges(2,2) + t216;
t101 = t208 * t103 + t211 * t104;
t99 = m(2) * t191 - t214 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t226;
t94 = mrSges(2,1) * g(3) + mrSges(2,3) * t191 + t214 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t101 - t238;
t93 = -mrSges(2,2) * g(3) - mrSges(2,3) * t190 + Ifges(2,5) * qJDD(1) - t214 * Ifges(2,6) - pkin(5) * t101 - t208 * t96 + t211 * t98;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t212 * t93 - t209 * t94 - pkin(4) * (t212 * t109 + t209 * t99), t93, t98, t219, t112; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t209 * t93 + t212 * t94 + pkin(4) * (-t209 * t109 + t212 * t99), t94, t96, t217 + (-t208 * t158 - t211 * t162) * qJD(1), t111; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t221, t221, t238, Ifges(4,5) * t180 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t181 - qJD(2) * t162 + t160 * t230 - t218, t224;];
m_new = t1;
