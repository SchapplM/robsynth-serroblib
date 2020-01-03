% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR10
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR10_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:12
% EndTime: 2019-12-31 17:11:15
% DurationCPUTime: 1.38s
% Computational Cost: add. (10189->241), mult. (20790->291), div. (0->0), fcn. (10237->6), ass. (0->91)
t191 = sin(qJ(1));
t194 = cos(qJ(1));
t178 = -t194 * g(1) - t191 * g(2);
t196 = qJD(1) ^ 2;
t153 = -t196 * pkin(1) + qJDD(1) * pkin(5) + t178;
t190 = sin(qJ(2));
t193 = cos(qJ(2));
t139 = -t193 * g(3) - t190 * t153;
t140 = -t190 * g(3) + t193 * t153;
t148 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t190 + Ifges(3,4) * t193) * qJD(1);
t166 = (mrSges(4,2) * t193 - mrSges(4,3) * t190) * qJD(1);
t213 = qJD(1) * qJD(2);
t211 = t193 * t213;
t168 = t190 * qJDD(1) + t211;
t212 = t190 * t213;
t169 = t193 * qJDD(1) - t212;
t215 = qJD(1) * t193;
t174 = -mrSges(4,1) * t215 - qJD(2) * mrSges(4,3);
t214 = t190 * qJD(1);
t176 = pkin(3) * t214 - qJD(2) * pkin(6);
t188 = t193 ^ 2;
t177 = t191 * g(1) - t194 * g(2);
t209 = -qJDD(1) * pkin(1) - t177;
t221 = -2 * qJD(3);
t200 = pkin(2) * t212 + t214 * t221 + (-t168 - t211) * qJ(3) + t209;
t118 = -t176 * t214 + (-pkin(3) * t188 - pkin(5)) * t196 + (-pkin(2) - pkin(6)) * t169 + t200;
t165 = (-pkin(2) * t193 - qJ(3) * t190) * qJD(1);
t195 = qJD(2) ^ 2;
t126 = -qJDD(2) * pkin(2) - t195 * qJ(3) + t165 * t214 + qJDD(3) - t139;
t121 = (-t190 * t193 * t196 - qJDD(2)) * pkin(6) + (t168 - t211) * pkin(3) + t126;
t189 = sin(qJ(4));
t192 = cos(qJ(4));
t116 = -t189 * t118 + t192 * t121;
t163 = -t189 * qJD(2) - t192 * t215;
t135 = t163 * qJD(4) + t192 * qJDD(2) - t189 * t169;
t164 = t192 * qJD(2) - t189 * t215;
t136 = -t163 * mrSges(5,1) + t164 * mrSges(5,2);
t180 = qJD(4) + t214;
t137 = -t180 * mrSges(5,2) + t163 * mrSges(5,3);
t162 = qJDD(4) + t168;
t112 = m(5) * t116 + t162 * mrSges(5,1) - t135 * mrSges(5,3) - t164 * t136 + t180 * t137;
t117 = t192 * t118 + t189 * t121;
t134 = -t164 * qJD(4) - t189 * qJDD(2) - t192 * t169;
t138 = t180 * mrSges(5,1) - t164 * mrSges(5,3);
t113 = m(5) * t117 - t162 * mrSges(5,2) + t134 * mrSges(5,3) + t163 * t136 - t180 * t138;
t102 = t192 * t112 + t189 * t113;
t204 = -t195 * pkin(2) + qJDD(2) * qJ(3) + t165 * t215 + t140;
t120 = -t188 * t196 * pkin(6) + t169 * pkin(3) + ((2 * qJD(3)) + t176) * qJD(2) + t204;
t127 = Ifges(5,5) * t164 + Ifges(5,6) * t163 + Ifges(5,3) * t180;
t129 = Ifges(5,1) * t164 + Ifges(5,4) * t163 + Ifges(5,5) * t180;
t105 = -mrSges(5,1) * t120 + mrSges(5,3) * t117 + Ifges(5,4) * t135 + Ifges(5,2) * t134 + Ifges(5,6) * t162 - t164 * t127 + t180 * t129;
t128 = Ifges(5,4) * t164 + Ifges(5,2) * t163 + Ifges(5,6) * t180;
t106 = mrSges(5,2) * t120 - mrSges(5,3) * t116 + Ifges(5,1) * t135 + Ifges(5,4) * t134 + Ifges(5,5) * t162 + t163 * t127 - t180 * t128;
t124 = qJD(2) * t221 - t204;
t150 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t190 - Ifges(4,6) * t193) * qJD(1);
t202 = -mrSges(4,2) * t126 + mrSges(4,3) * t124 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t168 + Ifges(4,5) * t169 + pkin(6) * t102 + t189 * t105 - t192 * t106 - t150 * t215;
t114 = -m(5) * t120 + t134 * mrSges(5,1) - t135 * mrSges(5,2) + t163 * t137 - t164 * t138;
t175 = mrSges(4,1) * t214 + qJD(2) * mrSges(4,2);
t203 = -m(4) * t124 + qJDD(2) * mrSges(4,3) + qJD(2) * t175 + t166 * t215 - t114;
t206 = -m(4) * t126 - t168 * mrSges(4,1) - t102;
t149 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t190 - Ifges(4,3) * t193) * qJD(1);
t216 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t190 + Ifges(3,2) * t193) * qJD(1) - t149;
t223 = (-t193 * t148 + t216 * t190) * qJD(1) + mrSges(3,1) * t139 - mrSges(3,2) * t140 + Ifges(3,5) * t168 + Ifges(3,6) * t169 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t174 - t166 * t214 + t206) + qJ(3) * (t169 * mrSges(4,1) + t203) - t202;
t220 = t196 * pkin(5);
t219 = Ifges(3,4) + Ifges(4,6);
t103 = -t189 * t112 + t192 * t113;
t151 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t190 - Ifges(4,5) * t193) * qJD(1);
t217 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t190 + Ifges(3,6) * t193) * qJD(1) + t151;
t167 = (-mrSges(3,1) * t193 + mrSges(3,2) * t190) * qJD(1);
t172 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t214;
t108 = t167 * t215 + m(3) * t140 - qJDD(2) * mrSges(3,2) - qJD(2) * t172 + (mrSges(3,3) + mrSges(4,1)) * t169 + t203;
t173 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t215;
t99 = m(3) * t139 - t168 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t173 - t174) * qJD(2) + (-t166 - t167) * t214 + t206;
t210 = t193 * t108 - t190 * t99;
t122 = -t169 * pkin(2) + t200 - t220;
t208 = -m(4) * t122 - t169 * mrSges(4,2) + t175 * t214 - t103;
t152 = t209 - t220;
t198 = -m(3) * t152 + t173 * t215 + t169 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t168 + (-t172 * t190 - t174 * t193) * qJD(1) + t208;
t100 = -t168 * mrSges(4,3) + t174 * t215 - t208;
t201 = -mrSges(4,1) * t124 + mrSges(4,2) * t122 - pkin(3) * t114 - pkin(6) * t103 - t192 * t105 - t189 * t106;
t91 = -mrSges(3,1) * t152 + mrSges(3,3) * t140 - pkin(2) * t100 + (Ifges(3,2) + Ifges(4,3)) * t169 + t219 * t168 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t148 - t150) * qJD(2) - t217 * t214 + t201;
t205 = mrSges(5,1) * t116 - mrSges(5,2) * t117 + Ifges(5,5) * t135 + Ifges(5,6) * t134 + Ifges(5,3) * t162 + t164 * t128 - t163 * t129;
t199 = mrSges(4,1) * t126 - mrSges(4,3) * t122 + pkin(3) * t102 + t205;
t93 = t217 * t215 + t199 + t219 * t169 + (Ifges(3,1) + Ifges(4,2)) * t168 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t216 * qJD(2) + mrSges(3,2) * t152 - mrSges(3,3) * t139 - qJ(3) * t100;
t207 = mrSges(2,1) * t177 - mrSges(2,2) * t178 + Ifges(2,3) * qJDD(1) + pkin(1) * t198 + pkin(5) * t210 + t190 * t93 + t193 * t91;
t97 = m(2) * t177 + qJDD(1) * mrSges(2,1) - t196 * mrSges(2,2) + t198;
t96 = t190 * t108 + t193 * t99;
t94 = m(2) * t178 - t196 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t210;
t89 = mrSges(2,1) * g(3) + mrSges(2,3) * t178 + t196 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t96 - t223;
t88 = -mrSges(2,2) * g(3) - mrSges(2,3) * t177 + Ifges(2,5) * qJDD(1) - t196 * Ifges(2,6) - pkin(5) * t96 - t190 * t91 + t193 * t93;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t194 * t88 - t191 * t89 - pkin(4) * (t191 * t94 + t194 * t97), t88, t93, -t149 * t214 - t202, t106; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t191 * t88 + t194 * t89 + pkin(4) * (-t191 * t97 + t194 * t94), t89, t91, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t168 - Ifges(4,6) * t169 - qJD(2) * t149 - t151 * t215 - t199, t105; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t207, t207, t223, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t168 - Ifges(4,3) * t169 + qJD(2) * t150 + t151 * t214 - t201, t205;];
m_new = t1;
