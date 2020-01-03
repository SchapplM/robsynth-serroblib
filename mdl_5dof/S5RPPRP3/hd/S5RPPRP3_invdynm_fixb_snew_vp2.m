% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:46
% EndTime: 2019-12-31 17:50:48
% DurationCPUTime: 1.12s
% Computational Cost: add. (10701->220), mult. (19122->263), div. (0->0), fcn. (8105->6), ass. (0->85)
t201 = sin(qJ(1));
t203 = cos(qJ(1));
t178 = t201 * g(1) - t203 * g(2);
t165 = qJDD(1) * pkin(1) + t178;
t179 = -t203 * g(1) - t201 * g(2);
t204 = qJD(1) ^ 2;
t168 = -t204 * pkin(1) + t179;
t198 = sin(pkin(7));
t199 = cos(pkin(7));
t135 = t198 * t165 + t199 * t168;
t218 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t135;
t236 = -pkin(2) - pkin(6);
t128 = t236 * t204 + t218;
t200 = sin(qJ(4));
t202 = cos(qJ(4));
t226 = qJD(1) * qJD(4);
t169 = -t200 * qJDD(1) - t202 * t226;
t170 = t202 * qJDD(1) - t200 * t226;
t228 = qJD(1) * t200;
t174 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t228;
t227 = qJD(1) * t202;
t177 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t227;
t175 = (qJD(4) * pkin(4)) - qJ(5) * t227;
t194 = t200 ^ 2;
t121 = t175 * t227 - t169 * pkin(4) + qJDD(5) + (-qJ(5) * t194 + t236) * t204 + t218;
t173 = -(qJD(4) * mrSges(6,2)) - mrSges(6,3) * t228;
t176 = (qJD(4) * mrSges(6,1)) - mrSges(6,3) * t227;
t219 = m(6) * t121 + t170 * mrSges(6,2) + t173 * t228 + t176 * t227;
t237 = -m(5) * t128 - t170 * mrSges(5,2) + (mrSges(5,1) + mrSges(6,1)) * t169 - t174 * t228 - t177 * t227 - t219;
t235 = pkin(4) * t204;
t234 = mrSges(3,1) - mrSges(4,2);
t232 = Ifges(3,5) - Ifges(4,4);
t231 = -Ifges(3,6) + Ifges(4,5);
t130 = t204 * pkin(2) - t218;
t209 = -m(4) * t130 + t204 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t237;
t105 = m(3) * t135 - t204 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t209;
t134 = t199 * t165 - t198 * t168;
t216 = -t204 * qJ(3) + qJDD(3) - t134;
t129 = t236 * qJDD(1) + t216;
t126 = t202 * t129;
t195 = -g(3) + qJDD(2);
t123 = -t200 * t195 + t126;
t166 = (mrSges(6,1) * t200 + mrSges(6,2) * t202) * qJD(1);
t220 = qJD(1) * (-t166 - (mrSges(5,1) * t200 + mrSges(5,2) * t202) * qJD(1));
t224 = -2 * qJD(1) * qJD(5);
t117 = t202 * t224 + (qJDD(4) * pkin(4)) - t170 * qJ(5) + t126 + (-qJ(5) * t226 - t202 * t235 - t195) * t200;
t223 = m(6) * t117 + (qJDD(4) * mrSges(6,1)) + qJD(4) * t173;
t108 = m(5) * t123 + (qJDD(4) * mrSges(5,1)) + qJD(4) * t174 + (-mrSges(5,3) - mrSges(6,3)) * t170 + t202 * t220 + t223;
t124 = t200 * t129 + t202 * t195;
t118 = t169 * qJ(5) - qJD(4) * t175 - t194 * t235 + t200 * t224 + t124;
t230 = m(6) * t118 + t169 * mrSges(6,3);
t110 = m(5) * t124 + t169 * mrSges(5,3) + ((-mrSges(5,2) - mrSges(6,2)) * qJDD(4)) + (-t176 - t177) * qJD(4) + t200 * t220 + t230;
t101 = t202 * t108 + t200 * t110;
t132 = -qJDD(1) * pkin(2) + t216;
t213 = -m(4) * t132 + t204 * mrSges(4,3) - t101;
t96 = m(3) * t134 - t204 * mrSges(3,2) + t234 * qJDD(1) + t213;
t93 = t198 * t105 + t199 * t96;
t144 = (Ifges(6,3) * qJD(4)) + (Ifges(6,5) * t202 - Ifges(6,6) * t200) * qJD(1);
t229 = -t144 - (Ifges(5,3) * qJD(4)) - (Ifges(5,5) * t202 - Ifges(5,6) * t200) * qJD(1);
t221 = t199 * t105 - t198 * t96;
t102 = -t200 * t108 + t202 * t110;
t100 = m(4) * t195 + t102;
t217 = mrSges(6,2) * t121 - mrSges(6,3) * t117 + Ifges(6,1) * t170 + Ifges(6,4) * t169 + Ifges(6,5) * qJDD(4);
t148 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t202 - Ifges(6,4) * t200) * qJD(1);
t215 = -mrSges(6,1) * t121 + mrSges(6,3) * t118 + Ifges(6,4) * t170 + Ifges(6,2) * t169 + Ifges(6,6) * qJDD(4) + qJD(4) * t148;
t146 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t202 - Ifges(6,2) * t200) * qJD(1);
t212 = mrSges(6,1) * t117 - mrSges(6,2) * t118 + Ifges(6,5) * t170 + Ifges(6,6) * t169 + (Ifges(6,3) * qJDD(4)) + t146 * t227 + t148 * t228;
t149 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t202 - Ifges(5,4) * t200) * qJD(1);
t94 = Ifges(5,4) * t170 + Ifges(5,2) * t169 + Ifges(5,6) * qJDD(4) + qJD(4) * t149 - mrSges(5,1) * t128 + mrSges(5,3) * t124 - pkin(4) * (-t169 * mrSges(6,1) + t219) + qJ(5) * (-(qJDD(4) * mrSges(6,2)) - qJD(4) * t176 + t230) + (-qJ(5) * t166 * t200 + t229 * t202) * qJD(1) + t215;
t112 = -t170 * mrSges(6,3) - t166 * t227 + t223;
t147 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t202 - Ifges(5,2) * t200) * qJD(1);
t98 = mrSges(5,2) * t128 - mrSges(5,3) * t123 + Ifges(5,1) * t170 + Ifges(5,4) * t169 + Ifges(5,5) * qJDD(4) - qJ(5) * t112 + (-t146 - t147) * qJD(4) + t229 * t228 + t217;
t211 = mrSges(4,2) * t132 - mrSges(4,3) * t130 + Ifges(4,1) * qJDD(1) - pkin(6) * t101 - t200 * t94 + t202 * t98;
t210 = -mrSges(4,1) * t130 - pkin(3) * t237 - pkin(6) * t102 - t200 * t98 - t202 * t94;
t208 = -mrSges(3,2) * t135 + qJ(3) * t209 + mrSges(3,1) * t134 + Ifges(3,3) * qJDD(1) + t211 + pkin(2) * (-qJDD(1) * mrSges(4,2) + t213);
t207 = mrSges(5,1) * t123 - mrSges(5,2) * t124 + Ifges(5,5) * t170 + Ifges(5,6) * t169 + (Ifges(5,3) * qJDD(4)) + pkin(4) * t112 + t147 * t227 + t149 * t228 + t212;
t206 = mrSges(2,1) * t178 - mrSges(2,2) * t179 + Ifges(2,3) * qJDD(1) + pkin(1) * t93 + t208;
t205 = mrSges(4,1) * t132 + pkin(3) * t101 + t207;
t91 = m(2) * t179 - t204 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t221;
t90 = m(2) * t178 + qJDD(1) * mrSges(2,1) - t204 * mrSges(2,2) + t93;
t89 = t205 + t231 * t204 + t232 * qJDD(1) + (mrSges(3,2) - mrSges(4,3)) * t195 - mrSges(3,3) * t134 - qJ(3) * t100;
t88 = mrSges(3,3) * t135 - pkin(2) * t100 - t231 * qJDD(1) - t234 * t195 + t232 * t204 + t210;
t87 = -mrSges(2,2) * g(3) - mrSges(2,3) * t178 + Ifges(2,5) * qJDD(1) - t204 * Ifges(2,6) - qJ(2) * t93 - t198 * t88 + t199 * t89;
t86 = Ifges(2,6) * qJDD(1) + t204 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t179 + t198 * t89 + t199 * t88 - pkin(1) * (m(3) * t195 + t100) + qJ(2) * t221;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t203 * t87 - t201 * t86 - pkin(5) * (t201 * t91 + t203 * t90), t87, t89, t211, t98, -qJD(4) * t146 - t144 * t228 + t217; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t201 * t87 + t203 * t86 + pkin(5) * (-t201 * t90 + t203 * t91), t86, t88, mrSges(4,3) * t195 + Ifges(4,4) * qJDD(1) - t204 * Ifges(4,5) - t205, t94, -t144 * t227 + t215; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t206, t206, t208, -mrSges(4,2) * t195 + t204 * Ifges(4,4) + Ifges(4,5) * qJDD(1) - t210, t207, t212;];
m_new = t1;
