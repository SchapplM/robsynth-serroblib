% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:53
% EndTime: 2019-12-05 16:39:58
% DurationCPUTime: 2.01s
% Computational Cost: add. (27185->210), mult. (36497->261), div. (0->0), fcn. (20254->8), ass. (0->89)
t189 = qJD(2) + qJD(3);
t199 = sin(qJ(4));
t202 = cos(qJ(4));
t167 = (-mrSges(6,1) * t202 + mrSges(6,2) * t199) * t189;
t188 = qJDD(2) + qJDD(3);
t224 = qJD(4) * t189;
t219 = t202 * t224;
t169 = t199 * t188 + t219;
t197 = sin(pkin(8));
t198 = cos(pkin(8));
t181 = t197 * g(1) - t198 * g(2);
t182 = -t198 * g(1) - t197 * g(2);
t201 = sin(qJ(2));
t204 = cos(qJ(2));
t144 = t204 * t181 - t201 * t182;
t141 = qJDD(2) * pkin(2) + t144;
t145 = t201 * t181 + t204 * t182;
t205 = qJD(2) ^ 2;
t142 = -t205 * pkin(2) + t145;
t200 = sin(qJ(3));
t203 = cos(qJ(3));
t137 = t200 * t141 + t203 * t142;
t187 = t189 ^ 2;
t134 = -t187 * pkin(3) + t188 * pkin(7) + t137;
t196 = -g(3) + qJDD(1);
t184 = t202 * t196;
t223 = qJD(5) * t189;
t231 = pkin(4) * t187;
t126 = qJDD(4) * pkin(4) + t184 + (-t169 + t219) * qJ(5) + (t202 * t231 - t134 - 0.2e1 * t223) * t199;
t228 = t189 * t202;
t178 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t228;
t222 = m(6) * t126 + qJDD(4) * mrSges(6,1) + qJD(4) * t178;
t229 = t189 * t199;
t121 = -t169 * mrSges(6,3) - t167 * t229 + t222;
t130 = -t199 * t134 + t184;
t131 = t202 * t134 + t199 * t196;
t153 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t199 + Ifges(5,2) * t202) * t189;
t154 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t199 + Ifges(6,4) * t202) * t189;
t155 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t199 + Ifges(5,4) * t202) * t189;
t170 = t202 * t188 - t199 * t224;
t175 = qJD(4) * pkin(4) - qJ(5) * t229;
t195 = t202 ^ 2;
t127 = t170 * qJ(5) - qJD(4) * t175 - t195 * t231 + 0.2e1 * t202 * t223 + t131;
t152 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t199 + Ifges(6,2) * t202) * t189;
t212 = -mrSges(6,1) * t126 + mrSges(6,2) * t127 - Ifges(6,5) * t169 - Ifges(6,6) * t170 - Ifges(6,3) * qJDD(4) - t152 * t229;
t233 = mrSges(5,1) * t130 - mrSges(5,2) * t131 + Ifges(5,5) * t169 + Ifges(5,6) * t170 + Ifges(5,3) * qJDD(4) + pkin(4) * t121 - (-t199 * t153 + (t154 + t155) * t202) * t189 - t212;
t230 = -mrSges(5,2) - mrSges(6,2);
t168 = (-mrSges(5,1) * t202 + mrSges(5,2) * t199) * t189;
t179 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t228;
t119 = m(5) * t130 + qJDD(4) * mrSges(5,1) + qJD(4) * t179 + (-t167 - t168) * t229 + (-mrSges(5,3) - mrSges(6,3)) * t169 + t222;
t221 = m(6) * t127 + t170 * mrSges(6,3) + t167 * t228;
t176 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t229;
t225 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t229 - t176;
t120 = m(5) * t131 + t170 * mrSges(5,3) + t225 * qJD(4) + t230 * qJDD(4) + t168 * t228 + t221;
t216 = -t199 * t119 + t202 * t120;
t110 = m(4) * t137 - t187 * mrSges(4,1) - t188 * mrSges(4,2) + t216;
t136 = t203 * t141 - t200 * t142;
t214 = -t188 * pkin(3) - t136;
t133 = -t187 * pkin(7) + t214;
t129 = t175 * t229 - t170 * pkin(4) + qJDD(5) + (-qJ(5) * t195 - pkin(7)) * t187 + t214;
t215 = -m(6) * t129 + t170 * mrSges(6,1) + t178 * t228;
t208 = -m(5) * t133 + t170 * mrSges(5,1) + t230 * t169 + t179 * t228 + t225 * t229 + t215;
t114 = m(4) * t136 + t188 * mrSges(4,1) - t187 * mrSges(4,2) + t208;
t103 = t200 * t110 + t203 * t114;
t100 = m(3) * t144 + qJDD(2) * mrSges(3,1) - t205 * mrSges(3,2) + t103;
t217 = t203 * t110 - t200 * t114;
t101 = m(3) * t145 - t205 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t217;
t95 = t204 * t100 + t201 * t101;
t112 = t202 * t119 + t199 * t120;
t220 = m(4) * t196 + t112;
t218 = -t201 * t100 + t204 * t101;
t213 = -mrSges(6,1) * t129 + mrSges(6,3) * t127 + Ifges(6,4) * t169 + Ifges(6,2) * t170 + Ifges(6,6) * qJDD(4) + qJD(4) * t154;
t150 = Ifges(6,3) * qJD(4) + (Ifges(6,5) * t199 + Ifges(6,6) * t202) * t189;
t211 = mrSges(6,2) * t129 - mrSges(6,3) * t126 + Ifges(6,1) * t169 + Ifges(6,4) * t170 + Ifges(6,5) * qJDD(4) + t150 * t228;
t151 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t199 + Ifges(5,6) * t202) * t189;
t105 = Ifges(5,4) * t169 + Ifges(5,2) * t170 + Ifges(5,6) * qJDD(4) + qJD(4) * t155 - mrSges(5,1) * t133 + mrSges(5,3) * t131 - pkin(4) * (t169 * mrSges(6,2) - t215) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t176 + t221) + (-pkin(4) * t176 - t150 - t151) * t229 + t213;
t107 = t151 * t228 + mrSges(5,2) * t133 - mrSges(5,3) * t130 + Ifges(5,1) * t169 + Ifges(5,4) * t170 + Ifges(5,5) * qJDD(4) - qJ(5) * t121 + (-t152 - t153) * qJD(4) + t211;
t210 = mrSges(4,1) * t136 - mrSges(4,2) * t137 + Ifges(4,3) * t188 + pkin(3) * t208 + pkin(7) * t216 + t202 * t105 + t199 * t107;
t209 = mrSges(3,1) * t144 - mrSges(3,2) * t145 + Ifges(3,3) * qJDD(2) + pkin(2) * t103 + t210;
t207 = mrSges(2,1) * t181 - mrSges(2,2) * t182 + pkin(1) * t95 + t209;
t96 = -mrSges(4,1) * t196 + mrSges(4,3) * t137 + t187 * Ifges(4,5) + Ifges(4,6) * t188 - pkin(3) * t112 - t233;
t93 = m(2) * t182 + t218;
t92 = m(2) * t181 + t95;
t91 = mrSges(4,2) * t196 - mrSges(4,3) * t136 + Ifges(4,5) * t188 - t187 * Ifges(4,6) - pkin(7) * t112 - t199 * t105 + t202 * t107;
t90 = mrSges(3,2) * t196 - mrSges(3,3) * t144 + Ifges(3,5) * qJDD(2) - t205 * Ifges(3,6) - pkin(6) * t103 - t200 * t96 + t203 * t91;
t89 = -mrSges(3,1) * t196 + mrSges(3,3) * t145 + t205 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t220 + pkin(6) * t217 + t200 * t91 + t203 * t96;
t88 = mrSges(2,2) * t196 - mrSges(2,3) * t181 - pkin(5) * t95 - t201 * t89 + t204 * t90;
t87 = -mrSges(2,1) * t196 + mrSges(2,3) * t182 + t201 * t90 + t204 * t89 - pkin(1) * (m(3) * t196 + t220) + pkin(5) * t218;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t198 * t88 - t197 * t87 - qJ(1) * (t197 * t93 + t198 * t92), t88, t90, t91, t107, -qJD(4) * t152 + t211; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t197 * t88 + t198 * t87 + qJ(1) * (-t197 * t92 + t198 * t93), t87, t89, t96, t105, -t150 * t229 + t213; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t207, t207, t209, t210, t233, -t154 * t228 - t212;];
m_new = t1;
