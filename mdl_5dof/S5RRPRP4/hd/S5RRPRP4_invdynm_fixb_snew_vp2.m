% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:39
% EndTime: 2019-12-31 19:52:41
% DurationCPUTime: 1.20s
% Computational Cost: add. (15552->221), mult. (18418->265), div. (0->0), fcn. (7745->6), ass. (0->85)
t190 = sin(qJ(1));
t193 = cos(qJ(1));
t169 = t190 * g(1) - t193 * g(2);
t162 = qJDD(1) * pkin(1) + t169;
t170 = -t193 * g(1) - t190 * g(2);
t195 = qJD(1) ^ 2;
t163 = -t195 * pkin(1) + t170;
t189 = sin(qJ(2));
t192 = cos(qJ(2));
t130 = t189 * t162 + t192 * t163;
t180 = qJDD(1) + qJDD(2);
t181 = (qJD(1) + qJD(2));
t224 = -t180 * qJ(3) - (2 * qJD(3) * t181) - t130;
t223 = (-pkin(2) - pkin(7));
t188 = sin(qJ(4));
t222 = t188 * g(3);
t221 = mrSges(3,1) - mrSges(4,2);
t220 = -mrSges(5,3) - mrSges(6,2);
t219 = Ifges(3,5) - Ifges(4,4);
t218 = (-Ifges(3,6) + Ifges(4,5));
t179 = t181 ^ 2;
t122 = (t223 * t179) - t224;
t191 = cos(qJ(4));
t214 = qJD(4) * t181;
t155 = t188 * t180 + t191 * t214;
t156 = t191 * t180 - t188 * t214;
t112 = t155 * pkin(4) - t156 * qJ(5) + (-0.2e1 * qJD(5) * t191 + (pkin(4) * t191 + qJ(5) * t188) * qJD(4)) * t181 + t122;
t216 = t181 * t191;
t166 = -(qJD(4) * mrSges(6,1)) + mrSges(6,2) * t216;
t217 = t181 * t188;
t167 = -mrSges(6,2) * t217 + (qJD(4) * mrSges(6,3));
t106 = m(6) * t112 + t155 * mrSges(6,1) - t156 * mrSges(6,3) - t166 * t216 + t167 * t217;
t164 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t217;
t165 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t216;
t102 = -m(5) * t122 - t155 * mrSges(5,1) - t156 * mrSges(5,2) - t164 * t217 - t165 * t216 - t106;
t125 = (t179 * pkin(2)) + t224;
t199 = -m(4) * t125 + (t179 * mrSges(4,2)) + t180 * mrSges(4,3) - t102;
t100 = m(3) * t130 - (t179 * mrSges(3,1)) - t180 * mrSges(3,2) + t199;
t129 = t192 * t162 - t189 * t163;
t206 = -t179 * qJ(3) + qJDD(3) - t129;
t127 = -t180 * pkin(2) + t206;
t124 = t223 * t180 + t206;
t120 = -t191 * g(3) + t188 * t124;
t153 = (mrSges(6,1) * t188 - mrSges(6,3) * t191) * t181;
t209 = t181 * (-t153 - (mrSges(5,1) * t188 + mrSges(5,2) * t191) * t181);
t152 = (pkin(4) * t188 - qJ(5) * t191) * t181;
t194 = qJD(4) ^ 2;
t115 = -t194 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t152 * t217 + t120;
t213 = m(6) * t115 + qJDD(4) * mrSges(6,3) + qJD(4) * t166;
t104 = m(5) * t120 - qJDD(4) * mrSges(5,2) - qJD(4) * t165 + t220 * t155 + t188 * t209 + t213;
t119 = t191 * t124 + t222;
t117 = -qJDD(4) * pkin(4) - t222 - t194 * qJ(5) + qJDD(5) + (t152 * t181 - t124) * t191;
t208 = -m(6) * t117 + qJDD(4) * mrSges(6,1) + qJD(4) * t167;
t105 = m(5) * t119 + qJDD(4) * mrSges(5,1) + qJD(4) * t164 + t220 * t156 + t191 * t209 + t208;
t96 = t188 * t104 + t191 * t105;
t204 = -m(4) * t127 + t179 * mrSges(4,3) - t96;
t93 = m(3) * t129 - t179 * mrSges(3,2) + t221 * t180 + t204;
t88 = t189 * t100 + t192 * t93;
t211 = t192 * t100 - t189 * t93;
t137 = (Ifges(6,2) * qJD(4)) + (Ifges(6,4) * t191 + Ifges(6,6) * t188) * t181;
t210 = t181 * (-(Ifges(5,3) * qJD(4)) - (Ifges(5,5) * t191 - Ifges(5,6) * t188) * t181 - t137);
t97 = t191 * t104 - t188 * t105;
t207 = -mrSges(6,1) * t112 + mrSges(6,2) * t115;
t135 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t191 + Ifges(6,3) * t188) * t181;
t205 = mrSges(6,2) * t117 - mrSges(6,3) * t112 + Ifges(6,1) * t156 + Ifges(6,4) * qJDD(4) + Ifges(6,5) * t155 + qJD(4) * t135;
t139 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t191 + Ifges(6,5) * t188) * t181;
t140 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t191 - Ifges(5,4) * t188) * t181;
t90 = -mrSges(5,1) * t122 + mrSges(5,3) * t120 - pkin(4) * t106 + t191 * t210 + (Ifges(5,4) - Ifges(6,5)) * t156 + (-Ifges(5,2) - Ifges(6,3)) * t155 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + (t139 + t140) * qJD(4) + t207;
t138 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t191 - Ifges(5,2) * t188) * t181;
t91 = mrSges(5,2) * t122 - mrSges(5,3) * t119 + Ifges(5,1) * t156 - Ifges(5,4) * t155 + Ifges(5,5) * qJDD(4) - qJ(5) * t106 - qJD(4) * t138 + t188 * t210 + t205;
t203 = mrSges(4,2) * t127 - mrSges(4,3) * t125 + Ifges(4,1) * t180 - pkin(7) * t96 - t188 * t90 + t191 * t91;
t202 = -mrSges(6,1) * t117 + mrSges(6,3) * t115 + Ifges(6,4) * t156 + Ifges(6,2) * qJDD(4) + Ifges(6,6) * t155 - t135 * t216 + t139 * t217;
t201 = -mrSges(4,1) * t125 - pkin(3) * t102 - pkin(7) * t97 - t188 * t91 - t191 * t90;
t200 = -mrSges(3,2) * t130 + qJ(3) * t199 + mrSges(3,1) * t129 + Ifges(3,3) * t180 + t203 + pkin(2) * (-t180 * mrSges(4,2) + t204);
t198 = mrSges(2,1) * t169 - mrSges(2,2) * t170 + Ifges(2,3) * qJDD(1) + pkin(1) * t88 + t200;
t197 = -mrSges(5,2) * t120 + qJ(5) * (-t155 * mrSges(6,2) - t153 * t217 + t213) + pkin(4) * (-t156 * mrSges(6,2) - t153 * t216 + t208) + mrSges(5,1) * t119 + t140 * t217 + t138 * t216 - Ifges(5,6) * t155 + Ifges(5,5) * t156 + Ifges(5,3) * qJDD(4) + t202;
t196 = mrSges(4,1) * t127 + pkin(3) * t96 + t197;
t95 = -m(4) * g(3) + t97;
t86 = m(2) * t170 - t195 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t211;
t85 = m(2) * t169 + qJDD(1) * mrSges(2,1) - t195 * mrSges(2,2) + t88;
t84 = (t218 * t179) + t219 * t180 + (-mrSges(3,2) + mrSges(4,3)) * g(3) + t196 - mrSges(3,3) * t129 - qJ(3) * t95;
t83 = mrSges(3,3) * t130 - pkin(2) * t95 + t221 * g(3) + t219 * t179 - t218 * t180 + t201;
t82 = -mrSges(2,2) * g(3) - mrSges(2,3) * t169 + Ifges(2,5) * qJDD(1) - t195 * Ifges(2,6) - pkin(6) * t88 - t189 * t83 + t192 * t84;
t81 = Ifges(2,6) * qJDD(1) + t195 * Ifges(2,5) + mrSges(2,3) * t170 + t189 * t84 + t192 * t83 - pkin(1) * t97 + pkin(6) * t211 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t193 * t82 - t190 * t81 - pkin(5) * (t190 * t86 + t193 * t85), t82, t84, t203, t91, -t137 * t217 + t205; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t190 * t82 + t193 * t81 + pkin(5) * (-t190 * t85 + t193 * t86), t81, t83, -mrSges(4,3) * g(3) + Ifges(4,4) * t180 - (t179 * Ifges(4,5)) - t196, t90, t202; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t198, t198, t200, mrSges(4,2) * g(3) + t179 * Ifges(4,4) + Ifges(4,5) * t180 - t201, t197, Ifges(6,5) * t156 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t155 - qJD(4) * t139 + t137 * t216 - t207;];
m_new = t1;
