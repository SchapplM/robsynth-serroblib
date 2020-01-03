% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:30
% EndTime: 2019-12-31 19:49:32
% DurationCPUTime: 2.10s
% Computational Cost: add. (31870->221), mult. (40655->275), div. (0->0), fcn. (19685->8), ass. (0->86)
t198 = sin(qJ(1));
t201 = cos(qJ(1));
t177 = t198 * g(1) - t201 * g(2);
t170 = qJDD(1) * pkin(1) + t177;
t178 = -t201 * g(1) - t198 * g(2);
t203 = qJD(1) ^ 2;
t171 = -t203 * pkin(1) + t178;
t197 = sin(qJ(2));
t200 = cos(qJ(2));
t140 = t200 * t170 - t197 * t171;
t186 = qJDD(1) + qJDD(2);
t137 = t186 * pkin(2) + t140;
t141 = t197 * t170 + t200 * t171;
t187 = qJD(1) + qJD(2);
t185 = t187 ^ 2;
t138 = -t185 * pkin(2) + t141;
t194 = sin(pkin(8));
t195 = cos(pkin(8));
t133 = t194 * t137 + t195 * t138;
t130 = -t185 * pkin(3) + t186 * pkin(7) + t133;
t196 = sin(qJ(4));
t193 = -g(3) + qJDD(3);
t199 = cos(qJ(4));
t221 = t199 * t193;
t126 = -t196 * t130 + t221;
t127 = t199 * t130 + t196 * t193;
t145 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t196 - Ifges(6,3) * t199) * t187;
t148 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t196 + Ifges(5,2) * t199) * t187;
t161 = (-mrSges(6,1) * t199 - mrSges(6,3) * t196) * t187;
t218 = qJD(4) * t187;
t163 = t196 * t186 + t199 * t218;
t164 = t199 * t186 - t196 * t218;
t160 = (-pkin(4) * t199 - qJ(5) * t196) * t187;
t202 = qJD(4) ^ 2;
t222 = t187 * t199;
t123 = -t202 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t160 * t222 + t127;
t125 = -qJDD(4) * pkin(4) - t202 * qJ(5) - t221 + qJDD(5) + (t160 * t187 + t130) * t196;
t210 = -mrSges(6,1) * t125 + mrSges(6,3) * t123 + Ifges(6,4) * t163 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t164;
t175 = mrSges(6,2) * t222 + qJD(4) * mrSges(6,3);
t212 = -m(6) * t125 + qJDD(4) * mrSges(6,1) + qJD(4) * t175;
t223 = t187 * t196;
t173 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t223;
t213 = m(6) * t123 + qJDD(4) * mrSges(6,3) + qJD(4) * t173 + t161 * t222;
t149 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t196 - Ifges(6,5) * t199) * t187;
t219 = t149 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t196 + Ifges(5,4) * t199) * t187;
t226 = -((t145 - t148) * t196 + t219 * t199) * t187 + mrSges(5,1) * t126 - mrSges(5,2) * t127 + Ifges(5,5) * t163 + Ifges(5,6) * t164 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t163 * mrSges(6,2) - t161 * t223 + t212) + qJ(5) * (t164 * mrSges(6,2) + t213) + t210;
t224 = mrSges(5,3) + mrSges(6,2);
t162 = (-mrSges(5,1) * t199 + mrSges(5,2) * t196) * t187;
t172 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t223;
t115 = m(5) * t127 - qJDD(4) * mrSges(5,2) - qJD(4) * t172 + t162 * t222 + t224 * t164 + t213;
t174 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t222;
t116 = m(5) * t126 + qJDD(4) * mrSges(5,1) + qJD(4) * t174 + (-t161 - t162) * t223 - t224 * t163 + t212;
t214 = t199 * t115 - t196 * t116;
t106 = m(4) * t133 - t185 * mrSges(4,1) - t186 * mrSges(4,2) + t214;
t132 = t195 * t137 - t194 * t138;
t129 = -t186 * pkin(3) - t185 * pkin(7) - t132;
t120 = -t164 * pkin(4) - t163 * qJ(5) + (-0.2e1 * qJD(5) * t196 + (pkin(4) * t196 - qJ(5) * t199) * qJD(4)) * t187 + t129;
t117 = m(6) * t120 - t164 * mrSges(6,1) - t163 * mrSges(6,3) - t173 * t223 - t175 * t222;
t206 = -m(5) * t129 + t164 * mrSges(5,1) - t163 * mrSges(5,2) - t172 * t223 + t174 * t222 - t117;
t110 = m(4) * t132 + t186 * mrSges(4,1) - t185 * mrSges(4,2) + t206;
t99 = t194 * t106 + t195 * t110;
t96 = m(3) * t140 + t186 * mrSges(3,1) - t185 * mrSges(3,2) + t99;
t215 = t195 * t106 - t194 * t110;
t97 = m(3) * t141 - t185 * mrSges(3,1) - t186 * mrSges(3,2) + t215;
t90 = t197 * t97 + t200 * t96;
t108 = t196 * t115 + t199 * t116;
t217 = m(4) * t193 + t108;
t216 = -t197 * t96 + t200 * t97;
t211 = -mrSges(6,1) * t120 + mrSges(6,2) * t123;
t146 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t196 + Ifges(5,6) * t199) * t187;
t147 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t196 - Ifges(6,6) * t199) * t187;
t102 = -mrSges(5,1) * t129 + mrSges(5,3) * t127 - pkin(4) * t117 + (-t146 - t147) * t223 + (Ifges(5,2) + Ifges(6,3)) * t164 + (Ifges(5,4) - Ifges(6,5)) * t163 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t219 * qJD(4) + t211;
t208 = mrSges(6,2) * t125 - mrSges(6,3) * t120 + Ifges(6,1) * t163 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t164 + qJD(4) * t145 + t147 * t222;
t103 = mrSges(5,2) * t129 - mrSges(5,3) * t126 + Ifges(5,1) * t163 + Ifges(5,4) * t164 + Ifges(5,5) * qJDD(4) - qJ(5) * t117 - qJD(4) * t148 + t146 * t222 + t208;
t209 = mrSges(4,1) * t132 - mrSges(4,2) * t133 + Ifges(4,3) * t186 + pkin(3) * t206 + pkin(7) * t214 + t199 * t102 + t196 * t103;
t207 = mrSges(3,1) * t140 - mrSges(3,2) * t141 + Ifges(3,3) * t186 + pkin(2) * t99 + t209;
t205 = mrSges(2,1) * t177 - mrSges(2,2) * t178 + Ifges(2,3) * qJDD(1) + pkin(1) * t90 + t207;
t92 = -mrSges(4,1) * t193 + mrSges(4,3) * t133 + t185 * Ifges(4,5) + Ifges(4,6) * t186 - pkin(3) * t108 - t226;
t91 = mrSges(4,2) * t193 - mrSges(4,3) * t132 + Ifges(4,5) * t186 - t185 * Ifges(4,6) - pkin(7) * t108 - t196 * t102 + t199 * t103;
t88 = m(2) * t178 - t203 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t216;
t87 = m(2) * t177 + qJDD(1) * mrSges(2,1) - t203 * mrSges(2,2) + t90;
t86 = -mrSges(3,2) * g(3) - mrSges(3,3) * t140 + Ifges(3,5) * t186 - t185 * Ifges(3,6) - qJ(3) * t99 - t194 * t92 + t195 * t91;
t85 = mrSges(3,1) * g(3) + mrSges(3,3) * t141 + t185 * Ifges(3,5) + Ifges(3,6) * t186 - pkin(2) * t217 + qJ(3) * t215 + t194 * t91 + t195 * t92;
t84 = -mrSges(2,2) * g(3) - mrSges(2,3) * t177 + Ifges(2,5) * qJDD(1) - t203 * Ifges(2,6) - pkin(6) * t90 - t197 * t85 + t200 * t86;
t83 = Ifges(2,6) * qJDD(1) + t203 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t178 + t197 * t86 + t200 * t85 - pkin(1) * (-m(3) * g(3) + t217) + pkin(6) * t216;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t201 * t84 - t198 * t83 - pkin(5) * (t198 * t88 + t201 * t87), t84, t86, t91, t103, t208; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t198 * t84 + t201 * t83 + pkin(5) * (-t198 * t87 + t201 * t88), t83, t85, t92, t102, (-t196 * t145 - t199 * t149) * t187 + t210; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t205, t205, t207, t209, t226, Ifges(6,5) * t163 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t164 - qJD(4) * t149 + t147 * t223 - t211;];
m_new = t1;
