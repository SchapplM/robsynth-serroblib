% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:42
% EndTime: 2019-12-05 16:41:46
% DurationCPUTime: 1.98s
% Computational Cost: add. (26592->210), mult. (35379->264), div. (0->0), fcn. (19685->8), ass. (0->84)
t190 = sin(pkin(8));
t191 = cos(pkin(8));
t175 = t190 * g(1) - t191 * g(2);
t176 = -t191 * g(1) - t190 * g(2);
t194 = sin(qJ(2));
t197 = cos(qJ(2));
t140 = t197 * t175 - t194 * t176;
t137 = qJDD(2) * pkin(2) + t140;
t141 = t194 * t175 + t197 * t176;
t199 = qJD(2) ^ 2;
t138 = -t199 * pkin(2) + t141;
t193 = sin(qJ(3));
t196 = cos(qJ(3));
t133 = t193 * t137 + t196 * t138;
t183 = qJD(2) + qJD(3);
t181 = t183 ^ 2;
t182 = qJDD(2) + qJDD(3);
t130 = -t181 * pkin(3) + t182 * pkin(7) + t133;
t192 = sin(qJ(4));
t189 = -g(3) + qJDD(1);
t195 = cos(qJ(4));
t217 = t195 * t189;
t126 = -t192 * t130 + t217;
t127 = t195 * t130 + t192 * t189;
t145 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t192 - Ifges(6,3) * t195) * t183;
t148 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t192 + Ifges(5,2) * t195) * t183;
t161 = (-mrSges(6,1) * t195 - mrSges(6,3) * t192) * t183;
t214 = qJD(4) * t183;
t163 = t192 * t182 + t195 * t214;
t164 = t195 * t182 - t192 * t214;
t160 = (-pkin(4) * t195 - qJ(5) * t192) * t183;
t198 = qJD(4) ^ 2;
t218 = t183 * t195;
t123 = -t198 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t160 * t218 + t127;
t125 = -qJDD(4) * pkin(4) - t198 * qJ(5) - t217 + qJDD(5) + (t160 * t183 + t130) * t192;
t206 = -mrSges(6,1) * t125 + mrSges(6,3) * t123 + Ifges(6,4) * t163 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t164;
t173 = mrSges(6,2) * t218 + qJD(4) * mrSges(6,3);
t208 = -m(6) * t125 + qJDD(4) * mrSges(6,1) + qJD(4) * t173;
t219 = t183 * t192;
t171 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t219;
t209 = m(6) * t123 + qJDD(4) * mrSges(6,3) + qJD(4) * t171 + t161 * t218;
t149 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t192 - Ifges(6,5) * t195) * t183;
t215 = t149 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t192 + Ifges(5,4) * t195) * t183;
t222 = -(t215 * t195 + (t145 - t148) * t192) * t183 + mrSges(5,1) * t126 - mrSges(5,2) * t127 + Ifges(5,5) * t163 + Ifges(5,6) * t164 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t163 * mrSges(6,2) - t161 * t219 + t208) + qJ(5) * (t164 * mrSges(6,2) + t209) + t206;
t220 = mrSges(5,3) + mrSges(6,2);
t162 = (-mrSges(5,1) * t195 + mrSges(5,2) * t192) * t183;
t170 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t219;
t115 = m(5) * t127 - qJDD(4) * mrSges(5,2) - qJD(4) * t170 + t162 * t218 + t220 * t164 + t209;
t172 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t218;
t116 = m(5) * t126 + qJDD(4) * mrSges(5,1) + qJD(4) * t172 + (-t161 - t162) * t219 - t220 * t163 + t208;
t210 = t195 * t115 - t192 * t116;
t106 = m(4) * t133 - t181 * mrSges(4,1) - t182 * mrSges(4,2) + t210;
t132 = t196 * t137 - t193 * t138;
t129 = -t182 * pkin(3) - t181 * pkin(7) - t132;
t120 = -t164 * pkin(4) - t163 * qJ(5) + (-0.2e1 * qJD(5) * t192 + (pkin(4) * t192 - qJ(5) * t195) * qJD(4)) * t183 + t129;
t117 = m(6) * t120 - t164 * mrSges(6,1) - t163 * mrSges(6,3) - t171 * t219 - t173 * t218;
t202 = -m(5) * t129 + t164 * mrSges(5,1) - t163 * mrSges(5,2) - t170 * t219 + t172 * t218 - t117;
t110 = m(4) * t132 + t182 * mrSges(4,1) - t181 * mrSges(4,2) + t202;
t99 = t193 * t106 + t196 * t110;
t96 = m(3) * t140 + qJDD(2) * mrSges(3,1) - t199 * mrSges(3,2) + t99;
t211 = t196 * t106 - t193 * t110;
t97 = m(3) * t141 - t199 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t211;
t90 = t194 * t97 + t197 * t96;
t108 = t192 * t115 + t195 * t116;
t213 = m(4) * t189 + t108;
t212 = -t194 * t96 + t197 * t97;
t207 = -mrSges(6,1) * t120 + mrSges(6,2) * t123;
t146 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t192 + Ifges(5,6) * t195) * t183;
t147 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t192 - Ifges(6,6) * t195) * t183;
t102 = -mrSges(5,1) * t129 + mrSges(5,3) * t127 - pkin(4) * t117 + (-t146 - t147) * t219 + (Ifges(5,2) + Ifges(6,3)) * t164 + (Ifges(5,4) - Ifges(6,5)) * t163 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t215 * qJD(4) + t207;
t204 = mrSges(6,2) * t125 - mrSges(6,3) * t120 + Ifges(6,1) * t163 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t164 + qJD(4) * t145 + t147 * t218;
t103 = mrSges(5,2) * t129 - mrSges(5,3) * t126 + Ifges(5,1) * t163 + Ifges(5,4) * t164 + Ifges(5,5) * qJDD(4) - qJ(5) * t117 - qJD(4) * t148 + t146 * t218 + t204;
t205 = mrSges(4,1) * t132 - mrSges(4,2) * t133 + Ifges(4,3) * t182 + pkin(3) * t202 + pkin(7) * t210 + t195 * t102 + t192 * t103;
t203 = mrSges(3,1) * t140 - mrSges(3,2) * t141 + Ifges(3,3) * qJDD(2) + pkin(2) * t99 + t205;
t201 = mrSges(2,1) * t175 - mrSges(2,2) * t176 + pkin(1) * t90 + t203;
t92 = -mrSges(4,1) * t189 + mrSges(4,3) * t133 + t181 * Ifges(4,5) + Ifges(4,6) * t182 - pkin(3) * t108 - t222;
t91 = mrSges(4,2) * t189 - mrSges(4,3) * t132 + Ifges(4,5) * t182 - t181 * Ifges(4,6) - pkin(7) * t108 - t192 * t102 + t195 * t103;
t88 = m(2) * t176 + t212;
t87 = m(2) * t175 + t90;
t86 = mrSges(3,2) * t189 - mrSges(3,3) * t140 + Ifges(3,5) * qJDD(2) - t199 * Ifges(3,6) - pkin(6) * t99 - t193 * t92 + t196 * t91;
t85 = -mrSges(3,1) * t189 + mrSges(3,3) * t141 + t199 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t213 + pkin(6) * t211 + t193 * t91 + t196 * t92;
t84 = mrSges(2,2) * t189 - mrSges(2,3) * t175 - pkin(5) * t90 - t194 * t85 + t197 * t86;
t83 = -mrSges(2,1) * t189 + mrSges(2,3) * t176 + t194 * t86 + t197 * t85 - pkin(1) * (m(3) * t189 + t213) + pkin(5) * t212;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t191 * t84 - t190 * t83 - qJ(1) * (t190 * t88 + t191 * t87), t84, t86, t91, t103, t204; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t190 * t84 + t191 * t83 + qJ(1) * (-t190 * t87 + t191 * t88), t83, t85, t92, t102, (-t192 * t145 - t195 * t149) * t183 + t206; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t201, t201, t203, t205, t222, Ifges(6,5) * t163 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t164 - qJD(4) * t149 + t147 * t219 - t207;];
m_new = t1;
