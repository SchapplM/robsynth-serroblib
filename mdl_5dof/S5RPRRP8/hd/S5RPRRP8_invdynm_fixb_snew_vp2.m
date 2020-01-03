% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:10
% EndTime: 2019-12-31 18:47:13
% DurationCPUTime: 1.41s
% Computational Cost: add. (17074->223), mult. (21740->267), div. (0->0), fcn. (7787->6), ass. (0->83)
t191 = qJD(1) ^ 2;
t185 = sin(qJ(1));
t188 = cos(qJ(1));
t163 = -t188 * g(1) - t185 * g(2);
t201 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t163;
t215 = -pkin(1) - pkin(2);
t127 = t215 * t191 + t201;
t162 = t185 * g(1) - t188 * g(2);
t200 = -t191 * qJ(2) + qJDD(2) - t162;
t131 = t215 * qJDD(1) + t200;
t184 = sin(qJ(3));
t187 = cos(qJ(3));
t121 = t187 * t127 + t184 * t131;
t170 = -qJD(1) + qJD(3);
t168 = t170 ^ 2;
t169 = -qJDD(1) + qJDD(3);
t118 = -(t168 * pkin(3)) + t169 * pkin(7) + t121;
t183 = sin(qJ(4));
t186 = cos(qJ(4));
t214 = t186 * g(3);
t114 = -t183 * t118 + t214;
t115 = t183 * g(3) + t186 * t118;
t132 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t183 - Ifges(6,3) * t186) * t170;
t135 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t183 + Ifges(5,2) * t186) * t170;
t151 = (-mrSges(6,1) * t186 - mrSges(6,3) * t183) * t170;
t207 = qJD(4) * t170;
t153 = t183 * t169 + t186 * t207;
t154 = t186 * t169 - t183 * t207;
t150 = (-pkin(4) * t186 - qJ(5) * t183) * t170;
t190 = qJD(4) ^ 2;
t210 = t170 * t186;
t111 = -t190 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t150 * t210 + t115;
t113 = -qJDD(4) * pkin(4) - t214 - t190 * qJ(5) + qJDD(5) + (t150 * t170 + t118) * t183;
t202 = -mrSges(6,1) * t113 + mrSges(6,3) * t111 + Ifges(6,4) * t153 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t154;
t160 = mrSges(6,2) * t210 + qJD(4) * mrSges(6,3);
t205 = -m(6) * t113 + qJDD(4) * mrSges(6,1) + qJD(4) * t160;
t211 = t170 * t183;
t158 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t211;
t206 = m(6) * t111 + qJDD(4) * mrSges(6,3) + qJD(4) * t158 + t151 * t210;
t136 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t183 - Ifges(6,5) * t186) * t170;
t208 = t136 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t183 + Ifges(5,4) * t186) * t170;
t217 = -(t208 * t186 + (t132 - t135) * t183) * t170 + mrSges(5,1) * t114 - mrSges(5,2) * t115 + Ifges(5,5) * t153 + Ifges(5,6) * t154 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t153 * mrSges(6,2) - t151 * t211 + t205) + qJ(5) * (t154 * mrSges(6,2) + t206) + t202;
t213 = mrSges(2,1) + mrSges(3,1);
t212 = mrSges(5,3) + mrSges(6,2);
t120 = -t184 * t127 + t187 * t131;
t117 = -t169 * pkin(3) - t168 * pkin(7) - t120;
t108 = -t154 * pkin(4) - t153 * qJ(5) + (-0.2e1 * qJD(5) * t183 + (pkin(4) * t183 - qJ(5) * t186) * qJD(4)) * t170 + t117;
t105 = m(6) * t108 - t154 * mrSges(6,1) - t153 * mrSges(6,3) - t158 * t211 - t160 * t210;
t157 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t211;
t159 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t210;
t101 = -m(5) * t117 + t154 * mrSges(5,1) - t153 * mrSges(5,2) - t157 * t211 + t159 * t210 - t105;
t100 = m(4) * t120 + t169 * mrSges(4,1) - (t168 * mrSges(4,2)) + t101;
t152 = (-mrSges(5,1) * t186 + mrSges(5,2) * t183) * t170;
t103 = m(5) * t115 - qJDD(4) * mrSges(5,2) - qJD(4) * t157 + t152 * t210 + t212 * t154 + t206;
t104 = m(5) * t114 + qJDD(4) * mrSges(5,1) + qJD(4) * t159 + (-t151 - t152) * t211 - t212 * t153 + t205;
t99 = t186 * t103 - t183 * t104;
t95 = m(4) * t121 - t168 * mrSges(4,1) - t169 * mrSges(4,2) + t99;
t91 = -t184 * t100 + t187 * t95;
t204 = -mrSges(6,1) * t108 + mrSges(6,2) * t111;
t90 = t187 * t100 + t184 * t95;
t98 = t183 * t103 + t186 * t104;
t138 = -t191 * pkin(1) + t201;
t203 = m(3) * t138 + qJDD(1) * mrSges(3,3) + t91;
t134 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t183 - Ifges(6,6) * t186) * t170;
t199 = mrSges(6,2) * t113 - mrSges(6,3) * t108 + Ifges(6,1) * t153 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t154 + qJD(4) * t132 + t134 * t210;
t149 = -qJDD(1) * pkin(1) + t200;
t198 = -m(3) * t149 + qJDD(1) * mrSges(3,1) + t191 * mrSges(3,3) - t90;
t133 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t183 + Ifges(5,6) * t186) * t170;
t92 = -mrSges(5,1) * t117 + mrSges(5,3) * t115 - pkin(4) * t105 + (-t133 - t134) * t211 + (Ifges(5,2) + Ifges(6,3)) * t154 + (Ifges(5,4) - Ifges(6,5)) * t153 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t208 * qJD(4) + t204;
t93 = mrSges(5,2) * t117 - mrSges(5,3) * t114 + Ifges(5,1) * t153 + Ifges(5,4) * t154 + Ifges(5,5) * qJDD(4) - qJ(5) * t105 - qJD(4) * t135 + t133 * t210 + t199;
t84 = mrSges(4,2) * g(3) - mrSges(4,3) * t120 + Ifges(4,5) * t169 - (t168 * Ifges(4,6)) - pkin(7) * t98 - t183 * t92 + t186 * t93;
t85 = -mrSges(4,1) * g(3) + mrSges(4,3) * t121 + t168 * Ifges(4,5) + Ifges(4,6) * t169 - pkin(3) * t98 - t217;
t197 = mrSges(3,2) * t149 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t191 * Ifges(3,6) - pkin(6) * t90 - t184 * t85 + t187 * t84;
t196 = mrSges(3,2) * t138 - pkin(2) * (-m(4) * g(3) - t98) - pkin(6) * t91 - t184 * t84 - t187 * t85;
t195 = mrSges(4,1) * t120 - mrSges(4,2) * t121 + Ifges(4,3) * t169 + pkin(3) * t101 + pkin(7) * t99 + t183 * t93 + t186 * t92;
t194 = -mrSges(3,1) * t149 + mrSges(3,3) * t138 + Ifges(3,2) * qJDD(1) - pkin(2) * t90 - t195;
t192 = -mrSges(2,2) * t163 + mrSges(2,1) * t162 + Ifges(2,3) * qJDD(1) + t194 + qJ(2) * (-t191 * mrSges(3,1) + t203) + pkin(1) * t198;
t96 = (-m(3) - m(4)) * g(3) - t98;
t87 = m(2) * t162 + qJDD(1) * mrSges(2,1) - t191 * mrSges(2,2) + t198;
t86 = m(2) * t163 - qJDD(1) * mrSges(2,2) - t213 * t191 + t203;
t82 = -mrSges(2,2) * g(3) - mrSges(2,3) * t162 + Ifges(2,5) * qJDD(1) - t191 * Ifges(2,6) - qJ(2) * t96 + t197;
t81 = mrSges(2,3) * t163 - pkin(1) * t96 + (Ifges(3,4) + Ifges(2,5)) * t191 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t213 * g(3) + t196;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t188 * t82 - t185 * t81 - pkin(5) * (t185 * t86 + t188 * t87), t82, t197, t84, t93, t199; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t185 * t82 + t188 * t81 + pkin(5) * (-t185 * t87 + t188 * t86), t81, t194, t85, t92, (-t183 * t132 - t186 * t136) * t170 + t202; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t192, t192, -mrSges(3,1) * g(3) - t191 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t196, t195, t217, Ifges(6,5) * t153 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t154 - qJD(4) * t136 + t134 * t211 - t204;];
m_new = t1;
