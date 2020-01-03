% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:07
% EndTime: 2020-01-03 11:36:13
% DurationCPUTime: 3.75s
% Computational Cost: add. (52562->203), mult. (73380->265), div. (0->0), fcn. (40823->10), ass. (0->101)
t228 = 2 * qJD(4);
t192 = sin(qJ(1));
t195 = cos(qJ(1));
t172 = -t195 * g(2) - t192 * g(3);
t165 = qJDD(1) * pkin(1) + t172;
t171 = -t192 * g(2) + t195 * g(3);
t196 = qJD(1) ^ 2;
t166 = -t196 * pkin(1) + t171;
t187 = sin(pkin(8));
t189 = cos(pkin(8));
t149 = t189 * t165 - t187 * t166;
t143 = qJDD(1) * pkin(2) + t149;
t150 = t187 * t165 + t189 * t166;
t144 = -t196 * pkin(2) + t150;
t191 = sin(qJ(3));
t194 = cos(qJ(3));
t139 = t191 * t143 + t194 * t144;
t182 = (qJD(1) + qJD(3));
t180 = t182 ^ 2;
t181 = qJDD(1) + qJDD(3);
t136 = -t180 * pkin(3) + t181 * qJ(4) + t139;
t227 = (t182 * t228) + t136;
t186 = sin(pkin(9));
t222 = t182 * t186;
t188 = cos(pkin(9));
t220 = t188 * t182;
t185 = -g(1) + qJDD(2);
t132 = t186 * t185 + t227 * t188;
t210 = -pkin(4) * t188 - pkin(7) * t186;
t161 = t210 * t182;
t130 = t161 * t220 + t132;
t138 = t194 * t143 - t191 * t144;
t203 = -t180 * qJ(4) + qJDD(4) - t138;
t133 = (-pkin(3) + t210) * t181 + t203;
t190 = sin(qJ(5));
t193 = cos(qJ(5));
t126 = -t190 * t130 + t193 * t133;
t169 = qJD(5) - t220;
t216 = t190 * t222;
t152 = -t169 * mrSges(6,2) - mrSges(6,3) * t216;
t154 = (mrSges(6,1) * t190 + mrSges(6,2) * t193) * t222;
t217 = qJD(5) * t182;
t156 = (t181 * t193 - t190 * t217) * t186;
t221 = t188 * t181;
t168 = qJDD(5) - t221;
t215 = t193 * t222;
t124 = m(6) * t126 + t168 * mrSges(6,1) - t156 * mrSges(6,3) + t169 * t152 - t154 * t215;
t127 = t193 * t130 + t190 * t133;
t153 = t169 * mrSges(6,1) - mrSges(6,3) * t215;
t155 = (-t181 * t190 - t193 * t217) * t186;
t125 = m(6) * t127 - t168 * mrSges(6,2) + t155 * mrSges(6,3) - t169 * t153 - t154 * t216;
t118 = -t190 * t124 + t193 * t125;
t219 = t188 * t185;
t129 = -t219 + (t136 + (t228 + t161) * t182) * t186;
t146 = Ifges(6,3) * t169 + (Ifges(6,5) * t193 - Ifges(6,6) * t190) * t222;
t148 = Ifges(6,5) * t169 + (Ifges(6,1) * t193 - Ifges(6,4) * t190) * t222;
t119 = -mrSges(6,1) * t129 + mrSges(6,3) * t127 + Ifges(6,4) * t156 + Ifges(6,2) * t155 + Ifges(6,6) * t168 - t146 * t215 + t169 * t148;
t147 = Ifges(6,6) * t169 + (Ifges(6,4) * t193 - Ifges(6,2) * t190) * t222;
t120 = mrSges(6,2) * t129 - mrSges(6,3) * t126 + Ifges(6,1) * t156 + Ifges(6,4) * t155 + Ifges(6,5) * t168 - t146 * t216 - t169 * t147;
t131 = -t227 * t186 + t219;
t204 = -m(6) * t129 + t155 * mrSges(6,1) - t156 * mrSges(6,2);
t206 = -t152 * t190 - t153 * t193;
t209 = Ifges(5,1) * t186 + Ifges(5,4) * t188;
t226 = -((Ifges(5,4) * t186 + Ifges(5,2) * t188) * t222 - t209 * t220) * t182 - mrSges(5,1) * t131 + mrSges(5,2) * t132 - pkin(4) * (t206 * t222 + t204) - pkin(7) * t118 - t193 * t119 - t190 * t120;
t225 = mrSges(5,2) * t186;
t157 = (-mrSges(5,1) * t188 + t225) * t182;
t223 = mrSges(5,3) * t181;
t115 = m(5) * t132 + (t157 * t182 + t223) * t188 + t118;
t122 = m(5) * t131 + (-t223 + (-t157 + t206) * t182) * t186 + t204;
t211 = t188 * t115 - t186 * t122;
t108 = m(4) * t139 - t180 * mrSges(4,1) - t181 * mrSges(4,2) + t211;
t117 = t193 * t124 + t190 * t125;
t135 = -t181 * pkin(3) + t203;
t201 = -m(5) * t135 + mrSges(5,1) * t221 - t117 + (t186 ^ 2 + t188 ^ 2) * mrSges(5,3) * t180;
t112 = m(4) * t138 - t180 * mrSges(4,2) + (mrSges(4,1) - t225) * t181 + t201;
t101 = t191 * t108 + t194 * t112;
t98 = m(3) * t149 + qJDD(1) * mrSges(3,1) - t196 * mrSges(3,2) + t101;
t212 = t194 * t108 - t191 * t112;
t99 = m(3) * t150 - t196 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t212;
t93 = t187 * t99 + t189 * t98;
t110 = t186 * t115 + t188 * t122;
t214 = m(4) * t185 + t110;
t213 = -t187 * t98 + t189 * t99;
t208 = Ifges(5,5) * t186 + Ifges(5,6) * t188;
t207 = t147 * t193 + t148 * t190;
t158 = t208 * t182;
t103 = mrSges(5,2) * t135 - mrSges(5,3) * t131 - pkin(7) * t117 - t190 * t119 + t193 * t120 + t158 * t220 + t209 * t181;
t200 = mrSges(6,1) * t126 - mrSges(6,2) * t127 + Ifges(6,5) * t156 + Ifges(6,6) * t155 + Ifges(6,3) * t168;
t105 = Ifges(5,2) * t221 - mrSges(5,1) * t135 + mrSges(5,3) * t132 - pkin(4) * t117 + (Ifges(5,4) * t181 + (-t158 - t207) * t182) * t186 - t200;
t202 = -mrSges(4,2) * t139 + qJ(4) * t211 + t186 * t103 + t188 * t105 + pkin(3) * (-t181 * t225 + t201) + mrSges(4,1) * t138 + Ifges(4,3) * t181;
t199 = mrSges(3,1) * t149 - mrSges(3,2) * t150 + Ifges(3,3) * qJDD(1) + pkin(2) * t101 + t202;
t197 = mrSges(2,1) * t172 - mrSges(2,2) * t171 + Ifges(2,3) * qJDD(1) + pkin(1) * t93 + t199;
t94 = -mrSges(4,1) * t185 + mrSges(4,3) * t139 + t180 * Ifges(4,5) - pkin(3) * t110 + (Ifges(4,6) - t208) * t181 + t226;
t91 = m(2) * t172 + qJDD(1) * mrSges(2,1) - t196 * mrSges(2,2) + t93;
t90 = m(2) * t171 - t196 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t213;
t89 = mrSges(4,2) * t185 - mrSges(4,3) * t138 + Ifges(4,5) * t181 - t180 * Ifges(4,6) - qJ(4) * t110 + t188 * t103 - t186 * t105;
t88 = mrSges(3,2) * t185 - mrSges(3,3) * t149 + Ifges(3,5) * qJDD(1) - t196 * Ifges(3,6) - pkin(6) * t101 - t191 * t94 + t194 * t89;
t87 = -mrSges(3,1) * t185 + mrSges(3,3) * t150 + t196 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t214 + pkin(6) * t212 + t191 * t89 + t194 * t94;
t86 = -mrSges(2,2) * g(1) - mrSges(2,3) * t172 + Ifges(2,5) * qJDD(1) - t196 * Ifges(2,6) - qJ(2) * t93 - t187 * t87 + t189 * t88;
t85 = Ifges(2,6) * qJDD(1) + t196 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t171 + t187 * t88 + t189 * t87 - pkin(1) * (m(3) * t185 + t214) + qJ(2) * t213;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t197, t86, t88, t89, t103, t120; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t192 * t86 + t195 * t85 - pkin(5) * (t192 * t91 - t195 * t90), t85, t87, t94, t105, t119; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t195 * t86 + t192 * t85 + pkin(5) * (t192 * t90 + t195 * t91), t197, t199, t202, t208 * t181 - t226, t207 * t222 + t200;];
m_new = t1;
