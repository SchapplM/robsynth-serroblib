% Calculate vector of cutting torques with Newton-Euler for
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:04
% EndTime: 2019-12-05 15:24:10
% DurationCPUTime: 2.96s
% Computational Cost: add. (36205->193), mult. (68404->247), div. (0->0), fcn. (43606->10), ass. (0->93)
t198 = qJD(2) ^ 2;
t190 = sin(pkin(7));
t193 = cos(pkin(7));
t174 = -t193 * g(1) - t190 * g(2);
t187 = -g(3) + qJDD(1);
t195 = sin(qJ(2));
t197 = cos(qJ(2));
t160 = -t195 * t174 + t197 * t187;
t156 = qJDD(2) * pkin(2) + t160;
t161 = t197 * t174 + t195 * t187;
t157 = -t198 * pkin(2) + t161;
t189 = sin(pkin(8));
t192 = cos(pkin(8));
t142 = t189 * t156 + t192 * t157;
t138 = -t198 * pkin(3) + qJDD(2) * qJ(4) + t142;
t188 = sin(pkin(9));
t173 = t190 * g(1) - t193 * g(2);
t172 = qJDD(3) - t173;
t191 = cos(pkin(9));
t220 = qJD(2) * qJD(4);
t222 = t191 * t172 - 0.2e1 * t188 * t220;
t225 = pkin(4) * t191;
t131 = (-pkin(6) * qJDD(2) + t198 * t225 - t138) * t188 + t222;
t134 = t188 * t172 + (t138 + 0.2e1 * t220) * t191;
t219 = qJDD(2) * t191;
t184 = t191 ^ 2;
t223 = t184 * t198;
t132 = -pkin(4) * t223 + pkin(6) * t219 + t134;
t194 = sin(qJ(5));
t196 = cos(qJ(5));
t129 = t196 * t131 - t194 * t132;
t207 = -t188 * t194 + t191 * t196;
t162 = t207 * qJD(2);
t208 = t188 * t196 + t191 * t194;
t163 = t208 * qJD(2);
t147 = -t162 * mrSges(6,1) + t163 * mrSges(6,2);
t152 = t162 * qJD(5) + t208 * qJDD(2);
t158 = -qJD(5) * mrSges(6,2) + t162 * mrSges(6,3);
t124 = m(6) * t129 + qJDD(5) * mrSges(6,1) - t152 * mrSges(6,3) + qJD(5) * t158 - t163 * t147;
t130 = t194 * t131 + t196 * t132;
t151 = -t163 * qJD(5) + t207 * qJDD(2);
t159 = qJD(5) * mrSges(6,1) - t163 * mrSges(6,3);
t125 = m(6) * t130 - qJDD(5) * mrSges(6,2) + t151 * mrSges(6,3) - qJD(5) * t159 + t162 * t147;
t116 = t196 * t124 + t194 * t125;
t133 = -t188 * t138 + t222;
t144 = Ifges(6,4) * t163 + Ifges(6,2) * t162 + Ifges(6,6) * qJD(5);
t145 = Ifges(6,1) * t163 + Ifges(6,4) * t162 + Ifges(6,5) * qJD(5);
t202 = -mrSges(6,1) * t129 + mrSges(6,2) * t130 - Ifges(6,5) * t152 - Ifges(6,6) * t151 - Ifges(6,3) * qJDD(5) - t163 * t144 + t162 * t145;
t212 = Ifges(5,4) * t188 + Ifges(5,2) * t191;
t213 = Ifges(5,1) * t188 + Ifges(5,4) * t191;
t226 = -mrSges(5,1) * t133 + mrSges(5,2) * t134 - pkin(4) * t116 - (t188 * t212 - t191 * t213) * t198 + t202;
t224 = mrSges(5,2) * t188;
t206 = mrSges(5,3) * qJDD(2) + t198 * (-mrSges(5,1) * t191 + t224);
t114 = m(5) * t133 - t206 * t188 + t116;
t215 = -t194 * t124 + t196 * t125;
t115 = m(5) * t134 + t206 * t191 + t215;
t216 = -t188 * t114 + t191 * t115;
t105 = m(4) * t142 - t198 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t216;
t141 = t192 * t156 - t189 * t157;
t210 = qJDD(4) - t141;
t137 = -qJDD(2) * pkin(3) - t198 * qJ(4) + t210;
t183 = t188 ^ 2;
t135 = (-pkin(3) - t225) * qJDD(2) + (-qJ(4) + (-t183 - t184) * pkin(6)) * t198 + t210;
t203 = m(6) * t135 - t151 * mrSges(6,1) + t152 * mrSges(6,2) - t162 * t158 + t163 * t159;
t201 = -m(5) * t137 + mrSges(5,1) * t219 - t203 + (t183 * t198 + t223) * mrSges(5,3);
t120 = m(4) * t141 - t198 * mrSges(4,2) + (mrSges(4,1) - t224) * qJDD(2) + t201;
t100 = t189 * t105 + t192 * t120;
t109 = t191 * t114 + t188 * t115;
t211 = Ifges(5,5) * t188 + Ifges(5,6) * t191;
t221 = t198 * t211;
t96 = m(3) * t160 + qJDD(2) * mrSges(3,1) - t198 * mrSges(3,2) + t100;
t217 = t192 * t105 - t189 * t120;
t97 = m(3) * t161 - t198 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t217;
t218 = -t195 * t96 + t197 * t97;
t214 = m(4) * t172 + t109;
t143 = Ifges(6,5) * t163 + Ifges(6,6) * t162 + Ifges(6,3) * qJD(5);
t117 = -mrSges(6,1) * t135 + mrSges(6,3) * t130 + Ifges(6,4) * t152 + Ifges(6,2) * t151 + Ifges(6,6) * qJDD(5) + qJD(5) * t145 - t163 * t143;
t118 = mrSges(6,2) * t135 - mrSges(6,3) * t129 + Ifges(6,1) * t152 + Ifges(6,4) * t151 + Ifges(6,5) * qJDD(5) - qJD(5) * t144 + t162 * t143;
t102 = mrSges(5,2) * t137 - mrSges(5,3) * t133 - pkin(6) * t116 + t213 * qJDD(2) - t194 * t117 + t196 * t118 + t191 * t221;
t99 = -mrSges(5,1) * t137 + mrSges(5,3) * t134 - pkin(4) * t203 + pkin(6) * t215 + t212 * qJDD(2) + t196 * t117 + t194 * t118 - t188 * t221;
t90 = mrSges(4,2) * t172 - mrSges(4,3) * t141 + Ifges(4,5) * qJDD(2) - t198 * Ifges(4,6) - qJ(4) * t109 + t191 * t102 - t188 * t99;
t94 = (Ifges(4,6) - t211) * qJDD(2) + t198 * Ifges(4,5) - mrSges(4,1) * t172 + mrSges(4,3) * t142 - pkin(3) * t109 + t226;
t86 = mrSges(3,1) * t173 + mrSges(3,3) * t161 + t198 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t214 + qJ(3) * t217 + t189 * t90 + t192 * t94;
t89 = -mrSges(3,2) * t173 - mrSges(3,3) * t160 + Ifges(3,5) * qJDD(2) - t198 * Ifges(3,6) - qJ(3) * t100 - t189 * t94 + t192 * t90;
t205 = -mrSges(2,2) * t174 + pkin(5) * t218 + t195 * t89 + t197 * t86 + pkin(1) * (m(3) * t173 - t214) + mrSges(2,1) * t173;
t204 = -mrSges(4,2) * t142 + qJ(4) * t216 + t188 * t102 + t191 * t99 + pkin(3) * (-qJDD(2) * t224 + t201) + mrSges(4,1) * t141 + Ifges(4,3) * qJDD(2);
t199 = mrSges(3,1) * t160 - mrSges(3,2) * t161 + Ifges(3,3) * qJDD(2) + pkin(2) * t100 + t204;
t106 = (m(2) + m(3)) * t173 - t214;
t93 = t195 * t97 + t197 * t96;
t91 = m(2) * t174 + t218;
t87 = -mrSges(2,1) * t187 + mrSges(2,3) * t174 - pkin(1) * t93 - t199;
t84 = mrSges(2,2) * t187 - mrSges(2,3) * t173 - pkin(5) * t93 - t195 * t86 + t197 * t89;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t193 * t84 - t190 * t87 - qJ(1) * (t193 * t106 + t190 * t91), t84, t89, t90, t102, t118; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t190 * t84 + t193 * t87 + qJ(1) * (-t190 * t106 + t193 * t91), t87, t86, t94, t99, t117; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t205, t205, t199, t204, t211 * qJDD(2) - t226, -t202;];
m_new = t1;
