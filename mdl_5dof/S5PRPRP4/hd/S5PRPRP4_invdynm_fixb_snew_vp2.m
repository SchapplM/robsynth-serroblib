% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRP4
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:08
% EndTime: 2019-12-05 15:35:11
% DurationCPUTime: 1.69s
% Computational Cost: add. (17715->210), mult. (30185->263), div. (0->0), fcn. (15841->8), ass. (0->82)
t188 = sin(pkin(7));
t190 = cos(pkin(7));
t172 = -t190 * g(1) - t188 * g(2);
t186 = -g(3) + qJDD(1);
t192 = sin(qJ(2));
t194 = cos(qJ(2));
t140 = -t192 * t172 + t194 * t186;
t138 = qJDD(2) * pkin(2) + t140;
t141 = t194 * t172 + t192 * t186;
t196 = qJD(2) ^ 2;
t139 = -t196 * pkin(2) + t141;
t187 = sin(pkin(8));
t189 = cos(pkin(8));
t134 = t187 * t138 + t189 * t139;
t131 = -t196 * pkin(3) + qJDD(2) * pkin(6) + t134;
t191 = sin(qJ(4));
t171 = t188 * g(1) - t190 * g(2);
t170 = qJDD(3) - t171;
t193 = cos(qJ(4));
t216 = t193 * t170;
t127 = -t191 * t131 + t216;
t128 = t193 * t131 + t191 * t170;
t145 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t191 - Ifges(6,3) * t193) * qJD(2);
t148 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t191 + Ifges(5,2) * t193) * qJD(2);
t163 = (-mrSges(6,1) * t193 - mrSges(6,3) * t191) * qJD(2);
t211 = qJD(2) * qJD(4);
t165 = t191 * qJDD(2) + t193 * t211;
t166 = t193 * qJDD(2) - t191 * t211;
t162 = (-pkin(4) * t193 - qJ(5) * t191) * qJD(2);
t195 = qJD(4) ^ 2;
t212 = qJD(2) * t193;
t124 = -t195 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t162 * t212 + t128;
t126 = -qJDD(4) * pkin(4) - t195 * qJ(5) - t216 + qJDD(5) + (qJD(2) * t162 + t131) * t191;
t203 = -mrSges(6,1) * t126 + mrSges(6,3) * t124 + Ifges(6,4) * t165 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t166;
t176 = mrSges(6,2) * t212 + qJD(4) * mrSges(6,3);
t205 = -m(6) * t126 + qJDD(4) * mrSges(6,1) + qJD(4) * t176;
t213 = qJD(2) * t191;
t174 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t213;
t207 = m(6) * t124 + qJDD(4) * mrSges(6,3) + qJD(4) * t174 + t163 * t212;
t149 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t191 - Ifges(6,5) * t193) * qJD(2);
t214 = t149 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t191 + Ifges(5,4) * t193) * qJD(2);
t219 = -((t145 - t148) * t191 + t193 * t214) * qJD(2) + mrSges(5,1) * t127 - mrSges(5,2) * t128 + Ifges(5,5) * t165 + Ifges(5,6) * t166 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t165 * mrSges(6,2) - t163 * t213 + t205) + qJ(5) * (t166 * mrSges(6,2) + t207) + t203;
t217 = mrSges(5,3) + mrSges(6,2);
t164 = (-mrSges(5,1) * t193 + mrSges(5,2) * t191) * qJD(2);
t173 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t213;
t116 = m(5) * t128 - qJDD(4) * mrSges(5,2) - qJD(4) * t173 + t164 * t212 + t166 * t217 + t207;
t175 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t212;
t117 = m(5) * t127 + qJDD(4) * mrSges(5,1) + qJD(4) * t175 - t217 * t165 + (-t163 - t164) * t213 + t205;
t208 = t193 * t116 - t191 * t117;
t103 = m(4) * t134 - t196 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t208;
t133 = t189 * t138 - t187 * t139;
t130 = -qJDD(2) * pkin(3) - t196 * pkin(6) - t133;
t121 = -t166 * pkin(4) - t165 * qJ(5) + (-0.2e1 * qJD(5) * t191 + (pkin(4) * t191 - qJ(5) * t193) * qJD(4)) * qJD(2) + t130;
t118 = m(6) * t121 - t166 * mrSges(6,1) - t165 * mrSges(6,3) - t174 * t213 - t176 * t212;
t198 = -m(5) * t130 + t166 * mrSges(5,1) - t165 * mrSges(5,2) - t173 * t213 + t175 * t212 - t118;
t110 = m(4) * t133 + qJDD(2) * mrSges(4,1) - t196 * mrSges(4,2) + t198;
t96 = t187 * t103 + t189 * t110;
t107 = t191 * t116 + t193 * t117;
t94 = m(3) * t140 + qJDD(2) * mrSges(3,1) - t196 * mrSges(3,2) + t96;
t209 = t189 * t103 - t187 * t110;
t95 = m(3) * t141 - t196 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t209;
t210 = -t192 * t94 + t194 * t95;
t206 = m(4) * t170 + t107;
t204 = -mrSges(6,1) * t121 + mrSges(6,2) * t124;
t146 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t191 + Ifges(5,6) * t193) * qJD(2);
t147 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t191 - Ifges(6,6) * t193) * qJD(2);
t200 = mrSges(6,2) * t126 - mrSges(6,3) * t121 + Ifges(6,1) * t165 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t166 + qJD(4) * t145 + t147 * t212;
t100 = mrSges(5,2) * t130 - mrSges(5,3) * t127 + Ifges(5,1) * t165 + Ifges(5,4) * t166 + Ifges(5,5) * qJDD(4) - qJ(5) * t118 - qJD(4) * t148 + t146 * t212 + t200;
t99 = -mrSges(5,1) * t130 + mrSges(5,3) * t128 - pkin(4) * t118 + (Ifges(5,2) + Ifges(6,3)) * t166 + (Ifges(5,4) - Ifges(6,5)) * t165 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t214 * qJD(4) + (-t146 - t147) * t213 + t204;
t91 = mrSges(4,2) * t170 - mrSges(4,3) * t133 + Ifges(4,5) * qJDD(2) - t196 * Ifges(4,6) - pkin(6) * t107 + t193 * t100 - t191 * t99;
t92 = -mrSges(4,1) * t170 + mrSges(4,3) * t134 + t196 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t107 - t219;
t85 = mrSges(3,1) * t171 + mrSges(3,3) * t141 + t196 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t206 + qJ(3) * t209 + t187 * t91 + t189 * t92;
t87 = -mrSges(3,2) * t171 - mrSges(3,3) * t140 + Ifges(3,5) * qJDD(2) - t196 * Ifges(3,6) - qJ(3) * t96 - t187 * t92 + t189 * t91;
t202 = -mrSges(2,2) * t172 + pkin(5) * t210 + t192 * t87 + t194 * t85 + pkin(1) * (m(3) * t171 - t206) + mrSges(2,1) * t171;
t201 = mrSges(4,1) * t133 - mrSges(4,2) * t134 + Ifges(4,3) * qJDD(2) + pkin(3) * t198 + pkin(6) * t208 + t191 * t100 + t193 * t99;
t199 = mrSges(3,1) * t140 - mrSges(3,2) * t141 + Ifges(3,3) * qJDD(2) + pkin(2) * t96 + t201;
t104 = (m(2) + m(3)) * t171 - t206;
t90 = t192 * t95 + t194 * t94;
t88 = m(2) * t172 + t210;
t83 = -mrSges(2,1) * t186 + mrSges(2,3) * t172 - pkin(1) * t90 - t199;
t82 = mrSges(2,2) * t186 - mrSges(2,3) * t171 - pkin(5) * t90 - t192 * t85 + t194 * t87;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t190 * t82 - t188 * t83 - qJ(1) * (t190 * t104 + t188 * t88), t82, t87, t91, t100, t200; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t188 * t82 + t190 * t83 + qJ(1) * (-t188 * t104 + t190 * t88), t83, t85, t92, t99, (-t191 * t145 - t193 * t149) * qJD(2) + t203; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t202, t202, t199, t201, t219, Ifges(6,5) * t165 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t166 - qJD(4) * t149 + t147 * t213 - t204;];
m_new = t1;
