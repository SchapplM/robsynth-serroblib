% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:36
% EndTime: 2019-12-05 15:08:40
% DurationCPUTime: 1.57s
% Computational Cost: add. (16107->199), mult. (27751->252), div. (0->0), fcn. (15841->8), ass. (0->80)
t177 = sin(pkin(7));
t179 = cos(pkin(7));
t164 = -t179 * g(1) - t177 * g(2);
t175 = -g(3) + qJDD(1);
t176 = sin(pkin(8));
t178 = cos(pkin(8));
t133 = -t176 * t164 + t178 * t175;
t134 = t178 * t164 + t176 * t175;
t181 = sin(qJ(3));
t183 = cos(qJ(3));
t129 = t181 * t133 + t183 * t134;
t185 = qJD(3) ^ 2;
t126 = -t185 * pkin(3) + qJDD(3) * pkin(6) + t129;
t180 = sin(qJ(4));
t163 = t177 * g(1) - t179 * g(2);
t162 = qJDD(2) - t163;
t182 = cos(qJ(4));
t205 = t182 * t162;
t122 = -t180 * t126 + t205;
t123 = t182 * t126 + t180 * t162;
t138 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t180 - Ifges(6,3) * t182) * qJD(3);
t141 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t180 + Ifges(5,2) * t182) * qJD(3);
t155 = (-mrSges(6,1) * t182 - mrSges(6,3) * t180) * qJD(3);
t200 = qJD(3) * qJD(4);
t157 = t180 * qJDD(3) + t182 * t200;
t158 = t182 * qJDD(3) - t180 * t200;
t154 = (-pkin(4) * t182 - qJ(5) * t180) * qJD(3);
t184 = qJD(4) ^ 2;
t201 = qJD(3) * t182;
t119 = -t184 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t154 * t201 + t123;
t121 = -qJDD(4) * pkin(4) - t184 * qJ(5) - t205 + qJDD(5) + (qJD(3) * t154 + t126) * t180;
t192 = -mrSges(6,1) * t121 + mrSges(6,3) * t119 + Ifges(6,4) * t157 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t158;
t168 = mrSges(6,2) * t201 + qJD(4) * mrSges(6,3);
t195 = -m(6) * t121 + qJDD(4) * mrSges(6,1) + qJD(4) * t168;
t202 = qJD(3) * t180;
t166 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t202;
t196 = m(6) * t119 + qJDD(4) * mrSges(6,3) + qJD(4) * t166 + t155 * t201;
t142 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t180 - Ifges(6,5) * t182) * qJD(3);
t203 = t142 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t180 + Ifges(5,4) * t182) * qJD(3);
t208 = -((t138 - t141) * t180 + t203 * t182) * qJD(3) + mrSges(5,1) * t122 - mrSges(5,2) * t123 + Ifges(5,5) * t157 + Ifges(5,6) * t158 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t157 * mrSges(6,2) - t155 * t202 + t195) + qJ(5) * (t158 * mrSges(6,2) + t196) + t192;
t206 = mrSges(5,3) + mrSges(6,2);
t128 = t183 * t133 - t181 * t134;
t125 = -qJDD(3) * pkin(3) - t185 * pkin(6) - t128;
t116 = -t158 * pkin(4) - t157 * qJ(5) + (-0.2e1 * qJD(5) * t180 + (pkin(4) * t180 - qJ(5) * t182) * qJD(4)) * qJD(3) + t125;
t113 = m(6) * t116 - t158 * mrSges(6,1) - t157 * mrSges(6,3) - t166 * t202 - t168 * t201;
t165 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t202;
t167 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t201;
t187 = -m(5) * t125 + t158 * mrSges(5,1) - t157 * mrSges(5,2) - t165 * t202 + t167 * t201 - t113;
t105 = m(4) * t128 + qJDD(3) * mrSges(4,1) - t185 * mrSges(4,2) + t187;
t156 = (-mrSges(5,1) * t182 + mrSges(5,2) * t180) * qJD(3);
t111 = m(5) * t123 - qJDD(4) * mrSges(5,2) - qJD(4) * t165 + t156 * t201 + t206 * t158 + t196;
t112 = m(5) * t122 + qJDD(4) * mrSges(5,1) + qJD(4) * t167 - t206 * t157 + (-t155 - t156) * t202 + t195;
t197 = t182 * t111 - t180 * t112;
t98 = m(4) * t129 - t185 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t197;
t91 = t183 * t105 + t181 * t98;
t102 = t180 * t111 + t182 * t112;
t89 = m(3) * t133 + t91;
t198 = -t181 * t105 + t183 * t98;
t90 = m(3) * t134 + t198;
t199 = -t176 * t89 + t178 * t90;
t194 = -mrSges(6,1) * t116 + mrSges(6,2) * t119;
t193 = (-m(3) - m(4)) * t162 - t102;
t139 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t180 + Ifges(5,6) * t182) * qJD(3);
t140 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t180 - Ifges(6,6) * t182) * qJD(3);
t94 = -mrSges(5,1) * t125 + mrSges(5,3) * t123 - pkin(4) * t113 + (Ifges(5,2) + Ifges(6,3)) * t158 + (Ifges(5,4) - Ifges(6,5)) * t157 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t203 * qJD(4) + (-t139 - t140) * t202 + t194;
t189 = mrSges(6,2) * t121 - mrSges(6,3) * t116 + Ifges(6,1) * t157 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t158 + qJD(4) * t138 + t140 * t201;
t95 = mrSges(5,2) * t125 - mrSges(5,3) * t122 + Ifges(5,1) * t157 + Ifges(5,4) * t158 + Ifges(5,5) * qJDD(4) - qJ(5) * t113 - qJD(4) * t141 + t139 * t201 + t189;
t86 = mrSges(4,2) * t162 - mrSges(4,3) * t128 + Ifges(4,5) * qJDD(3) - t185 * Ifges(4,6) - pkin(6) * t102 - t180 * t94 + t182 * t95;
t87 = -mrSges(4,1) * t162 + mrSges(4,3) * t129 + t185 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t102 - t208;
t80 = -mrSges(3,1) * t162 + mrSges(3,3) * t134 + t181 * t86 + t183 * t87 - pkin(2) * (m(4) * t162 + t102) + pkin(5) * t198;
t82 = mrSges(3,2) * t162 - mrSges(3,3) * t133 - pkin(5) * t91 - t181 * t87 + t183 * t86;
t191 = mrSges(2,1) * t163 - mrSges(2,2) * t164 + pkin(1) * t193 + qJ(2) * t199 + t176 * t82 + t178 * t80;
t190 = mrSges(4,1) * t128 - mrSges(4,2) * t129 + Ifges(4,3) * qJDD(3) + pkin(3) * t187 + pkin(6) * t197 + t180 * t95 + t182 * t94;
t188 = mrSges(3,1) * t133 - mrSges(3,2) * t134 + pkin(2) * t91 + t190;
t99 = m(2) * t163 + t193;
t85 = t176 * t90 + t178 * t89;
t83 = m(2) * t164 + t199;
t78 = -mrSges(2,1) * t175 + mrSges(2,3) * t164 - pkin(1) * t85 - t188;
t77 = mrSges(2,2) * t175 - mrSges(2,3) * t163 - qJ(2) * t85 - t176 * t80 + t178 * t82;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t179 * t77 - t177 * t78 - qJ(1) * (t177 * t83 + t179 * t99), t77, t82, t86, t95, t189; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t177 * t77 + t179 * t78 + qJ(1) * (-t177 * t99 + t179 * t83), t78, t80, t87, t94, (-t180 * t138 - t182 * t142) * qJD(3) + t192; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t191, t191, t188, t190, t208, Ifges(6,5) * t157 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t158 - qJD(4) * t142 + t140 * t202 - t194;];
m_new = t1;
