% Calculate vector of inverse dynamics joint torques for
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:18
% EndTime: 2019-12-05 16:17:25
% DurationCPUTime: 2.20s
% Computational Cost: add. (1706->244), mult. (2491->348), div. (0->0), fcn. (1389->10), ass. (0->129)
t110 = qJD(2) + qJD(3);
t112 = cos(pkin(9));
t82 = -t110 * t112 + qJD(5);
t218 = Ifges(6,5) * t82;
t113 = sin(qJ(5));
t217 = -t113 / 0.2e1;
t216 = -mrSges(5,3) + mrSges(4,2);
t111 = sin(pkin(9));
t107 = t111 ^ 2;
t215 = (t112 ^ 2 + t107) * mrSges(5,3);
t115 = cos(qJ(5));
t190 = Ifges(6,4) * t115;
t121 = (-t113 * Ifges(6,2) + t190) * t111;
t191 = Ifges(6,4) * t113;
t122 = (Ifges(6,1) * t115 - t191) * t111;
t116 = cos(qJ(3));
t189 = pkin(2) * qJD(2);
t155 = t116 * t189;
t137 = qJD(4) - t155;
t78 = -t112 * pkin(4) - t111 * pkin(7) - pkin(3);
t37 = t110 * t78 + t137;
t114 = sin(qJ(3));
t156 = t114 * t189;
t73 = qJ(4) * t110 + t156;
t59 = qJD(1) * t111 + t112 * t73;
t10 = t113 * t37 + t115 * t59;
t9 = -t113 * t59 + t115 * t37;
t135 = t10 * t115 - t113 * t9;
t58 = -t112 * qJD(1) + t111 * t73;
t214 = -t135 * mrSges(6,3) + (t110 * t122 + t218) * t217 + (-t58 * mrSges(6,2) - t218 / 0.2e1) * t113 + (t58 * mrSges(6,1) - Ifges(6,6) * t82 - t110 * t121 / 0.2e1) * t115;
t106 = qJDD(2) + qJDD(3);
t166 = qJD(4) * t110;
t188 = pkin(2) * qJD(3);
t148 = qJD(2) * t188;
t175 = pkin(2) * qJDD(2);
t71 = t114 * t175 + t116 * t148;
t38 = qJ(4) * t106 + t166 + t71;
t25 = qJDD(1) * t111 + t112 * t38;
t179 = t112 * t25;
t24 = -t112 * qJDD(1) + t111 * t38;
t212 = t111 * t24 + t179;
t134 = -t112 * mrSges(5,1) + mrSges(5,2) * t111;
t210 = t115 * (-Ifges(6,1) * t113 - t190) / 0.2e1 + (-Ifges(6,2) * t115 - t191) * t217;
t133 = mrSges(6,1) * t113 + mrSges(6,2) * t115;
t172 = t110 * t111;
t53 = t133 * t172;
t184 = t111 * t53;
t209 = -t110 * t215 - t184;
t165 = qJD(4) * t112;
t167 = t112 * t116;
t168 = t112 * t115;
t55 = qJ(4) * t168 + t113 * t78;
t208 = -qJD(5) * t55 - t113 * t165 - (-t113 * t167 + t114 * t115) * t189;
t169 = t112 * t113;
t54 = -qJ(4) * t169 + t115 * t78;
t207 = qJD(5) * t54 + t115 * t165 - (t113 * t114 + t115 * t167) * t189;
t205 = t112 * t59;
t109 = pkin(8) + qJ(2);
t104 = qJ(3) + t109;
t94 = sin(t104);
t95 = cos(t104);
t49 = t115 * t95 + t169 * t94;
t50 = t113 * t95 - t168 * t94;
t202 = -t50 * mrSges(6,1) - t49 * mrSges(6,2) + t216 * t95 + (-m(6) * t78 + t111 * mrSges(6,3) + mrSges(4,1) - t134) * t94;
t178 = t112 * t95;
t182 = t111 * t95;
t51 = t115 * t94 - t169 * t95;
t52 = t113 * t94 + t168 * t95;
t201 = -t95 * mrSges(4,1) - mrSges(5,1) * t178 - t52 * mrSges(6,1) - t51 * mrSges(6,2) + t216 * t94 + (mrSges(5,2) - mrSges(6,3)) * t182;
t164 = qJD(5) * t110;
t56 = (t106 * t115 - t113 * t164) * t111;
t57 = (-t106 * t113 - t115 * t164) * t111;
t14 = -mrSges(6,1) * t57 + mrSges(6,2) * t56;
t200 = t106 * t215 + t111 * t14;
t199 = pkin(3) * t94;
t102 = sin(t109);
t198 = pkin(2) * t102;
t103 = cos(t109);
t93 = pkin(2) * t103;
t197 = pkin(2) * t114;
t196 = pkin(2) * t116;
t194 = t95 * pkin(3) + t94 * qJ(4);
t183 = t111 * t58;
t180 = t112 * t24;
t174 = t106 * t111;
t173 = t106 * t112;
t171 = t111 * t113;
t170 = t111 * t115;
t162 = m(2) + m(3) + m(4);
t81 = qJDD(5) - t173;
t160 = Ifges(6,5) * t56 + Ifges(6,6) * t57 + Ifges(6,3) * t81;
t158 = mrSges(6,3) * t171;
t157 = mrSges(6,3) * t170;
t154 = t114 * t188;
t153 = t116 * t188;
t85 = t95 * qJ(4);
t145 = t85 - t198;
t143 = pkin(4) * t178 + pkin(7) * t182 + t194;
t142 = t110 * mrSges(4,1) * t197;
t141 = t110 * t155;
t140 = t111 * t155;
t65 = -mrSges(5,1) * t173 + mrSges(5,2) * t174;
t139 = -g(1) * t94 + g(2) * t95;
t70 = -t114 * t148 + t116 * t175;
t132 = t183 + t205;
t47 = -mrSges(6,2) * t82 - t110 * t158;
t48 = mrSges(6,1) * t82 - t110 * t157;
t131 = t113 * t48 - t115 * t47;
t128 = qJDD(4) - t70;
t69 = t78 - t196;
t96 = qJ(4) + t197;
t32 = t113 * t69 + t168 * t96;
t31 = t115 * t69 - t169 * t96;
t89 = qJD(4) + t153;
t127 = (t24 * t96 + t58 * t89) * t111;
t21 = t106 * t78 + t128;
t3 = qJD(5) * t9 + t113 * t21 + t115 * t25;
t4 = -qJD(5) * t10 - t113 * t25 + t115 * t21;
t60 = -t106 * pkin(3) + t128;
t67 = t133 * t111;
t117 = t81 * (-Ifges(6,3) * t112 + (Ifges(6,5) * t115 - Ifges(6,6) * t113) * t111) / 0.2e1 + Ifges(4,3) * t106 + t70 * mrSges(4,1) - t71 * mrSges(4,2) + t24 * t67 + t57 * (-Ifges(6,6) * t112 + t121) / 0.2e1 + t56 * (-Ifges(6,5) * t112 + t122) / 0.2e1 + t60 * t134 + t4 * (-mrSges(6,1) * t112 - t157) + t3 * (mrSges(6,2) * t112 - t158) - t112 * t160 / 0.2e1 + (Ifges(6,1) * t56 + Ifges(6,4) * t57 + Ifges(6,5) * t81) * t170 / 0.2e1 - (Ifges(6,4) * t56 + Ifges(6,2) * t57 + Ifges(6,6) * t81) * t171 / 0.2e1 + (Ifges(5,4) * t111 + Ifges(5,2) * t112) * t173 + (Ifges(5,1) * t111 + Ifges(5,4) * t112) * t174 + t210 * t107 * t164 + t212 * mrSges(5,3) + t214 * qJD(5) * t111;
t97 = -pkin(3) - t196;
t72 = -t110 * pkin(3) + t137;
t66 = t134 * t110;
t23 = -mrSges(6,2) * t81 + mrSges(6,3) * t57;
t22 = mrSges(6,1) * t81 - mrSges(6,3) * t56;
t12 = -qJD(5) * t32 + t115 * t154 - t169 * t89;
t11 = qJD(5) * t31 + t113 * t154 + t168 * t89;
t1 = [-t112 * t14 + (-t113 * t22 + t115 * t23 + (-t113 * t47 - t115 * t48) * qJD(5)) * t111 + m(5) * (t111 * t25 - t180) + m(6) * (-t180 + (-t113 * t4 + t115 * t3 + (-t10 * t113 - t115 * t9) * qJD(5)) * t111) + t162 * qJDD(1) + (-m(5) - m(6) - t162) * g(3); t97 * t65 + t11 * t47 + t12 * t48 + t31 * t22 + t32 * t23 + t117 + Ifges(3,3) * qJDD(2) + m(4) * (t114 * t71 + t116 * t70) * pkin(2) + m(6) * (t10 * t11 + t12 * t9 + t3 * t32 + t31 * t4 + t127) + t106 * mrSges(4,1) * t196 + t66 * t154 - qJD(3) * t142 + m(5) * (t72 * t154 + t60 * t97 + t127) + (m(5) * t179 + t200) * t96 + (-t106 * t197 - t110 * t153) * mrSges(4,2) + (-t103 * mrSges(3,1) + t102 * mrSges(3,2) - m(6) * (t93 + t143) - m(5) * (t93 + t194) - m(4) * t93 + t201) * g(2) + (t102 * mrSges(3,1) + t103 * mrSges(3,2) - m(6) * t145 + m(4) * t198 - m(5) * (t145 - t199) + t202) * g(1) + (m(5) * t205 - t209) * t89; mrSges(4,2) * t141 - pkin(3) * t65 + qJD(2) * t142 + qJD(4) * t184 - t53 * t140 - t66 * t156 + t54 * t22 + t55 * t23 + t117 + t208 * t48 + t207 * t47 + t200 * qJ(4) + (-t58 * t140 + t3 * t55 + t4 * t54 + (qJ(4) * t24 + qJD(4) * t58) * t111 + t208 * t9 + t207 * t10) * m(6) + (-(t114 * t72 + t116 * t132) * t189 - pkin(3) * t60 + t132 * qJD(4) + t212 * qJ(4)) * m(5) + (-m(5) * t194 - m(6) * t143 + t201) * g(2) + (-m(6) * t85 - m(5) * (t85 - t199) + t202) * g(1) + (-t141 + t166) * t215; t113 * t23 + t115 * t22 - t131 * qJD(5) + (t135 * qJD(5) + t113 * t3 + t115 * t4 + t139) * m(6) + (t139 + t60) * m(5) + t65 + (t131 * t112 - m(5) * t132 - m(6) * (t10 * t168 - t169 * t9 + t183) + t209) * t110; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t9 * t47 + t10 * t48 - g(1) * (mrSges(6,1) * t51 - mrSges(6,2) * t52) - g(2) * (-mrSges(6,1) * t49 + mrSges(6,2) * t50) + g(3) * t67 + (-t210 * t172 - t214) * t172 + t160;];
tau = t1;
