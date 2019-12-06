% Calculate vector of inverse dynamics joint torques for
% S5PRRRP3
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:31
% EndTime: 2019-12-05 16:43:45
% DurationCPUTime: 5.16s
% Computational Cost: add. (2267->329), mult. (4942->417), div. (0->0), fcn. (3121->8), ass. (0->144)
t222 = Ifges(5,4) + Ifges(6,4);
t137 = sin(qJ(3));
t139 = cos(qJ(3));
t110 = -mrSges(4,1) * t139 + mrSges(4,2) * t137;
t135 = qJ(3) + qJ(4);
t129 = sin(t135);
t130 = cos(t135);
t184 = mrSges(5,2) + mrSges(6,2);
t185 = mrSges(5,1) + mrSges(6,1);
t151 = t129 * t184 - t185 * t130;
t230 = -t110 - t151;
t223 = Ifges(5,1) + Ifges(6,1);
t221 = Ifges(5,5) + Ifges(6,5);
t220 = Ifges(5,2) + Ifges(6,2);
t219 = Ifges(5,6) + Ifges(6,6);
t136 = sin(qJ(4));
t138 = cos(qJ(4));
t95 = -t136 * t137 + t138 * t139;
t89 = t95 * qJD(2);
t229 = t222 * t89;
t96 = t136 * t139 + t137 * t138;
t90 = t96 * qJD(2);
t228 = t222 * t90;
t134 = qJD(3) + qJD(4);
t227 = t219 * t134 + t220 * t89 + t228;
t226 = t221 * t134 + t223 * t90 + t229;
t167 = qJD(2) * qJD(3);
t103 = qJDD(2) * t139 - t137 * t167;
t133 = pkin(8) + qJ(2);
t124 = sin(t133);
t125 = cos(t133);
t215 = g(1) * t125 + g(2) * t124;
t140 = -pkin(7) - pkin(6);
t113 = t140 * t139;
t175 = qJD(1) * t137;
t86 = -qJD(2) * t113 + t175;
t75 = t136 * t86;
t112 = t140 * t137;
t128 = t139 * qJD(1);
t85 = qJD(2) * t112 + t128;
t78 = qJD(3) * pkin(3) + t85;
t40 = t138 * t78 - t75;
t81 = t90 * qJ(5);
t12 = t40 - t81;
t11 = pkin(4) * t134 + t12;
t119 = pkin(3) * t139 + pkin(2);
t111 = t119 * qJD(2);
t62 = -pkin(4) * t89 + qJD(5) - t111;
t224 = t111 * mrSges(5,2) - t62 * mrSges(6,2) + t40 * mrSges(5,3) + t11 * mrSges(6,3);
t197 = t137 / 0.2e1;
t69 = -mrSges(6,2) * t134 + mrSges(6,3) * t89;
t70 = -mrSges(5,2) * t134 + mrSges(5,3) * t89;
t218 = t69 + t70;
t71 = mrSges(6,1) * t134 - t90 * mrSges(6,3);
t72 = mrSges(5,1) * t134 - t90 * mrSges(5,3);
t217 = t71 + t72;
t176 = qJDD(2) * pkin(2);
t168 = qJD(1) * qJD(3);
t65 = t103 * pkin(6) + t137 * qJDD(1) + t139 * t168;
t104 = qJDD(2) * t137 + t139 * t167;
t66 = -pkin(6) * t104 + t139 * qJDD(1) - t137 * t168;
t216 = -t137 * t66 + t139 * t65;
t180 = qJ(5) * t89;
t77 = t138 * t86;
t41 = t136 * t78 + t77;
t13 = t41 + t180;
t156 = t41 * mrSges(5,3) + t13 * mrSges(6,3);
t117 = pkin(4) * t130;
t214 = mrSges(3,1) + m(6) * (t117 + t119) + m(5) * t119 + m(4) * pkin(2) + t230;
t213 = mrSges(3,2) + m(6) * (-qJ(5) + t140) - mrSges(6,3) + m(5) * t140 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t205 = t90 / 0.2e1;
t202 = m(2) + m(3);
t201 = pkin(4) * t90;
t195 = pkin(3) * t137;
t192 = g(3) * t139;
t131 = qJDD(3) + qJDD(4);
t148 = t96 * qJD(4);
t37 = -qJD(2) * t148 + t103 * t138 - t104 * t136;
t25 = -mrSges(6,2) * t131 + mrSges(6,3) * t37;
t26 = -mrSges(5,2) * t131 + mrSges(5,3) * t37;
t183 = t25 + t26;
t45 = t138 * t85 - t75;
t64 = t136 * t112 - t138 * t113;
t182 = Ifges(4,4) * t137;
t181 = Ifges(4,4) * t139;
t178 = t139 * Ifges(4,2);
t174 = qJD(3) * t137;
t173 = qJD(3) * t139;
t172 = qJD(4) * t136;
t171 = qJD(4) * t138;
t170 = t137 * qJD(2);
t169 = t139 * qJD(2);
t165 = pkin(3) * t174;
t163 = qJD(3) * t140;
t147 = t95 * qJD(4);
t36 = qJD(2) * t147 + t103 * t136 + t104 * t138;
t160 = -t37 * mrSges(6,1) + t36 * mrSges(6,2);
t159 = t184 * t130;
t44 = -t136 * t85 - t77;
t63 = t138 * t112 + t113 * t136;
t155 = mrSges(4,1) * t137 + mrSges(4,2) * t139;
t154 = t178 + t182;
t153 = Ifges(4,5) * t139 - Ifges(4,6) * t137;
t82 = -pkin(3) * t103 - t176;
t105 = -pkin(6) * t170 + t128;
t106 = pkin(6) * t169 + t175;
t152 = t105 * t139 + t106 * t137;
t150 = pkin(2) * t155;
t43 = qJDD(3) * pkin(3) - pkin(7) * t104 + t66;
t55 = pkin(7) * t103 + t65;
t5 = t136 * t43 + t138 * t55 + t78 * t171 - t172 * t86;
t101 = t137 * t163;
t102 = t139 * t163;
t15 = t138 * t101 + t136 * t102 + t112 * t171 + t113 * t172;
t149 = t137 * (Ifges(4,1) * t139 - t182);
t6 = -qJD(4) * t41 - t136 * t55 + t138 * t43;
t16 = -qJD(4) * t64 - t101 * t136 + t138 * t102;
t2 = pkin(4) * t131 - qJ(5) * t36 - qJD(5) * t90 + t6;
t3 = qJ(5) * t37 + qJD(5) * t89 + t5;
t141 = -t5 * mrSges(5,2) - t3 * mrSges(6,2) + (-t62 * t90 + t2) * mrSges(6,1) + (t111 * t90 + t6) * mrSges(5,1) + t224 * t89 + t219 * t37 + t221 * t36 - (t223 * t89 - t228) * t90 / 0.2e1 + t227 * t205 - (-t219 * t90 + t221 * t89) * t134 / 0.2e1 + (Ifges(6,3) + Ifges(5,3)) * t131 - (-t220 * t90 + t226 + t229) * t89 / 0.2e1;
t121 = Ifges(4,4) * t169;
t109 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t169;
t108 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t170;
t88 = Ifges(4,1) * t170 + Ifges(4,5) * qJD(3) + t121;
t87 = Ifges(4,6) * qJD(3) + qJD(2) * t154;
t80 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t104;
t79 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t103;
t74 = -pkin(4) * t95 - t119;
t68 = pkin(3) * t170 + t201;
t61 = -qJD(3) * t96 - t148;
t60 = qJD(3) * t95 + t147;
t59 = -mrSges(5,1) * t89 + mrSges(5,2) * t90;
t58 = -mrSges(6,1) * t89 + mrSges(6,2) * t90;
t54 = -pkin(4) * t61 + t165;
t53 = qJ(5) * t95 + t64;
t52 = -qJ(5) * t96 + t63;
t24 = mrSges(5,1) * t131 - mrSges(5,3) * t36;
t23 = mrSges(6,1) * t131 - mrSges(6,3) * t36;
t18 = -t81 + t45;
t17 = t44 - t180;
t10 = -pkin(4) * t37 + qJDD(5) + t82;
t8 = -qJ(5) * t60 - qJD(5) * t96 + t16;
t7 = qJ(5) * t61 + qJD(5) * t95 + t15;
t1 = [t137 * t79 + t139 * t80 + t183 * t96 + (t23 + t24) * t95 + t217 * t61 + t218 * t60 + t202 * qJDD(1) + (-t108 * t137 + t109 * t139) * qJD(3) + m(4) * (t137 * t65 + t139 * t66 + (-t105 * t137 + t106 * t139) * qJD(3)) + m(5) * (t40 * t61 + t41 * t60 + t5 * t96 + t6 * t95) + m(6) * (t11 * t61 + t13 * t60 + t2 * t95 + t3 * t96) + (-m(4) - m(5) - m(6) - t202) * g(3); (t139 * (-Ifges(4,2) * t137 + t181) + t149) * t167 / 0.2e1 + m(5) * (-t111 * t165 - t119 * t82 + t15 * t41 + t16 * t40 + t5 * t64 + t6 * t63) + (t82 * mrSges(5,2) + t10 * mrSges(6,2) - t6 * mrSges(5,3) - t2 * mrSges(6,3) + t221 * t131 + t222 * t37 + t223 * t36) * t96 + (-t82 * mrSges(5,1) - t10 * mrSges(6,1) + t5 * mrSges(5,3) + t3 * mrSges(6,3) + t219 * t131 + t220 * t37 + t222 * t36) * t95 + m(6) * (t10 * t74 + t11 * t8 + t13 * t7 + t2 * t52 + t3 * t53 + t54 * t62) + qJD(3) ^ 2 * t153 / 0.2e1 + t103 * t154 / 0.2e1 + (t214 * t124 + t213 * t125) * g(1) + (t213 * t124 - t214 * t125) * g(2) + t139 * (Ifges(4,4) * t104 + Ifges(4,2) * t103) / 0.2e1 + (t219 * t61 + t221 * t60) * t134 / 0.2e1 + (t220 * t61 + t222 * t60) * t89 / 0.2e1 + (t222 * t61 + t223 * t60) * t205 - t119 * (-mrSges(5,1) * t37 + mrSges(5,2) * t36) + t226 * t60 / 0.2e1 + t227 * t61 / 0.2e1 + t63 * t24 + t64 * t26 + t7 * t69 + t15 * t70 + t8 * t71 + t16 * t72 + t52 * t23 + t53 * t25 + t54 * t58 - t224 * t60 + (Ifges(4,1) * t104 + Ifges(4,4) * t103) * t197 + (m(4) * t176 + mrSges(4,1) * t103 - mrSges(4,2) * t104) * pkin(2) + (t111 * mrSges(5,1) - t62 * mrSges(6,1) + t156) * t61 + (-t105 * t173 - t106 * t174 + t216) * mrSges(4,3) + (-t108 * t173 - t109 * t174 - t137 * t80 + t139 * t79 + m(4) * (-qJD(3) * t152 + t216)) * pkin(6) + t59 * t165 + (0.2e1 * Ifges(4,5) * t197 + Ifges(4,6) * t139) * qJDD(3) + t74 * t160 - t150 * t167 + t88 * t173 / 0.2e1 - t87 * t174 / 0.2e1 - t110 * t176 + Ifges(3,3) * qJDD(2) + t104 * (t137 * Ifges(4,1) + t181) / 0.2e1; -m(6) * (t11 * t17 + t13 * t18 + t62 * t68) + t156 * t90 - m(5) * (t40 * t44 + t41 * t45) + Ifges(4,5) * t104 + t106 * t108 - t105 * t109 + Ifges(4,6) * t103 - t65 * mrSges(4,2) + t66 * mrSges(4,1) - t68 * t58 - t18 * t69 - t45 * t70 - t17 * t71 - t44 * t72 + (-m(6) * t117 - t230) * g(3) + t141 + (-qJD(3) * t153 / 0.2e1 + t87 * t197 + (t178 * t197 - t149 / 0.2e1 + t150) * qJD(2) + (m(5) * t111 - t59) * t195 + t152 * mrSges(4,3) - (t121 + t88) * t139 / 0.2e1) * qJD(2) + (-m(6) * t192 + (t215 * t137 + t41 * t171 - t40 * t172 - t192) * m(5) + (t24 + (m(6) * t13 + t218) * qJD(4) + t6 * m(5)) * t138 + (m(6) * t3 + t183 + (-m(6) * t11 - t217) * qJD(4) + t5 * m(5)) * t136) * pkin(3) + Ifges(4,3) * qJDD(3) + t215 * (-m(6) * (-pkin(4) * t129 - t195) + t129 * t185 + t155 + t159) + (t2 * m(6) + t23) * (pkin(3) * t138 + pkin(4)); t151 * g(3) + (-pkin(4) * t58 + t156) * t90 - t12 * t69 - t40 * t70 + t13 * t71 + t41 * t72 + pkin(4) * t23 + t141 + (-g(3) * t117 - t62 * t201 - (-t11 + t12) * t13 + t2 * pkin(4)) * m(6) + t215 * (t159 + (m(6) * pkin(4) + t185) * t129); -t89 * t69 + t90 * t71 + (-g(1) * t124 + g(2) * t125 + t11 * t90 - t13 * t89 + t10) * m(6) + t160;];
tau = t1;
