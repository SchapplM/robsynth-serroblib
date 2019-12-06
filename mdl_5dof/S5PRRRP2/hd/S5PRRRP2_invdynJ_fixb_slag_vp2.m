% Calculate vector of inverse dynamics joint torques for
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:37
% EndTime: 2019-12-05 16:41:44
% DurationCPUTime: 2.27s
% Computational Cost: add. (1429->251), mult. (2140->317), div. (0->0), fcn. (969->8), ass. (0->118)
t229 = Ifges(6,4) + Ifges(5,5);
t228 = Ifges(6,6) - Ifges(5,6);
t232 = mrSges(5,1) + mrSges(6,1);
t107 = sin(qJ(4));
t109 = cos(qJ(4));
t134 = t109 * mrSges(6,1) + t107 * mrSges(6,3);
t77 = -t109 * mrSges(5,1) + t107 * mrSges(5,2);
t231 = -t134 + t77;
t106 = qJD(2) + qJD(3);
t173 = t106 * t109;
t162 = mrSges(5,3) * t173;
t154 = mrSges(6,2) * t173;
t67 = qJD(4) * mrSges(6,3) + t154;
t216 = -qJD(4) * mrSges(5,2) + t162 + t67;
t174 = t106 * t107;
t155 = mrSges(5,3) * t174;
t163 = mrSges(6,2) * t174;
t217 = t232 * qJD(4) - t155 - t163;
t222 = -t217 * t107 + t216 * t109;
t230 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t186 = Ifges(6,5) * t109;
t132 = Ifges(6,1) * t107 - t186;
t79 = Ifges(5,4) * t173;
t227 = Ifges(5,1) * t174 + t229 * qJD(4) + t106 * t132 + t79;
t133 = t107 * mrSges(6,1) - t109 * mrSges(6,3);
t135 = mrSges(5,1) * t107 + mrSges(5,2) * t109;
t110 = cos(qJ(3));
t185 = pkin(2) * qJD(2);
t160 = t110 * t185;
t172 = t107 * qJ(5);
t128 = t109 * pkin(4) + t172;
t75 = -pkin(3) - t128;
t23 = t106 * t75 - t160;
t69 = -t106 * pkin(3) - t160;
t226 = t23 * t133 + t69 * t135;
t225 = t228 * t107 + t229 * t109;
t169 = qJD(4) * t109;
t104 = qJDD(2) + qJDD(3);
t108 = sin(qJ(3));
t184 = pkin(2) * qJD(3);
t153 = qJD(2) * t184;
t175 = pkin(2) * qJDD(2);
t61 = t108 * t175 + t110 * t153;
t48 = pkin(7) * t104 + t61;
t164 = qJD(1) * t169 + t107 * qJDD(1) + t109 * t48;
t170 = qJD(4) * t107;
t161 = t108 * t185;
t68 = pkin(7) * t106 + t161;
t8 = -t170 * t68 + t164;
t42 = qJD(1) * t107 + t109 * t68;
t9 = -qJD(4) * t42 + qJDD(1) * t109 - t107 * t48;
t224 = -t9 * t107 + t109 * t8;
t182 = t107 * t68;
t5 = qJDD(4) * qJ(5) + (qJD(5) - t182) * qJD(4) + t164;
t6 = -qJDD(4) * pkin(4) + qJDD(5) - t9;
t223 = t107 * t6 + t109 * t5;
t60 = -t108 * t153 + t110 * t175;
t47 = -t104 * pkin(3) - t60;
t54 = -t109 * t104 + t106 * t170;
t55 = t104 * t107 + t106 * t169;
t221 = m(5) * t47 + mrSges(5,1) * t54 + mrSges(5,2) * t55;
t35 = -qJDD(4) * mrSges(6,1) + t55 * mrSges(6,2);
t218 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t55 - t35;
t36 = -mrSges(6,2) * t54 + qJDD(4) * mrSges(6,3);
t219 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t54 + t36;
t41 = qJD(1) * t109 - t182;
t28 = -qJD(4) * pkin(4) + qJD(5) - t41;
t30 = qJD(4) * qJ(5) + t42;
t220 = m(5) * ((-t107 * t42 - t109 * t41) * qJD(4) + t224) + m(6) * ((-t107 * t30 + t109 * t28) * qJD(4) + t223) - t218 * t107 + t219 * t109 - t170 * t216 - t169 * t217;
t105 = pkin(8) + qJ(2);
t103 = qJ(3) + t105;
t92 = sin(t103);
t93 = cos(t103);
t215 = g(1) * t93 + g(2) * t92;
t214 = t230 * t93 + (-m(6) * t75 + mrSges(4,1) - t231) * t92;
t178 = t109 * t93;
t213 = -t232 * t178 + t230 * t92 + (-mrSges(4,1) + (mrSges(5,2) - mrSges(6,3)) * t107) * t93;
t125 = -t107 * t41 + t109 * t42;
t197 = pkin(2) * t110;
t198 = pkin(2) * t108;
t206 = m(5) * pkin(2) * (t108 * t69 + t110 * t125) + (-mrSges(4,1) * t198 - mrSges(4,2) * t197) * t106;
t203 = pkin(3) * t92;
t200 = t107 / 0.2e1;
t100 = sin(t105);
t199 = pkin(2) * t100;
t101 = cos(t105);
t91 = pkin(2) * t101;
t190 = t93 * pkin(3) + t92 * pkin(7);
t189 = Ifges(5,4) * t107;
t188 = Ifges(5,4) * t109;
t187 = Ifges(6,5) * t107;
t171 = qJD(4) * t106;
t168 = qJD(5) * t107;
t167 = m(2) + m(3) + m(4);
t159 = t108 * t184;
t86 = t93 * pkin(7);
t147 = t86 - t199;
t146 = -t171 / 0.2e1;
t143 = pkin(4) * t178 + t93 * t172 + t190;
t131 = Ifges(5,2) * t109 + t189;
t127 = pkin(4) * t107 - qJ(5) * t109;
t126 = t107 * t28 + t109 * t30;
t119 = t107 * (Ifges(5,1) * t109 - t189);
t118 = t109 * (Ifges(6,3) * t107 + t186);
t117 = t126 * t110;
t50 = pkin(4) * t170 - qJ(5) * t169 - t168;
t2 = t54 * pkin(4) - t55 * qJ(5) - t106 * t168 + t47;
t78 = Ifges(6,5) * t174;
t43 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t173 + t78;
t44 = Ifges(5,6) * qJD(4) + t106 * t131;
t111 = -t44 * t170 / 0.2e1 + t119 * t171 / 0.2e1 + t109 * (Ifges(5,4) * t55 + Ifges(5,6) * qJDD(4)) / 0.2e1 - t109 * (Ifges(6,5) * t55 + Ifges(6,6) * qJDD(4)) / 0.2e1 + Ifges(4,3) * t104 + t60 * mrSges(4,1) - t61 * mrSges(4,2) + t47 * t77 - t2 * t134 + t118 * t146 + (Ifges(5,1) * t107 + t132 + t188) * t55 / 0.2e1 + ((Ifges(5,1) + Ifges(6,1)) * t55 + t229 * qJDD(4)) * t200 + (t229 * t107 - t228 * t109) * qJDD(4) / 0.2e1 + (t106 * (Ifges(6,1) * t109 + t187) + t43) * t170 / 0.2e1 + (-t41 * t169 - t42 * t170 + t224) * mrSges(5,3) + (t187 / 0.2e1 - t131 / 0.2e1 + (Ifges(6,5) - Ifges(5,4)) * t200 + (-Ifges(6,3) - Ifges(5,2) / 0.2e1) * t109) * t54 + (t106 * (-Ifges(5,2) * t107 + t188) + t227) * t169 / 0.2e1 + (t169 * t28 - t170 * t30 + t223) * mrSges(6,2) + (t226 + t225 * qJD(4) / 0.2e1) * qJD(4);
t59 = t75 - t197;
t53 = t127 * t106;
t52 = t77 * t106;
t51 = t134 * t106;
t31 = t50 + t159;
t15 = mrSges(6,1) * t54 - mrSges(6,3) * t55;
t1 = [t218 * t109 + t219 * t107 + t222 * qJD(4) + m(5) * (qJD(4) * t125 + t107 * t8 + t109 * t9) + m(6) * (qJD(4) * t126 + t107 * t5 - t109 * t6) + t167 * qJDD(1) + (-m(5) - m(6) - t167) * g(3); m(4) * (t108 * t61 + t110 * t60) * pkin(2) + m(6) * (t117 * t184 + t2 * t59 + t23 * t31) + t111 + t52 * t159 - t31 * t51 + t59 * t15 + Ifges(3,3) * qJDD(2) + (mrSges(4,1) * t197 - mrSges(4,2) * t198) * t104 + t206 * qJD(3) + (-m(4) * t91 - m(5) * (t91 + t190) - m(6) * (t91 + t143) - mrSges(3,1) * t101 + mrSges(3,2) * t100 + t213) * g(2) + (m(4) * t199 - m(5) * (t147 - t203) - m(6) * t147 + mrSges(3,1) * t100 + mrSges(3,2) * t101 + t214) * g(1) + t221 * (-pkin(3) - t197) + t222 * t110 * t184 + t220 * (pkin(7) + t198); t111 + t75 * t15 - t50 * t51 + (-t52 + t51) * t161 + (t2 * t75 + t23 * t50 - (t108 * t23 + t117) * t185) * m(6) - t206 * qJD(2) + (-m(5) * t190 - m(6) * t143 + t213) * g(2) + (-m(6) * t86 - m(5) * (t86 - t203) + t214) * g(1) - t222 * t160 - t221 * pkin(3) + t220 * pkin(7); t225 * t146 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + t229 * t55 + t228 * t54 + t44 * t174 / 0.2e1 - t28 * t154 + t30 * t163 - (Ifges(6,1) * t173 + t43 + t78) * t174 / 0.2e1 - (-Ifges(5,2) * t174 + t227 + t79) * t173 / 0.2e1 + qJD(5) * t67 + t53 * t51 + qJ(5) * t36 - pkin(4) * t35 + t5 * mrSges(6,3) - t6 * mrSges(6,1) - t8 * mrSges(5,2) + t9 * mrSges(5,1) + (-t6 * pkin(4) - g(3) * t128 + t5 * qJ(5) + t30 * qJD(5) - t23 * t53) * m(6) + t231 * g(3) + (-m(6) * t30 + t162 - t216) * t41 + (-m(6) * t28 + t155 + t217) * t42 + (m(6) * t127 + t133 + t135) * t215 + ((-t119 / 0.2e1 + t118 / 0.2e1) * t106 - t226) * t106; -t51 * t174 - qJD(4) * t67 + (g(3) * t109 - t30 * qJD(4) - t107 * t215 + t23 * t174 + t6) * m(6) + t35;];
tau = t1;
