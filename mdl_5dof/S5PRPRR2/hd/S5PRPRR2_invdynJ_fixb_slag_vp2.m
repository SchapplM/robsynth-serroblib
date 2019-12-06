% Calculate vector of inverse dynamics joint torques for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:47
% EndTime: 2019-12-05 15:44:53
% DurationCPUTime: 2.22s
% Computational Cost: add. (1909->237), mult. (3534->321), div. (0->0), fcn. (2447->14), ass. (0->121)
t108 = cos(qJ(5));
t160 = t108 * mrSges(6,1);
t185 = m(6) * pkin(4);
t100 = qJ(2) + pkin(9);
t96 = qJ(4) + t100;
t89 = sin(t96);
t90 = cos(t96);
t114 = mrSges(5,2) * t90 + (mrSges(5,1) + t160 + t185) * t89;
t105 = sin(qJ(5));
t170 = mrSges(6,2) * t105;
t200 = -t89 * t170 + t90 * (-m(6) * pkin(7) - mrSges(6,3));
t184 = t105 / 0.2e1;
t107 = sin(qJ(2));
t110 = cos(qJ(2));
t133 = mrSges(3,1) * t107 + mrSges(3,2) * t110;
t180 = pkin(2) * t107;
t93 = sin(t100);
t94 = cos(t100);
t199 = m(4) * t180 + mrSges(4,1) * t93 + mrSges(4,2) * t94 + t133 + (-m(6) - m(5)) * (-pkin(3) * t93 - t180) + t114;
t149 = qJD(5) * t105;
t98 = qJDD(2) + qJDD(4);
t99 = qJD(2) + qJD(4);
t60 = t108 * t98 - t149 * t99;
t37 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t60;
t148 = qJD(5) * t108;
t61 = t105 * t98 + t148 * t99;
t38 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t61;
t126 = -t105 * t38 + t108 * t37;
t162 = t105 * t99;
t70 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t162;
t157 = t108 * t99;
t71 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t157;
t198 = -t70 * t148 - t71 * t149 + t126;
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t101 = sin(pkin(9));
t103 = cos(pkin(9));
t150 = qJD(1) * t107;
t81 = qJD(2) * pkin(2) + qJD(1) * t110;
t39 = -t101 * t150 + t103 * t81;
t36 = qJD(2) * pkin(3) + t39;
t40 = t101 * t81 + t103 * t150;
t17 = -t106 * t40 + t109 * t36;
t197 = t160 - t170;
t131 = mrSges(6,1) * t105 + mrSges(6,2) * t108;
t15 = -pkin(4) * t99 - t17;
t195 = t15 * t131 + qJD(5) * (Ifges(6,5) * t108 - Ifges(6,6) * t105) / 0.2e1;
t193 = m(5) + m(4);
t55 = t197 * t99;
t192 = t99 * mrSges(5,1) + t55;
t191 = (-mrSges(5,1) - t197) * t90 + (mrSges(5,2) - mrSges(6,3)) * t89;
t125 = -t105 * t70 + t108 * t71;
t118 = t99 * mrSges(5,2) - t125;
t18 = t106 * t36 + t109 * t40;
t16 = pkin(7) * t99 + t18;
t13 = qJD(3) * t108 - t105 * t16;
t14 = qJD(3) * t105 + t108 * t16;
t155 = qJDD(2) * pkin(2);
t147 = qJD(1) * qJD(2);
t140 = t107 * t147;
t72 = t110 * qJDD(1) - t140;
t63 = t72 + t155;
t139 = t110 * t147;
t73 = qJDD(1) * t107 + t139;
t27 = -t101 * t73 + t103 * t63;
t25 = qJDD(2) * pkin(3) + t27;
t28 = t101 * t63 + t103 * t73;
t8 = t17 * qJD(4) + t106 * t25 + t109 * t28;
t5 = pkin(7) * t98 + t8;
t3 = -qJD(5) * t14 + qJDD(3) * t108 - t105 * t5;
t176 = t3 * t105;
t188 = -t13 * t148 - t14 * t149 - t176;
t9 = -qJD(4) * t18 - t106 * t28 + t109 * t25;
t187 = m(5) * t17 - m(6) * t15 + t192;
t127 = -t105 * t13 + t108 * t14;
t186 = m(5) * t18 + m(6) * t127 - t118;
t181 = pkin(2) * t101;
t2 = qJD(5) * t13 + qJDD(3) * t105 + t108 * t5;
t177 = t108 * t2;
t97 = t110 * pkin(2);
t175 = t98 * mrSges(5,2);
t172 = t90 * pkin(4) + t89 * pkin(7);
t91 = pkin(2) * t103 + pkin(3);
t53 = t106 * t91 + t109 * t181;
t171 = pkin(3) * t94 + t97;
t169 = Ifges(6,4) * t105;
t168 = Ifges(6,4) * t108;
t167 = Ifges(6,2) * t108;
t102 = sin(pkin(8));
t154 = t102 * t105;
t153 = t102 * t108;
t104 = cos(pkin(8));
t152 = t104 * t105;
t151 = t104 * t108;
t137 = t200 * t102;
t136 = t200 * t104;
t135 = -g(1) * t102 + g(2) * t104;
t134 = mrSges(3,1) * t110 - mrSges(3,2) * t107;
t130 = t167 + t169;
t128 = t105 * t14 + t108 * t13;
t66 = -t101 * t107 + t103 * t110;
t67 = t101 * t110 + t103 * t107;
t29 = t106 * t67 - t109 * t66;
t30 = t106 * t66 + t109 * t67;
t52 = -t106 * t181 + t109 * t91;
t120 = t105 * (Ifges(6,1) * t108 - t169);
t115 = -qJD(5) * t128 - t176;
t113 = t115 + t177;
t45 = Ifges(6,6) * qJD(5) + t130 * t99;
t82 = Ifges(6,4) * t157;
t46 = Ifges(6,1) * t162 + Ifges(6,5) * qJD(5) + t82;
t6 = -pkin(4) * t98 - t9;
t112 = t9 * mrSges(5,1) - t8 * mrSges(5,2) + mrSges(6,3) * t177 - t6 * t197 + (Ifges(6,1) * t61 + Ifges(6,4) * t60) * t184 + t108 * (Ifges(6,4) * t61 + Ifges(6,2) * t60) / 0.2e1 + t60 * t130 / 0.2e1 + t61 * (Ifges(6,1) * t105 + t168) / 0.2e1 - t45 * t149 / 0.2e1 + Ifges(5,3) * t98 + (t46 + t99 * (-Ifges(6,2) * t105 + t168)) * t148 / 0.2e1 + (0.2e1 * Ifges(6,5) * t184 + Ifges(6,6) * t108) * qJDD(5) + (t120 * t99 / 0.2e1 + t195) * qJD(5);
t59 = t66 * qJD(2);
t58 = t66 * qJD(1);
t57 = t67 * qJD(2);
t56 = t67 * qJD(1);
t26 = -mrSges(6,1) * t60 + mrSges(6,2) * t61;
t11 = qJD(4) * t30 + t106 * t59 + t109 * t57;
t10 = -qJD(4) * t29 - t106 * t57 + t109 * t59;
t1 = [m(2) * qJDD(1) - t11 * t55 + t29 * t26 + (-t11 * t99 - t29 * t98) * mrSges(5,1) - t118 * t10 + (mrSges(4,1) * t66 - mrSges(4,2) * t67 + t134) * qJDD(2) + (-t175 + (-t105 * t71 - t108 * t70) * qJD(5) + t126) * t30 + (-m(2) - m(3) - m(6) - t193) * g(3) + m(3) * (t107 * t73 + t110 * t72) + m(6) * (t10 * t127 + t11 * t15 + t113 * t30 + t29 * t6) + m(4) * (t27 * t66 + t28 * t67 - t39 * t57 + t40 * t59) + m(5) * (t10 * t18 - t11 * t17 - t29 * t9 + t30 * t8) + (-mrSges(4,1) * t57 - mrSges(4,2) * t59 - qJD(2) * t133) * qJD(2); -t53 * t175 + t52 * t98 * mrSges(5,1) + t112 + m(5) * (t52 * t9 + t53 * t8) + (m(6) * t6 + t26) * (-pkin(4) - t52) + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t139 - t73) * mrSges(3,2) + (t140 + t72) * mrSges(3,1) + ((t101 * t28 + t103 * t27) * pkin(2) + t39 * t56 - t40 * t58) * m(4) + t187 * (t106 * t58 + t109 * t56) - t186 * (-t106 * t56 + t109 * t58) + t188 * mrSges(6,3) + (qJD(2) * t58 - t101 * t155 - t28) * mrSges(4,2) + (qJD(2) * t56 + t103 * t155 + t27) * mrSges(4,1) + (m(6) * t113 + t198) * (pkin(7) + t53) + (-t134 - m(4) * t97 - mrSges(4,1) * t94 + mrSges(4,2) * t93 - m(5) * t171 - m(6) * (t171 + t172) + t191) * g(3) + (t186 * t52 - t187 * t53) * qJD(4) + (t199 * t102 + t137) * g(2) + (t199 * t104 + t136) * g(1); t105 * t37 + t108 * t38 + t125 * qJD(5) + (qJD(5) * t127 + t105 * t2 + t108 * t3 + t135) * m(6) + t193 * (qJDD(3) + t135); (-m(6) * t172 + t191) * g(3) + t115 * mrSges(6,3) + (t104 * t114 + t136) * g(1) + (t102 * t114 + t137) * g(2) + t192 * t18 + t118 * t17 + (m(6) * (t177 + t188) + t198) * pkin(7) - pkin(4) * t26 + t112 - t6 * t185 - m(6) * (t127 * t17 + t15 * t18); Ifges(6,5) * t61 + Ifges(6,6) * t60 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t13 * t71 + t14 * t70 - g(1) * ((-t152 * t90 + t153) * mrSges(6,1) + (-t151 * t90 - t154) * mrSges(6,2)) - g(2) * ((-t154 * t90 - t151) * mrSges(6,1) + (-t153 * t90 + t152) * mrSges(6,2)) + g(3) * t131 * t89 + (t45 * t184 + (-t120 / 0.2e1 + t167 * t184) * t99 + t128 * mrSges(6,3) - (t46 + t82) * t108 / 0.2e1 - t195) * t99;];
tau = t1;
