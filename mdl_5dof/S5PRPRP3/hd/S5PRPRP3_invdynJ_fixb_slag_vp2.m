% Calculate vector of inverse dynamics joint torques for
% S5PRPRP3
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:40
% EndTime: 2019-12-05 15:32:52
% DurationCPUTime: 4.06s
% Computational Cost: add. (1156->272), mult. (2431->337), div. (0->0), fcn. (1482->10), ass. (0->142)
t198 = Ifges(5,2) + Ifges(6,2);
t200 = Ifges(5,4) + Ifges(6,4);
t84 = sin(qJ(4));
t86 = cos(qJ(4));
t107 = -mrSges(5,1) * t86 + mrSges(5,2) * t84;
t214 = -mrSges(4,1) + t107;
t213 = Ifges(5,1) + Ifges(6,1);
t199 = Ifges(6,5) + Ifges(5,5);
t197 = Ifges(6,6) + Ifges(5,6);
t165 = Ifges(6,4) * t84;
t167 = Ifges(5,4) * t84;
t212 = t198 * t86 + t165 + t167;
t211 = m(5) + m(6) + m(4);
t105 = -mrSges(6,1) * t86 + mrSges(6,2) * t84;
t68 = pkin(4) * t86 + pkin(3);
t210 = m(5) * pkin(3) + m(6) * t68 - t105 - t214;
t144 = qJD(2) * t86;
t208 = t200 * t144;
t140 = qJDD(2) * pkin(2);
t138 = qJD(1) * qJD(2);
t85 = sin(qJ(2));
t123 = t85 * t138;
t87 = cos(qJ(2));
t55 = t87 * qJDD(1) - t123;
t44 = t55 + t140;
t122 = t87 * t138;
t56 = qJDD(1) * t85 + t122;
t79 = sin(pkin(8));
t81 = cos(pkin(8));
t13 = t79 * t44 + t81 * t56;
t11 = qJDD(2) * pkin(6) + t13;
t207 = qJD(3) * qJD(4) + t11;
t206 = t214 * qJD(2);
t142 = qJD(4) * t84;
t149 = qJD(1) * t85;
t148 = qJD(1) * t87;
t62 = qJD(2) * pkin(2) + t148;
t24 = t149 * t81 + t79 * t62;
t20 = qJD(2) * pkin(6) + t24;
t3 = t84 * qJDD(3) - t142 * t20 + t207 * t86;
t143 = qJD(3) * t84;
t15 = t20 * t86 + t143;
t74 = t86 * qJDD(3);
t4 = -t15 * qJD(4) - t11 * t84 + t74;
t111 = t3 * t86 - t4 * t84;
t76 = t86 * qJD(3);
t14 = -t20 * t84 + t76;
t141 = qJD(4) * t86;
t205 = -t14 * t141 - t15 * t142 + t111;
t80 = sin(pkin(7));
t82 = cos(pkin(7));
t204 = g(1) * t82 + g(2) * t80;
t78 = qJ(2) + pkin(8);
t71 = sin(t78);
t203 = g(3) * t71;
t202 = -mrSges(5,2) - mrSges(6,2);
t196 = qJD(2) * t212 + qJD(4) * t197;
t145 = qJD(2) * t84;
t195 = t199 * qJD(4) + t145 * t213 + t208;
t106 = mrSges(5,1) * t84 + mrSges(5,2) * t86;
t63 = t79 * t149;
t23 = t62 * t81 - t63;
t16 = -qJD(2) * t68 + qJD(5) - t23;
t168 = mrSges(6,2) * t86;
t19 = -qJD(2) * pkin(3) - t23;
t194 = t19 * t106 + t16 * (mrSges(6,1) * t84 + t168);
t193 = m(6) * pkin(4) + mrSges(6,1);
t192 = -t197 * t84 + t199 * t86;
t125 = mrSges(6,3) * t145;
t58 = qJD(4) * mrSges(6,1) - t125;
t124 = mrSges(6,3) * t144;
t60 = -qJD(4) * mrSges(6,2) + t124;
t190 = t84 * t58 - t86 * t60;
t189 = t200 * t86;
t188 = -mrSges(5,1) - t193;
t50 = t105 * qJD(2);
t187 = -m(6) * t16 - t50;
t136 = qJD(2) * qJD(5);
t137 = qJD(2) * qJD(4);
t54 = qJDD(2) * t84 + t137 * t86;
t1 = -t20 * t141 + qJDD(4) * pkin(4) - qJ(5) * t54 + t74 + (-t136 - t207) * t84;
t53 = qJDD(2) * t86 - t137 * t84;
t2 = qJ(5) * t53 + t136 * t86 + t3;
t114 = qJ(5) * qJD(2) + t20;
t8 = -t114 * t84 + t76;
t6 = qJD(4) * pkin(4) + t8;
t9 = t114 * t86 + t143;
t185 = -t1 * t84 - t6 * t141 - t9 * t142 + t2 * t86;
t184 = qJD(2) ^ 2;
t178 = pkin(2) * t79;
t177 = pkin(2) * t81;
t77 = t87 * pkin(2);
t48 = t79 * t85 - t81 * t87;
t37 = t48 * qJD(2);
t163 = t37 * t84;
t162 = t37 * t86;
t161 = t80 * t84;
t160 = t80 * t86;
t159 = t82 * t84;
t158 = t82 * t86;
t30 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t53;
t31 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t53;
t155 = t30 + t31;
t32 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t54;
t33 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t54;
t154 = t32 + t33;
t129 = mrSges(5,3) * t145;
t59 = qJD(4) * mrSges(5,1) - t129;
t153 = -t58 - t59;
t128 = mrSges(5,3) * t144;
t61 = -qJD(4) * mrSges(5,2) + t128;
t152 = t60 + t61;
t65 = pkin(6) + t178;
t150 = qJ(5) + t65;
t146 = qJD(2) * mrSges(4,2);
t139 = qJDD(2) * mrSges(4,2);
t132 = pkin(4) * t142;
t17 = -t53 * mrSges(6,1) + t54 * mrSges(6,2);
t12 = t44 * t81 - t79 * t56;
t115 = qJD(4) * t150;
t112 = -g(1) * t80 + g(2) * t82;
t110 = -t6 * t84 + t86 * t9;
t109 = mrSges(3,1) * t87 - mrSges(3,2) * t85;
t108 = mrSges(3,1) * t85 + mrSges(3,2) * t87;
t100 = -t14 * t84 + t15 * t86;
t99 = t79 * t87 + t81 * t85;
t10 = -qJDD(2) * pkin(3) - t12;
t96 = t84 * (Ifges(5,1) * t86 - t167);
t95 = t84 * (Ifges(6,1) * t86 - t165);
t92 = t152 * t86 + t153 * t84;
t83 = -qJ(5) - pkin(6);
t72 = cos(t78);
t66 = -pkin(3) - t177;
t57 = -t68 - t177;
t46 = t150 * t86;
t45 = t150 * t84;
t35 = t99 * qJD(2);
t22 = -qJD(5) * t84 - t115 * t86;
t21 = qJD(5) * t86 - t115 * t84;
t18 = -mrSges(5,1) * t53 + mrSges(5,2) * t54;
t5 = -pkin(4) * t53 + qJDD(5) + t10;
t7 = [m(2) * qJDD(1) - t108 * t184 + (t17 + t18) * t48 + (t50 + t206) * t35 + (-mrSges(4,1) * t48 + t109) * qJDD(2) - (t92 - t146) * t37 + (-m(2) - m(3) - t211) * g(3) + m(6) * (t16 * t35 - t162 * t9 + t163 * t6 + t48 * t5) + m(5) * (t10 * t48 + t14 * t163 - t15 * t162 + t19 * t35) + m(4) * (-t12 * t48 - t23 * t35 - t24 * t37) + m(3) * (t55 * t87 + t56 * t85) + (-t139 + t155 * t86 - t154 * t84 + (-t152 * t84 + t153 * t86) * qJD(4) + m(6) * t185 + m(5) * t205 + m(4) * t13) * t99; (-t203 + t205) * mrSges(5,3) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t213 * t84 * t54 + t204 * (t108 + t211 * pkin(2) * t85 + (-m(5) * pkin(6) + m(6) * t83 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3)) * t72 + t210 * t71) + (m(5) * t66 + t107) * t10 - t139 * t178 + m(4) * (t12 * t81 + t13 * t79) * pkin(2) + t189 * t54 / 0.2e1 + (t198 * t53 + t200 * t54) * t86 / 0.2e1 + (-t59 * t141 - t61 * t142 + m(5) * ((-t14 * t86 - t15 * t84) * qJD(4) + t111) + t86 * t31 - t84 * t33) * t65 + (m(4) * t23 - m(5) * t19 + t187 - t206) * t99 * qJD(1) + (t185 - t203) * mrSges(6,3) + t195 * t141 / 0.2e1 - t196 * t142 / 0.2e1 + (t197 * t86 + t199 * t84) * qJDD(4) + t194 * qJD(4) + (-m(4) * t24 - m(5) * t100 - m(6) * t110 + t59 * t84 - t61 * t86 + t146 + t190) * (t148 * t81 - t63) + t192 * qJD(4) ^ 2 / 0.2e1 + m(6) * (-t1 * t45 + t16 * t132 + t2 * t46 + t21 * t9 + t22 * t6 + t5 * t57) + (-m(6) * (-t71 * t83 + t77) - m(5) * (pkin(6) * t71 + t77) - m(4) * t77 + mrSges(4,2) * t71 - t109 - t210 * t72) * g(3) + (-t56 + t122) * mrSges(3,2) + t57 * t17 + t22 * t58 + t21 * t60 + t66 * t18 + t50 * t132 + t5 * t105 + (t55 + t123) * mrSges(3,1) + ((-t198 * t84 + t189) * t86 + t96 + t95) * t137 / 0.2e1 - t45 * t32 + t46 * t30 - t13 * mrSges(4,2) + (t200 * t84 + t212) * t53 / 0.2e1 + (t140 * t81 + t12) * mrSges(4,1); t154 * t86 + t155 * t84 + t92 * qJD(4) + (qJD(4) * t110 + t1 * t86 + t2 * t84 + t112) * m(6) + (qJD(4) * t100 + t3 * t84 + t4 * t86 + t112) * m(5) + (qJDD(3) + t112) * m(4); t4 * mrSges(5,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) + t6 * t124 - t8 * t60 + (-m(6) * (-t6 + t8) + t58 + t125) * t9 + t199 * t54 + t197 * t53 + (-t96 / 0.2e1 - t95 / 0.2e1) * t184 + (t193 * t84 + t106 + t168) * t203 + (t59 + t129) * t15 + (-t61 + t128) * t14 + t196 * t145 / 0.2e1 - t192 * t137 / 0.2e1 + t193 * t1 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) - t194 * qJD(2) + (t202 * (-t160 * t72 + t159) + t188 * (-t161 * t72 - t158)) * g(2) + (t202 * (-t158 * t72 - t161) + t188 * (-t159 * t72 + t160)) * g(1) - (-t145 * t198 + t195 + t208) * t144 / 0.2e1 + (t145 * t187 + t32) * pkin(4); t190 * qJD(2) + (g(3) * t72 - t110 * qJD(2) - t204 * t71 + t5) * m(6) + t17;];
tau = t7;
