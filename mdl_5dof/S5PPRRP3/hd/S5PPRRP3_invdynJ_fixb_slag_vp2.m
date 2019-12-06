% Calculate vector of inverse dynamics joint torques for
% S5PPRRP3
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:10:41
% DurationCPUTime: 3.20s
% Computational Cost: add. (1032->266), mult. (2472->359), div. (0->0), fcn. (1646->8), ass. (0->121)
t195 = mrSges(5,1) + mrSges(6,1);
t181 = Ifges(6,4) + Ifges(5,5);
t180 = Ifges(6,6) - Ifges(5,6);
t182 = -m(5) - m(6);
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t100 = mrSges(5,1) * t80 - mrSges(5,2) * t78;
t99 = t80 * mrSges(6,1) + t78 * mrSges(6,3);
t192 = -mrSges(4,1) - t100 - t99;
t94 = pkin(4) * t80 + qJ(5) * t78;
t194 = -m(6) * t94 + t182 * pkin(3) + t192;
t76 = cos(pkin(8));
t135 = qJD(1) * t76;
t74 = sin(pkin(8));
t136 = qJD(1) * t74;
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t45 = qJD(2) * t79 + t136 * t81;
t36 = qJD(3) * pkin(6) + t45;
t19 = -t135 * t80 - t36 * t78;
t191 = qJD(5) - t19;
t124 = qJDD(1) * t76;
t134 = qJD(3) * t74;
t187 = -qJD(1) * t134 + qJDD(2);
t188 = qJD(2) * qJD(3) + qJDD(1) * t74;
t17 = t187 * t79 + t188 * t81;
t8 = qJDD(3) * pkin(6) + t17;
t104 = -t124 * t78 + t80 * t8;
t3 = qJD(4) * t19 + t104;
t91 = t78 * t135 - t36 * t80;
t4 = qJD(4) * t91 - t124 * t80 - t78 * t8;
t102 = t3 * t80 - t4 * t78;
t131 = qJD(4) * t80;
t132 = qJD(4) * t78;
t190 = -t19 * t131 + t132 * t91 + t102;
t1 = qJDD(4) * qJ(5) + (qJD(5) + t19) * qJD(4) + t104;
t2 = -qJDD(4) * pkin(4) + qJDD(5) - t4;
t103 = t1 * t80 + t2 * t78;
t10 = qJD(4) * qJ(5) - t91;
t130 = t10 * qJD(4);
t7 = -qJD(4) * pkin(4) + t191;
t189 = -t78 * t130 + t7 * t131 + t103;
t44 = qJD(2) * t81 - t79 * t136;
t51 = -pkin(3) - t94;
t21 = qJD(3) * t51 - t44;
t35 = -qJD(3) * pkin(3) - t44;
t186 = t35 * (mrSges(5,1) * t78 + mrSges(5,2) * t80) + t21 * (mrSges(6,1) * t78 - mrSges(6,3) * t80) + (t180 * t78 + t181 * t80) * qJD(4) / 0.2e1;
t126 = t80 * qJD(3);
t56 = mrSges(6,2) * t126 + qJD(4) * mrSges(6,3);
t138 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t126 + t56;
t127 = t78 * qJD(3);
t139 = (-mrSges(6,2) - mrSges(5,3)) * t127 + t195 * qJD(4);
t185 = t138 * t80 - t139 * t78;
t122 = qJD(3) * qJD(4);
t50 = qJDD(3) * t78 + t122 * t80;
t28 = -qJDD(4) * mrSges(6,1) + t50 * mrSges(6,2);
t140 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t50 + t28;
t49 = -t80 * qJDD(3) + t122 * t78;
t29 = -mrSges(6,2) * t49 + qJDD(4) * mrSges(6,3);
t141 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t49 + t29;
t184 = t140 * t78 + t141 * t80;
t149 = t74 * t81;
t77 = cos(pkin(7));
t146 = t77 * t79;
t75 = sin(pkin(7));
t147 = t75 * t81;
t32 = t147 * t76 - t146;
t145 = t77 * t81;
t148 = t75 * t79;
t34 = t145 * t76 + t148;
t183 = g(1) * t34 + g(2) * t32 + t149 * g(3);
t70 = Ifges(5,4) * t126;
t144 = t78 * Ifges(6,1);
t155 = Ifges(6,5) * t80;
t98 = t144 - t155;
t179 = Ifges(5,1) * t127 + qJD(3) * t98 + t181 * qJD(4) + t70;
t46 = t99 * qJD(3);
t178 = t100 * qJD(3) + t46;
t174 = -t44 * qJD(3) + t17;
t18 = t187 * t81 - t188 * t79;
t173 = t45 * qJD(3) + t18;
t166 = t78 / 0.2e1;
t165 = -t80 / 0.2e1;
t158 = Ifges(5,4) * t78;
t157 = Ifges(5,4) * t80;
t156 = Ifges(6,5) * t78;
t151 = t74 * t78;
t150 = t74 * t80;
t143 = t80 * Ifges(5,2);
t133 = qJD(3) * t81;
t120 = t80 * t149;
t116 = t79 * t134;
t112 = m(3) + m(4) - t182;
t106 = m(6) * pkin(4) + t195;
t105 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t101 = t10 * t78 - t7 * t80;
t97 = t143 + t158;
t93 = pkin(4) * t78 - qJ(5) * t80;
t92 = t19 * t80 - t78 * t91;
t37 = t149 * t78 + t76 * t80;
t87 = t78 * (Ifges(5,1) * t80 - t158);
t86 = t80 * (Ifges(6,3) * t78 + t155);
t22 = mrSges(6,1) * t49 - mrSges(6,3) * t50;
t23 = mrSges(5,1) * t49 + mrSges(5,2) * t50;
t82 = qJD(3) ^ 2;
t85 = qJDD(3) * mrSges(4,1) - mrSges(4,2) * t82 - t22 - t23;
t9 = -qJDD(3) * pkin(3) - t18;
t84 = -mrSges(4,1) * t82 - qJDD(3) * mrSges(4,2) - qJD(3) * t178;
t71 = t76 ^ 2 * qJDD(1);
t69 = Ifges(6,5) * t127;
t48 = t93 * qJD(3);
t40 = Ifges(5,6) * qJD(4) + qJD(3) * t97;
t39 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t126 + t69;
t38 = -t76 * t78 + t120;
t30 = qJD(4) * t93 - qJD(5) * t78;
t16 = -qJD(4) * t37 - t116 * t80;
t15 = -qJD(4) * t120 + (qJD(4) * t76 + t116) * t78;
t13 = -t150 * t77 + t34 * t78;
t11 = -t150 * t75 + t32 * t78;
t5 = pkin(4) * t49 - qJ(5) * t50 - qJD(5) * t127 + t9;
t6 = [m(2) * qJDD(1) + t141 * t38 + t140 * t37 + t138 * t16 + t139 * t15 + (-m(2) - t112) * g(3) + (-t79 * t85 + t81 * t84) * t74 + m(4) * (t71 + (t17 * t81 - t18 * t79 + (-t44 * t81 - t45 * t79) * qJD(3)) * t74) + m(6) * (t1 * t38 + t10 * t16 - t15 * t7 + t2 * t37 + (t133 * t21 + t5 * t79) * t74) + m(5) * (t15 * t19 - t16 * t91 + t3 * t38 - t37 * t4 + (t133 * t35 + t79 * t9) * t74) + m(3) * (qJDD(1) * t74 ^ 2 + t71); m(3) * qJDD(2) + (t185 * qJD(3) + m(5) * (-t126 * t91 - t127 * t19 - t9) + m(6) * (t10 * t126 + t127 * t7 - t5) + m(4) * t173 + t85) * t81 + ((-t138 * t78 - t139 * t80) * qJD(4) + m(5) * (qJD(3) * t35 + t190) + m(6) * (qJD(3) * t21 + t189) + m(4) * t174 + t84 + t184) * t79 + (-g(1) * t75 + g(2) * t77) * t112; ((t19 * t78 + t80 * t91) * m(5) + (-t10 * t80 - t7 * t78) * m(6) - t185) * t44 + (mrSges(4,2) * t32 + t194 * (-t148 * t76 - t145)) * g(2) + (mrSges(4,2) * t34 + t194 * (-t146 * t76 + t147)) * g(1) + ((Ifges(5,1) + Ifges(6,1)) * t50 + (-Ifges(5,4) + Ifges(6,5)) * t49 + t181 * qJDD(4)) * t166 + t178 * t45 + t179 * t131 / 0.2e1 + (-t180 * t80 + t181 * t78) * qJDD(4) / 0.2e1 + t173 * mrSges(4,1) - t174 * mrSges(4,2) + (-t40 / 0.2e1 + t39 / 0.2e1) * t132 + (-t138 * t132 - t139 * t131 + m(5) * (-qJD(4) * t92 + t102) + m(6) * (-qJD(4) * t101 + t103) + t183 * t182 + t184) * pkin(6) + t186 * qJD(4) + (t80 * (-Ifges(5,2) * t78 + t157) + t78 * (Ifges(6,1) * t80 + t156) + t87) * t122 / 0.2e1 + (t78 * Ifges(5,1) + t157 + t98) * t50 / 0.2e1 + (mrSges(4,2) * t81 + (m(5) * pkin(3) - m(6) * t51 - t192) * t79) * t74 * g(3) + (-t183 + t189) * mrSges(6,2) + (-t183 + t190) * mrSges(5,3) - t9 * t100 - t49 * t97 / 0.2e1 - t5 * t99 + (-pkin(3) * t9 - t35 * t45) * m(5) + (t5 * t51 + (-t45 + t30) * t21) * m(6) - t86 * t122 / 0.2e1 + t80 * (Ifges(5,4) * t50 - Ifges(5,2) * t49 + Ifges(5,6) * qJDD(4)) / 0.2e1 + t51 * t22 - t30 * t46 - pkin(3) * t23 + (Ifges(6,5) * t50 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t49) * t165 + Ifges(4,3) * qJDD(3) + t49 * (-t80 * Ifges(6,3) + t156) / 0.2e1; t4 * mrSges(5,1) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + t1 * mrSges(6,3) - pkin(4) * t28 + qJ(5) * t29 + qJD(5) * t56 + t48 * t46 + t181 * t50 + t180 * t49 - t139 * t91 - t138 * t19 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + (t105 * t38 + t106 * t37) * g(3) + (t105 * (t151 * t75 + t32 * t80) + t106 * t11) * g(2) + (t105 * (t151 * t77 + t34 * t80) + t106 * t13) * g(1) + (t40 * t166 + (t143 * t166 + t86 / 0.2e1 - t87 / 0.2e1) * qJD(3) + t92 * mrSges(5,3) + t101 * mrSges(6,2) - (t69 + t39) * t78 / 0.2e1 + (qJD(3) * t144 + t179 + t70) * t165 - t186) * qJD(3) + (-pkin(4) * t2 + qJ(5) * t1 + t191 * t10 - t21 * t48 + t7 * t91) * m(6); -t46 * t127 - qJD(4) * t56 + (-g(1) * t13 - g(2) * t11 - g(3) * t37 + t127 * t21 - t130 + t2) * m(6) + t28;];
tau = t6;
