% Calculate vector of inverse dynamics joint torques for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:39
% EndTime: 2019-12-05 15:06:45
% DurationCPUTime: 2.85s
% Computational Cost: add. (990->236), mult. (2281->303), div. (0->0), fcn. (1473->10), ass. (0->123)
t185 = Ifges(5,2) + Ifges(6,2);
t75 = sin(qJ(4));
t77 = cos(qJ(4));
t96 = -mrSges(5,1) * t77 + mrSges(5,2) * t75;
t203 = -mrSges(4,1) + t96;
t187 = Ifges(5,4) + Ifges(6,4);
t202 = Ifges(5,1) + Ifges(6,1);
t186 = Ifges(6,5) + Ifges(5,5);
t184 = Ifges(6,6) + Ifges(5,6);
t154 = Ifges(6,4) * t75;
t156 = Ifges(5,4) * t75;
t201 = t185 * t77 + t154 + t156;
t61 = pkin(4) * t77 + pkin(3);
t94 = -mrSges(6,1) * t77 + mrSges(6,2) * t75;
t200 = m(5) * pkin(3) + m(6) * t61 - t203 - t94;
t74 = -qJ(5) - pkin(6);
t199 = -m(5) * pkin(6) + m(6) * t74 + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t128 = qJD(1) * qJD(3);
t76 = sin(qJ(3));
t78 = cos(qJ(3));
t192 = qJDD(1) * t78 - t76 * t128;
t193 = qJDD(1) * t76 + t78 * t128;
t70 = sin(pkin(8));
t72 = cos(pkin(8));
t10 = t192 * t70 + t193 * t72;
t7 = qJDD(3) * pkin(6) + t10;
t197 = qJD(2) * qJD(4) + t7;
t131 = t77 * qJD(3);
t196 = t187 * t131;
t195 = t203 * qJD(3);
t133 = qJD(4) * t77;
t134 = qJD(4) * t75;
t88 = t70 * t78 + t76 * t72;
t32 = t88 * qJD(1);
t20 = qJD(3) * pkin(6) + t32;
t68 = t77 * qJD(2);
t14 = -t20 * t75 + t68;
t137 = qJD(2) * t75;
t15 = t20 * t77 + t137;
t3 = t75 * qJDD(2) - t134 * t20 + t197 * t77;
t67 = t77 * qJDD(2);
t4 = -t15 * qJD(4) - t7 * t75 + t67;
t99 = t3 * t77 - t4 * t75;
t194 = -t14 * t133 - t15 * t134 + t99;
t191 = -t70 * t76 + t78 * t72;
t71 = sin(pkin(7));
t73 = cos(pkin(7));
t190 = g(1) * t73 + g(2) * t71;
t189 = -mrSges(5,2) - mrSges(6,2);
t183 = t201 * qJD(3) + t184 * qJD(4);
t132 = t75 * qJD(3);
t182 = t186 * qJD(4) + t202 * t132 + t196;
t157 = mrSges(6,2) * t77;
t31 = t191 * qJD(1);
t18 = -qJD(3) * t61 + qJD(5) - t31;
t19 = -qJD(3) * pkin(3) - t31;
t95 = mrSges(5,1) * t75 + mrSges(5,2) * t77;
t181 = t19 * t95 + t18 * (mrSges(6,1) * t75 + t157);
t180 = -t184 * t75 + t186 * t77;
t179 = m(6) * pkin(4) + mrSges(6,1);
t118 = mrSges(6,3) * t132;
t46 = qJD(4) * mrSges(6,1) - t118;
t114 = mrSges(6,3) * t131;
t48 = -qJD(4) * mrSges(6,2) + t114;
t177 = t75 * t46 - t77 * t48;
t176 = t187 * t77;
t174 = -mrSges(5,1) - t179;
t42 = t94 * qJD(3);
t173 = -m(6) * t18 - t42;
t125 = qJD(3) * qJD(5);
t126 = qJD(3) * qJD(4);
t45 = qJDD(3) * t75 + t126 * t77;
t1 = -t20 * t133 + qJDD(4) * pkin(4) - qJ(5) * t45 + t67 + (-t125 - t197) * t75;
t101 = qJ(5) * qJD(3) + t20;
t13 = t101 * t77 + t137;
t44 = qJDD(3) * t77 - t126 * t75;
t2 = qJ(5) * t44 + t125 * t77 + t3;
t12 = -t101 * t75 + t68;
t9 = qJD(4) * pkin(4) + t12;
t172 = -t1 * t75 - t13 * t134 - t9 * t133 + t2 * t77;
t165 = m(3) + m(4);
t33 = t191 * qJD(3);
t152 = t33 * t75;
t151 = t33 * t77;
t149 = t71 * t75;
t148 = t71 * t77;
t147 = t73 * t75;
t146 = t73 * t77;
t25 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t44;
t26 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t44;
t142 = t25 + t26;
t27 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t45;
t28 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t45;
t141 = t27 + t28;
t120 = mrSges(5,3) * t132;
t47 = qJD(4) * mrSges(5,1) - t120;
t140 = -t46 - t47;
t119 = mrSges(5,3) * t131;
t49 = -qJD(4) * mrSges(5,2) + t119;
t139 = t48 + t49;
t135 = qJD(3) * mrSges(4,2);
t122 = pkin(4) * t134;
t113 = m(5) + m(6) + t165;
t16 = -t44 * mrSges(6,1) + t45 * mrSges(6,2);
t104 = qJD(4) * t74;
t98 = t13 * t77 - t75 * t9;
t89 = -t14 * t75 + t15 * t77;
t85 = t75 * (Ifges(5,1) * t77 - t156);
t84 = t75 * (Ifges(6,1) * t77 - t154);
t11 = t192 * t72 - t193 * t70;
t81 = t139 * t77 + t140 * t75;
t8 = -qJDD(3) * pkin(3) - t11;
t69 = pkin(8) + qJ(3);
t65 = cos(t69);
t64 = sin(t69);
t51 = t74 * t77;
t50 = t74 * t75;
t34 = t88 * qJD(3);
t30 = -qJD(5) * t75 + t104 * t77;
t29 = qJD(5) * t77 + t104 * t75;
t17 = -mrSges(5,1) * t44 + mrSges(5,2) * t45;
t5 = -t44 * pkin(4) + qJDD(5) + t8;
t6 = [-(-qJDD(3) * mrSges(4,1) + t16 + t17) * t191 + (t42 + t195) * t34 + (t81 - t135) * t33 + (-m(2) - t113) * g(3) + m(6) * (t13 * t151 - t152 * t9 + t18 * t34 - t191 * t5) + m(5) * (-t14 * t152 + t15 * t151 + t19 * t34 - t191 * t8) + m(4) * (t11 * t191 - t31 * t34 + t32 * t33) + (-qJDD(3) * mrSges(4,2) + t142 * t77 - t141 * t75 + (-t139 * t75 + t140 * t77) * qJD(4) + m(6) * t172 + m(5) * t194 + m(4) * t10) * t88 + (m(2) + m(3) * (t70 ^ 2 + t72 ^ 2)) * qJDD(1); t141 * t77 + t142 * t75 + t165 * qJDD(2) + t81 * qJD(4) + m(5) * (qJD(4) * t89 + t3 * t75 + t4 * t77) + m(6) * (qJD(4) * t98 + t1 * t77 + t2 * t75) + (-g(1) * t71 + g(2) * t73) * t113; t176 * t45 / 0.2e1 + m(6) * (t1 * t50 + t122 * t18 + t13 * t29 - t2 * t51 + t30 * t9 - t5 * t61) + (t185 * t44 + t187 * t45) * t77 / 0.2e1 + ((-t185 * t75 + t176) * t77 + t85 + t84) * t126 / 0.2e1 + (t187 * t75 + t201) * t44 / 0.2e1 + (t199 * g(3) + t190 * t200) * t64 + (-t200 * g(3) + t190 * t199) * t65 + (-t47 * t133 - t49 * t134 + t77 * t26 - t75 * t28 + m(5) * ((-t14 * t77 - t15 * t75) * qJD(4) + t99)) * pkin(6) + (t184 * t77 + t186 * t75) * qJDD(4) + t180 * qJD(4) ^ 2 / 0.2e1 + t181 * qJD(4) + t182 * t133 / 0.2e1 - t183 * t134 / 0.2e1 + (-m(5) * t89 - m(6) * t98 + t47 * t75 - t49 * t77 + t135 + t177) * t31 + (-m(5) * t8 - t17) * pkin(3) + t172 * mrSges(6,3) + t42 * t122 + t50 * t27 - t51 * t25 - t61 * t16 + t30 * t46 + t29 * t48 - t10 * mrSges(4,2) + t11 * mrSges(4,1) + t5 * t94 + t8 * t96 + t202 * t75 * t45 + t194 * mrSges(5,3) + (-m(5) * t19 + t173 - t195) * t32 + Ifges(4,3) * qJDD(3); t4 * mrSges(5,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) + t9 * t114 - t12 * t48 + t186 * t45 + t184 * t44 + (t179 * t75 + t157 + t95) * g(3) * t64 + (t120 + t47) * t15 + (t119 - t49) * t14 + (-m(6) * (t12 - t9) + t118 + t46) * t13 + t183 * t132 / 0.2e1 - t180 * t126 / 0.2e1 + t179 * t1 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t189 * (-t148 * t65 + t147) + t174 * (-t149 * t65 - t146)) * g(2) + (t189 * (-t146 * t65 - t149) + t174 * (-t147 * t65 + t148)) * g(1) - (-t132 * t185 + t182 + t196) * t131 / 0.2e1 + (-t181 + (-t84 / 0.2e1 - t85 / 0.2e1) * qJD(3)) * qJD(3) + (t132 * t173 + t27) * pkin(4); t177 * qJD(3) + (g(3) * t65 - t98 * qJD(3) - t190 * t64 + t5) * m(6) + t16;];
tau = t6;
