% Calculate vector of inverse dynamics joint torques for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:07
% EndTime: 2020-01-03 11:45:13
% DurationCPUTime: 2.67s
% Computational Cost: add. (1896->309), mult. (3271->377), div. (0->0), fcn. (1667->12), ass. (0->138)
t225 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t140 = cos(qJ(4));
t223 = -mrSges(6,1) - mrSges(5,1);
t224 = t140 * t223 - mrSges(4,1);
t222 = Ifges(5,1) + Ifges(6,1);
t220 = Ifges(6,5) + Ifges(5,5);
t219 = Ifges(5,2) + Ifges(6,2);
t218 = Ifges(6,6) + Ifges(5,6);
t131 = qJDD(1) + qJDD(3);
t138 = sin(qJ(3));
t141 = cos(qJ(3));
t134 = sin(pkin(8));
t201 = pkin(1) * t134;
t135 = cos(pkin(8));
t117 = pkin(1) * t135 + pkin(2);
t92 = t117 * qJD(1);
t215 = qJD(3) * t92 + qJDD(1) * t201;
t170 = qJD(1) * t201;
t216 = -qJD(3) * t170 + t117 * qJDD(1);
t17 = t216 * t138 + t215 * t141;
t13 = pkin(7) * t131 + t17;
t217 = qJD(2) * qJD(4) + t13;
t214 = (Ifges(5,4) + Ifges(6,4)) * t140;
t213 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t211 = m(3) * pkin(1);
t196 = mrSges(5,2) + mrSges(6,2);
t137 = sin(qJ(4));
t175 = qJD(2) * t137;
t132 = qJD(1) + qJD(3);
t47 = t138 * t92 + t141 * t170;
t37 = pkin(7) * t132 + t47;
t31 = t140 * t37 + t175;
t179 = t31 * qJD(4);
t124 = t140 * qJDD(2);
t6 = -t13 * t137 + t124 - t179;
t210 = -t6 - t179;
t65 = t117 * t141 - t138 * t201;
t127 = t140 * qJD(2);
t162 = qJ(5) * t132 + t37;
t24 = -t137 * t162 + t127;
t19 = qJD(4) * pkin(4) + t24;
t30 = -t137 * t37 + t127;
t209 = -t30 * mrSges(5,3) - t19 * mrSges(6,3);
t207 = m(5) * pkin(7);
t204 = m(3) + m(4);
t200 = pkin(4) * t140;
t174 = qJD(4) * t137;
t5 = t137 * qJDD(2) + t217 * t140 - t174 * t37;
t199 = t140 * t5;
t136 = -qJ(5) - pkin(7);
t133 = qJ(1) + pkin(8);
t125 = qJ(3) + t133;
t115 = sin(t125);
t116 = cos(t125);
t119 = pkin(3) + t200;
t195 = t115 * t119 + t116 * t136;
t181 = t132 * t137;
t82 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t181;
t83 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t181;
t194 = t82 + t83;
t180 = t132 * t140;
t84 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t180;
t85 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t180;
t193 = t84 + t85;
t192 = Ifges(5,4) * t137;
t190 = Ifges(6,4) * t137;
t66 = t138 * t117 + t141 * t201;
t63 = pkin(7) + t66;
t187 = -qJ(5) - t63;
t185 = qJD(4) * t63;
t178 = t116 * pkin(3) + t115 * pkin(7);
t121 = sin(t133);
t139 = sin(qJ(1));
t177 = t139 * pkin(1) + pkin(2) * t121;
t122 = cos(t133);
t142 = cos(qJ(1));
t176 = t142 * pkin(1) + pkin(2) * t122;
t173 = qJD(4) * t140;
t126 = t140 * qJD(5);
t169 = pkin(4) * t174;
t168 = -mrSges(2,1) - t211;
t166 = t132 * t174;
t71 = t131 * t140 - t166;
t72 = t131 * t137 + t132 * t173;
t28 = -t71 * mrSges(6,1) + t72 * mrSges(6,2);
t163 = qJD(4) * t136;
t46 = -t138 * t170 + t141 * t92;
t161 = -t115 * t136 + t116 * t119;
t160 = qJD(4) * t187;
t62 = -pkin(3) - t65;
t25 = t140 * t162 + t175;
t154 = t31 * mrSges(5,3) + t25 * mrSges(6,3);
t153 = g(2) * t116 + g(3) * t115;
t95 = -mrSges(5,1) * t140 + mrSges(5,2) * t137;
t152 = -mrSges(6,1) * t140 + mrSges(6,2) * t137;
t151 = Ifges(5,2) * t140 + t192;
t150 = Ifges(6,2) * t140 + t190;
t149 = -t137 * t19 + t140 * t25;
t148 = -t137 * t30 + t140 * t31;
t18 = -t215 * t138 + t216 * t141;
t147 = t213 * t115 + t224 * t116;
t14 = -pkin(3) * t131 - t18;
t57 = t66 * qJD(3);
t146 = (t207 - t213) * t116 + t224 * t115;
t145 = t196 * t153;
t27 = -t119 * t132 + qJD(5) - t46;
t3 = qJ(5) * t71 + t126 * t132 + t5;
t36 = -pkin(3) * t132 - t46;
t58 = Ifges(6,6) * qJD(4) + t132 * t150;
t59 = Ifges(5,6) * qJD(4) + t132 * t151;
t101 = Ifges(6,4) * t180;
t60 = Ifges(6,1) * t181 + Ifges(6,5) * qJD(4) + t101;
t102 = Ifges(5,4) * t180;
t61 = Ifges(5,1) * t181 + Ifges(5,5) * qJD(4) + t102;
t8 = -pkin(4) * t71 + qJDD(5) + t14;
t144 = t3 * t140 * mrSges(6,3) + t18 * mrSges(4,1) + mrSges(5,3) * t199 + Ifges(4,3) * t131 + t14 * t95 + t8 * t152 + (-t218 * t137 + t220 * t140) * qJD(4) ^ 2 / 0.2e1 - (t59 + t58) * t174 / 0.2e1 + (t222 * t140 - t190 - t192) * t166 / 0.2e1 + (t27 * (mrSges(6,1) * t137 + mrSges(6,2) * t140) + t36 * (mrSges(5,1) * t137 + mrSges(5,2) * t140)) * qJD(4) + (t151 / 0.2e1 + t150 / 0.2e1 + t137 * t225 + t219 * t140 / 0.2e1) * t71 + (t222 * t137 + t214 / 0.2e1 + t140 * t225) * t72 + (t61 + t60 + (-t219 * t137 + t214) * t132) * t173 / 0.2e1 + (t220 * t137 + t218 * t140) * qJDD(4);
t128 = t140 * qJ(5);
t109 = t115 * pkin(3);
t96 = pkin(7) * t140 + t128;
t94 = t136 * t137;
t70 = t95 * t132;
t69 = t152 * t132;
t68 = -qJD(5) * t137 + t140 * t163;
t67 = t137 * t163 + t126;
t56 = t65 * qJD(3);
t51 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t72;
t50 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t72;
t49 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t71;
t48 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t71;
t45 = t62 - t200;
t40 = t57 + t169;
t39 = t140 * t63 + t128;
t38 = t187 * t137;
t29 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t11 = (-qJD(5) - t56) * t137 + t140 * t160;
t10 = t137 * t160 + t140 * t56 + t126;
t2 = -t37 * t173 + qJDD(4) * pkin(4) - qJ(5) * t72 + t124 + (-qJD(5) * t132 - t217) * t137;
t1 = [(t63 * t49 + m(5) * (-t185 * t30 + t31 * t56 + t5 * t63) + t56 * t85 - t83 * t185 + t209 * qJD(4)) * t140 + (-t85 * t185 + (-t25 * qJD(4) - t2) * mrSges(6,3) + t210 * mrSges(5,3) + t145 + (-m(5) * t30 - t83) * t56 + (m(5) * t210 - t51) * t63) * t137 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t135 * mrSges(3,1) - 0.2e1 * t134 * mrSges(3,2) + (t134 ^ 2 + t135 ^ 2) * t211) * pkin(1)) * qJDD(1) + (-t66 * t131 - t56 * t132 - t17) * mrSges(4,2) + (-m(6) * (t161 + t176) - m(5) * (t176 + t178) + mrSges(2,2) * t139 - mrSges(3,1) * t122 + mrSges(3,2) * t121 - m(4) * t176 + t168 * t142 + t147) * g(2) + (-m(6) * (t177 + t195) - m(5) * (t109 + t177) - mrSges(2,2) * t142 - mrSges(3,1) * t121 - mrSges(3,2) * t122 - m(4) * t177 + t168 * t139 + t146) * g(3) + m(5) * (t14 * t62 + t36 * t57) + (t65 * t131 - t57 * t132) * mrSges(4,1) + t144 + t11 * t82 + t10 * t84 + t40 * t69 + t57 * t70 + t38 * t50 + t62 * t29 + t45 * t28 + t39 * t48 + m(4) * (t17 * t66 + t18 * t65 - t46 * t57 + t47 * t56) + m(6) * (t10 * t25 + t11 * t19 + t2 * t38 + t27 * t40 + t3 * t39 + t45 * t8); (t50 + t51) * t140 + (t48 + t49) * t137 + t204 * qJDD(2) + (-t194 * t137 + t193 * t140) * qJD(4) + m(5) * (qJD(4) * t148 + t137 * t5 + t140 * t6) + m(6) * (qJD(4) * t149 + t137 * t3 + t140 * t2) + (-m(5) - m(6) - t204) * g(1); (-t6 * mrSges(5,3) - t2 * mrSges(6,3) - pkin(7) * t51 + t194 * t46 + (pkin(4) * t69 - pkin(7) * t85 - t154) * qJD(4) + t145) * t137 + (pkin(7) * t49 - t193 * t46 + (-pkin(7) * t83 + t209) * qJD(4)) * t140 + t146 * g(3) + (t132 * t46 - t17) * mrSges(4,2) + t147 * g(2) + t144 - t119 * t28 + (t132 * mrSges(4,1) - t69 - t70) * t47 + t68 * t82 + t67 * t84 + t94 * t50 + t96 * t48 - pkin(3) * t29 + (-t6 * t137 + t199 + (-t137 * t31 - t140 * t30) * qJD(4)) * t207 + (-t161 * g(2) - t195 * g(3) - t119 * t8 - t149 * t46 + t19 * t68 + t2 * t94 + t25 * t67 + t3 * t96 + (t169 - t47) * t27) * m(6) + (-pkin(3) * t14 - g(2) * t178 - g(3) * t109 - t148 * t46 - t36 * t47) * m(5); t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) - t24 * t84 - t30 * t85 + t31 * t83 + t220 * t72 + t218 * t71 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t152 + t95) * g(1) + (t50 + (-g(1) * t140 + t2) * m(6)) * pkin(4) + (-m(6) * (-t19 + t24) + t82) * t25 + ((-t102 / 0.2e1 - t101 / 0.2e1 - t60 / 0.2e1 - t61 / 0.2e1 - t36 * mrSges(5,2) - t27 * mrSges(6,2) + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(4) - t209) * t140 + (t58 / 0.2e1 + t59 / 0.2e1 - t36 * mrSges(5,1) - t27 * mrSges(6,1) + t225 * t181 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t27 - t69) * pkin(4) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1) * t180 + t154) * t137) * t132 + (g(2) * t115 - g(3) * t116) * (t137 * (m(6) * pkin(4) - t223) + t140 * t196); (t137 * t82 - t140 * t84) * t132 + (-t132 * t149 + t153 + t8) * m(6) + t28;];
tau = t1;
