% Calculate vector of inverse dynamics joint torques for
% S5RPRRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:15
% DurationCPUTime: 3.42s
% Computational Cost: add. (2110->273), mult. (3168->325), div. (0->0), fcn. (1446->6), ass. (0->128)
t102 = sin(qJ(4));
t104 = cos(qJ(4));
t119 = pkin(4) * t104 + qJ(5) * t102;
t124 = mrSges(5,1) * t104 - mrSges(5,2) * t102;
t71 = t104 * mrSges(6,1) + t102 * mrSges(6,3);
t221 = m(6) * t119 + t124 + t71;
t255 = mrSges(4,1) + t221;
t249 = -mrSges(4,2) + mrSges(6,2) + mrSges(5,3);
t103 = sin(qJ(3));
t105 = cos(qJ(3));
t194 = sin(qJ(1));
t195 = cos(qJ(1));
t49 = -t103 * t194 - t105 * t195;
t50 = t103 * t195 - t105 * t194;
t254 = -t249 * t49 - t255 * t50;
t253 = -t249 * t50 + t255 * t49;
t99 = -qJD(1) + qJD(3);
t160 = t104 * t99;
t145 = mrSges(6,2) * t160;
t65 = qJD(4) * mrSges(6,3) + t145;
t177 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t160 + t65;
t166 = t102 * t99;
t146 = mrSges(6,2) * t166;
t178 = -mrSges(5,3) * t166 - t146 + (mrSges(5,1) + mrSges(6,1)) * qJD(4);
t154 = qJ(2) * qJD(1);
t106 = -pkin(1) - pkin(2);
t80 = qJD(1) * t106 + qJD(2);
t35 = t103 * t80 + t105 * t154;
t22 = pkin(7) * t99 + t35;
t170 = t102 * t22;
t19 = -qJD(4) * pkin(4) + qJD(5) + t170;
t164 = t104 * t22;
t20 = qJD(4) * qJ(5) + t164;
t252 = t178 * t102 - t177 * t104 + t99 * mrSges(4,2) - m(6) * (t102 * t19 + t104 * t20);
t248 = Ifges(6,4) + Ifges(5,5);
t247 = Ifges(6,6) - Ifges(5,6);
t122 = t102 * mrSges(6,1) - t104 * mrSges(6,3);
t123 = mrSges(5,1) * t102 + mrSges(5,2) * t104;
t34 = -t103 * t154 + t105 * t80;
t70 = -pkin(3) - t119;
t15 = t70 * t99 - t34;
t21 = -pkin(3) * t99 - t34;
t251 = t122 * t15 + t21 * t123;
t250 = g(1) * t194 - g(2) * t195;
t245 = -qJD(3) * t154 + qJDD(1) * t106 + qJDD(2);
t149 = qJD(1) * qJD(2);
t83 = qJDD(1) * qJ(2) + t149;
t244 = qJD(3) * t80 + t83;
t9 = t245 * t103 + t244 * t105;
t98 = -qJDD(1) + qJDD(3);
t7 = pkin(7) * t98 + t9;
t6 = t104 * t7;
t2 = qJDD(4) * qJ(5) + t6 + (qJD(5) - t170) * qJD(4);
t152 = qJD(4) * t104;
t5 = -t102 * t7 - t152 * t22;
t3 = -qJDD(4) * pkin(4) + qJDD(5) - t5;
t126 = t102 * t3 + t104 * t2;
t153 = qJD(4) * t102;
t243 = t19 * t152 - t20 * t153 + t126;
t223 = (mrSges(4,1) + t124) * t99;
t38 = t71 * t99;
t130 = t38 + t223;
t235 = m(5) * t21;
t241 = -m(6) * t15 + t130 - t235;
t42 = t102 * t98 + t152 * t99;
t25 = -qJDD(4) * mrSges(6,1) + t42 * mrSges(6,2);
t41 = -t104 * t98 + t153 * t99;
t26 = -mrSges(6,2) * t41 + qJDD(4) * mrSges(6,3);
t240 = (-qJDD(4) * mrSges(5,1) + mrSges(5,3) * t42 + t25) * t102 + (-qJDD(4) * mrSges(5,2) - mrSges(5,3) * t41 + t26) * t104;
t117 = -t102 * t20 + t104 * t19;
t4 = -t153 * t22 + t6;
t125 = -t102 * t5 + t104 * t4;
t239 = m(5) * t125 + m(6) * (qJD(4) * t117 + t126) - t177 * t153 - t178 * t152 + t240;
t238 = -m(6) - m(5);
t237 = m(4) * t34;
t236 = m(4) * t35;
t131 = (t102 ^ 2 + t104 ^ 2) * t22;
t233 = m(5) * t131;
t228 = mrSges(2,1) + mrSges(3,1);
t227 = mrSges(2,2) - mrSges(3,3);
t172 = Ifges(6,5) * t104;
t77 = Ifges(6,1) * t102 - t172;
t82 = Ifges(5,4) * t160;
t226 = Ifges(5,1) * t166 + qJD(4) * t248 + t77 * t99 + t82;
t225 = (t102 * t247 + t104 * t248) * qJD(4);
t66 = -t103 * qJ(2) + t105 * t106;
t219 = -t50 * pkin(3) - t49 * pkin(7);
t218 = t49 * pkin(3) - pkin(7) * t50;
t217 = g(1) * t49 + g(2) * t50;
t17 = mrSges(5,1) * t41 + mrSges(5,2) * t42;
t10 = -t244 * t103 + t245 * t105;
t8 = -pkin(3) * t98 - t10;
t215 = m(5) * t8 + t17;
t118 = pkin(4) * t102 - qJ(5) * t104;
t151 = qJD(5) * t102;
t37 = qJD(4) * t118 - t151;
t208 = -t233 + t252;
t207 = -t236 + t252;
t196 = t104 / 0.2e1;
t183 = t98 * mrSges(4,1);
t182 = t98 * mrSges(4,2);
t67 = t105 * qJ(2) + t103 * t106;
t176 = t195 * pkin(1) + t194 * qJ(2);
t175 = Ifges(5,4) * t102;
t174 = Ifges(5,4) * t104;
t173 = Ifges(6,5) * t102;
t150 = qJDD(1) * mrSges(3,1);
t147 = t195 * pkin(2) + t176;
t137 = t166 / 0.2e1;
t136 = -t160 / 0.2e1;
t129 = -pkin(1) * t194 + t195 * qJ(2);
t121 = Ifges(5,1) * t104 - t175;
t76 = Ifges(5,2) * t104 + t175;
t120 = Ifges(6,3) * t102 + t172;
t114 = -pkin(2) * t194 + t129;
t33 = qJD(2) * t103 + qJD(3) * t67;
t1 = pkin(4) * t41 - qJ(5) * t42 - t151 * t99 + t8;
t81 = Ifges(6,5) * t166;
t28 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t160 + t81;
t29 = Ifges(5,6) * qJD(4) + t76 * t99;
t108 = Ifges(4,3) * t98 + (Ifges(5,4) * t42 - Ifges(5,2) * t41) * t196 - t1 * t71 - t8 * t124 - t9 * mrSges(4,2) + t10 * mrSges(4,1) + (t173 / 0.2e1 - t76 / 0.2e1) * t41 + (Ifges(5,1) * t102 + t174 + t77) * t42 / 0.2e1 + ((Ifges(5,1) + Ifges(6,1)) * t42 + (-Ifges(5,4) + Ifges(6,5)) * t41) * t102 / 0.2e1 + (-t29 / 0.2e1 + t28 / 0.2e1) * t153 + t226 * t152 / 0.2e1 + t125 * mrSges(5,3) + t243 * mrSges(6,2) + ((-Ifges(5,2) * t102 + t174) * t160 + t225) * qJD(4) / 0.2e1 + (Ifges(5,6) * t196 + t248 * t102) * qJDD(4) + (t120 * t136 + (Ifges(6,1) * t104 + t121 + t173) * t137 + t251) * qJD(4) - (Ifges(6,5) * t42 + 0.2e1 * Ifges(6,3) * t41 + (Ifges(6,6) + t247) * qJDD(4)) * t104 / 0.2e1;
t107 = qJD(1) ^ 2;
t91 = -qJDD(1) * pkin(1) + qJDD(2);
t40 = t118 * t99;
t27 = -t66 - t70;
t18 = t33 - t37;
t16 = mrSges(6,1) * t41 - mrSges(6,3) * t42;
t11 = [-t108 - t67 * t182 + m(3) * (-pkin(1) * t91 + (t83 + t149) * qJ(2)) + t66 * t183 - t91 * mrSges(3,1) - t18 * t38 + t27 * t16 + pkin(1) * t150 + m(6) * (t1 * t27 + t15 * t18) + m(4) * (t10 * t66 + t67 * t9) + t215 * (pkin(3) - t66) + (-t223 + t235 - t237) * t33 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1) + 0.2e1 * t83 * mrSges(3,3) + (-t208 + t236) * (qJD(2) * t105 + t66 * qJD(3)) + (-m(3) * t176 - m(4) * t147 - t228 * t195 + t227 * t194 + t238 * (t147 - t218) + t253) * g(2) + (-m(3) * t129 - m(4) * t114 + t227 * t195 + t228 * t194 + t238 * (t114 - t219) + t254) * g(1) + t239 * (-pkin(7) + t67); -t150 - t107 * mrSges(3,3) + (-t107 * qJ(2) - t250 + t91) * m(3) + (t183 - t16 - t17 + m(5) * (qJD(3) * t131 - t8) - m(6) * t1 + m(4) * t10 - t207 * qJD(3)) * t105 + (-t182 - t130 * qJD(3) + (-t102 * t177 - t104 * t178) * qJD(4) + m(5) * (qJD(3) * t21 + t125) + m(6) * (qJD(3) * t15 + t243) + m(4) * (-qJD(3) * t34 + t9) + t240) * t103 + ((t237 + t241) * t103 + (-t233 + t207) * t105) * qJD(1) + t250 * (-m(4) + t238); t108 + t70 * t16 - t37 * t38 + m(6) * (t1 * t70 + t15 * t37) - t215 * pkin(3) + (t238 * t218 - t253) * g(2) + (t238 * t219 - t254) * g(1) + t241 * t35 + t208 * t34 + t239 * pkin(7); -t19 * t145 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + t248 * t42 + t247 * t41 - (Ifges(6,1) * t160 + t28 + t81) * t166 / 0.2e1 + t221 * g(3) + (-pkin(4) * t3 + qJ(5) * t2 + qJD(5) * t20 - t117 * t22 - t15 * t40) * m(6) + qJD(5) * t65 + t40 * t38 - pkin(4) * t25 + qJ(5) * t26 + t2 * mrSges(6,3) - t3 * mrSges(6,1) - t4 * mrSges(5,2) + t5 * mrSges(5,1) + t20 * t146 + t177 * t170 + t178 * t164 + (-Ifges(5,2) * t166 + t226 + t82) * t136 + t29 * t137 + (-m(6) * t118 - t122 - t123) * t217 + (-t225 / 0.2e1 + (-t102 * t121 / 0.2e1 + t120 * t196) * t99 - t251) * t99; -t38 * t166 - qJD(4) * t65 + (-g(3) * t104 - t20 * qJD(4) + t217 * t102 + t15 * t166 + t3) * m(6) + t25;];
tau = t11;
