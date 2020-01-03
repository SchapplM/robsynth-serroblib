% Calculate vector of inverse dynamics joint torques for
% S5RPPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:57:51
% DurationCPUTime: 5.43s
% Computational Cost: add. (3094->376), mult. (6831->516), div. (0->0), fcn. (4698->14), ass. (0->166)
t111 = sin(pkin(8));
t92 = pkin(1) * t111 + qJ(3);
t84 = qJD(1) * qJD(3) + qJDD(1) * t92;
t110 = sin(pkin(9));
t112 = cos(pkin(9));
t116 = sin(qJ(4));
t119 = cos(qJ(4));
t86 = t110 * t116 - t119 * t112;
t78 = t86 * qJD(1);
t224 = qJD(5) + t78;
t223 = t110 ^ 2 + t112 ^ 2;
t100 = t112 * qJD(2);
t177 = pkin(6) * qJD(1);
t91 = t92 * qJD(1);
t58 = t100 + (-t91 - t177) * t110;
t71 = t110 * qJD(2) + t112 * t91;
t59 = t112 * t177 + t71;
t30 = -t116 * t59 + t119 * t58;
t87 = t110 * t119 + t112 * t116;
t81 = t87 * qJD(4);
t51 = -qJD(1) * t81 - qJDD(1) * t86;
t115 = sin(qJ(5));
t118 = cos(qJ(5));
t80 = t86 * qJD(4);
t50 = -qJD(1) * t80 + qJDD(1) * t87;
t79 = t87 * qJD(1);
t60 = qJD(4) * t118 - t115 * t79;
t26 = qJD(5) * t60 + qJDD(4) * t115 + t118 * t50;
t203 = t26 / 0.2e1;
t61 = qJD(4) * t115 + t118 * t79;
t27 = -qJD(5) * t61 + qJDD(4) * t118 - t115 * t50;
t202 = t27 / 0.2e1;
t47 = qJDD(5) - t51;
t201 = t47 / 0.2e1;
t222 = -m(4) - m(3);
t221 = -m(6) - m(5);
t7 = -mrSges(6,1) * t27 + mrSges(6,2) * t26;
t220 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t50 + t7;
t31 = t116 * t58 + t119 * t59;
t29 = qJD(4) * pkin(7) + t31;
t113 = cos(pkin(8));
t189 = pkin(1) * t113;
t95 = pkin(3) * t112 + pkin(2);
t90 = -t95 - t189;
t77 = qJD(1) * t90 + qJD(3);
t35 = pkin(4) * t78 - pkin(7) * t79 + t77;
t12 = -t115 * t29 + t118 * t35;
t219 = t12 * mrSges(6,1);
t13 = t115 * t35 + t118 * t29;
t218 = t13 * mrSges(6,2);
t182 = t79 * mrSges(5,3);
t180 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t60 + mrSges(6,2) * t61 + t182;
t157 = qJD(5) * t118;
t129 = -t115 * t80 + t87 * t157;
t14 = mrSges(6,1) * t47 - mrSges(6,3) * t26;
t15 = -mrSges(6,2) * t47 + mrSges(6,3) * t27;
t215 = -t115 * t14 + t118 * t15;
t98 = t112 * qJDD(2);
t62 = -t110 * t84 + t98;
t63 = t110 * qJDD(2) + t112 * t84;
t214 = -t110 * t62 + t112 * t63;
t74 = qJDD(1) * t90 + qJDD(3);
t18 = -pkin(4) * t51 - pkin(7) * t50 + t74;
t55 = t98 + (-pkin(6) * qJDD(1) - t84) * t110;
t153 = qJDD(1) * t112;
t56 = pkin(6) * t153 + t63;
t10 = t30 * qJD(4) + t116 * t55 + t119 * t56;
t8 = qJDD(4) * pkin(7) + t10;
t1 = qJD(5) * t12 + t115 * t18 + t118 * t8;
t2 = -qJD(5) * t13 - t115 * t8 + t118 * t18;
t213 = t1 * t118 - t115 * t2;
t11 = -qJD(4) * t31 - t116 * t56 + t119 * t55;
t108 = pkin(9) + qJ(4);
t101 = sin(t108);
t103 = cos(t108);
t140 = mrSges(5,1) * t103 - mrSges(5,2) * t101;
t141 = -mrSges(4,1) * t112 + mrSges(4,2) * t110;
t209 = m(4) * pkin(2) + t101 * mrSges(6,3) + mrSges(3,1) + t140 - t141;
t28 = -qJD(4) * pkin(4) - t30;
t208 = -m(6) * t28 - t180;
t207 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t206 = -m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t205 = Ifges(6,1) * t203 + Ifges(6,4) * t202 + Ifges(6,5) * t201;
t200 = -t60 / 0.2e1;
t199 = -t61 / 0.2e1;
t198 = t61 / 0.2e1;
t197 = -t224 / 0.2e1;
t194 = t79 / 0.2e1;
t193 = t118 / 0.2e1;
t192 = pkin(6) + t92;
t191 = Ifges(5,4) * t79;
t190 = Ifges(6,4) * t61;
t117 = sin(qJ(1));
t188 = pkin(1) * t117;
t120 = cos(qJ(1));
t105 = t120 * pkin(1);
t183 = t30 * mrSges(5,3);
t179 = Ifges(6,4) * t115;
t178 = Ifges(6,4) * t118;
t172 = t115 * t78;
t170 = t115 * t87;
t167 = t118 * t78;
t166 = t118 * t87;
t109 = qJ(1) + pkin(8);
t102 = sin(t109);
t163 = t102 * t115;
t162 = t102 * t118;
t104 = cos(t109);
t161 = t104 * t115;
t160 = t104 * t118;
t158 = qJD(5) * t115;
t156 = m(4) - t221;
t154 = qJDD(1) * t110;
t152 = Ifges(6,5) * t26 + Ifges(6,6) * t27 + Ifges(6,3) * t47;
t151 = m(6) * pkin(7) + mrSges(6,3);
t57 = Ifges(6,4) * t60;
t23 = Ifges(6,1) * t61 + Ifges(6,5) * t224 + t57;
t148 = t23 * t193;
t96 = -pkin(2) - t189;
t146 = -t51 * mrSges(5,1) + t50 * mrSges(5,2);
t145 = -t158 / 0.2e1;
t143 = -mrSges(4,1) * t153 + mrSges(4,2) * t154;
t142 = pkin(4) * t103 + pkin(7) * t101;
t138 = -mrSges(6,1) * t118 + mrSges(6,2) * t115;
t137 = mrSges(6,1) * t115 + mrSges(6,2) * t118;
t136 = Ifges(6,1) * t118 - t179;
t135 = -Ifges(6,2) * t115 + t178;
t134 = Ifges(6,5) * t118 - Ifges(6,6) * t115;
t132 = -t110 * (-t110 * t91 + t100) + t112 * t71;
t131 = -t12 * t115 + t13 * t118;
t38 = pkin(4) * t86 - pkin(7) * t87 + t90;
t82 = t192 * t110;
t83 = t192 * t112;
t44 = -t116 * t82 + t119 * t83;
t19 = -t115 * t44 + t118 * t38;
t20 = t115 * t38 + t118 * t44;
t43 = t116 * t83 + t119 * t82;
t128 = t118 * t80 + t158 * t87;
t127 = t28 * t137;
t125 = m(6) * pkin(4) - t138;
t124 = (-t13 * t115 - t12 * t118) * qJD(5) + t213;
t114 = -pkin(6) - qJ(3);
t89 = qJDD(1) * t96 + qJDD(3);
t73 = Ifges(5,4) * t78;
t68 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t78;
t67 = t103 * t160 + t163;
t66 = -t103 * t161 + t162;
t65 = -t103 * t162 + t161;
t64 = t103 * t163 + t160;
t49 = pkin(4) * t81 + pkin(7) * t80;
t48 = pkin(4) * t79 + pkin(7) * t78;
t42 = t79 * Ifges(5,1) + Ifges(5,5) * qJD(4) - t73;
t41 = -t78 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t191;
t40 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t51;
t37 = mrSges(6,1) * t224 - mrSges(6,3) * t61;
t36 = -mrSges(6,2) * t224 + mrSges(6,3) * t60;
t32 = -qJD(3) * t86 - qJD(4) * t43;
t22 = Ifges(6,2) * t60 + Ifges(6,6) * t224 + t190;
t21 = t61 * Ifges(6,5) + t60 * Ifges(6,6) + Ifges(6,3) * t224;
t17 = t115 * t48 + t118 * t30;
t16 = -t115 * t30 + t118 * t48;
t9 = -qJDD(4) * pkin(4) - t11;
t5 = t26 * Ifges(6,4) + t27 * Ifges(6,2) + t47 * Ifges(6,6);
t4 = -qJD(5) * t20 - t115 * t32 + t118 * t49;
t3 = qJD(5) * t19 + t115 * t49 + t118 * t32;
t6 = [-t80 * t148 + t90 * t146 + t89 * t141 + t96 * t143 + t28 * (mrSges(6,1) * t129 - mrSges(6,2) * t128) + t60 * (-Ifges(6,4) * t128 - Ifges(6,2) * t129 + Ifges(6,6) * t81) / 0.2e1 + m(5) * (t10 * t44 + t31 * t32 + t74 * t90) + m(6) * (t1 * t20 + t12 * t4 + t13 * t3 + t19 * t2) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t113 * mrSges(3,1) - 0.2e1 * t111 * mrSges(3,2) + m(3) * (t111 ^ 2 + t113 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-t1 * t170 + t12 * t128 - t129 * t13 - t166 * t2) * mrSges(6,3) + (-Ifges(5,1) * t80 - Ifges(5,4) * t81) * t194 + qJD(4) * (-Ifges(5,5) * t80 - Ifges(5,6) * t81) / 0.2e1 - t78 * (-Ifges(5,4) * t80 - Ifges(5,2) * t81) / 0.2e1 + t77 * (mrSges(5,1) * t81 - mrSges(5,2) * t80) + t81 * t219 + (-Ifges(6,1) * t128 - Ifges(6,4) * t129 + Ifges(6,5) * t81) * t198 + t166 * t205 + (Ifges(4,4) * t110 + Ifges(4,2) * t112) * t153 + (Ifges(4,1) * t110 + Ifges(4,4) * t112) * t154 + (t74 * mrSges(5,2) - t11 * mrSges(5,3) + Ifges(5,1) * t50 + Ifges(5,4) * t51 + Ifges(5,5) * qJDD(4) + t134 * t201 + t135 * t202 + t136 * t203 + t137 * t9 + t145 * t23) * t87 + (-m(5) * t11 + m(6) * t9 + t220) * t43 + (-mrSges(2,1) * t120 - t67 * mrSges(6,1) + mrSges(2,2) * t117 - t66 * mrSges(6,2) + t221 * (-t102 * t114 + t104 * t95 + t105) + t222 * t105 + t206 * t102 + (-m(6) * t142 - t209) * t104) * g(2) + (mrSges(2,1) * t117 - t65 * mrSges(6,1) + mrSges(2,2) * t120 - t64 * mrSges(6,2) - t222 * t188 + t221 * (-t104 * t114 - t188) + t206 * t104 + (m(5) * t95 - m(6) * (-t142 - t95) + t209) * t102) * g(1) - t81 * t218 + m(4) * (t132 * qJD(3) + t214 * t92 + t89 * t96) - t129 * t22 / 0.2e1 + (-m(5) * t30 - t208) * (qJD(3) * t87 + qJD(4) * t44) + (t152 / 0.2e1 + Ifges(6,3) * t201 + Ifges(6,6) * t202 + Ifges(6,5) * t203 - Ifges(5,4) * t50 - Ifges(5,2) * t51 - Ifges(5,6) * qJDD(4) + t74 * mrSges(5,1) - t10 * mrSges(5,3) + t207) * t86 + t224 * (-Ifges(6,5) * t128 - Ifges(6,6) * t129 + Ifges(6,3) * t81) / 0.2e1 + (t223 * t84 + t214) * mrSges(4,3) + t80 * t183 - t5 * t170 / 0.2e1 - t80 * t42 / 0.2e1 + t81 * t21 / 0.2e1 - t81 * t41 / 0.2e1 + t32 * t68 + t44 * t40 + t3 * t36 + t4 * t37 - t31 * t81 * mrSges(5,3) + t19 * t14 + t20 * t15; m(3) * qJDD(2) + t220 * t86 + t180 * t81 - (-t115 * t37 + t118 * t36 + t68) * t80 + (t40 + (-t115 * t36 - t118 * t37) * qJD(5) + t215) * t87 + m(4) * (t110 * t63 + t112 * t62) + m(5) * (t10 * t87 - t11 * t86 - t30 * t81 - t31 * t80) + m(6) * (t124 * t87 - t131 * t80 + t28 * t81 + t86 * t9) + (-m(3) - t156) * g(3); t78 * t68 - t180 * t79 - t223 * qJD(1) ^ 2 * mrSges(4,3) + (t224 * t36 + t14) * t118 + (-t224 * t37 + t15) * t115 + t143 + t146 + (-g(1) * t102 + g(2) * t104) * t156 + (t1 * t115 + t2 * t118 + t131 * t224 - t28 * t79) * m(6) + (t30 * t79 + t31 * t78 + t74) * m(5) + (-qJD(1) * t132 + t89) * m(4); t9 * t138 + t78 * t127 + ((mrSges(5,2) - t151) * t103 + (mrSges(5,1) + t125) * t101) * (g(1) * t104 + g(2) * t102) + (-Ifges(5,2) * t79 + t42 - t73) * t78 / 0.2e1 - (-Ifges(5,1) * t78 - t191 + t21) * t79 / 0.2e1 + (Ifges(6,5) * t79 - t136 * t78) * t199 + (Ifges(6,6) * t79 - t135 * t78) * t200 + (Ifges(6,3) * t79 - t134 * t78) * t197 - t77 * (mrSges(5,1) * t79 - mrSges(5,2) * t78) - qJD(4) * (-Ifges(5,5) * t78 - Ifges(5,6) * t79) / 0.2e1 + (t145 - t172 / 0.2e1) * t22 + t79 * t218 + (Ifges(6,5) * t115 + Ifges(6,6) * t118) * t201 + (Ifges(6,2) * t118 + t179) * t202 + (Ifges(6,1) * t115 + t178) * t203 + t115 * t205 + t5 * t193 + t41 * t194 - t78 * t183 - t79 * t219 + (-pkin(4) * t9 - t12 * t16 - t13 * t17) * m(6) + ((-t158 - t172) * t13 + (-t157 - t167) * t12 + t213) * mrSges(6,3) + (m(6) * t124 - t157 * t37 - t158 * t36 + t215) * pkin(7) + (t182 + t208) * t31 + (t134 * t224 + t135 * t60 + t136 * t61) * qJD(5) / 0.2e1 + t23 * t167 / 0.2e1 + (t148 + t127) * qJD(5) + (-t101 * t151 - t103 * t125 - t140) * g(3) - t30 * t68 + Ifges(5,5) * t50 + Ifges(5,6) * t51 - t17 * t36 - t16 * t37 - pkin(4) * t7 - t10 * mrSges(5,2) + t11 * mrSges(5,1) + Ifges(5,3) * qJDD(4); -t28 * (mrSges(6,1) * t61 + mrSges(6,2) * t60) + (Ifges(6,1) * t60 - t190) * t199 + t22 * t198 + (Ifges(6,5) * t60 - Ifges(6,6) * t61) * t197 - t12 * t36 + t13 * t37 - g(1) * (mrSges(6,1) * t66 - mrSges(6,2) * t67) - g(2) * (-mrSges(6,1) * t64 + mrSges(6,2) * t65) + g(3) * t137 * t101 + (t12 * t60 + t13 * t61) * mrSges(6,3) + t152 + (-Ifges(6,2) * t61 + t23 + t57) * t200 + t207;];
tau = t6;
