% Calculate vector of inverse dynamics joint torques for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:13
% EndTime: 2022-01-23 09:12:22
% DurationCPUTime: 5.08s
% Computational Cost: add. (1469->292), mult. (3110->390), div. (0->0), fcn. (1845->10), ass. (0->150)
t215 = Ifges(5,4) + Ifges(6,4);
t228 = Ifges(5,1) + Ifges(6,1);
t227 = Ifges(5,2) + Ifges(6,2);
t97 = cos(qJ(4));
t231 = t215 * t97;
t95 = sin(qJ(4));
t230 = t215 * t95;
t90 = sin(pkin(8));
t229 = -t90 / 0.2e1;
t214 = Ifges(5,5) + Ifges(6,5);
t213 = Ifges(5,6) + Ifges(6,6);
t226 = Ifges(5,3) + Ifges(6,3);
t225 = (-t227 * t95 + t231) * t90;
t224 = (t228 * t97 - t230) * t90;
t87 = t90 ^ 2;
t92 = cos(pkin(8));
t221 = mrSges(4,3) * (t92 ^ 2 + t87);
t91 = sin(pkin(7));
t77 = pkin(1) * t91 + qJ(3);
t65 = qJD(1) * qJD(3) + qJDD(1) * t77;
t74 = -qJD(1) * t92 + qJD(4);
t220 = t97 * (t225 * qJD(1) + t213 * t74) + t95 * (t224 * qJD(1) + t214 * t74);
t68 = t77 * qJD(1);
t83 = t92 * qJD(2);
t44 = t68 * t90 - t83;
t181 = t44 * t90;
t155 = t95 * qJD(1);
t21 = qJD(5) - t83 + (pkin(4) * t155 + t68) * t90;
t184 = t21 * t90;
t219 = -(mrSges(6,1) * t97 - mrSges(6,2) * t95) * t184 - (mrSges(5,1) * t97 - mrSges(5,2) * t95) * t181 + (-t213 * t97 - t214 * t95) * t74 * t229;
t196 = m(6) * pkin(4);
t218 = -mrSges(5,1) - mrSges(6,1);
t217 = mrSges(5,2) + mrSges(6,2);
t216 = -mrSges(4,3) + mrSges(3,2);
t120 = pkin(3) * t92 + pkin(6) * t90;
t114 = -pkin(2) - t120;
t93 = cos(pkin(7));
t191 = pkin(1) * t93;
t63 = t114 - t191;
t172 = t92 * t97;
t66 = t77 * t172;
t20 = t95 * t63 + t66;
t211 = mrSges(6,1) + t196;
t150 = qJD(1) * qJD(4);
t60 = (qJDD(1) * t97 - t150 * t95) * t90;
t61 = (-qJDD(1) * t95 - t150 * t97) * t90;
t153 = qJDD(1) * t92;
t73 = qJDD(4) - t153;
t208 = t213 * t61 + t214 * t60 + t226 * t73;
t36 = -t92 * qJDD(2) + t65 * t90;
t37 = qJDD(2) * t90 + t65 * t92;
t207 = t36 * t90 + t37 * t92;
t152 = m(4) + m(5) + m(6);
t205 = mrSges(5,1) + t211;
t118 = mrSges(4,1) * t92 - mrSges(4,2) * t90;
t201 = mrSges(3,1) + t118 + (mrSges(5,3) + mrSges(6,3)) * t90;
t160 = qJD(1) * t90;
t142 = t97 * t160;
t125 = mrSges(6,3) * t142;
t55 = mrSges(6,1) * t74 - t125;
t127 = mrSges(5,3) * t142;
t56 = mrSges(5,1) * t74 - t127;
t163 = t55 + t56;
t143 = t90 * t155;
t126 = mrSges(6,3) * t143;
t53 = -mrSges(6,2) * t74 - t126;
t128 = mrSges(5,3) * t143;
t54 = -mrSges(5,2) * t74 - t128;
t164 = t53 + t54;
t200 = t163 * t95 - t164 * t97;
t199 = ((-t228 * t95 - t231) * t97 / 0.2e1 - (-t227 * t97 - t230) * t95 / 0.2e1) * t87;
t156 = qJD(5) * t90;
t132 = qJD(1) * t156;
t157 = qJD(4) * t97;
t35 = qJDD(1) * t63 + qJDD(3);
t38 = qJD(1) * t63 + qJD(3);
t148 = t38 * t157 + t95 * t35 + t97 * t37;
t45 = qJD(2) * t90 + t68 * t92;
t2 = qJ(5) * t61 + (-qJD(4) * t45 - t132) * t95 + t148;
t158 = qJD(4) * t95;
t3 = -t158 * t45 + t148;
t11 = t38 * t95 + t45 * t97;
t4 = -qJD(4) * t11 + t97 * t35 - t37 * t95;
t198 = -t4 * mrSges(5,1) + t3 * mrSges(5,2) + t2 * mrSges(6,2);
t197 = qJD(1) ^ 2;
t195 = t60 / 0.2e1;
t194 = t61 / 0.2e1;
t193 = t73 / 0.2e1;
t96 = sin(qJ(1));
t190 = pkin(1) * t96;
t98 = cos(qJ(1));
t86 = t98 * pkin(1);
t189 = mrSges(6,2) * t97;
t89 = qJ(1) + pkin(7);
t84 = sin(t89);
t179 = t84 * t95;
t85 = cos(t89);
t178 = t85 * t95;
t175 = t90 * t95;
t174 = t90 * t97;
t173 = t92 * t95;
t22 = mrSges(6,1) * t73 - mrSges(6,3) * t60;
t23 = mrSges(5,1) * t73 - mrSges(5,3) * t60;
t167 = t22 + t23;
t24 = -mrSges(6,2) * t73 + mrSges(6,3) * t61;
t25 = -mrSges(5,2) * t73 + mrSges(5,3) * t61;
t166 = t24 + t25;
t159 = qJD(3) * t92;
t165 = t63 * t157 + t97 * t159;
t162 = qJ(5) * t90;
t154 = qJDD(1) * t90;
t149 = t77 * t173;
t144 = t97 * t162;
t141 = t95 * t159;
t140 = t90 * t158;
t139 = t90 * t157;
t80 = -pkin(2) - t191;
t135 = qJ(5) * t160;
t17 = -t61 * mrSges(6,1) + t60 * mrSges(6,2);
t10 = t97 * t38 - t45 * t95;
t119 = -mrSges(4,1) * t153 + mrSges(4,2) * t154;
t117 = mrSges(5,1) * t95 + mrSges(5,2) * t97;
t116 = t45 * t92 + t181;
t115 = (pkin(4) * t97 + pkin(3)) * t92 - t90 * (-qJ(5) - pkin(6));
t41 = -t173 * t85 + t84 * t97;
t39 = t173 * t84 + t85 * t97;
t113 = (mrSges(6,1) * t95 + t189) * t90;
t8 = -t135 * t97 + t10;
t67 = qJDD(1) * t80 + qJDD(3);
t64 = (pkin(4) * t157 + qJD(3)) * t90;
t62 = t117 * t90;
t59 = (pkin(4) * t95 + t77) * t90;
t58 = t117 * t160;
t57 = qJD(1) * t113;
t52 = t97 * t63;
t42 = t172 * t85 + t179;
t40 = -t172 * t84 + t178;
t19 = t52 - t149;
t18 = -mrSges(5,1) * t61 + mrSges(5,2) * t60;
t16 = -t162 * t95 + t20;
t15 = -qJD(4) * t20 - t141;
t14 = -qJD(4) * t149 + t165;
t13 = -t144 + t52 + (-t77 * t95 - pkin(4)) * t92;
t12 = -pkin(4) * t61 + qJDD(5) + t36;
t9 = -t135 * t95 + t11;
t7 = -t141 - t97 * t156 + (-t66 + (-t63 + t162) * t95) * qJD(4);
t6 = -t95 * t156 + (-t144 - t149) * qJD(4) + t165;
t5 = pkin(4) * t74 + t8;
t1 = pkin(4) * t73 - qJ(5) * t60 - t132 * t97 + t4;
t26 = [(-t179 * t196 - m(3) * t86 - mrSges(2,1) * t98 + mrSges(2,2) * t96 + t216 * t84 + t218 * t42 - t217 * t41 - t152 * (t85 * pkin(2) + t84 * qJ(3) + t86) + (-m(5) * t120 - m(6) * t115 - t201) * t85) * g(2) + (-t178 * t196 + m(3) * t190 + mrSges(2,1) * t96 + mrSges(2,2) * t98 + t216 * t85 + t218 * t40 - t217 * t39 - t152 * (t85 * qJ(3) - t190) + (-m(5) * t114 - m(6) * (-pkin(2) - t115) + m(4) * pkin(2) + t201) * t84) * g(1) + (m(5) * (qJD(3) * t44 + t36 * t77) + t77 * t18 + qJD(3) * t58 + Ifges(4,4) * t153 + Ifges(4,1) * t154 + (-t213 * t95 + t214 * t97) * t193) * t90 + t207 * mrSges(4,3) + m(4) * (t116 * qJD(3) + t207 * t77 + t67 * t80) - t67 * t118 + t80 * t119 + m(5) * (t10 * t15 + t11 * t14 + t19 * t4 + t20 * t3) + t199 * t150 + (-t1 * t174 - t139 * t9 + t140 * t5 - t175 * t2) * mrSges(6,3) + (t10 * t140 - t11 * t139 - t174 * t4 - t175 * t3) * mrSges(5,3) + m(6) * (t1 * t13 + t12 * t59 + t16 * t2 + t21 * t64 + t5 * t7 + t6 * t9) + t59 * t17 + t36 * t62 + t64 * t57 + t6 * t53 + t14 * t54 + t7 * t55 + t15 * t56 + t13 * t22 + t19 * t23 + t16 * t24 + t20 * t25 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t93 * mrSges(3,1) - 0.2e1 * t91 * mrSges(3,2) + m(3) * (t91 ^ 2 + t93 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (t220 * t229 - t219) * qJD(4) + (t214 * t73 + t215 * t61 + t228 * t60) * t174 / 0.2e1 + (-t1 * mrSges(6,1) + Ifges(4,4) * t154 + Ifges(4,2) * t153 - t214 * t195 - t213 * t194 - t226 * t193 + t198 - t208 / 0.2e1) * t92 - (t213 * t73 + t215 * t60 + t227 * t61) * t175 / 0.2e1 + t225 * t194 + t224 * t195 + t65 * t221 + t12 * t113; m(3) * qJDD(2) + (-m(3) - t152) * g(3) + (-m(5) * t36 - m(6) * t12 - t17 - t18) * t92 + (t166 * t97 - t167 * t95 + (-t163 * t97 - t164 * t95) * qJD(4) + m(5) * (-t10 * t157 - t11 * t158 + t3 * t97 - t4 * t95) + m(6) * (-t1 * t95 - t157 * t5 - t158 * t9 + t2 * t97)) * t90 + m(4) * (-t36 * t92 + t37 * t90); t167 * t97 + t166 * t95 - t197 * t221 - t200 * qJD(4) + m(6) * (t1 * t97 + t2 * t95 + (-t5 * t95 + t9 * t97) * qJD(4)) + m(5) * (t3 * t95 + t4 * t97 + (-t10 * t95 + t11 * t97) * qJD(4)) + m(4) * t67 + ((-t57 - t58) * t90 + t200 * t92 - m(5) * (-t10 * t173 + t11 * t172 + t181) - m(6) * (t172 * t9 - t173 * t5 + t184) - m(4) * t116) * qJD(1) + t119 + (-g(1) * t84 + g(2) * t85) * t152; -t199 * t197 - t5 * t126 + (t127 + t56) * t11 + (-t128 - t54) * t10 + t211 * t1 + (-t205 * t41 + t217 * t42) * g(1) + (t205 * t39 - t217 * t40) * g(2) + (-m(6) * (-t5 + t8) + t125 + t55) * t9 - t8 * t53 + (t62 - (-t211 * t95 - t189) * t90) * g(3) - t198 + t208 + t220 * t160 / 0.2e1 + t219 * qJD(1) + ((-m(6) * t21 - t57) * t142 + t22) * pkin(4); (g(3) * t92 + t12) * m(6) + ((-g(1) * t85 - g(2) * t84) * m(6) + (t97 * t55 + t95 * t53 - m(6) * (-t5 * t97 - t9 * t95)) * qJD(1)) * t90 + t17;];
tau = t26;
