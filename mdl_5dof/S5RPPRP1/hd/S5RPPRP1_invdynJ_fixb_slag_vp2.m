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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:35:36
% EndTime: 2019-12-05 17:35:46
% DurationCPUTime: 4.99s
% Computational Cost: add. (1469->289), mult. (3110->388), div. (0->0), fcn. (1845->10), ass. (0->147)
t214 = Ifges(5,4) + Ifges(6,4);
t227 = Ifges(5,1) + Ifges(6,1);
t226 = Ifges(5,2) + Ifges(6,2);
t96 = cos(qJ(4));
t229 = t214 * t96;
t94 = sin(qJ(4));
t228 = t214 * t94;
t213 = Ifges(5,5) + Ifges(6,5);
t212 = Ifges(5,6) + Ifges(6,6);
t225 = Ifges(5,3) + Ifges(6,3);
t89 = sin(pkin(8));
t224 = (-t226 * t94 + t229) * t89;
t223 = (t227 * t96 - t228) * t89;
t86 = t89 ^ 2;
t91 = cos(pkin(8));
t220 = mrSges(4,3) * (t91 ^ 2 + t86);
t90 = sin(pkin(7));
t78 = pkin(1) * t90 + qJ(3);
t65 = qJD(1) * qJD(3) + qJDD(1) * t78;
t76 = -qJD(1) * t91 + qJD(4);
t219 = t96 * (t224 * qJD(1) + t212 * t76) + t94 * (t223 * qJD(1) + t213 * t76);
t68 = t78 * qJD(1);
t83 = t91 * qJD(2);
t44 = t68 * t89 - t83;
t179 = t44 * t89;
t155 = t94 * qJD(1);
t21 = qJD(5) - t83 + (pkin(4) * t155 + t68) * t89;
t182 = t21 * t89;
t218 = (mrSges(5,1) * t96 - mrSges(5,2) * t94) * t179 + (mrSges(6,1) * t96 - mrSges(6,2) * t94) * t182 + (-t212 * t96 - t213 * t94) * t76 * t89 / 0.2e1;
t195 = m(6) * pkin(4);
t217 = mrSges(5,1) + mrSges(6,1);
t216 = -mrSges(5,2) - mrSges(6,2);
t215 = mrSges(4,3) - mrSges(3,2);
t116 = -pkin(3) * t91 - pkin(6) * t89 - pkin(2);
t92 = cos(pkin(7));
t190 = pkin(1) * t92;
t63 = t116 - t190;
t172 = t91 * t96;
t66 = t78 * t172;
t20 = t94 * t63 + t66;
t210 = mrSges(6,1) + t195;
t150 = qJD(1) * qJD(4);
t60 = (qJDD(1) * t96 - t150 * t94) * t89;
t61 = (-qJDD(1) * t94 - t150 * t96) * t89;
t153 = qJDD(1) * t91;
t75 = qJDD(4) - t153;
t207 = t212 * t61 + t213 * t60 + t225 * t75;
t36 = -t91 * qJDD(2) + t65 * t89;
t37 = qJDD(2) * t89 + t65 * t91;
t206 = t36 * t89 + t37 * t91;
t152 = -m(4) - m(5) - m(6);
t204 = mrSges(5,1) + t210;
t120 = -mrSges(4,1) * t91 + mrSges(4,2) * t89;
t200 = -m(5) * t116 + m(4) * pkin(2) - t120 - m(6) * (-(pkin(4) * t96 + pkin(3)) * t91 - pkin(2)) + mrSges(3,1) + (mrSges(5,3) - m(6) * (-qJ(5) - pkin(6)) + mrSges(6,3)) * t89;
t161 = qJD(1) * t89;
t143 = t96 * t161;
t126 = mrSges(6,3) * t143;
t55 = mrSges(6,1) * t76 - t126;
t128 = mrSges(5,3) * t143;
t56 = mrSges(5,1) * t76 - t128;
t163 = t55 + t56;
t144 = t89 * t155;
t127 = mrSges(6,3) * t144;
t53 = -mrSges(6,2) * t76 - t127;
t129 = mrSges(5,3) * t144;
t54 = -mrSges(5,2) * t76 - t129;
t164 = t53 + t54;
t199 = t163 * t94 - t164 * t96;
t198 = ((-t227 * t94 - t229) * t96 / 0.2e1 - (-t226 * t96 - t228) * t94 / 0.2e1) * t86;
t157 = qJD(5) * t89;
t133 = qJD(1) * t157;
t158 = qJD(4) * t96;
t35 = qJDD(1) * t63 + qJDD(3);
t38 = qJD(1) * t63 + qJD(3);
t148 = t38 * t158 + t94 * t35 + t96 * t37;
t45 = qJD(2) * t89 + t68 * t91;
t2 = qJ(5) * t61 + (-qJD(4) * t45 - t133) * t94 + t148;
t159 = qJD(4) * t94;
t3 = -t159 * t45 + t148;
t11 = t38 * t94 + t45 * t96;
t4 = -qJD(4) * t11 + t96 * t35 - t37 * t94;
t197 = -t4 * mrSges(5,1) + t3 * mrSges(5,2) + t2 * mrSges(6,2);
t196 = qJD(1) ^ 2;
t194 = t60 / 0.2e1;
t193 = t61 / 0.2e1;
t192 = t75 / 0.2e1;
t95 = sin(qJ(1));
t189 = pkin(1) * t95;
t97 = cos(qJ(1));
t188 = pkin(1) * t97;
t187 = mrSges(6,2) * t96;
t88 = qJ(1) + pkin(7);
t84 = sin(t88);
t177 = t84 * t94;
t85 = cos(t88);
t176 = t85 * t94;
t175 = t89 * t94;
t174 = t89 * t96;
t173 = t91 * t94;
t22 = mrSges(6,1) * t75 - mrSges(6,3) * t60;
t23 = mrSges(5,1) * t75 - mrSges(5,3) * t60;
t167 = t22 + t23;
t24 = -mrSges(6,2) * t75 + mrSges(6,3) * t61;
t25 = -mrSges(5,2) * t75 + mrSges(5,3) * t61;
t166 = t24 + t25;
t160 = qJD(3) * t91;
t165 = t63 * t158 + t96 * t160;
t162 = qJ(5) * t89;
t154 = qJDD(1) * t89;
t149 = t78 * t173;
t145 = t96 * t162;
t142 = t94 * t160;
t141 = t89 * t159;
t140 = t89 * t158;
t80 = -pkin(2) - t190;
t136 = qJ(5) * t161;
t17 = -t61 * mrSges(6,1) + t60 * mrSges(6,2);
t10 = t96 * t38 - t45 * t94;
t121 = -mrSges(4,1) * t153 + mrSges(4,2) * t154;
t119 = mrSges(5,1) * t94 + mrSges(5,2) * t96;
t117 = t45 * t91 + t179;
t41 = t173 * t85 - t84 * t96;
t39 = t173 * t84 + t85 * t96;
t115 = (mrSges(6,1) * t94 + t187) * t89;
t8 = -t136 * t96 + t10;
t67 = qJDD(1) * t80 + qJDD(3);
t64 = (pkin(4) * t158 + qJD(3)) * t89;
t62 = t119 * t89;
t59 = (pkin(4) * t94 + t78) * t89;
t58 = t119 * t161;
t57 = qJD(1) * t115;
t52 = t96 * t63;
t42 = -t172 * t85 - t177;
t40 = t172 * t84 - t176;
t19 = t52 - t149;
t18 = -mrSges(5,1) * t61 + mrSges(5,2) * t60;
t16 = -t162 * t94 + t20;
t15 = -qJD(4) * t20 - t142;
t14 = -qJD(4) * t149 + t165;
t13 = -t145 + t52 + (-t78 * t94 - pkin(4)) * t91;
t12 = -pkin(4) * t61 + qJDD(5) + t36;
t9 = -t136 * t94 + t11;
t7 = -t142 - t96 * t157 + (-t66 + (-t63 + t162) * t94) * qJD(4);
t6 = -t94 * t157 + (-t145 - t149) * qJD(4) + t165;
t5 = pkin(4) * t76 + t8;
t1 = pkin(4) * t75 - qJ(5) * t60 - t133 * t96 + t4;
t26 = [(-t1 * mrSges(6,1) + Ifges(4,4) * t154 + Ifges(4,2) * t153 - t213 * t194 - t212 * t193 - t225 * t192 + t197 - t207 / 0.2e1) * t91 - (t212 * t75 + t214 * t60 + t226 * t61) * t175 / 0.2e1 + (t213 * t75 + t214 * t61 + t227 * t60) * t174 / 0.2e1 + t224 * t193 + t223 * t194 + m(5) * (t10 * t15 + t11 * t14 + t19 * t4 + t20 * t3) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t92 * mrSges(3,1) - 0.2e1 * t90 * mrSges(3,2) + m(3) * (t90 ^ 2 + t92 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t65 * t220 + (-t1 * t174 - t140 * t9 + t141 * t5 - t175 * t2) * mrSges(6,3) + (t218 - t219 * t89 / 0.2e1) * qJD(4) + t198 * t150 + (t10 * t141 - t11 * t140 - t174 * t4 - t175 * t3) * mrSges(5,3) + m(6) * (t1 * t13 + t12 * t59 + t16 * t2 + t21 * t64 + t5 * t7 + t6 * t9) + t36 * t62 + t64 * t57 + t6 * t53 + t14 * t54 + t7 * t55 + t15 * t56 + t59 * t17 + t13 * t22 + t19 * t23 + t16 * t24 + t20 * t25 + t12 * t115 + (-t176 * t195 + m(3) * t189 + mrSges(2,1) * t95 + mrSges(2,2) * t97 - t215 * t85 + t217 * t40 + t216 * t39 + t152 * (t85 * qJ(3) - t189) + t200 * t84) * g(3) + (t177 * t195 + m(3) * t188 + mrSges(2,1) * t97 - mrSges(2,2) * t95 + t215 * t84 - t217 * t42 + t216 * t41 + t152 * (-qJ(3) * t84 - t188) + t200 * t85) * g(2) + (m(5) * (qJD(3) * t44 + t36 * t78) + t78 * t18 + qJD(3) * t58 + Ifges(4,4) * t153 + Ifges(4,1) * t154 + (-t212 * t94 + t213 * t96) * t192) * t89 + t206 * mrSges(4,3) + m(4) * (t117 * qJD(3) + t206 * t78 + t67 * t80) + t67 * t120 + t80 * t121; m(3) * qJDD(2) + (-m(3) + t152) * g(1) + (-m(5) * t36 - m(6) * t12 - t17 - t18) * t91 + (t166 * t96 - t167 * t94 + (-t163 * t96 - t164 * t94) * qJD(4) + m(5) * (-t10 * t158 - t11 * t159 + t3 * t96 - t4 * t94) + m(6) * (-t1 * t94 - t158 * t5 - t159 * t9 + t2 * t96)) * t89 + (-t36 * t91 + t37 * t89) * m(4); t167 * t96 + t166 * t94 - t196 * t220 - t199 * qJD(4) + m(6) * (t1 * t96 + t2 * t94 + (-t5 * t94 + t9 * t96) * qJD(4)) + m(5) * (t3 * t94 + t4 * t96 + (-t10 * t94 + t11 * t96) * qJD(4)) + m(4) * t67 + ((-t57 - t58) * t89 + t199 * t91 - m(5) * (-t10 * t173 + t11 * t172 + t179) - m(6) * (t172 * t9 - t173 * t5 + t182) - m(4) * t117) * qJD(1) + t121 + (g(2) * t85 + g(3) * t84) * t152; -t198 * t196 + (t56 + t128) * t11 + (-t54 - t129) * t10 + t210 * t1 + (-(-t210 * t94 - t187) * t89 + t62) * g(1) + (-t204 * t39 + t216 * t40) * g(2) + (t204 * t41 + t216 * t42) * g(3) - t8 * t53 - t197 + (t55 + t126 - m(6) * (-t5 + t8)) * t9 - t5 * t127 + t207 + t219 * t161 / 0.2e1 - t218 * qJD(1) + ((-m(6) * t21 - t57) * t143 + t22) * pkin(4); (g(1) * t91 + t12) * m(6) + ((g(2) * t84 - g(3) * t85) * m(6) + (t94 * t53 + t96 * t55 - m(6) * (-t5 * t96 - t9 * t94)) * qJD(1)) * t89 + t17;];
tau = t26;
