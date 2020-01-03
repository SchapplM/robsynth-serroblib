% Calculate vector of inverse dynamics joint torques for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:44
% EndTime: 2020-01-03 11:55:49
% DurationCPUTime: 2.16s
% Computational Cost: add. (2407->283), mult. (3901->370), div. (0->0), fcn. (2356->16), ass. (0->139)
t155 = sin(pkin(8));
t160 = sin(qJ(2));
t195 = qJD(1) * pkin(1);
t183 = t160 * t195;
t108 = t155 * t183;
t157 = cos(pkin(8));
t163 = cos(qJ(2));
t182 = t163 * t195;
t80 = t157 * t182 - t108;
t209 = qJD(4) - t80;
t154 = sin(pkin(9));
t156 = cos(pkin(9));
t174 = -mrSges(5,1) * t156 + t154 * mrSges(5,2);
t151 = pkin(9) + qJ(5);
t139 = sin(t151);
t140 = cos(t151);
t213 = mrSges(6,1) * t140 - t139 * mrSges(6,2);
t212 = mrSges(5,3) + mrSges(6,3);
t159 = sin(qJ(5));
t162 = cos(qJ(5));
t121 = pkin(2) * t155 + qJ(4);
t90 = (-pkin(7) - t121) * t154;
t144 = t156 * pkin(7);
t91 = t121 * t156 + t144;
t54 = t159 * t90 + t162 * t91;
t97 = t154 * t162 + t156 * t159;
t211 = -t54 * qJD(5) - t209 * t97;
t53 = -t159 * t91 + t162 * t90;
t96 = -t154 * t159 + t156 * t162;
t210 = t53 * qJD(5) + t209 * t96;
t153 = qJ(1) + qJ(2);
t141 = pkin(8) + t153;
t125 = sin(t141);
t126 = cos(t141);
t208 = -t125 * pkin(3) + qJ(4) * t126;
t207 = -mrSges(4,1) - t213 + t174;
t186 = t154 ^ 2 + t156 ^ 2;
t206 = t186 * mrSges(5,3);
t152 = qJD(1) + qJD(2);
t71 = t97 * t152;
t204 = t71 / 0.2e1;
t203 = Ifges(6,4) * t71;
t202 = pkin(1) * t163;
t142 = sin(t153);
t129 = pkin(2) * t142;
t143 = cos(t153);
t130 = pkin(2) * t143;
t201 = pkin(2) * t157;
t200 = pkin(4) * t156;
t199 = g(3) * t126;
t70 = t96 * t152;
t198 = mrSges(6,1) * t70 - mrSges(6,2) * t71 - t174 * t152;
t148 = qJDD(1) + qJDD(2);
t185 = qJD(2) * t160;
t181 = pkin(1) * t185;
t94 = -qJD(1) * t181 + qJDD(1) * t202;
t75 = pkin(2) * t148 + t94;
t184 = qJD(2) * t163;
t95 = (qJD(1) * t184 + qJDD(1) * t160) * pkin(1);
t48 = t155 * t75 + t157 * t95;
t32 = qJ(4) * t148 + qJD(4) * t152 + t48;
t24 = t154 * qJDD(3) + t156 * t32;
t192 = t156 * t24;
t136 = t156 * qJDD(3);
t23 = -t154 * t32 + t136;
t191 = t23 * t154;
t102 = pkin(2) * t152 + t182;
t64 = t155 * t102 + t157 * t183;
t58 = qJ(4) * t152 + t64;
t51 = t154 * qJD(3) + t156 * t58;
t190 = t148 * t154;
t189 = t148 * t156;
t188 = t157 * t160;
t131 = pkin(2) + t202;
t84 = pkin(1) * t188 + t155 * t131;
t161 = sin(qJ(1));
t145 = t161 * pkin(1);
t187 = t129 + t145;
t127 = pkin(3) + t200;
t158 = -pkin(7) - qJ(4);
t180 = t125 * t127 + t126 * t158 + t129;
t179 = t126 * pkin(3) + t125 * qJ(4) + t130;
t88 = t96 * qJD(5);
t44 = t97 * t148 + t152 * t88;
t89 = t97 * qJD(5);
t45 = t96 * t148 - t152 * t89;
t11 = -t45 * mrSges(6,1) + t44 * mrSges(6,2);
t47 = -t155 * t95 + t157 * t75;
t63 = t102 * t157 - t108;
t83 = -t155 * t160 * pkin(1) + t131 * t157;
t77 = -pkin(3) - t83;
t82 = -mrSges(5,1) * t189 + mrSges(5,2) * t190;
t178 = -t125 * t158 + t126 * t127 + t130;
t177 = g(2) * t126 + g(3) * t125;
t176 = qJDD(4) - t47;
t175 = qJD(4) - t63;
t138 = t156 * qJD(3);
t50 = -t154 * t58 + t138;
t171 = -t154 * t50 + t156 * t51;
t40 = t138 + (-pkin(7) * t152 - t58) * t154;
t41 = t152 * t144 + t51;
t9 = -t159 * t41 + t162 * t40;
t10 = t159 * t40 + t162 * t41;
t76 = qJ(4) + t84;
t59 = (-pkin(7) - t76) * t154;
t60 = t156 * t76 + t144;
t21 = -t159 * t60 + t162 * t59;
t22 = t159 * t59 + t162 * t60;
t81 = t157 * pkin(1) * t184 - t155 * t181;
t170 = pkin(1) * (t155 * t163 + t188);
t167 = -t142 * mrSges(3,1) - t143 * mrSges(3,2) - t126 * mrSges(4,2) + t207 * t125;
t166 = -t143 * mrSges(3,1) + t142 * mrSges(3,2) + (mrSges(4,2) - t212) * t125 + t207 * t126;
t16 = t136 + (-pkin(7) * t148 - t32) * t154;
t17 = pkin(7) * t189 + t24;
t2 = t9 * qJD(5) + t159 * t16 + t162 * t17;
t3 = -t10 * qJD(5) - t159 * t17 + t16 * t162;
t30 = -t127 * t148 + t176;
t37 = Ifges(6,2) * t70 + Ifges(6,6) * qJD(5) + t203;
t69 = Ifges(6,4) * t70;
t38 = Ifges(6,1) * t71 + Ifges(6,5) * qJD(5) + t69;
t39 = -pkin(3) * t148 + t176;
t52 = -t127 * t152 + t175;
t165 = (Ifges(6,1) * t88 - Ifges(6,4) * t89) * t204 + mrSges(5,3) * t192 + (Ifges(5,4) * t154 + Ifges(5,2) * t156) * t189 + (Ifges(5,1) * t154 + Ifges(5,4) * t156) * t190 + t94 * mrSges(3,1) + t88 * t38 / 0.2e1 + qJD(5) * (Ifges(6,5) * t88 - Ifges(6,6) * t89) / 0.2e1 + t70 * (Ifges(6,4) * t88 - Ifges(6,2) * t89) / 0.2e1 - t89 * t37 / 0.2e1 + t52 * (mrSges(6,1) * t89 + mrSges(6,2) * t88) - t95 * mrSges(3,2) + t47 * mrSges(4,1) + t39 * t174 + (Ifges(4,3) + Ifges(3,3)) * t148 + (-t10 * t89 - t9 * t88) * mrSges(6,3) + (mrSges(6,2) * t30 - mrSges(6,3) * t3 + Ifges(6,1) * t44 + Ifges(6,4) * t45 + Ifges(6,5) * qJDD(5)) * t97 + (-mrSges(6,1) * t30 + mrSges(6,3) * t2 + Ifges(6,4) * t44 + Ifges(6,2) * t45 + Ifges(6,6) * qJDD(5)) * t96;
t164 = cos(qJ(1));
t146 = t164 * pkin(1);
t128 = -pkin(3) - t201;
t105 = -t127 - t201;
t79 = qJD(2) * t170;
t78 = qJD(1) * t170;
t72 = qJD(4) + t81;
t67 = t77 - t200;
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t71;
t61 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t70;
t57 = -pkin(3) * t152 + t175;
t34 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t45;
t33 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t44;
t6 = -t22 * qJD(5) - t97 * t72;
t5 = t21 * qJD(5) + t96 * t72;
t1 = [(-t152 * mrSges(4,1) - t198) * t79 + (-t148 * t84 - t152 * t81 - t48) * mrSges(4,2) + ((-t148 * t160 - t152 * t184) * mrSges(3,2) + (t148 * t163 - t152 * t185) * mrSges(3,1) + (-g(2) * t164 - g(3) * t161 + t160 * t95 + t163 * t94) * m(3)) * pkin(1) + m(5) * (t39 * t77 + t57 * t79 + (t24 * t76 + t51 * t72) * t156 + (-t23 * t76 - t50 * t72) * t154) + t165 + t83 * t148 * mrSges(4,1) + t77 * t82 + t5 * t61 + t6 * t62 + t67 * t11 + t21 * t33 + t22 * t34 + (-t191 + t199 + (t148 * t76 + t152 * t72) * t186) * mrSges(5,3) + (-t164 * mrSges(2,1) + t161 * mrSges(2,2) - m(5) * (t146 + t179) - m(6) * (t146 + t178) - m(4) * (t130 + t146) + t166) * g(2) + m(6) * (t10 * t5 + t2 * t22 + t21 * t3 + t30 * t67 + t52 * t79 + t6 * t9) + m(4) * (t47 * t83 + t48 * t84 - t63 * t79 + t64 * t81) + (-m(5) * (t187 - t208) - t161 * mrSges(2,1) - t164 * mrSges(2,2) - m(4) * t187 - m(6) * (t145 + t180) + t126 * mrSges(6,3) + t167) * g(3) + Ifges(2,3) * qJDD(1); t165 - mrSges(5,3) * t191 + t128 * t82 + t105 * t11 + t53 * t33 + t54 * t34 - t48 * mrSges(4,2) + (t78 * mrSges(4,1) + t80 * mrSges(4,2) + (mrSges(3,1) * t160 + mrSges(3,2) * t163) * t195 + t209 * t206) * t152 + t166 * g(2) + (t212 * t126 + t167) * g(3) + ((mrSges(4,1) * t157 - mrSges(4,2) * t155) * pkin(2) + t121 * t206) * t148 + t198 * t78 + t210 * t61 + t211 * t62 + (-t178 * g(2) - t180 * g(3) + t210 * t10 + t105 * t30 + t2 * t54 + t211 * t9 + t3 * t53 - t52 * t78) * m(6) + (t128 * t39 + (-t191 + t192) * t121 - t57 * t78 - t179 * g(2) + (-t129 + t208) * g(3) + t209 * t171) * m(5) + ((t155 * t48 + t157 * t47) * pkin(2) - t130 * g(2) - t129 * g(3) + t63 * t78 - t64 * t80) * m(4); m(4) * qJDD(3) + t96 * t33 + t97 * t34 + t88 * t61 - t89 * t62 + m(5) * (t154 * t24 + t156 * t23) + m(6) * (t10 * t88 + t2 * t97 + t3 * t96 - t89 * t9) + (-m(4) - m(5) - m(6)) * g(1); -t152 ^ 2 * t206 - t70 * t61 + t71 * t62 + t11 + t82 + (-t10 * t70 + t9 * t71 + t177 + t30) * m(6) + (-t171 * t152 + t177 + t39) * m(5); Ifges(6,5) * t44 + Ifges(6,6) * t45 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t52 * (mrSges(6,1) * t71 + mrSges(6,2) * t70) - t71 * (Ifges(6,1) * t70 - t203) / 0.2e1 + t37 * t204 - qJD(5) * (Ifges(6,5) * t70 - Ifges(6,6) * t71) / 0.2e1 - t9 * t61 + t10 * t62 - g(1) * t213 + (t10 * t71 + t70 * t9) * mrSges(6,3) - (-Ifges(6,2) * t71 + t38 + t69) * t70 / 0.2e1 + (g(2) * t125 - t199) * (mrSges(6,1) * t139 + mrSges(6,2) * t140);];
tau = t1;
