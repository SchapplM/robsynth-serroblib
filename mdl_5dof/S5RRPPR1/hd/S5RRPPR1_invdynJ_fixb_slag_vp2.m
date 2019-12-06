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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:18:10
% EndTime: 2019-12-05 18:18:15
% DurationCPUTime: 1.79s
% Computational Cost: add. (2407->255), mult. (3901->343), div. (0->0), fcn. (2356->16), ass. (0->130)
t143 = sin(pkin(8));
t148 = sin(qJ(2));
t183 = qJD(1) * pkin(1);
t172 = t148 * t183;
t106 = t143 * t172;
t145 = cos(pkin(8));
t151 = cos(qJ(2));
t171 = t151 * t183;
t80 = t145 * t171 - t106;
t199 = qJD(4) - t80;
t142 = sin(pkin(9));
t144 = cos(pkin(9));
t163 = -t144 * mrSges(5,1) + mrSges(5,2) * t142;
t139 = pkin(9) + qJ(5);
t129 = sin(t139);
t130 = cos(t139);
t202 = t130 * mrSges(6,1) - mrSges(6,2) * t129;
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t115 = pkin(2) * t143 + qJ(4);
t90 = (-pkin(7) - t115) * t142;
t134 = t144 * pkin(7);
t91 = t115 * t144 + t134;
t53 = -t147 * t91 + t150 * t90;
t96 = -t142 * t147 + t144 * t150;
t201 = qJD(5) * t53 + t199 * t96;
t54 = t147 * t90 + t150 * t91;
t97 = t142 * t150 + t144 * t147;
t200 = -qJD(5) * t54 - t199 * t97;
t126 = t144 * qJDD(3);
t136 = qJDD(1) + qJDD(2);
t140 = qJD(1) + qJD(2);
t184 = pkin(1) * qJD(2);
t170 = t148 * t184;
t190 = pkin(1) * t151;
t94 = -qJD(1) * t170 + qJDD(1) * t190;
t75 = pkin(2) * t136 + t94;
t174 = qJD(2) * t151;
t95 = (qJD(1) * t174 + qJDD(1) * t148) * pkin(1);
t48 = t143 * t75 + t145 * t95;
t32 = qJ(4) * t136 + qJD(4) * t140 + t48;
t23 = -t142 * t32 + t126;
t24 = t142 * qJDD(3) + t144 * t32;
t198 = -t142 * t23 + t144 * t24;
t188 = pkin(4) * t144;
t119 = pkin(3) + t188;
t197 = m(5) * pkin(3) + m(6) * t119 + mrSges(4,1) - t163 + t202;
t196 = -m(5) * qJ(4) - mrSges(6,3) - mrSges(5,3) + mrSges(4,2) + m(6) * (-pkin(7) - qJ(4));
t169 = mrSges(5,3) * (t142 ^ 2 + t144 ^ 2);
t71 = t97 * t140;
t192 = t71 / 0.2e1;
t191 = Ifges(6,4) * t71;
t189 = pkin(2) * t145;
t70 = t96 * t140;
t187 = -mrSges(6,1) * t70 + mrSges(6,2) * t71 + t163 * t140;
t100 = pkin(2) * t140 + t171;
t64 = t143 * t100 + t145 * t172;
t58 = qJ(4) * t140 + t64;
t51 = t142 * qJD(3) + t144 * t58;
t178 = t136 * t142;
t177 = t136 * t144;
t176 = t145 * t148;
t141 = qJ(1) + qJ(2);
t121 = pkin(2) + t190;
t84 = pkin(1) * t176 + t143 * t121;
t173 = m(4) + m(5) + m(6);
t88 = t96 * qJD(5);
t44 = t136 * t97 + t140 * t88;
t89 = t97 * qJD(5);
t45 = t136 * t96 - t140 * t89;
t11 = -t45 * mrSges(6,1) + t44 * mrSges(6,2);
t47 = -t143 * t95 + t145 * t75;
t63 = t100 * t145 - t106;
t83 = -t143 * t148 * pkin(1) + t121 * t145;
t77 = -pkin(3) - t83;
t82 = -mrSges(5,1) * t177 + mrSges(5,2) * t178;
t131 = pkin(8) + t141;
t117 = sin(t131);
t118 = cos(t131);
t168 = -g(2) * t118 - g(3) * t117;
t167 = qJDD(4) - t47;
t166 = qJD(4) - t63;
t165 = pkin(2) * t173 + mrSges(3,1);
t164 = mrSges(3,1) * t148 + mrSges(3,2) * t151;
t128 = t144 * qJD(3);
t50 = -t142 * t58 + t128;
t161 = -t142 * t50 + t144 * t51;
t40 = t128 + (-pkin(7) * t140 - t58) * t142;
t41 = t134 * t140 + t51;
t9 = -t147 * t41 + t150 * t40;
t10 = t147 * t40 + t150 * t41;
t76 = qJ(4) + t84;
t59 = (-pkin(7) - t76) * t142;
t60 = t144 * t76 + t134;
t21 = -t147 * t60 + t150 * t59;
t22 = t147 * t59 + t150 * t60;
t81 = t145 * pkin(1) * t174 - t143 * t170;
t160 = mrSges(2,1) + (m(3) + t173) * pkin(1);
t159 = pkin(1) * (t143 * t151 + t176);
t132 = sin(t141);
t133 = cos(t141);
t155 = -t132 * mrSges(3,2) - t196 * t117 + t197 * t118 + t165 * t133;
t154 = t133 * mrSges(3,2) + t197 * t117 + t196 * t118 + t165 * t132;
t16 = t126 + (-pkin(7) * t136 - t32) * t142;
t17 = pkin(7) * t177 + t24;
t2 = qJD(5) * t9 + t147 * t16 + t150 * t17;
t3 = -qJD(5) * t10 - t147 * t17 + t150 * t16;
t30 = -t119 * t136 + t167;
t37 = Ifges(6,2) * t70 + Ifges(6,6) * qJD(5) + t191;
t69 = Ifges(6,4) * t70;
t38 = Ifges(6,1) * t71 + Ifges(6,5) * qJD(5) + t69;
t39 = -pkin(3) * t136 + t167;
t52 = -t119 * t140 + t166;
t153 = (Ifges(5,4) * t142 + Ifges(5,2) * t144) * t177 + (Ifges(5,1) * t142 + Ifges(5,4) * t144) * t178 - t89 * t37 / 0.2e1 + t94 * mrSges(3,1) - t95 * mrSges(3,2) + t88 * t38 / 0.2e1 + t47 * mrSges(4,1) - t48 * mrSges(4,2) + t39 * t163 + (Ifges(6,1) * t88 - Ifges(6,4) * t89) * t192 + t70 * (Ifges(6,4) * t88 - Ifges(6,2) * t89) / 0.2e1 + t52 * (mrSges(6,1) * t89 + mrSges(6,2) * t88) + qJD(5) * (Ifges(6,5) * t88 - Ifges(6,6) * t89) / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t136 + (-t10 * t89 - t9 * t88) * mrSges(6,3) + t198 * mrSges(5,3) + (t30 * mrSges(6,2) - t3 * mrSges(6,3) + Ifges(6,1) * t44 + Ifges(6,4) * t45 + Ifges(6,5) * qJDD(5)) * t97 + (-t30 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t44 + Ifges(6,2) * t45 + Ifges(6,6) * qJDD(5)) * t96;
t152 = cos(qJ(1));
t149 = sin(qJ(1));
t120 = -pkin(3) - t189;
t103 = -t119 - t189;
t79 = qJD(2) * t159;
t78 = qJD(1) * t159;
t72 = qJD(4) + t81;
t67 = t77 - t188;
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t71;
t61 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t70;
t57 = -pkin(3) * t140 + t166;
t34 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t45;
t33 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t44;
t6 = -qJD(5) * t22 - t72 * t97;
t5 = qJD(5) * t21 + t72 * t96;
t1 = [m(3) * (t148 * t95 + t151 * t94) * pkin(1) + (t152 * mrSges(2,2) + t149 * t160 + t154) * g(3) + m(5) * (t39 * t77 + t57 * t79 + (t24 * t76 + t51 * t72) * t144 + (-t23 * t76 - t50 * t72) * t142) + (mrSges(4,1) * t83 - mrSges(4,2) * t84 + (mrSges(3,1) * t151 - mrSges(3,2) * t148) * pkin(1) + t76 * t169) * t136 + (-t79 * mrSges(4,1) - t81 * mrSges(4,2) - t164 * t184 + t169 * t72) * t140 + t187 * t79 + t153 + t77 * t82 + t5 * t61 + t6 * t62 + t67 * t11 + t21 * t33 + t22 * t34 + (-t149 * mrSges(2,2) + t152 * t160 + t155) * g(2) + Ifges(2,3) * qJDD(1) + m(4) * (t47 * t83 + t48 * t84 - t63 * t79 + t64 * t81) + m(6) * (t10 * t5 + t2 * t22 + t21 * t3 + t30 * t67 + t52 * t79 + t6 * t9); -t187 * t78 + t200 * t62 + t201 * t61 + t153 + t120 * t82 + t103 * t11 + t53 * t33 + t54 * t34 + (t78 * mrSges(4,1) + t80 * mrSges(4,2) + t164 * t183 + t199 * t169) * t140 + t155 * g(2) + ((mrSges(4,1) * t145 - mrSges(4,2) * t143) * pkin(2) + t115 * t169) * t136 + t154 * g(3) + (t201 * t10 + t103 * t30 + t2 * t54 + t200 * t9 + t3 * t53 - t52 * t78) * m(6) + (t198 * t115 + t120 * t39 + t199 * t161 - t57 * t78) * m(5) + ((t143 * t48 + t145 * t47) * pkin(2) + t63 * t78 - t64 * t80) * m(4); m(4) * qJDD(3) + t96 * t33 + t97 * t34 + t88 * t61 - t89 * t62 + m(5) * (t142 * t24 + t144 * t23) + m(6) * (t10 * t88 + t2 * t97 + t3 * t96 - t89 * t9) - t173 * g(1); -t140 ^ 2 * t169 - t70 * t61 + t71 * t62 + t11 + t82 + (-t10 * t70 + t9 * t71 + t168 + t30) * m(6) + (-t140 * t161 + t168 + t39) * m(5); Ifges(6,5) * t44 + Ifges(6,6) * t45 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t52 * (mrSges(6,1) * t71 + mrSges(6,2) * t70) - t71 * (Ifges(6,1) * t70 - t191) / 0.2e1 + t37 * t192 - qJD(5) * (Ifges(6,5) * t70 - Ifges(6,6) * t71) / 0.2e1 - t9 * t61 + t10 * t62 - g(1) * t202 + (t10 * t71 + t70 * t9) * mrSges(6,3) - (-Ifges(6,2) * t71 + t38 + t69) * t70 / 0.2e1 + (-g(2) * t117 + g(3) * t118) * (mrSges(6,1) * t129 + mrSges(6,2) * t130);];
tau = t1;
