% Calculate vector of inverse dynamics joint torques for
% S5RPRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:21
% EndTime: 2019-12-05 17:49:25
% DurationCPUTime: 1.68s
% Computational Cost: add. (2279->239), mult. (3954->318), div. (0->0), fcn. (2390->16), ass. (0->124)
t139 = sin(qJ(3));
t142 = cos(qJ(3));
t134 = sin(pkin(8));
t181 = pkin(1) * t134;
t166 = qJD(1) * t181;
t136 = cos(pkin(8));
t114 = pkin(1) * t136 + pkin(2);
t98 = t114 * qJD(1);
t66 = -t139 * t166 + t142 * t98;
t158 = qJD(4) - t66;
t133 = sin(pkin(9));
t135 = cos(pkin(9));
t155 = -t135 * mrSges(5,1) + mrSges(5,2) * t133;
t194 = mrSges(4,1) - t155;
t130 = pkin(9) + qJ(5);
t120 = sin(t130);
t122 = cos(t130);
t193 = t122 * mrSges(6,1) - mrSges(6,2) * t120;
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t137 = -pkin(7) - qJ(4);
t94 = t137 * t133;
t125 = t135 * pkin(7);
t95 = qJ(4) * t135 + t125;
t58 = t138 * t94 + t141 * t95;
t87 = t133 * t141 + t135 * t138;
t192 = -t58 * qJD(5) - t158 * t87;
t57 = -t138 * t95 + t141 * t94;
t86 = -t133 * t138 + t135 * t141;
t191 = t57 * qJD(5) + t158 * t86;
t165 = qJD(3) * t181;
t190 = -qJD(1) * t165 + t114 * qJDD(1);
t180 = pkin(4) * t135;
t113 = pkin(3) + t180;
t189 = m(5) * pkin(3) + m(6) * t113 + t193 + t194;
t188 = m(5) * qJ(4) - m(6) * t137 - mrSges(4,2) + mrSges(6,3);
t131 = qJD(1) + qJD(3);
t169 = t133 ^ 2 + t135 ^ 2;
t185 = t169 * t131;
t80 = t139 * t114 + t142 * t181;
t71 = t87 * t131;
t183 = t71 / 0.2e1;
t182 = Ifges(6,4) * t71;
t132 = qJ(1) + pkin(8);
t124 = qJ(3) + t132;
t111 = sin(t124);
t179 = g(2) * t111;
t117 = t135 * qJDD(2);
t127 = qJDD(1) + qJDD(3);
t163 = qJDD(1) * t181;
t168 = qJD(3) * t142;
t39 = t190 * t139 + t142 * t163 + t98 * t168;
t30 = qJ(4) * t127 + qJD(4) * t131 + t39;
t19 = -t133 * t30 + t117;
t175 = t133 * t19;
t20 = t133 * qJDD(2) + t135 * t30;
t173 = t135 * t20;
t172 = t139 * t98;
t67 = t142 * t166 + t172;
t56 = qJ(4) * t131 + t67;
t51 = t133 * qJD(2) + t135 * t56;
t171 = t127 * t133;
t170 = t127 * t135;
t167 = m(4) + m(5) + m(6);
t164 = m(3) + t167;
t82 = t86 * qJD(5);
t47 = t87 * t127 + t131 * t82;
t83 = t87 * qJD(5);
t48 = t86 * t127 - t131 * t83;
t12 = -t48 * mrSges(6,1) + t47 * mrSges(6,2);
t162 = t169 * t127;
t79 = t114 * t142 - t139 * t181;
t70 = t86 * t131;
t160 = mrSges(6,1) * t70 - mrSges(6,2) * t71 + t194 * t131;
t77 = -pkin(3) - t79;
t78 = -mrSges(5,1) * t170 + mrSges(5,2) * t171;
t112 = cos(t124);
t159 = -g(2) * t112 - g(3) * t111;
t157 = -t175 + t179;
t156 = t167 * pkin(2) + mrSges(3,1);
t119 = t135 * qJD(2);
t50 = -t133 * t56 + t119;
t153 = -t133 * t50 + t135 * t51;
t43 = t119 + (-pkin(7) * t131 - t56) * t133;
t44 = t131 * t125 + t51;
t9 = -t138 * t44 + t141 * t43;
t10 = t138 * t43 + t141 * t44;
t76 = qJ(4) + t80;
t59 = (-pkin(7) - t76) * t133;
t60 = t135 * t76 + t125;
t25 = -t138 * t60 + t141 * t59;
t26 = t138 * t59 + t141 * t60;
t40 = -qJD(3) * t172 - t139 * t163 + t190 * t142;
t72 = t114 * t168 - t139 * t165;
t152 = t164 * pkin(1) + mrSges(2,1);
t151 = qJDD(4) - t40;
t147 = t188 * t111 + t189 * t112;
t146 = t189 * t111 + (-mrSges(5,3) - t188) * t112;
t13 = t117 + (-pkin(7) * t127 - t30) * t133;
t14 = pkin(7) * t170 + t20;
t2 = t9 * qJD(5) + t13 * t138 + t14 * t141;
t24 = -t113 * t127 + t151;
t3 = -t10 * qJD(5) + t13 * t141 - t138 * t14;
t33 = -pkin(3) * t127 + t151;
t37 = Ifges(6,2) * t70 + Ifges(6,6) * qJD(5) + t182;
t68 = Ifges(6,4) * t70;
t38 = Ifges(6,1) * t71 + Ifges(6,5) * qJD(5) + t68;
t52 = -t113 * t131 + t158;
t145 = Ifges(4,3) * t127 + t52 * (mrSges(6,1) * t83 + mrSges(6,2) * t82) + mrSges(5,3) * t173 + t82 * t38 / 0.2e1 - t83 * t37 / 0.2e1 + t33 * t155 + t70 * (Ifges(6,4) * t82 - Ifges(6,2) * t83) / 0.2e1 + (Ifges(6,1) * t82 - Ifges(6,4) * t83) * t183 + t40 * mrSges(4,1) + qJD(5) * (Ifges(6,5) * t82 - Ifges(6,6) * t83) / 0.2e1 + (Ifges(5,1) * t133 + Ifges(5,4) * t135) * t171 + (Ifges(5,4) * t133 + Ifges(5,2) * t135) * t170 + (-t10 * t83 - t9 * t82) * mrSges(6,3) + (t24 * mrSges(6,2) - t3 * mrSges(6,3) + Ifges(6,1) * t47 + Ifges(6,4) * t48 + Ifges(6,5) * qJDD(5)) * t87 + (-t24 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t47 + Ifges(6,2) * t48 + Ifges(6,6) * qJDD(5)) * t86;
t143 = cos(qJ(1));
t140 = sin(qJ(1));
t123 = cos(t132);
t121 = sin(t132);
t73 = t80 * qJD(3);
t69 = qJD(4) + t72;
t65 = t77 - t180;
t62 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t71;
t61 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t70;
t55 = -pkin(3) * t131 + t158;
t35 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t48;
t34 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t47;
t5 = -t26 * qJD(5) - t87 * t69;
t4 = t25 * qJD(5) + t86 * t69;
t1 = [(Ifges(3,3) + Ifges(2,3) + (0.2e1 * t136 * mrSges(3,1) - 0.2e1 * t134 * mrSges(3,2) + m(3) * (t134 ^ 2 + t136 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + m(5) * (t33 * t77 + t55 * t73 + (t20 * t76 + t51 * t69) * t135 + (-t19 * t76 - t50 * t69) * t133) + m(6) * (t10 * t4 + t2 * t26 + t24 * t65 + t25 * t3 + t5 * t9 + t52 * t73) + m(4) * (t39 * t80 + t40 * t79 - t66 * t73 + t67 * t72) + t145 + (t76 * t162 + t69 * t185 + t157) * mrSges(5,3) + (-mrSges(2,2) * t140 - mrSges(3,2) * t121 + t156 * t123 + t152 * t143 + t147) * g(2) + t79 * t127 * mrSges(4,1) + (mrSges(2,2) * t143 + mrSges(3,2) * t123 + t156 * t121 + t152 * t140 + t146) * g(3) - t160 * t73 + t77 * t78 + t4 * t61 + t5 * t62 + t65 * t12 + t25 * t34 + t26 * t35 + (-t80 * t127 - t72 * t131 - t39) * mrSges(4,2); t86 * t34 + t87 * t35 + t82 * t61 - t83 * t62 + (m(3) + m(4)) * qJDD(2) + m(5) * (t133 * t20 + t135 * t19) + m(6) * (t10 * t82 + t2 * t87 + t3 * t86 - t83 * t9) - t164 * g(1); t192 * t62 + t191 * t61 + (qJ(4) * t162 + t158 * t185 + t157) * mrSges(5,3) + t145 + t147 * g(2) + t146 * g(3) + (t131 * t66 - t39) * mrSges(4,2) + t160 * t67 - t113 * t12 - pkin(3) * t78 + t57 * t34 + t58 * t35 + (t191 * t10 - t113 * t24 + t192 * t9 + t2 * t58 + t3 * t57 - t52 * t67) * m(6) + (-pkin(3) * t33 + (t173 - t175) * qJ(4) - t55 * t67 + t158 * t153) * m(5); -t169 * t131 ^ 2 * mrSges(5,3) - t70 * t61 + t71 * t62 + t12 + t78 + (-t10 * t70 + t9 * t71 + t159 + t24) * m(6) + (-t153 * t131 + t159 + t33) * m(5); Ifges(6,5) * t47 + Ifges(6,6) * t48 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t52 * (mrSges(6,1) * t71 + mrSges(6,2) * t70) - t71 * (Ifges(6,1) * t70 - t182) / 0.2e1 + t37 * t183 - qJD(5) * (Ifges(6,5) * t70 - Ifges(6,6) * t71) / 0.2e1 - t9 * t61 + t10 * t62 - g(1) * t193 + (t10 * t71 + t70 * t9) * mrSges(6,3) - (-Ifges(6,2) * t71 + t38 + t68) * t70 / 0.2e1 + (g(3) * t112 - t179) * (mrSges(6,1) * t120 + mrSges(6,2) * t122);];
tau = t1;
