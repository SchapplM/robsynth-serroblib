% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:37:07
% EndTime: 2018-11-23 15:37:09
% DurationCPUTime: 1.77s
% Computational Cost: add. (1638->248), mult. (3324->339), div. (0->0), fcn. (1621->6), ass. (0->120)
t53 = sin(pkin(9)) * pkin(1) + qJ(3);
t111 = qJD(1) * t53;
t152 = m(4) * t111;
t49 = -pkin(7) + t53;
t151 = m(6) * t49 - mrSges(6,3);
t59 = sin(qJ(6));
t122 = Ifges(7,6) * t59;
t61 = cos(qJ(6));
t123 = Ifges(7,5) * t61;
t105 = qJD(5) * t59;
t62 = cos(qJ(5));
t109 = qJD(1) * t62;
t42 = t109 * t61 + t105;
t126 = Ifges(7,4) * t42;
t104 = qJD(5) * t61;
t41 = -t109 * t59 + t104;
t60 = sin(qJ(5));
t110 = qJD(1) * t60;
t52 = qJD(6) + t110;
t13 = Ifges(7,2) * t41 + Ifges(7,6) * t52 + t126;
t131 = t61 / 0.2e1;
t132 = -t59 / 0.2e1;
t135 = t42 / 0.2e1;
t37 = Ifges(7,4) * t41;
t14 = Ifges(7,1) * t42 + Ifges(7,5) * t52 + t37;
t45 = qJD(4) + t111;
t38 = -qJD(1) * pkin(7) + t45;
t31 = -qJD(2) * t60 + t38 * t62;
t22 = -qJD(5) * pkin(5) - t31;
t124 = Ifges(7,4) * t61;
t77 = -Ifges(7,2) * t59 + t124;
t125 = Ifges(7,4) * t59;
t79 = Ifges(7,1) * t61 - t125;
t80 = mrSges(7,1) * t59 + mrSges(7,2) * t61;
t32 = qJD(2) * t62 + t38 * t60;
t23 = qJD(5) * pkin(8) + t32;
t50 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t34 = pkin(5) * t60 - pkin(8) * t62 + t50;
t26 = qJD(1) * t34 - qJD(3);
t7 = -t23 * t59 + t26 * t61;
t8 = t23 * t61 + t26 * t59;
t83 = t8 * t59 + t7 * t61;
t150 = -t83 * mrSges(7,3) + t13 * t132 + t131 * t14 + (-t122 + t123) * t52 / 0.2e1 + t79 * t135 + t77 * t41 / 0.2e1 + t22 * t80;
t113 = Ifges(6,5) * qJD(5);
t127 = Ifges(6,4) * t60;
t148 = t113 / 0.2e1 + (t62 * Ifges(6,1) - t127) * qJD(1) / 0.2e1 - t31 * mrSges(6,3) + t150;
t108 = qJD(5) * mrSges(6,1);
t94 = mrSges(6,3) * t109;
t114 = mrSges(7,1) * t41 - mrSges(7,2) * t42 + t108 - t94;
t147 = m(6) * t31 - m(7) * t22 + t114;
t91 = -Ifges(6,6) * qJD(5) / 0.2e1;
t106 = qJD(5) * t31;
t97 = qJD(1) * qJD(3);
t19 = t60 * t97 + t106;
t87 = pkin(5) * t62 + pkin(8) * t60;
t40 = qJD(5) * t87 + qJD(4);
t33 = t40 * qJD(1);
t1 = qJD(6) * t7 + t19 * t61 + t33 * t59;
t2 = -qJD(6) * t8 - t19 * t59 + t33 * t61;
t85 = t1 * t61 - t2 * t59;
t100 = qJD(6) * t62;
t96 = qJD(5) * qJD(6);
t27 = t61 * t96 + (-t100 * t59 - t104 * t60) * qJD(1);
t28 = -t59 * t96 + (-t100 * t61 + t105 * t60) * qJD(1);
t146 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t27 + Ifges(7,6) * t28;
t29 = -mrSges(7,2) * t52 + mrSges(7,3) * t41;
t30 = mrSges(7,1) * t52 - mrSges(7,3) * t42;
t72 = -t59 * t29 - t61 * t30;
t145 = -m(7) * t83 + t72;
t99 = t32 * qJD(5);
t20 = -t62 * t97 + t99;
t9 = -mrSges(7,1) * t28 + mrSges(7,2) * t27;
t144 = -m(6) * (-t20 + t99) + m(7) * (-t104 * t8 + t105 * t7 + t20) + t9;
t141 = m(5) / 0.2e1;
t140 = t27 / 0.2e1;
t139 = t28 / 0.2e1;
t138 = -t41 / 0.2e1;
t136 = -t42 / 0.2e1;
t134 = -t52 / 0.2e1;
t117 = t59 * t60;
t116 = t60 * t61;
t115 = -mrSges(4,3) - mrSges(5,2);
t107 = qJD(5) * mrSges(6,2);
t103 = qJD(5) * t62;
t102 = qJD(6) * t59;
t101 = qJD(6) * t61;
t98 = t50 * qJD(1);
t95 = mrSges(6,3) * t110;
t93 = qJD(1) * t103;
t92 = -t113 / 0.2e1;
t39 = -qJD(3) + t98;
t90 = t39 + t98;
t86 = -qJD(6) * t49 * t60 + t40;
t84 = -t1 * t59 - t2 * t61;
t82 = t7 * t59 - t8 * t61;
t81 = t61 * mrSges(7,1) - t59 * mrSges(7,2);
t78 = Ifges(7,1) * t59 + t124;
t76 = Ifges(7,2) * t61 + t125;
t74 = Ifges(7,5) * t59 + Ifges(7,6) * t61;
t15 = mrSges(7,1) * t93 - mrSges(7,3) * t27;
t16 = -mrSges(7,2) * t93 + mrSges(7,3) * t28;
t73 = -t59 * t15 + t61 * t16;
t47 = -t95 - t107;
t71 = t61 * t29 - t59 * t30 + t47;
t67 = qJD(3) * t60 + qJD(6) * t34 + t103 * t49;
t66 = t39 * mrSges(6,1) + t7 * mrSges(7,1) + t52 * Ifges(7,3) + t42 * Ifges(7,5) + t41 * Ifges(7,6) + t91 - (Ifges(6,4) * t62 - t60 * Ifges(6,2)) * qJD(1) / 0.2e1 - t8 * mrSges(7,2);
t65 = m(6) * (t19 - t106) + m(7) * (qJD(5) * t22 - t101 * t7 - t102 * t8 + t85) + t72 * qJD(6) + t73;
t63 = qJD(1) ^ 2;
t51 = Ifges(7,3) * t93;
t44 = t87 * qJD(1);
t43 = qJD(1) * (mrSges(6,1) * t60 + mrSges(6,2) * t62);
t18 = t116 * t49 + t34 * t59;
t17 = -t117 * t49 + t34 * t61;
t11 = t31 * t61 + t44 * t59;
t10 = -t31 * t59 + t44 * t61;
t6 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t93;
t5 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t93;
t4 = -t59 * t67 + t61 * t86;
t3 = t59 * t86 + t61 * t67;
t12 = [qJD(4) * t43 + t17 * t15 + t18 * t16 + t3 * t29 + t4 * t30 + m(7) * (t1 * t18 + t17 * t2 + t3 * t8 + t4 * t7) + (qJD(1) * qJD(4) * mrSges(6,1) + t51 / 0.2e1 + (0.3e1 / 0.2e1 * Ifges(6,4) * t110 + t92 - t90 * mrSges(6,2) - t147 * t49 - t148) * qJD(5) + t146 + t151 * t19) * t60 + (t79 * t140 + t77 * t139 + t6 * t131 + t5 * t132 - t49 * t9 + t84 * mrSges(7,3) + (t151 * t32 + t49 * t47 + t66 + t91) * qJD(5) + (t74 * t134 + t78 * t136 + t76 * t138 + t22 * t81 - t61 * t13 / 0.2e1 + t14 * t132 + t82 * mrSges(7,3)) * qJD(6) + (qJD(4) * mrSges(6,2) + ((-0.3e1 / 0.2e1 * Ifges(6,4) + t123 / 0.2e1 - t122 / 0.2e1) * t62 + t50 * mrSges(6,1) + (0.3e1 / 0.2e1 * Ifges(6,2) + Ifges(7,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(6,1)) * t60) * qJD(5)) * qJD(1) + (mrSges(6,3) + t80 + (-m(7) - m(6)) * t49) * t20) * t62 + 0.2e1 * (qJD(1) * mrSges(5,3) + (m(6) / 0.2e1 + t141) * t90) * qJD(4) + (-0.2e1 * t115 * qJD(1) + 0.2e1 * (t45 + t111) * t141 + 0.2e1 * t152 + (m(6) * t32 + t47) * t60 + t147 * t62) * qJD(3); ((-t71 - t95) * qJD(5) + t144) * t60 + ((-t94 - t114) * qJD(5) + t65) * t62; -t29 * t101 + t30 * t102 + m(7) * (t82 * qJD(6) + t84) - t59 * t16 - t61 * t15 + t115 * t63 + (-t152 + (-t45 - qJD(4)) * m(5) + (-t108 - t114) * t62 + (-t71 + t107) * t60 - m(7) * (t116 * t8 - t117 * t7 - t22 * t62) + (-t31 * t62 - t32 * t60 - qJD(4)) * m(6)) * qJD(1); -t63 * mrSges(5,3) + (qJD(5) * t71 - t144) * t62 + (-qJD(5) * t114 + t65) * t60 + (-m(6) * t39 - t43 + (-t39 + qJD(3)) * m(5) + t145) * qJD(1); t78 * t140 + t76 * t139 + t5 * t131 + t59 * t6 / 0.2e1 - t31 * t47 - t11 * t29 - t10 * t30 - t19 * mrSges(6,2) - pkin(5) * t9 + t114 * t32 + (-mrSges(6,1) - t81) * t20 + t85 * mrSges(7,3) + t150 * qJD(6) + ((t91 + qJD(5) * t74 / 0.2e1 + t32 * mrSges(6,3) + Ifges(6,4) * t109 / 0.2e1 - t66) * t62 + (t39 * mrSges(6,2) + t92 + (-t127 / 0.2e1 + (Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1) * t62) * qJD(1) + t148) * t60) * qJD(1) + (-pkin(5) * t20 - t10 * t7 - t11 * t8 - t22 * t32) * m(7) + (m(7) * t85 + qJD(6) * t145 + t73) * pkin(8); t51 - t22 * (mrSges(7,1) * t42 + mrSges(7,2) * t41) + (Ifges(7,1) * t41 - t126) * t136 + t13 * t135 + (Ifges(7,5) * t41 - Ifges(7,6) * t42) * t134 - t7 * t29 + t8 * t30 + (t41 * t7 + t42 * t8) * mrSges(7,3) + (-Ifges(7,2) * t42 + t14 + t37) * t138 + t146;];
tauc  = t12(:);
