% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR16_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:38:49
% DurationCPUTime: 2.01s
% Computational Cost: add. (1437->285), mult. (3060->374), div. (0->0), fcn. (1402->4), ass. (0->138)
t162 = -qJD(1) / 0.2e1;
t77 = sin(qJ(3));
t115 = t77 * qJD(1);
t78 = cos(qJ(5));
t119 = qJD(3) * t78;
t76 = sin(qJ(5));
t51 = t115 * t76 + t119;
t140 = Ifges(6,4) * t51;
t121 = qJD(3) * t76;
t50 = t115 * t78 - t121;
t79 = cos(qJ(3));
t114 = t79 * qJD(1);
t66 = qJD(5) + t114;
t13 = Ifges(6,2) * t50 + Ifges(6,6) * t66 + t140;
t47 = Ifges(6,4) * t50;
t14 = Ifges(6,1) * t51 + Ifges(6,5) * t66 + t47;
t146 = t78 / 0.2e1;
t147 = -t76 / 0.2e1;
t112 = qJD(3) * qJ(4);
t81 = -pkin(1) - pkin(6);
t63 = qJD(1) * t81 + qJD(2);
t55 = t77 * t63;
t36 = -pkin(4) * t115 + t55;
t32 = t36 + t112;
t98 = mrSges(6,1) * t78 - mrSges(6,2) * t76;
t102 = pkin(4) * qJD(1) - t63;
t37 = t102 * t79;
t165 = qJD(4) + t37;
t80 = -pkin(3) - pkin(7);
t23 = qJD(3) * t80 + t165;
t126 = pkin(3) * t115 + qJD(1) * qJ(2);
t123 = qJ(4) * t79;
t91 = pkin(7) * t77 - t123;
t31 = qJD(1) * t91 + t126;
t7 = t23 * t78 - t31 * t76;
t8 = t23 * t76 + t31 * t78;
t99 = t7 * t76 - t8 * t78;
t167 = t99 * mrSges(6,3) - t13 * t146 + t14 * t147 + t32 * t98;
t148 = t66 / 0.2e1;
t150 = t51 / 0.2e1;
t152 = t50 / 0.2e1;
t160 = qJD(3) / 0.2e1;
t161 = -qJD(3) / 0.2e1;
t113 = qJ(4) * qJD(1);
t40 = -t113 * t79 + t126;
t45 = -t55 - t112;
t134 = Ifges(6,6) * t78;
t137 = Ifges(6,5) * t76;
t92 = t134 + t137;
t139 = Ifges(6,4) * t76;
t93 = Ifges(6,2) * t78 + t139;
t138 = Ifges(6,4) * t78;
t95 = Ifges(6,1) * t76 + t138;
t166 = -t45 * mrSges(5,1) + t40 * mrSges(5,2) - Ifges(5,5) * t160 - Ifges(4,6) * t161 - t92 * t148 - t95 * t150 - t93 * t152 + t167 + ((-Ifges(4,4) - Ifges(5,6)) * t79 + (Ifges(4,2) + Ifges(5,3)) * t77) * t162;
t164 = qJD(1) * (qJ(2) * (m(3) + m(4)) + mrSges(3,3));
t163 = (-mrSges(4,1) + mrSges(5,2)) * qJD(3);
t159 = t77 * pkin(3) + qJ(2);
t105 = qJD(3) * t115;
t71 = pkin(3) * t114;
t74 = qJD(1) * qJD(2);
t110 = qJ(4) * t105 + qJD(3) * t71 + t74;
t88 = (qJD(3) * pkin(7) - qJD(4)) * t79;
t19 = qJD(1) * t88 + t110;
t120 = qJD(3) * t77;
t33 = t102 * t120;
t1 = qJD(5) * t7 + t19 * t78 - t33 * t76;
t2 = -qJD(5) * t8 - t19 * t76 - t33 * t78;
t100 = t1 * t76 + t2 * t78;
t111 = qJD(3) * qJD(5);
t116 = qJD(5) * t78;
t118 = qJD(3) * t79;
t25 = -t76 * t111 + (t116 * t77 + t118 * t76) * qJD(1);
t117 = qJD(5) * t76;
t26 = -t78 * t111 + (-t117 * t77 + t118 * t78) * qJD(1);
t158 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t25 + Ifges(6,6) * t26;
t155 = t25 / 0.2e1;
t154 = t26 / 0.2e1;
t153 = -t50 / 0.2e1;
t151 = -t51 / 0.2e1;
t149 = -t66 / 0.2e1;
t142 = pkin(4) - t81;
t141 = Ifges(4,4) * t77;
t136 = Ifges(6,5) * t78;
t135 = Ifges(6,6) * t76;
t131 = t63 * t79;
t122 = qJD(3) * mrSges(4,2);
t59 = -mrSges(4,3) * t115 - t122;
t61 = mrSges(5,1) * t115 - qJD(3) * mrSges(5,3);
t129 = t59 - t61;
t20 = -mrSges(6,1) * t50 + mrSges(6,2) * t51;
t128 = -t61 + t20;
t127 = (mrSges(5,1) + mrSges(4,3)) * t114 + t163;
t52 = t77 * t113 + t71;
t125 = qJ(2) * mrSges(4,1);
t124 = qJ(2) * mrSges(4,2);
t109 = -0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * Ifges(5,6);
t108 = Ifges(5,4) / 0.2e1 - Ifges(4,5) / 0.2e1;
t107 = Ifges(5,5) / 0.2e1 - Ifges(4,6) / 0.2e1;
t106 = pkin(3) * t118 + t77 * t112 + qJD(2);
t39 = -qJD(3) * pkin(3) + qJD(4) - t131;
t104 = -t39 + t131;
t103 = qJD(3) * t142;
t97 = mrSges(6,1) * t76 + mrSges(6,2) * t78;
t96 = Ifges(6,1) * t78 - t139;
t94 = -Ifges(6,2) * t76 + t138;
t15 = -mrSges(6,1) * t105 - mrSges(6,3) * t25;
t16 = mrSges(6,2) * t105 + mrSges(6,3) * t26;
t90 = t78 * t15 + t76 * t16;
t29 = -mrSges(6,2) * t66 + mrSges(6,3) * t50;
t30 = mrSges(6,1) * t66 - mrSges(6,3) * t51;
t89 = t78 * t29 - t76 * t30;
t46 = t91 + t159;
t58 = t142 * t79;
t18 = t46 * t78 + t58 * t76;
t17 = -t46 * t76 + t58 * t78;
t85 = -m(6) * t99 + t89;
t69 = Ifges(5,6) * t115;
t84 = t40 * mrSges(5,3) + t8 * mrSges(6,2) - t66 * Ifges(6,3) - t51 * Ifges(6,5) - t50 * Ifges(6,6) + Ifges(4,5) * t161 + (Ifges(4,1) * t79 - t141) * t162 + Ifges(5,4) * t160 - Ifges(5,2) * t114 / 0.2e1 + t69 / 0.2e1 - t7 * mrSges(6,1);
t57 = t142 * t77;
t56 = -t123 + t159;
t54 = qJD(1) * (mrSges(4,1) * t77 + mrSges(4,2) * t79);
t53 = (-t77 * mrSges(5,2) - t79 * mrSges(5,3)) * qJD(1);
t49 = t79 * t103;
t48 = t77 * t103;
t38 = (-qJD(4) - t131) * qJD(3);
t35 = pkin(7) * t114 + t52;
t34 = -qJD(4) * t79 + t106;
t28 = (-t37 + qJD(4)) * qJD(3);
t27 = t88 + t106;
t24 = -qJD(4) * t114 + t110;
t11 = t35 * t78 + t36 * t76;
t10 = -t35 * t76 + t36 * t78;
t9 = -mrSges(6,1) * t26 + mrSges(6,2) * t25;
t6 = t25 * Ifges(6,1) + t26 * Ifges(6,4) - Ifges(6,5) * t105;
t5 = t25 * Ifges(6,4) + t26 * Ifges(6,2) - Ifges(6,6) * t105;
t4 = -qJD(5) * t18 - t27 * t76 - t48 * t78;
t3 = qJD(5) * t17 + t27 * t78 - t48 * t76;
t12 = [t17 * t15 + t18 * t16 - t49 * t20 + t3 * t29 + t4 * t30 + t34 * t53 - t57 * t9 + m(6) * (t1 * t18 + t17 * t2 - t28 * t57 + t3 * t8 - t32 * t49 + t4 * t7) + m(5) * (t24 * t56 + t34 * t40) + (t54 + 0.2e1 * t164) * qJD(2) + (-t24 * mrSges(5,3) + mrSges(4,2) * t74 + ((-m(5) * t45 + t129) * t81 + (-t56 * mrSges(5,2) + t109 * t79 + 0.2e1 * t125) * qJD(1) + t107 * qJD(3) - t166) * qJD(3) + t158) * t79 + (-t28 * t98 + t76 * t6 / 0.2e1 + t95 * t155 + t93 * t154 - t24 * mrSges(5,2) + t5 * t146 + mrSges(4,1) * t74 + (-m(5) * t81 + mrSges(5,1)) * t38 + (t1 * t78 - t2 * t76) * mrSges(6,3) + ((-t135 + t136) * t148 + t94 * t152 + t96 * t150 + t32 * t97 + t14 * t146 + t13 * t147 + (-t7 * t78 - t76 * t8) * mrSges(6,3)) * qJD(5) + (t108 * qJD(3) + t104 * mrSges(5,1) + (-m(5) * t104 + t127) * t81 + ((-t137 / 0.2e1 - t134 / 0.2e1 - t109) * t77 - 0.2e1 * t124 + t56 * mrSges(5,3) + (0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - Ifges(6,3)) * t79) * qJD(1) + t84) * qJD(3)) * t77; (t9 + (t76 * t29 + t78 * t30 + t127) * qJD(3) + m(5) * (qJD(3) * t39 - t38) + m(6) * (t119 * t7 + t121 * t8 + t28)) * t77 + (m(6) * (-t116 * t8 + t117 * t7 - t100) - t29 * t116 + t30 * t117 + (t59 + m(5) * (-t45 - t55) + m(6) * t32 + t128) * qJD(3) - t90) * t79 + (-m(5) * t40 - t164 - t53 - t54 - t85) * qJD(1); t96 * t155 + t94 * t154 + t28 * t97 + t6 * t146 + t5 * t147 - t52 * t53 + t37 * t20 - t38 * mrSges(5,3) - t11 * t29 - t10 * t30 + qJ(4) * t9 + t90 * t80 + t128 * qJD(4) - t100 * mrSges(6,3) + ((-t122 - t129) * t79 + (t163 - t127) * t77) * t63 + (t149 * t92 + t151 * t95 + t153 * t93 + t80 * t85 + t167) * qJD(5) + ((t39 * mrSges(5,1) - t69 / 0.2e1 + (t124 - t141 / 0.2e1) * qJD(1) + (-t136 / 0.2e1 + t135 / 0.2e1 + pkin(3) * mrSges(5,1) + t108) * qJD(3) - t84) * t77 + ((-t125 + (Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1) * t79) * qJD(1) + (-qJ(4) * mrSges(5,1) + t107) * qJD(3) + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t115 + t166) * t79) * qJD(1) + (qJ(4) * t28 - t10 * t7 + t100 * t80 - t11 * t8 + t165 * t32) * m(6) + (-qJ(4) * t38 - qJD(4) * t45 - t40 * t52 + (-pkin(3) * t120 - t39 * t77 + t45 * t79) * t63) * m(5); t89 * qJD(5) + (m(5) * t55 - t128) * qJD(3) + (-mrSges(5,1) * t120 + (t53 + t89) * t79) * qJD(1) - m(5) * (-qJD(3) * t45 - t114 * t40) + t90 + (-qJD(3) * t32 - t66 * t99 + t100) * m(6); -Ifges(6,3) * t105 - t32 * (mrSges(6,1) * t51 + mrSges(6,2) * t50) + (Ifges(6,1) * t50 - t140) * t151 + t13 * t150 + (Ifges(6,5) * t50 - Ifges(6,6) * t51) * t149 - t7 * t29 + t8 * t30 + (t50 * t7 + t51 * t8) * mrSges(6,3) + (-Ifges(6,2) * t51 + t14 + t47) * t153 + t158;];
tauc = t12(:);
