% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:28
% DurationCPUTime: 1.90s
% Computational Cost: add. (1731->244), mult. (4766->356), div. (0->0), fcn. (3273->6), ass. (0->124)
t110 = qJD(2) + qJD(4);
t113 = sin(qJ(4));
t115 = cos(qJ(4));
t111 = sin(pkin(7));
t112 = cos(pkin(7));
t114 = sin(qJ(2));
t116 = cos(qJ(2));
t93 = -t111 * t114 + t112 * t116;
t82 = t93 * qJD(1);
t127 = qJD(1) * t116;
t128 = qJD(1) * t114;
t83 = -t111 * t127 - t112 * t128;
t120 = t113 * t83 + t115 * t82;
t38 = Ifges(5,4) * t120;
t44 = t113 * t82 - t115 * t83;
t13 = Ifges(5,1) * t44 + Ifges(5,5) * t110 + t38;
t139 = Ifges(5,4) * t44;
t94 = t111 * t116 + t112 * t114;
t84 = t94 * qJD(2);
t73 = qJD(1) * t84;
t85 = t93 * qJD(2);
t74 = qJD(1) * t85;
t17 = qJD(4) * t120 - t113 * t73 + t115 * t74;
t18 = -qJD(4) * t44 - t113 * t74 - t115 * t73;
t137 = -qJ(3) - pkin(5);
t119 = qJD(2) * t137;
t80 = qJD(3) * t116 + t114 * t119;
t63 = t80 * qJD(1);
t81 = -t114 * qJD(3) + t116 * t119;
t64 = t81 * qJD(1);
t28 = -t111 * t63 + t112 * t64;
t22 = -pkin(6) * t74 + t28;
t29 = t111 * t64 + t112 * t63;
t23 = -pkin(6) * t73 + t29;
t142 = pkin(6) * t83;
t104 = t137 * t116;
t99 = qJD(1) * t104;
t86 = t111 * t99;
t103 = t137 * t114;
t98 = qJD(1) * t103;
t92 = qJD(2) * pkin(2) + t98;
t46 = t112 * t92 + t86;
t26 = qJD(2) * pkin(3) + t142 + t46;
t143 = pkin(6) * t82;
t135 = t112 * t99;
t47 = t111 * t92 - t135;
t27 = t47 + t143;
t6 = -t113 * t27 + t115 * t26;
t2 = qJD(4) * t6 + t113 * t22 + t115 * t23;
t7 = t113 * t26 + t115 * t27;
t3 = -qJD(4) * t7 - t113 * t23 + t115 * t22;
t108 = -pkin(2) * t116 - pkin(1);
t129 = qJD(1) * t108;
t100 = qJD(3) + t129;
t52 = -t82 * pkin(3) + t100;
t163 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t17 + Ifges(5,6) * t18 - (Ifges(5,5) * t120 - Ifges(5,6) * t44) * t110 / 0.2e1 - (-Ifges(5,2) * t44 + t13 + t38) * t120 / 0.2e1 - t52 * (mrSges(5,1) * t44 + mrSges(5,2) * t120) - (Ifges(5,1) * t120 - t139) * t44 / 0.2e1;
t162 = t120 * t6 + t44 * t7;
t12 = Ifges(5,2) * t120 + Ifges(5,6) * t110 + t139;
t160 = t12 / 0.2e1;
t50 = -t111 * t98 + t135;
t30 = t50 - t143;
t51 = t112 * t98 + t86;
t31 = t51 + t142;
t107 = pkin(2) * t112 + pkin(3);
t138 = pkin(2) * t111;
t78 = t107 * t115 - t113 * t138;
t159 = t78 * qJD(4) - t113 * t30 - t115 * t31;
t79 = t107 * t113 + t115 * t138;
t158 = -t79 * qJD(4) + t113 * t31 - t115 * t30;
t152 = t120 / 0.2e1;
t150 = t44 / 0.2e1;
t148 = -t83 / 0.2e1;
t147 = -t84 / 0.2e1;
t146 = t85 / 0.2e1;
t145 = pkin(1) * mrSges(3,1);
t144 = pkin(1) * mrSges(3,2);
t140 = Ifges(4,4) * t83;
t35 = t111 * t81 + t112 * t80;
t55 = t111 * t103 - t112 * t104;
t136 = Ifges(3,4) * t114;
t133 = Ifges(3,5) * qJD(2);
t132 = Ifges(3,6) * qJD(2);
t131 = qJD(2) * mrSges(3,1);
t130 = qJD(2) * mrSges(3,2);
t126 = qJD(2) * t114;
t125 = pkin(2) * t126;
t124 = t133 / 0.2e1;
t123 = -t132 / 0.2e1;
t122 = t73 * mrSges(4,1) + t74 * mrSges(4,2);
t121 = -t18 * mrSges(5,1) + t17 * mrSges(5,2);
t34 = -t111 * t80 + t112 * t81;
t54 = t112 * t103 + t104 * t111;
t118 = qJD(1) * t125;
t36 = -pkin(6) * t94 + t54;
t37 = pkin(6) * t93 + t55;
t10 = -t113 * t37 + t115 * t36;
t11 = t113 * t36 + t115 * t37;
t48 = -t113 * t94 + t115 * t93;
t49 = t113 * t93 + t115 * t94;
t109 = Ifges(3,4) * t127;
t102 = mrSges(3,3) * t127 - t130;
t101 = -mrSges(3,3) * t128 + t131;
t91 = Ifges(3,1) * t128 + t109 + t133;
t90 = t132 + (t116 * Ifges(3,2) + t136) * qJD(1);
t77 = Ifges(4,4) * t82;
t65 = -t93 * pkin(3) + t108;
t62 = qJD(2) * mrSges(4,1) + t83 * mrSges(4,3);
t61 = -qJD(2) * mrSges(4,2) + t82 * mrSges(4,3);
t60 = pkin(3) * t84 + t125;
t59 = pkin(2) * t128 - pkin(3) * t83;
t53 = pkin(3) * t73 + t118;
t45 = -mrSges(4,1) * t82 - mrSges(4,2) * t83;
t40 = -Ifges(4,1) * t83 + Ifges(4,5) * qJD(2) + t77;
t39 = Ifges(4,2) * t82 + Ifges(4,6) * qJD(2) - t140;
t33 = mrSges(5,1) * t110 - mrSges(5,3) * t44;
t32 = -mrSges(5,2) * t110 + mrSges(5,3) * t120;
t25 = -pkin(6) * t84 + t35;
t24 = -pkin(6) * t85 + t34;
t21 = -qJD(4) * t49 - t113 * t85 - t115 * t84;
t20 = qJD(4) * t48 - t113 * t84 + t115 * t85;
t19 = -mrSges(5,1) * t120 + mrSges(5,2) * t44;
t5 = -qJD(4) * t11 - t113 * t25 + t115 * t24;
t4 = qJD(4) * t10 + t113 * t24 + t115 * t25;
t1 = [(t147 * t82 - t73 * t93) * Ifges(4,2) + (t146 * t82 - t148 * t84 - t73 * t94 + t74 * t93) * Ifges(4,4) + (-t28 * t94 + t29 * t93 - t46 * t85 - t47 * t84 - t54 * t74 - t55 * t73) * mrSges(4,3) + m(5) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + t52 * t60 + t53 * t65) + m(4) * (t28 * t54 + t29 * t55 + t34 * t46 + t35 * t47) + t110 * (Ifges(5,5) * t20 + Ifges(5,6) * t21) / 0.2e1 + t60 * t19 + t35 * t61 + t34 * t62 + t52 * (-mrSges(5,1) * t21 + mrSges(5,2) * t20) + t53 * (-mrSges(5,1) * t48 + mrSges(5,2) * t49) + t4 * t32 + t5 * t33 + t20 * t13 / 0.2e1 + (-t90 / 0.2e1 - pkin(5) * t102 + t123 + (-0.2e1 * t145 - 0.3e1 / 0.2e1 * t136 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t116) * qJD(1) + (m(4) * (t100 + t129) + t45 + qJD(1) * (-mrSges(4,1) * t93 + mrSges(4,2) * t94)) * pkin(2)) * t126 + (Ifges(4,5) * t146 + Ifges(4,6) * t147 + (t91 / 0.2e1 - pkin(5) * t101 + t124 + (-0.2e1 * t144 + 0.3e1 / 0.2e1 * Ifges(3,4) * t116) * qJD(1)) * t116) * qJD(2) + (t148 * t85 + t74 * t94) * Ifges(4,1) + (t150 * t20 + t17 * t49) * Ifges(5,1) + (t152 * t21 + t48 * t18) * Ifges(5,2) + (t150 * t21 + t152 * t20 + t48 * t17 + t18 * t49) * Ifges(5,4) + t39 * t147 + t40 * t146 + t100 * (mrSges(4,1) * t84 + mrSges(4,2) * t85) + (-t10 * t17 + t11 * t18 + t2 * t48 - t6 * t20 + t7 * t21 - t3 * t49) * mrSges(5,3) + t65 * t121 + t108 * t122 + t21 * t160; -m(4) * (t46 * t50 + t47 * t51) + t44 * t160 + (t46 * t82 - t47 * t83 + (-t111 * t73 - t112 * t74) * pkin(2)) * mrSges(4,3) + t158 * t33 + (t158 * t6 + t159 * t7 + t2 * t79 + t3 * t78 - t52 * t59) * m(5) + t159 * t32 + (-t17 * t78 + t18 * t79 + t162) * mrSges(5,3) - t100 * (-t83 * mrSges(4,1) + t82 * mrSges(4,2)) - qJD(2) * (Ifges(4,5) * t82 + Ifges(4,6) * t83) / 0.2e1 - Ifges(4,6) * t73 + Ifges(4,5) * t74 - t51 * t61 - t50 * t62 - t59 * t19 + t28 * mrSges(4,1) - t29 * mrSges(4,2) + t163 + ((t124 - t109 / 0.2e1 - t91 / 0.2e1 + qJD(1) * t144 + (t101 - t131) * pkin(5)) * t116 + (t123 + t90 / 0.2e1 + (t145 + t136 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t116) * qJD(1) + (t102 + t130) * pkin(5) + (-m(4) * t100 - t45) * pkin(2)) * t114) * qJD(1) + t83 * (Ifges(4,1) * t82 + t140) / 0.2e1 + t39 * t148 + m(4) * (t111 * t29 + t112 * t28) * pkin(2) - (Ifges(4,2) * t83 + t40 + t77) * t82 / 0.2e1; -t120 * t32 + t44 * t33 - t82 * t61 - t83 * t62 + t121 + t122 + (-t120 * t7 + t44 * t6 + t53) * m(5) + (-t46 * t83 - t47 * t82 + t118) * m(4); t162 * mrSges(5,3) + t12 * t150 - t6 * t32 + t7 * t33 + t163;];
tauc = t1(:);
