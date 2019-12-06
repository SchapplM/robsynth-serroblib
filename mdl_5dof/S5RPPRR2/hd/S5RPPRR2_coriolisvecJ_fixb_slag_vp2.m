% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:34
% EndTime: 2019-12-05 17:39:40
% DurationCPUTime: 1.97s
% Computational Cost: add. (2377->231), mult. (5353->327), div. (0->0), fcn. (3676->6), ass. (0->110)
t101 = sin(pkin(8));
t105 = sin(qJ(4));
t102 = cos(pkin(8));
t107 = cos(qJ(4));
t128 = t102 * t107;
t112 = t101 * t105 - t128;
t104 = sin(qJ(5));
t106 = cos(qJ(5));
t80 = -t107 * t101 - t105 * t102;
t51 = -t104 * t112 - t106 * t80;
t125 = qJD(4) * t107;
t126 = qJD(4) * t105;
t75 = -t101 * t126 + t102 * t125;
t76 = t80 * qJD(4);
t109 = -qJD(5) * t51 - t104 * t75 + t106 * t76;
t152 = t104 * t80 - t106 * t112;
t111 = t80 * qJD(3);
t103 = -pkin(1) - qJ(3);
t87 = qJD(1) * t103 + qJD(2);
t120 = -pkin(6) * qJD(1) + t87;
t68 = t120 * t101;
t69 = t120 * t102;
t28 = qJD(1) * t111 + t69 * t125 - t126 * t68;
t67 = t112 * qJD(4) * qJD(1);
t19 = pkin(7) * t67 + t28;
t162 = t112 * qJD(3);
t40 = t105 * t69 + t107 * t68;
t29 = t162 * qJD(1) - qJD(4) * t40;
t66 = qJD(1) * t76;
t20 = -pkin(7) * t66 + t29;
t73 = t80 * qJD(1);
t32 = pkin(7) * t73 + t40;
t133 = t104 * t32;
t39 = -t105 * t68 + t107 * t69;
t127 = qJD(1) * t101;
t74 = qJD(1) * t128 - t105 * t127;
t31 = -pkin(7) * t74 + t39;
t30 = qJD(4) * pkin(4) + t31;
t6 = t106 * t30 - t133;
t2 = qJD(5) * t6 + t104 * t20 + t106 * t19;
t22 = qJD(5) * t152 + t104 * t76 + t106 * t75;
t130 = t106 * t32;
t7 = t104 * t30 + t130;
t3 = -qJD(5) * t7 - t104 * t19 + t106 * t20;
t165 = t109 * t6 + t152 * t3 + t2 * t51 + t22 * t7;
t119 = -t104 * t74 + t106 * t73;
t41 = Ifges(6,4) * t119;
t47 = t104 * t73 + t106 * t74;
t99 = qJD(4) + qJD(5);
t13 = Ifges(6,1) * t47 + Ifges(6,5) * t99 + t41;
t139 = Ifges(6,4) * t47;
t17 = qJD(5) * t119 + t104 * t67 + t106 * t66;
t18 = -qJD(5) * t47 - t104 * t66 + t106 * t67;
t100 = qJD(1) * qJ(2);
t95 = qJD(3) + t100;
t85 = pkin(3) * t127 + t95;
t55 = -pkin(4) * t73 + t85;
t164 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t18 - (Ifges(6,5) * t119 - Ifges(6,6) * t47) * t99 / 0.2e1 + (t119 * t6 + t47 * t7) * mrSges(6,3) - (-Ifges(6,2) * t47 + t13 + t41) * t119 / 0.2e1 - t55 * (mrSges(6,1) * t47 + mrSges(6,2) * t119) - (Ifges(6,1) * t119 - t139) * t47 / 0.2e1;
t163 = t18 * t51;
t12 = Ifges(6,2) * t119 + Ifges(6,6) * t99 + t139;
t159 = t12 / 0.2e1;
t154 = t152 * t17;
t21 = -mrSges(6,1) * t119 + mrSges(6,2) * t47;
t153 = -m(6) * t55 - t21;
t134 = t101 ^ 2 + t102 ^ 2;
t117 = qJD(1) * t134;
t149 = t119 / 0.2e1;
t147 = t47 / 0.2e1;
t145 = -t75 / 0.2e1;
t144 = t76 / 0.2e1;
t141 = m(3) * qJ(2);
t140 = Ifges(5,4) * t74;
t138 = t66 * t112;
t137 = t80 * t67;
t136 = -pkin(6) + t103;
t83 = t136 * t101;
t84 = t136 * t102;
t54 = t105 * t84 + t107 * t83;
t115 = mrSges(4,1) * t101 + mrSges(4,2) * t102;
t135 = -mrSges(5,1) * t73 + mrSges(5,2) * t74 + qJD(1) * t115;
t91 = t101 * pkin(3) + qJ(2);
t124 = qJD(1) * qJD(2);
t122 = -t67 * mrSges(5,1) + t66 * mrSges(5,2);
t121 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t53 = -t105 * t83 + t107 * t84;
t37 = pkin(7) * t112 + t53;
t38 = pkin(7) * t80 + t54;
t10 = -t104 * t38 + t106 * t37;
t11 = t104 * t37 + t106 * t38;
t110 = -t112 * t29 - t28 * t80 + t39 * t76 + t40 * t75;
t34 = -qJD(4) * t54 + t162;
t33 = t84 * t125 - t126 * t83 + t111;
t108 = qJD(1) ^ 2;
t70 = Ifges(5,4) * t73;
t63 = pkin(4) * t75 + qJD(2);
t61 = -pkin(4) * t80 + t91;
t60 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t74;
t59 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t73;
t56 = -pkin(4) * t67 + t124;
t43 = Ifges(5,1) * t74 + Ifges(5,5) * qJD(4) + t70;
t42 = Ifges(5,2) * t73 + Ifges(5,6) * qJD(4) + t140;
t36 = mrSges(6,1) * t99 - mrSges(6,3) * t47;
t35 = -mrSges(6,2) * t99 + mrSges(6,3) * t119;
t27 = -pkin(7) * t76 + t34;
t26 = -pkin(7) * t75 + t33;
t9 = t106 * t31 - t133;
t8 = -t104 * t31 - t130;
t5 = -qJD(5) * t11 - t104 * t26 + t106 * t27;
t4 = qJD(5) * t10 + t104 * t27 + t106 * t26;
t1 = [(t144 * t74 - t138) * Ifges(5,1) + t109 * t13 / 0.2e1 + m(5) * (t28 * t54 + t29 * t53 + t33 * t40 + t34 * t39) + t99 * (Ifges(6,5) * t109 - Ifges(6,6) * t22) / 0.2e1 + t55 * (mrSges(6,1) * t22 + mrSges(6,2) * t109) + (t109 * t149 - t147 * t22 + t152 * t18 - t17 * t51) * Ifges(6,4) + (-t149 * t22 - t163) * Ifges(6,2) - t22 * t159 + ((-mrSges(5,1) * t80 - mrSges(5,2) * t112 + (2 * mrSges(3,3)) + t115 + 0.2e1 * t141) * qJD(1) + t135 + m(4) * (t95 + t100) + m(5) * (qJD(1) * t91 + t85)) * qJD(2) + (-t112 * t67 + t144 * t73 + t145 * t74 + t80 * t66) * Ifges(5,4) + t56 * (mrSges(6,1) * t51 + mrSges(6,2) * t152) + t85 * (mrSges(5,1) * t75 + mrSges(5,2) * t76) + qJD(4) * (Ifges(5,5) * t76 - Ifges(5,6) * t75) / 0.2e1 + t33 * t59 + t34 * t60 + t63 * t21 + t4 * t35 + t5 * t36 + (-t10 * t17 + t11 * t18 - t165) * mrSges(6,3) + m(6) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + t55 * t63 + t56 * t61) + (t109 * t147 + t154) * Ifges(6,1) + (0.2e1 * mrSges(4,3) * t117 + m(4) * (-t103 * t117 - t134 * t87)) * qJD(3) + (t145 * t73 + t137) * Ifges(5,2) + t43 * t144 + t42 * t145 + t61 * t121 + t91 * t122 + (-t53 * t66 + t54 * t67 - t110) * mrSges(5,3); t22 * t35 + t109 * t36 + t75 * t59 + t76 * t60 + (-mrSges(3,3) - t141) * t108 + (-t154 + t163) * mrSges(6,3) + (-t137 + t138) * mrSges(5,3) + m(5) * t110 + m(6) * t165 + (-m(5) * t85 + (-qJD(3) * t134 - t95) * m(4) - t135 + t153) * qJD(1); -t119 * t35 + t47 * t36 - t73 * t59 + t74 * t60 + (m(4) + m(5)) * t124 - t134 * t108 * mrSges(4,3) - m(5) * (-t39 * t74 + t40 * t73) + m(4) * t87 * t117 + t121 + t122 + (-t119 * t7 + t47 * t6 + t56) * m(6); (t39 * t73 + t40 * t74) * mrSges(5,3) - (-Ifges(5,2) * t74 + t43 + t70) * t73 / 0.2e1 + t47 * t159 - m(6) * (t6 * t8 + t7 * t9) - t85 * (mrSges(5,1) * t74 + mrSges(5,2) * t73) - qJD(4) * (Ifges(5,5) * t73 - Ifges(5,6) * t74) / 0.2e1 + t74 * t42 / 0.2e1 + Ifges(5,5) * t66 + Ifges(5,6) * t67 - t39 * t59 + t40 * t60 - t9 * t35 - t8 * t36 - t28 * mrSges(5,2) + t29 * mrSges(5,1) + ((t104 * t18 - t106 * t17) * mrSges(6,3) + m(6) * (t104 * t2 + t106 * t3) + t153 * t74 + (-t104 * t36 + t106 * t35 + m(6) * (-t104 * t6 + t106 * t7)) * qJD(5)) * pkin(4) - t74 * (Ifges(5,1) * t73 - t140) / 0.2e1 + t164; t12 * t147 - t6 * t35 + t7 * t36 + t164;];
tauc = t1(:);
