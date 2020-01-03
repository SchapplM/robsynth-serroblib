% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:10
% EndTime: 2019-12-31 17:11:14
% DurationCPUTime: 1.75s
% Computational Cost: add. (1009->239), mult. (2525->332), div. (0->0), fcn. (1216->4), ass. (0->119)
t142 = qJD(1) / 0.2e1;
t148 = -mrSges(3,1) + mrSges(4,2);
t73 = cos(qJ(4));
t107 = qJD(2) * t73;
t74 = cos(qJ(2));
t109 = qJD(1) * t74;
t71 = sin(qJ(4));
t44 = -t71 * t109 + t107;
t122 = Ifges(5,4) * t44;
t43 = -qJD(2) * t71 - t73 * t109;
t72 = sin(qJ(2));
t110 = qJD(1) * t72;
t64 = qJD(4) + t110;
t11 = Ifges(5,2) * t43 + Ifges(5,6) * t64 + t122;
t115 = Ifges(5,6) * t73;
t119 = Ifges(5,5) * t71;
t42 = Ifges(5,4) * t43;
t12 = Ifges(5,1) * t44 + Ifges(5,5) * t64 + t42;
t128 = -t73 / 0.2e1;
t129 = -t71 / 0.2e1;
t130 = -t64 / 0.2e1;
t131 = -t44 / 0.2e1;
t132 = -t43 / 0.2e1;
t103 = qJD(2) * qJ(3);
t68 = pkin(5) * t109;
t55 = -t68 - t103;
t69 = pkin(3) * t109;
t38 = -t55 + t69;
t121 = Ifges(5,4) * t71;
t86 = Ifges(5,2) * t73 + t121;
t120 = Ifges(5,4) * t73;
t88 = Ifges(5,1) * t71 + t120;
t91 = mrSges(5,1) * t73 - mrSges(5,2) * t71;
t75 = -pkin(2) - pkin(6);
t95 = -qJ(3) * t72 - pkin(1);
t41 = t75 * t74 + t95;
t29 = t41 * qJD(1);
t66 = pkin(5) * t110;
t47 = -pkin(3) * t110 - t66;
t30 = t75 * qJD(2) + qJD(3) - t47;
t8 = -t29 * t71 + t30 * t73;
t9 = t29 * t73 + t30 * t71;
t92 = t8 * t71 - t9 * t73;
t147 = t92 * mrSges(5,3) + t11 * t128 + t12 * t129 + (t115 + t119) * t130 + t88 * t131 + t86 * t132 + t38 * t91;
t117 = Ifges(4,6) * t74;
t140 = qJD(2) / 0.2e1;
t141 = -qJD(2) / 0.2e1;
t143 = -Ifges(3,4) * t109 / 0.2e1;
t144 = -Ifges(3,1) / 0.2e1;
t53 = -pkin(2) * t74 + t95;
t39 = t53 * qJD(1);
t52 = -qJD(2) * pkin(2) + qJD(3) + t66;
t146 = (m(4) * t52 + (mrSges(4,1) + mrSges(3,3)) * t110 + t148 * qJD(2)) * pkin(5) - t39 * mrSges(4,3) - t9 * mrSges(5,2) + t64 * Ifges(5,3) + t44 * Ifges(5,5) + t43 * Ifges(5,6) - t110 * t144 - Ifges(3,5) * t141 - t143 - Ifges(4,4) * t140 - (-t72 * Ifges(4,2) - t117) * t142 + t52 * mrSges(4,1) + t8 * mrSges(5,1);
t57 = -mrSges(4,1) * t109 - qJD(2) * mrSges(4,3);
t145 = (m(4) * t55 + qJD(2) * mrSges(3,2) - mrSges(3,3) * t109 + t57) * pkin(5) + t55 * mrSges(4,1) - Ifges(3,6) * t140 - Ifges(4,5) * t141 - t39 * mrSges(4,2) - t147 + ((-Ifges(3,2) - Ifges(4,3)) * t74 + (-Ifges(3,4) - Ifges(4,6)) * t72) * t142;
t127 = pkin(3) + pkin(5);
t139 = t127 * t72;
t67 = pkin(2) * t110;
t63 = qJD(2) * t67;
t105 = qJD(3) * t72;
t84 = pkin(6) * t72 - qJ(3) * t74;
t78 = t84 * qJD(2) - t105;
t20 = t78 * qJD(1) + t63;
t60 = t127 * t74;
t50 = qJD(2) * t60;
t40 = qJD(1) * t50;
t1 = t8 * qJD(4) + t20 * t73 + t40 * t71;
t2 = -t9 * qJD(4) - t20 * t71 + t40 * t73;
t138 = t1 * t71 + t2 * t73;
t102 = qJD(2) * qJD(4);
t104 = qJD(4) * t74;
t108 = qJD(2) * t72;
t23 = -t71 * t102 + (-t73 * t104 + t71 * t108) * qJD(1);
t24 = -t73 * t102 + (t71 * t104 + t72 * t107) * qJD(1);
t137 = t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t23 + Ifges(5,6) * t24;
t133 = t11 / 0.2e1;
t126 = pkin(1) * mrSges(3,1);
t125 = pkin(1) * mrSges(3,2);
t118 = Ifges(5,5) * t73;
t116 = Ifges(5,6) * t71;
t19 = -mrSges(5,1) * t43 + mrSges(5,2) * t44;
t111 = -t57 + t19;
t106 = qJD(2) * t74;
t100 = -0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(4,6);
t99 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t98 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t97 = m(4) * pkin(5) + mrSges(4,1);
t96 = qJD(1) * t106;
t90 = mrSges(5,1) * t71 + mrSges(5,2) * t73;
t89 = Ifges(5,1) * t73 - t121;
t87 = -Ifges(5,2) * t71 + t120;
t15 = mrSges(5,1) * t96 - mrSges(5,3) * t23;
t16 = -mrSges(5,2) * t96 + mrSges(5,3) * t24;
t83 = t73 * t15 + t71 * t16;
t25 = -mrSges(5,2) * t64 + mrSges(5,3) * t43;
t26 = mrSges(5,1) * t64 - mrSges(5,3) * t44;
t82 = t73 * t25 - t71 * t26;
t18 = t139 * t71 + t41 * t73;
t17 = t139 * t73 - t41 * t71;
t79 = -t74 * t103 - t105;
t70 = pkin(2) * t108;
t62 = Ifges(5,3) * t96;
t51 = (-qJD(3) + t66) * qJD(2);
t49 = qJD(2) * t139;
t48 = t68 + t69;
t45 = (t74 * mrSges(4,2) - t72 * mrSges(4,3)) * qJD(1);
t33 = t70 + t79;
t32 = t84 * qJD(1) + t67;
t31 = (-qJD(1) * t139 + qJD(3)) * qJD(2);
t28 = t79 * qJD(1) + t63;
t27 = t70 + t78;
t14 = t32 * t73 + t48 * t71;
t13 = -t32 * t71 + t48 * t73;
t7 = -mrSges(5,1) * t24 + mrSges(5,2) * t23;
t6 = t23 * Ifges(5,1) + t24 * Ifges(5,4) + Ifges(5,5) * t96;
t5 = t23 * Ifges(5,4) + t24 * Ifges(5,2) + Ifges(5,6) * t96;
t4 = -t18 * qJD(4) - t27 * t71 + t50 * t73;
t3 = t17 * qJD(4) + t27 * t73 + t50 * t71;
t10 = [t17 * t15 + t18 * t16 - t49 * t19 + t3 * t25 + t4 * t26 + t33 * t45 + t60 * t7 + m(4) * (t28 * t53 + t33 * t39) + m(5) * (t1 * t18 + t17 * t2 + t3 * t9 + t31 * t60 - t38 * t49 + t4 * t8) + (-t28 * mrSges(4,3) + t62 / 0.2e1 + (t98 * qJD(2) + (-t53 * mrSges(4,2) + t100 * t72 - 0.2e1 * t126) * qJD(1) + t145) * qJD(2) + t137) * t72 + (t31 * t91 - t23 * t88 / 0.2e1 - t24 * t86 / 0.2e1 + t5 * t128 + t6 * t129 + t28 * mrSges(4,2) - t97 * t51 + (-t1 * t73 + t2 * t71) * mrSges(5,3) + (t64 * (t116 - t118) / 0.2e1 + t89 * t131 - t38 * t90 + t87 * t132 + t12 * t128 + t71 * t133 + (t71 * t9 + t73 * t8) * mrSges(5,3)) * qJD(4) + (t99 * qJD(2) + ((-t119 / 0.2e1 - t115 / 0.2e1 - t100) * t74 - 0.2e1 * t125 - t53 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + t97 * pkin(5)) * t72) * qJD(1) + t146) * qJD(2)) * t74; t31 * t90 + t23 * t89 / 0.2e1 + t24 * t87 / 0.2e1 + t73 * t6 / 0.2e1 + t5 * t129 - t51 * mrSges(4,3) - t47 * t19 - t14 * t25 - t13 * t26 + qJ(3) * t7 + t83 * t75 + t111 * qJD(3) - t138 * mrSges(5,3) + m(4) * (-qJ(3) * t51 - qJD(3) * t55) + m(5) * (t31 * qJ(3) + t38 * qJD(3) + t138 * t75) - m(5) * (t13 * t8 + t14 * t9 + t38 * t47) + ((-m(5) * t92 + t82) * t75 + t147) * qJD(4) + (((Ifges(4,3) / 0.2e1 + t144 - Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t109 + (-qJ(3) * mrSges(4,1) + pkin(5) * mrSges(3,2) + t98) * qJD(2) + (t126 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t72) * qJD(1) - t145) * t72 + (t143 + (t125 - t117 / 0.2e1) * qJD(1) + (t118 / 0.2e1 - t116 / 0.2e1 - pkin(2) * mrSges(4,1) + (-m(4) * pkin(2) + t148) * pkin(5) + t99) * qJD(2) - t146) * t74) * qJD(1) + (-m(4) * t39 - t45) * (-qJ(3) * t109 + t67); t82 * qJD(4) - t111 * qJD(2) + (t97 * t106 + (t45 + t82) * t72) * qJD(1) - m(4) * (-qJD(2) * t55 - t39 * t110) + t83 + (-qJD(2) * t38 - t64 * t92 + t138) * m(5); t62 - t38 * (mrSges(5,1) * t44 + mrSges(5,2) * t43) + (Ifges(5,1) * t43 - t122) * t131 + t44 * t133 + (Ifges(5,5) * t43 - Ifges(5,6) * t44) * t130 - t8 * t25 + t9 * t26 + (t43 * t8 + t44 * t9) * mrSges(5,3) + (-Ifges(5,2) * t44 + t12 + t42) * t132 + t137;];
tauc = t10(:);
