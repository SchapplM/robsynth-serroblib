% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:43
% EndTime: 2019-12-05 16:08:54
% DurationCPUTime: 3.03s
% Computational Cost: add. (1199->261), mult. (3263->357), div. (0->0), fcn. (2019->6), ass. (0->134)
t169 = -mrSges(6,1) - mrSges(5,1);
t167 = Ifges(5,1) + Ifges(6,1);
t168 = Ifges(6,4) + Ifges(5,5);
t134 = mrSges(6,2) + mrSges(5,3);
t133 = -Ifges(5,4) + Ifges(6,5);
t165 = -Ifges(5,6) + Ifges(6,6);
t93 = cos(qJ(2));
t125 = qJD(1) * t93;
t129 = cos(pkin(8));
t132 = -qJ(4) - pkin(6);
t105 = qJD(3) * t132;
t92 = cos(qJ(3));
t119 = qJD(4) * t92;
t90 = sin(qJ(3));
t56 = t105 * t90 + t119;
t107 = t129 * t90;
t89 = sin(pkin(8));
t67 = t89 * t92 + t107;
t97 = -qJD(4) * t90 + t105 * t92;
t164 = -t67 * t125 - t129 * t97 + t56 * t89;
t58 = t67 * qJD(3);
t48 = qJD(2) * t58;
t138 = t89 * t90;
t99 = t129 * t92 - t138;
t60 = t99 * qJD(3);
t49 = qJD(2) * t60;
t18 = t48 * mrSges(6,1) - t49 * mrSges(6,3);
t19 = t48 * mrSges(5,1) + t49 * mrSges(5,2);
t163 = t18 + t19;
t104 = qJD(2) * t129;
t124 = qJD(2) * t90;
t57 = -t104 * t92 + t124 * t89;
t144 = Ifges(6,5) * t57;
t55 = Ifges(5,4) * t57;
t59 = t67 * qJD(2);
t162 = t168 * qJD(3) + t167 * t59 + t144 - t55;
t27 = mrSges(6,1) * t57 - mrSges(6,3) * t59;
t28 = mrSges(5,1) * t57 + mrSges(5,2) * t59;
t161 = -t27 - t28;
t142 = t57 * mrSges(5,3);
t149 = mrSges(6,2) * t57;
t41 = qJD(3) * mrSges(6,3) - t149;
t131 = -qJD(3) * mrSges(5,2) - t142 + t41;
t141 = t59 * mrSges(5,3);
t148 = mrSges(6,2) * t59;
t130 = -t169 * qJD(3) - t141 - t148;
t96 = t93 * t99;
t160 = -qJD(1) * t96 + t129 * t56 + t89 * t97;
t116 = qJD(2) * qJD(3);
t109 = t90 * t116;
t73 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t124;
t122 = qJD(2) * t92;
t74 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t122;
t159 = t92 * t73 + t90 * t74;
t87 = t90 ^ 2;
t88 = t92 ^ 2;
t94 = qJD(2) ^ 2;
t158 = -t57 / 0.2e1;
t157 = t57 / 0.2e1;
t155 = t59 / 0.2e1;
t153 = pkin(3) * t89;
t118 = qJ(4) * qJD(3);
t117 = qJD(1) * qJD(2);
t110 = t93 * t117;
t121 = qJD(3) * t90;
t91 = sin(qJ(2));
t126 = qJD(1) * t91;
t77 = qJD(2) * pkin(6) + t126;
t34 = t92 * t110 - t121 * t77;
t29 = (-t118 * t90 + t119) * qJD(2) + t34;
t113 = qJD(3) * t77 * t92;
t95 = -t113 + (-t92 * t118 + (-qJD(4) - t125) * t90) * qJD(2);
t2 = -t129 * t95 + t29 * t89;
t75 = t132 * t92;
t31 = -t107 * t132 - t75 * t89;
t152 = t2 * t31;
t50 = t67 * t91;
t151 = t2 * t50;
t150 = t2 * t67;
t147 = Ifges(4,4) * t90;
t146 = Ifges(5,4) * t59;
t143 = t49 * mrSges(6,2);
t78 = -qJD(2) * pkin(2) - t125;
t140 = t78 * t91;
t103 = qJ(4) * qJD(2) + t77;
t53 = t103 * t92;
t139 = t89 * t53;
t135 = -qJD(3) / 0.2e1;
t3 = t129 * t29 + t89 * t95;
t36 = t129 * t53;
t52 = t103 * t90;
t42 = qJD(3) * pkin(3) - t52;
t10 = t89 * t42 + t36;
t84 = t91 * t117;
t65 = pkin(3) * t109 + t84;
t128 = Ifges(4,5) * qJD(3);
t127 = Ifges(4,6) * qJD(3);
t123 = qJD(2) * t91;
t120 = qJD(3) * t91;
t115 = pkin(3) * t124;
t114 = pkin(3) * t121;
t85 = -pkin(3) * t92 - pkin(2);
t112 = (t87 + t88) * t77;
t111 = t129 * pkin(3);
t35 = -t110 * t90 - t113;
t102 = t34 * t92 - t35 * t90;
t100 = t90 * t73 - t92 * t74;
t9 = t129 * t42 - t139;
t98 = (mrSges(4,1) * t90 + mrSges(4,2) * t92) * qJD(2);
t61 = qJD(2) * t85 + qJD(4) - t125;
t86 = Ifges(4,4) * t122;
t83 = -t111 - pkin(4);
t81 = qJ(5) + t153;
t64 = qJD(3) * t98;
t63 = Ifges(4,1) * t124 + t128 + t86;
t62 = t127 + (t92 * Ifges(4,2) + t147) * qJD(2);
t54 = Ifges(6,5) * t59;
t51 = t99 * t91;
t32 = -t129 * t75 + t132 * t138;
t30 = -pkin(4) * t99 - qJ(5) * t67 + t85;
t24 = -t57 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t146;
t23 = Ifges(6,6) * qJD(3) + t57 * Ifges(6,3) + t54;
t15 = pkin(4) * t59 + qJ(5) * t57 + t115;
t14 = qJD(2) * t96 - t120 * t67;
t13 = (-t104 * t90 - t122 * t89) * t93 - t99 * t120;
t12 = -t129 * t52 - t139;
t11 = -t52 * t89 + t36;
t8 = qJD(3) * qJ(5) + t10;
t7 = pkin(4) * t57 - qJ(5) * t59 + t61;
t6 = -qJD(3) * pkin(4) + qJD(5) - t9;
t5 = pkin(4) * t58 - qJ(5) * t60 - qJD(5) * t67 + t114;
t4 = pkin(4) * t48 - qJ(5) * t49 - qJD(5) * t59 + t65;
t1 = qJD(3) * qJD(5) + t3;
t16 = [t131 * t14 + t130 * t13 + (-t94 * mrSges(3,2) - t100 * qJD(2) - t163 - t64) * t93 + (-t94 * mrSges(3,1) - t159 * qJD(3) + (qJD(2) * (-mrSges(4,1) * t92 + mrSges(4,2) * t90) - t161) * qJD(2)) * t91 + m(4) * (t102 * t91 + (t140 + (t112 - t126) * t93) * qJD(2)) + m(5) * (t10 * t14 + t123 * t61 + t13 * t9 + t3 * t51 - t65 * t93 + t151) + m(6) * (t1 * t51 + t123 * t7 - t13 * t6 + t14 * t8 - t4 * t93 + t151) + t134 * (-t48 * t51 + t49 * t50); (t100 * t93 + t161 * t91) * qJD(1) + t102 * mrSges(4,3) + t85 * t19 + t65 * (-mrSges(5,1) * t99 + mrSges(5,2) * t67) - pkin(2) * t64 + t4 * (-mrSges(6,1) * t99 - mrSges(6,3) * t67) + t5 * t27 + t30 * t18 + (-0.3e1 / 0.2e1 * t87 + 0.3e1 / 0.2e1 * t88) * Ifges(4,4) * t116 + ((t78 * mrSges(4,2) + t63 / 0.2e1 - pkin(6) * t73 + t128 / 0.2e1) * t92 + (t78 * mrSges(4,1) - t62 / 0.2e1 - pkin(6) * t74 + pkin(3) * t28 - t127 / 0.2e1 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t122) * t90) * qJD(3) + (-t133 * t99 + t134 * t31 + t167 * t67) * t49 + (t133 * t67 - (Ifges(5,2) + Ifges(6,3)) * t99 - t134 * t32) * t48 + (t3 * t99 + t150) * mrSges(5,3) + (t1 * t99 + t150) * mrSges(6,2) + t160 * t131 - t164 * t130 + (t1 * t32 + t30 * t4 + t152 + t160 * t8 + (-t126 + t5) * t7 + t164 * t6) * m(6) + (t3 * t32 + t65 * t85 + t152 - t164 * t9 + (t114 - t126) * t61 + t160 * t10) * m(5) + (-(t112 * t93 + t140) * qJD(1) - pkin(2) * t84 + t102 * pkin(6)) * m(4) + (t7 * mrSges(6,1) + t61 * mrSges(5,1) - t24 / 0.2e1 + t23 / 0.2e1 + Ifges(6,3) * t157 - Ifges(5,2) * t158 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * qJD(3) - t10 * mrSges(5,3) - t8 * mrSges(6,2) + t133 * t155) * t58 + (-t7 * mrSges(6,3) + t61 * mrSges(5,2) + Ifges(6,5) * t157 + Ifges(5,4) * t158 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * qJD(3) - t9 * mrSges(5,3) + t6 * mrSges(6,2) + t167 * t155 + t162 / 0.2e1) * t60; (-mrSges(5,3) * t111 + t168) * t49 + (t165 * t59 - t168 * t57) * t135 + t169 * t2 - (-Ifges(4,2) * t124 + t63 + t86) * t122 / 0.2e1 - t78 * t98 + t116 * Ifges(4,5) * t92 / 0.2e1 - t61 * (mrSges(5,1) * t59 - mrSges(5,2) * t57) - t7 * (mrSges(6,1) * t59 + mrSges(6,3) * t57) - t34 * mrSges(4,2) + t35 * mrSges(4,1) + qJD(5) * t41 - t15 * t27 + t1 * mrSges(6,3) - t3 * mrSges(5,2) - t90 * t94 * (Ifges(4,1) * t92 - t147) / 0.2e1 + ((-t129 * t2 + t3 * t89) * pkin(3) - t10 * t12 + t11 * t9 - t115 * t61) * m(5) + (t1 * t81 - t11 * t6 - t15 * t7 + t2 * t83 + (-t12 + qJD(5)) * t8) * m(6) - Ifges(4,6) * t109 / 0.2e1 + t10 * t141 + t83 * t143 + t8 * t148 + t6 * t149 + t24 * t155 + (Ifges(6,3) * t59 - t144) * t158 - t28 * t115 + t62 * t124 / 0.2e1 - t9 * t142 + t159 * t77 + t130 * t11 - t131 * t12 + (-Ifges(5,2) * t59 + t162 - t55) * t157 + (-t81 * mrSges(6,2) - mrSges(5,3) * t153 + t165) * t48 - (-t167 * t57 - t146 + t23 + t54) * t59 / 0.2e1; t130 * t59 + t131 * t57 + (t57 * t8 - t59 * t6 + t4) * m(6) + (t10 * t57 + t9 * t59 + t65) * m(5) + t163; t143 - qJD(3) * t41 + t59 * t27 + 0.2e1 * (t2 / 0.2e1 + t8 * t135 + t7 * t155) * m(6);];
tauc = t16(:);
