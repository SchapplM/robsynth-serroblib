% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:18
% DurationCPUTime: 3.13s
% Computational Cost: add. (2330->288), mult. (6293->378), div. (0->0), fcn. (4396->6), ass. (0->132)
t103 = -qJD(3) + qJD(5);
t178 = t103 / 0.2e1;
t109 = sin(qJ(3));
t146 = cos(qJ(3));
t106 = sin(pkin(8));
t134 = pkin(6) + qJ(2);
t95 = t134 * t106;
t90 = qJD(1) * t95;
t107 = cos(pkin(8));
t96 = t134 * t107;
t91 = qJD(1) * t96;
t59 = -t109 * t91 - t146 * t90;
t158 = -t59 + qJD(4);
t165 = Ifges(4,1) + Ifges(5,1);
t164 = Ifges(5,4) + Ifges(4,5);
t108 = sin(qJ(5));
t110 = cos(qJ(5));
t125 = t146 * t107;
t115 = -t109 * t106 + t125;
t80 = t115 * qJD(1);
t89 = t106 * t146 + t109 * t107;
t81 = t89 * qJD(1);
t119 = -t108 * t80 + t110 * t81;
t50 = t108 * t81 + t110 * t80;
t43 = Ifges(6,4) * t50;
t177 = Ifges(6,2) * t119 + t43;
t176 = -pkin(7) * t81 + t158;
t124 = qJD(3) * t146;
t127 = qJD(1) * qJD(2);
t98 = qJD(2) * t125;
t28 = -t90 * t124 + qJD(1) * t98 + (-qJD(3) * t91 - t106 * t127) * t109;
t27 = qJD(3) * qJD(4) + t28;
t83 = t89 * qJD(3);
t70 = qJD(1) * t83;
t17 = pkin(7) * t70 + t27;
t114 = t89 * qJD(2);
t60 = -t109 * t90 + t146 * t91;
t29 = qJD(1) * t114 + qJD(3) * t60;
t82 = t115 * qJD(3);
t69 = qJD(1) * t82;
t20 = -t69 * pkin(7) + t29;
t111 = -pkin(3) - pkin(4);
t26 = qJD(3) * t111 + t176;
t105 = qJD(3) * qJ(4);
t35 = -pkin(7) * t80 + t60;
t30 = t105 + t35;
t5 = -t108 * t30 + t110 * t26;
t1 = qJD(5) * t5 + t108 * t20 + t110 * t17;
t13 = -qJD(5) * t50 + t108 * t70 + t110 * t69;
t14 = -qJD(5) * t119 - t108 * t69 + t110 * t70;
t6 = t108 * t26 + t110 * t30;
t2 = -qJD(5) * t6 - t108 * t17 + t110 * t20;
t126 = pkin(2) * t107 + pkin(1);
t92 = -qJD(1) * t126 + qJD(2);
t173 = -qJ(4) * t81 + t92;
t22 = -t111 * t80 - t173;
t175 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t14 + (Ifges(6,5) * t50 + Ifges(6,6) * t119) * t178 + t22 * (-mrSges(6,1) * t119 + mrSges(6,2) * t50);
t155 = -t50 / 0.2e1;
t174 = -mrSges(5,1) - mrSges(4,1);
t172 = t119 * t6 - t5 * t50;
t143 = Ifges(6,4) * t119;
t170 = -Ifges(6,1) * t50 - t143;
t154 = -t119 / 0.2e1;
t93 = -qJ(4) * t108 + t110 * t111;
t168 = qJD(5) * t93 - t108 * t35 + t176 * t110;
t94 = qJ(4) * t110 + t108 * t111;
t167 = -qJD(5) * t94 - t176 * t108 - t110 * t35;
t135 = -Ifges(4,4) + Ifges(5,5);
t163 = Ifges(5,6) - Ifges(4,6);
t142 = Ifges(5,5) * t80;
t74 = Ifges(4,4) * t80;
t162 = t164 * qJD(3) + t165 * t81 - t142 + t74;
t157 = (m(3) * qJ(2) + mrSges(3,3)) * (t106 ^ 2 + t107 ^ 2);
t156 = t50 / 0.2e1;
t153 = t119 / 0.2e1;
t152 = t80 / 0.2e1;
t151 = -t80 / 0.2e1;
t149 = t81 / 0.2e1;
t145 = mrSges(4,3) * t80;
t144 = Ifges(4,4) * t81;
t61 = t109 * t96 + t146 * t95;
t141 = t29 * t61;
t140 = t29 * t89;
t139 = t81 * mrSges(4,3);
t136 = mrSges(5,2) + mrSges(4,3);
t66 = t80 * mrSges(5,2) + qJD(3) * mrSges(5,3);
t133 = -qJD(3) * mrSges(4,2) + t145 + t66;
t132 = -mrSges(5,2) * t81 - t174 * qJD(3) - t139;
t62 = -t109 * t95 + t146 * t96;
t131 = qJ(4) * t80;
t122 = -t14 * mrSges(6,1) + t13 * mrSges(6,2);
t36 = -mrSges(6,2) * t103 - mrSges(6,3) * t50;
t37 = mrSges(6,1) * t103 - mrSges(6,3) * t119;
t120 = -t108 * t37 + t110 * t36;
t41 = -pkin(7) * t89 + t61;
t42 = -pkin(7) * t115 + t62;
t9 = -t108 * t42 + t110 * t41;
t10 = t108 * t41 + t110 * t42;
t57 = -t108 * t89 - t110 * t115;
t58 = -t108 * t115 + t110 * t89;
t118 = qJ(4) * t69 + qJD(4) * t81;
t117 = qJ(4) * t82 + qJD(4) * t89;
t116 = qJ(4) * t89 + t126;
t39 = -t95 * t124 + t98 + (-qJD(2) * t106 - qJD(3) * t96) * t109;
t40 = qJD(3) * t62 + t114;
t73 = Ifges(5,5) * t81;
t68 = t70 * mrSges(5,1);
t67 = t69 * mrSges(4,2);
t56 = -pkin(3) * t115 - t116;
t55 = -mrSges(5,1) * t80 - mrSges(5,3) * t81;
t54 = pkin(3) * t81 - t131;
t53 = t105 + t60;
t52 = -qJD(3) * pkin(3) + t158;
t45 = t80 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t144;
t44 = Ifges(5,6) * qJD(3) - t80 * Ifges(5,3) + t73;
t38 = -pkin(3) * t80 + t173;
t33 = -t111 * t115 + t116;
t32 = pkin(3) * t83 - t117;
t31 = t111 * t81 + t131;
t25 = pkin(3) * t70 - t118;
t24 = -t82 * pkin(7) + t40;
t23 = pkin(7) * t83 + t39;
t21 = t111 * t83 + t117;
t19 = -qJD(5) * t58 - t108 * t82 + t110 * t83;
t18 = qJD(5) * t57 + t108 * t83 + t110 * t82;
t16 = t111 * t70 + t118;
t15 = mrSges(6,1) * t50 + mrSges(6,2) * t119;
t12 = Ifges(6,1) * t119 + Ifges(6,5) * t103 - t43;
t11 = -Ifges(6,2) * t50 + Ifges(6,6) * t103 + t143;
t4 = -qJD(5) * t10 - t108 * t23 + t110 * t24;
t3 = qJD(5) * t9 + t108 * t24 + t110 * t23;
t7 = [0.2e1 * t157 * t127 + (-mrSges(5,3) * t56 - t115 * t135 + t136 * t61 + t165 * t89) * t69 + (t115 * t28 + t140) * mrSges(4,3) + (t115 * t27 + t140) * mrSges(5,2) + t25 * (-mrSges(5,1) * t115 - mrSges(5,3) * t89) + (-t126 * mrSges(4,1) + t135 * t89 - (Ifges(4,2) + Ifges(5,3)) * t115 - t136 * t62) * t70 + (t92 * mrSges(4,1) + t38 * mrSges(5,1) + t44 / 0.2e1 - t45 / 0.2e1 - t60 * mrSges(4,3) - t53 * mrSges(5,2) + Ifges(5,3) * t151 - Ifges(4,2) * t152) * t83 + (t92 * mrSges(4,2) + t52 * mrSges(5,2) - t59 * mrSges(4,3) - t38 * mrSges(5,3) + Ifges(4,4) * t152 + Ifges(5,5) * t151) * t82 + (t135 * t83 + t165 * t82) * t149 + t162 * t82 / 0.2e1 + (t163 * t83 + t164 * t82) * qJD(3) / 0.2e1 + m(6) * (t1 * t10 + t16 * t33 + t2 * t9 + t21 * t22 + t3 * t6 + t4 * t5) + (t1 * t57 + t10 * t14 - t13 * t9 - t18 * t5 + t19 * t6 - t2 * t58) * mrSges(6,3) + t56 * t68 + t32 * t55 + t16 * (-mrSges(6,1) * t57 + mrSges(6,2) * t58) + t3 * t36 + t4 * t37 + t21 * t15 + t22 * (-mrSges(6,1) * t19 + mrSges(6,2) * t18) + t18 * t12 / 0.2e1 + t19 * t11 / 0.2e1 + (Ifges(6,5) * t18 + Ifges(6,6) * t19) * t178 + m(5) * (t25 * t56 + t27 * t62 + t32 * t38 + t39 * t53 + t40 * t52 + t141) + m(4) * (t28 * t62 + t39 * t60 - t40 * t59 + t141) - t132 * t40 + t133 * t39 - t126 * t67 + t33 * t122 + (t13 * t58 + t153 * t18) * Ifges(6,1) + (t57 * t14 + t155 * t19) * Ifges(6,2) + (t57 * t13 + t14 * t58 + t153 * t19 + t155 * t18) * Ifges(6,4); t70 * mrSges(4,1) - t69 * mrSges(5,3) - t50 * t36 - t119 * t37 + t67 + t68 + t132 * t81 - t133 * t80 - m(4) * (-t59 * t81 + t60 * t80) - t122 - t157 * qJD(1) ^ 2 + (-t119 * t5 - t50 * t6 - t16) * m(6) + (-t52 * t81 - t53 * t80 + t25) * m(5); -(t165 * t80 - t144 + t44 + t73) * t81 / 0.2e1 - (t163 * t81 + t164 * t80) * qJD(3) / 0.2e1 + (-pkin(3) * t69 - qJ(4) * t70 - t52 * t80 + t53 * t81) * mrSges(5,2) - t92 * (mrSges(4,1) * t81 + mrSges(4,2) * t80) - t38 * (mrSges(5,1) * t81 - mrSges(5,3) * t80) + (Ifges(5,3) * t81 + t142) * t152 + (-t133 + t145) * t59 + t174 * t29 + (t11 - t170) * t154 + (-t13 * t93 + t14 * t94 - t172) * mrSges(6,3) + t167 * t37 + (t1 * t94 + t167 * t5 + t168 * t6 + t2 * t93 - t22 * t31) * m(6) + t168 * t36 + t163 * t70 + t164 * t69 + (-pkin(3) * t29 + qJ(4) * t27 + t158 * t53 - t38 * t54 - t52 * t60) * m(5) + qJD(4) * t66 - t54 * t55 - t28 * mrSges(4,2) - t31 * t15 + t27 * mrSges(5,3) - t175 + (t132 + t139) * t60 + (-Ifges(4,2) * t81 + t162 + t74) * t151 + t177 * t156 + t12 * t155 + t45 * t149; t69 * mrSges(5,2) + (-t15 + t55) * t81 + t120 * qJD(5) + (t108 * t14 - t110 * t13) * mrSges(6,3) + (-t120 - t66) * qJD(3) + (t1 * t108 + t110 * t2 - t22 * t81 + t103 * (-t108 * t5 + t110 * t6)) * m(6) + (-qJD(3) * t53 + t38 * t81 + t29) * m(5); t170 * t154 + t11 * t153 - t5 * t36 + t6 * t37 + t172 * mrSges(6,3) + (t12 - t177) * t156 + t175;];
tauc = t7(:);
