% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:03
% EndTime: 2019-12-05 17:47:11
% DurationCPUTime: 2.59s
% Computational Cost: add. (2554->275), mult. (5682->395), div. (0->0), fcn. (3713->6), ass. (0->135)
t132 = sin(qJ(5));
t134 = cos(qJ(5));
t130 = sin(pkin(8));
t131 = cos(pkin(8));
t133 = sin(qJ(3));
t135 = cos(qJ(3));
t104 = -t130 * t135 - t131 * t133;
t141 = t130 * t133 - t131 * t135;
t53 = t104 * t134 + t132 * t141;
t154 = qJD(3) * t135;
t155 = qJD(3) * t133;
t97 = t130 * t155 - t131 * t154;
t98 = t104 * qJD(3);
t138 = qJD(5) * t53 + t132 * t97 + t134 * t98;
t188 = t104 * t132 - t134 * t141;
t136 = -pkin(1) - pkin(6);
t114 = qJD(1) * t136 + qJD(2);
t152 = qJD(4) * t135;
t63 = -t114 * t155 + (qJ(4) * t155 - t152) * qJD(1);
t153 = qJD(4) * t133;
t64 = t114 * t154 + (-qJ(4) * t154 - t153) * qJD(1);
t30 = -t130 * t64 + t131 * t63;
t86 = qJD(1) * t98;
t20 = -pkin(7) * t86 + t30;
t31 = t130 * t63 + t131 * t64;
t85 = t141 * qJD(3) * qJD(1);
t21 = pkin(7) * t85 + t31;
t156 = qJD(1) * t135;
t157 = qJD(1) * t133;
t95 = t130 * t157 - t131 * t156;
t176 = pkin(7) * t95;
t89 = -qJ(4) * t157 + t114 * t133;
t74 = t130 * t89;
t90 = -qJ(4) * t156 + t135 * t114;
t78 = qJD(3) * pkin(3) + t90;
t36 = t131 * t78 - t74;
t26 = qJD(3) * pkin(4) + t176 + t36;
t96 = t104 * qJD(1);
t175 = pkin(7) * t96;
t165 = t131 * t89;
t37 = t130 * t78 + t165;
t27 = t37 + t175;
t6 = -t132 * t27 + t134 * t26;
t2 = qJD(5) * t6 + t132 * t20 + t134 * t21;
t22 = qJD(5) * t188 + t132 * t98 - t134 * t97;
t7 = t132 * t26 + t134 * t27;
t3 = -qJD(5) * t7 - t132 * t21 + t134 * t20;
t202 = t138 * t6 + t188 * t3 - t2 * t53 + t22 * t7;
t127 = qJD(3) + qJD(5);
t148 = t132 * t95 + t134 * t96;
t44 = Ifges(6,4) * t148;
t50 = t132 * t96 - t134 * t95;
t13 = Ifges(6,1) * t50 + Ifges(6,5) * t127 + t44;
t17 = qJD(5) * t148 + t132 * t85 + t134 * t86;
t170 = Ifges(6,4) * t50;
t18 = -qJD(5) * t50 - t132 * t86 + t134 * t85;
t109 = pkin(3) * t157 + qJD(1) * qJ(2) + qJD(4);
t62 = -pkin(4) * t96 + t109;
t201 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t18 - (Ifges(6,5) * t148 - Ifges(6,6) * t50) * t127 / 0.2e1 - (-Ifges(6,2) * t50 + t13 + t44) * t148 / 0.2e1 - t62 * (mrSges(6,1) * t50 + mrSges(6,2) * t148) - (Ifges(6,1) * t148 - t170) * t50 / 0.2e1;
t200 = t18 * t53;
t115 = pkin(3) * t154 + qJD(2);
t198 = t148 * t6 + t50 * t7;
t12 = Ifges(6,2) * t148 + Ifges(6,6) * t127 + t170;
t196 = t12 / 0.2e1;
t42 = -t130 * t90 - t165;
t32 = t42 - t175;
t43 = t131 * t90 - t74;
t33 = t43 + t176;
t121 = pkin(3) * t131 + pkin(4);
t169 = pkin(3) * t130;
t92 = t121 * t134 - t132 * t169;
t195 = t92 * qJD(5) - t132 * t32 - t134 * t33;
t93 = t121 * t132 + t134 * t169;
t194 = -t93 * qJD(5) + t132 * t33 - t134 * t32;
t189 = t188 * t17;
t187 = qJ(2) * (m(3) + m(4));
t184 = t148 / 0.2e1;
t182 = t50 / 0.2e1;
t181 = -t95 / 0.2e1;
t179 = t97 / 0.2e1;
t178 = t98 / 0.2e1;
t173 = t133 / 0.2e1;
t172 = -t135 / 0.2e1;
t171 = Ifges(5,4) * t95;
t158 = qJ(4) - t136;
t87 = t155 * t158 - t152;
t111 = t158 * t135;
t88 = -qJD(3) * t111 - t153;
t39 = t130 * t87 + t131 * t88;
t168 = Ifges(4,4) * t133;
t167 = Ifges(4,4) * t135;
t166 = t141 * t86;
t163 = t85 * t104;
t162 = Ifges(4,5) * qJD(3);
t161 = Ifges(4,6) * qJD(3);
t122 = t133 * pkin(3) + qJ(2);
t110 = t158 * t133;
t61 = -t131 * t110 - t130 * t111;
t108 = t115 * qJD(1);
t51 = -mrSges(5,1) * t96 - mrSges(5,2) * t95;
t151 = -m(5) * t109 - t51;
t150 = -t85 * mrSges(5,1) + t86 * mrSges(5,2);
t149 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t38 = -t130 * t88 + t131 * t87;
t60 = t110 * t130 - t131 * t111;
t145 = mrSges(4,1) * t133 + mrSges(4,2) * t135;
t40 = pkin(7) * t141 + t60;
t41 = pkin(7) * t104 + t61;
t10 = -t132 * t41 + t134 * t40;
t11 = t132 * t40 + t134 * t41;
t112 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t157;
t113 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t156;
t142 = t135 * t112 - t133 * t113;
t140 = qJ(2) * (mrSges(4,1) * t135 - mrSges(4,2) * t133);
t139 = -t104 * t31 - t141 * t30 + t36 * t98 - t37 * t97;
t106 = qJD(1) * t145;
t100 = t162 + (t135 * Ifges(4,1) - t168) * qJD(1);
t99 = t161 + (-t133 * Ifges(4,2) + t167) * qJD(1);
t91 = Ifges(5,4) * t96;
t79 = -pkin(4) * t104 + t122;
t73 = qJD(3) * mrSges(5,1) + t95 * mrSges(5,3);
t72 = -qJD(3) * mrSges(5,2) + t96 * mrSges(5,3);
t68 = pkin(3) * t156 - pkin(4) * t95;
t65 = -pkin(4) * t97 + t115;
t56 = -pkin(4) * t85 + t108;
t46 = -t95 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t91;
t45 = t96 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t171;
t35 = mrSges(6,1) * t127 - mrSges(6,3) * t50;
t34 = -mrSges(6,2) * t127 + mrSges(6,3) * t148;
t29 = pkin(7) * t97 + t39;
t28 = -pkin(7) * t98 + t38;
t19 = -mrSges(6,1) * t148 + mrSges(6,2) * t50;
t5 = -qJD(5) * t11 - t132 * t29 + t134 * t28;
t4 = qJD(5) * t10 + t132 * t28 + t134 * t29;
t1 = [(-t10 * t17 + t11 * t18 - t202) * mrSges(6,3) + (t104 * t86 - t141 * t85 + t178 * t96 + t181 * t97) * Ifges(5,4) + t108 * (-mrSges(5,1) * t104 - mrSges(5,2) * t141) + t56 * (-mrSges(6,1) * t53 + mrSges(6,2) * t188) + t138 * t13 / 0.2e1 + (t181 * t98 - t166) * Ifges(5,1) + (t179 * t96 + t163) * Ifges(5,2) + (Ifges(5,5) * t178 + Ifges(5,6) * t179 + (-t99 / 0.2e1 + t136 * t112 - t161 / 0.2e1) * t135 + (-t100 / 0.2e1 - t136 * t113 - t162 / 0.2e1) * t133) * qJD(3) + qJD(2) * t106 + t109 * (-mrSges(5,1) * t97 + mrSges(5,2) * t98) + t115 * t51 + t38 * t73 + t65 * t19 + t39 * t72 + t4 * t34 + t5 * t35 + m(5) * (t108 * t122 + t109 * t115 + t30 * t60 + t31 * t61 + t36 * t38 + t37 * t39) + m(6) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + t56 * t79 + t62 * t65) + t79 * t149 + t122 * t150 + (t138 * t184 + t17 * t53 + t18 * t188 - t182 * t22) * Ifges(6,4) + t127 * (Ifges(6,5) * t138 - Ifges(6,6) * t22) / 0.2e1 + t62 * (mrSges(6,1) * t22 + mrSges(6,2) * t138) - t22 * t196 + (-t184 * t22 + t200) * Ifges(6,2) + (-t60 * t86 + t61 * t85 - t139) * mrSges(5,3) + (((2 * mrSges(3,3)) + t145 + 0.2e1 * t187) * qJD(2) + (0.2e1 * t140 + (-0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2)) * t133 * t135 + (-0.3e1 / 0.2e1 * t135 ^ 2 + 0.3e1 / 0.2e1 * t133 ^ 2) * Ifges(4,4)) * qJD(3)) * qJD(1) + t46 * t178 + t45 * t179 + (t138 * t182 + t189) * Ifges(6,1); t22 * t34 + t138 * t35 - t97 * t72 + t98 * t73 + t142 * qJD(3) + (-t189 - t200) * mrSges(6,3) + (-t163 + t166) * mrSges(5,3) + m(5) * t139 + m(6) * t202 + (-m(6) * t62 - t106 + t151 - t19 + (-mrSges(3,3) - t187) * qJD(1)) * qJD(1); (t36 * t96 - t37 * t95 + (t130 * t85 - t131 * t86) * pkin(3)) * mrSges(5,3) + ((t130 * t31 + t131 * t30) * pkin(3) - t36 * t42 - t37 * t43) * m(5) + t50 * t196 + (-t17 * t92 + t18 * t93 + t198) * mrSges(6,3) - (Ifges(5,2) * t95 + t46 + t91) * t96 / 0.2e1 - t109 * (-t95 * mrSges(5,1) + t96 * mrSges(5,2)) - qJD(3) * (Ifges(5,5) * t96 + Ifges(5,6) * t95) / 0.2e1 + Ifges(5,5) * t86 + Ifges(5,6) * t85 - t68 * t19 - t43 * t72 - t42 * t73 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t95 * (Ifges(5,1) * t96 + t171) / 0.2e1 + (t135 * t99 / 0.2e1 + t100 * t173 + ((-Ifges(4,1) * t133 - t167) * t172 + (-Ifges(4,2) * t135 - t168) * t173 - t140) * qJD(1) + t151 * t135 * pkin(3) + (-Ifges(4,5) * t133 / 0.2e1 + Ifges(4,6) * t172) * qJD(3)) * qJD(1) + (-qJD(3) * t145 - t142) * t114 + t201 + t194 * t35 + (t194 * t6 + t195 * t7 + t2 * t93 + t3 * t92 - t62 * t68) * m(6) + t195 * t34 + t45 * t181; -t148 * t34 + t50 * t35 - t96 * t72 - t95 * t73 + t149 + t150 + (-t148 * t7 + t50 * t6 + t56) * m(6) + (-t36 * t95 - t37 * t96 + t108) * m(5); t198 * mrSges(6,3) + t12 * t182 - t6 * t34 + t7 * t35 + t201;];
tauc = t1(:);
