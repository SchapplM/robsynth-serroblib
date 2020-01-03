% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:29:04
% DurationCPUTime: 4.70s
% Computational Cost: add. (2537->341), mult. (6678->452), div. (0->0), fcn. (4432->6), ass. (0->163)
t124 = -qJD(2) + qJD(5);
t222 = t124 / 0.2e1;
t129 = sin(qJ(2));
t168 = -qJ(3) - pkin(6);
t115 = t168 * t129;
t110 = qJD(1) * t115;
t131 = cos(qJ(2));
t116 = t168 * t131;
t111 = qJD(1) * t116;
t127 = sin(pkin(8));
t160 = t127 * t111;
t163 = cos(pkin(8));
t63 = t163 * t110 + t160;
t221 = -t63 + qJD(4);
t220 = -mrSges(5,1) - mrSges(4,1);
t205 = Ifges(4,1) + Ifges(5,1);
t219 = Ifges(5,4) + Ifges(4,5);
t128 = sin(qJ(5));
t130 = cos(qJ(5));
t145 = t163 * t131;
t159 = qJD(1) * t129;
t93 = -qJD(1) * t145 + t127 * t159;
t106 = t127 * t131 + t129 * t163;
t95 = t106 * qJD(1);
t139 = t128 * t93 + t130 * t95;
t52 = t128 * t95 - t130 * t93;
t44 = Ifges(6,4) * t52;
t218 = Ifges(6,2) * t139 + t44;
t184 = t95 * pkin(7);
t217 = -t184 + t221;
t147 = qJD(2) * t168;
t92 = -t129 * qJD(3) + t131 * t147;
t132 = qJD(1) * t92;
t91 = qJD(3) * t131 + t129 * t147;
t75 = t91 * qJD(1);
t33 = t127 * t132 + t163 * t75;
t30 = qJD(2) * qJD(4) + t33;
t94 = t106 * qJD(2);
t85 = qJD(1) * t94;
t20 = pkin(7) * t85 + t30;
t32 = t127 * t75 - t163 * t132;
t135 = -t127 * t129 + t145;
t96 = t135 * qJD(2);
t86 = qJD(1) * t96;
t21 = -pkin(7) * t86 + t32;
t104 = qJD(2) * pkin(2) + t110;
t58 = t104 * t163 + t160;
t136 = qJD(4) - t58;
t186 = -pkin(3) - pkin(4);
t24 = qJD(2) * t186 + t136 - t184;
t185 = pkin(7) * t93;
t146 = t163 * t111;
t59 = t127 * t104 - t146;
t54 = qJD(2) * qJ(4) + t59;
t29 = t54 + t185;
t5 = -t128 * t29 + t130 * t24;
t1 = qJD(5) * t5 + t128 * t21 + t130 * t20;
t14 = -qJD(5) * t52 + t128 * t85 + t130 * t86;
t15 = -qJD(5) * t139 - t128 * t86 + t130 * t85;
t6 = t128 * t24 + t130 * t29;
t2 = -qJD(5) * t6 - t128 * t20 + t130 * t21;
t154 = pkin(2) * t131 + pkin(1);
t214 = qJD(1) * t154;
t112 = qJD(3) - t214;
t215 = -t95 * qJ(4) + t112;
t23 = t186 * t93 - t215;
t216 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t14 + Ifges(6,6) * t15 + (Ifges(6,5) * t52 + Ifges(6,6) * t139) * t222 + t23 * (-mrSges(6,1) * t139 + mrSges(6,2) * t52);
t194 = -t52 / 0.2e1;
t213 = t139 * t6 - t5 * t52;
t178 = Ifges(6,4) * t139;
t211 = -Ifges(6,1) * t52 - t178;
t193 = -t139 / 0.2e1;
t62 = t110 * t127 - t146;
t34 = t62 + t185;
t153 = t163 * pkin(2);
t121 = -t153 - pkin(3);
t118 = -pkin(4) + t121;
t176 = pkin(2) * t127;
t119 = qJ(4) + t176;
t76 = t118 * t130 - t119 * t128;
t209 = qJD(5) * t76 - t128 * t34 + t217 * t130;
t77 = t118 * t128 + t119 * t130;
t208 = -qJD(5) * t77 - t217 * t128 - t130 * t34;
t206 = mrSges(5,2) + mrSges(4,3);
t169 = -Ifges(4,4) + Ifges(5,5);
t203 = -Ifges(4,6) + Ifges(5,6);
t177 = Ifges(5,5) * t93;
t90 = Ifges(4,4) * t93;
t202 = t219 * qJD(2) + t205 * t95 + t177 - t90;
t171 = t93 * mrSges(4,3);
t182 = mrSges(5,2) * t93;
t74 = qJD(2) * mrSges(5,3) - t182;
t167 = -qJD(2) * mrSges(4,2) - t171 + t74;
t180 = mrSges(4,3) * t95;
t181 = mrSges(5,2) * t95;
t166 = -t220 * qJD(2) - t180 - t181;
t158 = qJD(1) * t131;
t161 = Ifges(3,6) * qJD(2);
t165 = Ifges(3,4) * t129;
t198 = t161 / 0.2e1 + (t131 * Ifges(3,2) + t165) * qJD(1) / 0.2e1 + pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t158);
t197 = -0.2e1 * pkin(1);
t195 = t52 / 0.2e1;
t192 = t139 / 0.2e1;
t191 = -t93 / 0.2e1;
t190 = t93 / 0.2e1;
t188 = t95 / 0.2e1;
t179 = Ifges(4,4) * t95;
t175 = pkin(6) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t159);
t66 = -t163 * t115 - t116 * t127;
t173 = t32 * t66;
t172 = t86 * mrSges(5,2);
t41 = t127 * t92 + t163 * t91;
t162 = Ifges(3,5) * qJD(2);
t67 = t127 * t115 - t163 * t116;
t157 = qJD(2) * t129;
t156 = qJD(1) * qJD(2);
t155 = pkin(2) * t159;
t152 = t85 * mrSges(4,1) + t86 * mrSges(4,2);
t151 = t85 * mrSges(5,1) - t86 * mrSges(5,3);
t150 = t129 * t156;
t149 = t131 * t156;
t40 = t127 * t91 - t163 * t92;
t144 = pkin(2) * t150;
t143 = -t15 * mrSges(6,1) + t14 * mrSges(6,2);
t36 = -mrSges(6,2) * t124 - t52 * mrSges(6,3);
t37 = mrSges(6,1) * t124 - mrSges(6,3) * t139;
t140 = -t128 * t37 + t130 * t36;
t42 = -pkin(7) * t106 + t66;
t43 = -pkin(7) * t135 + t67;
t9 = -t128 * t43 + t130 * t42;
t10 = t128 * t42 + t130 * t43;
t60 = -t106 * t128 - t130 * t135;
t61 = t106 * t130 - t128 * t135;
t138 = t106 * qJ(4) + t154;
t137 = -qJ(4) * t93 - t155;
t22 = t85 * pkin(3) - t86 * qJ(4) - t95 * qJD(4) + t144;
t134 = -pkin(2) * t157 + qJ(4) * t96 + qJD(4) * t106;
t123 = Ifges(3,4) * t158;
t103 = Ifges(3,1) * t159 + t123 + t162;
t89 = Ifges(5,5) * t95;
t57 = -pkin(3) * t135 - t138;
t56 = mrSges(4,1) * t93 + mrSges(4,2) * t95;
t55 = mrSges(5,1) * t93 - mrSges(5,3) * t95;
t47 = -t93 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t179;
t46 = Ifges(5,6) * qJD(2) + t93 * Ifges(5,3) + t89;
t45 = -qJD(2) * pkin(3) + t136;
t39 = pkin(3) * t95 - t137;
t38 = t93 * pkin(3) + t215;
t31 = -t135 * t186 + t138;
t28 = pkin(3) * t94 - t134;
t27 = pkin(7) * t94 + t41;
t26 = -pkin(7) * t96 + t40;
t25 = t186 * t95 + t137;
t19 = t186 * t94 + t134;
t18 = -qJD(5) * t61 - t128 * t96 + t130 * t94;
t17 = qJD(5) * t60 + t128 * t94 + t130 * t96;
t16 = mrSges(6,1) * t52 + mrSges(6,2) * t139;
t13 = pkin(4) * t85 + t22;
t12 = Ifges(6,1) * t139 + Ifges(6,5) * t124 - t44;
t11 = -Ifges(6,2) * t52 + Ifges(6,6) * t124 + t178;
t4 = -qJD(5) * t10 - t128 * t27 + t130 * t26;
t3 = qJD(5) * t9 + t128 * t26 + t130 * t27;
t7 = [((Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * qJD(2) - t58 * mrSges(4,3) + t45 * mrSges(5,2) + t112 * mrSges(4,2) - t38 * mrSges(5,3) + Ifges(5,5) * t190 + Ifges(4,4) * t191 + t202 / 0.2e1 + t205 * t188) * t96 + (t14 * t61 + t17 * t192) * Ifges(6,1) + (t60 * t15 + t18 * t194) * Ifges(6,2) + (t60 * t14 + t15 * t61 + t17 * t194 + t18 * t192) * Ifges(6,4) + m(4) * (t33 * t67 - t40 * t58 + t41 * t59 + t173) + m(5) * (t22 * t57 + t28 * t38 + t30 * t67 + t40 * t45 + t41 * t54 + t173) - t166 * t40 + t167 * t41 + ((Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(2) - t59 * mrSges(4,3) - t54 * mrSges(5,2) + t112 * mrSges(4,1) + t38 * mrSges(5,1) + t46 / 0.2e1 - t47 / 0.2e1 + Ifges(5,3) * t190 - Ifges(4,2) * t191 + t169 * t188) * t94 + t57 * t151 + t31 * t143 + (-t161 / 0.2e1 + (mrSges(3,1) * t197 - 0.3e1 / 0.2e1 * t165 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t131) * qJD(1) + (m(4) * (t112 - t214) + t56 + qJD(1) * (-mrSges(4,1) * t135 + mrSges(4,2) * t106)) * pkin(2) - t198) * t157 + (-mrSges(5,3) * t22 + t169 * t85 + t205 * t86 + t206 * t32) * t106 + t206 * (t66 * t86 - t67 * t85) + (t103 / 0.2e1 - t175 + t162 / 0.2e1 + (mrSges(3,2) * t197 + 0.3e1 / 0.2e1 * Ifges(3,4) * t131) * qJD(1)) * t131 * qJD(2) + t28 * t55 - t13 * (-mrSges(6,1) * t60 + mrSges(6,2) * t61) + t3 * t36 + t4 * t37 + t23 * (-mrSges(6,1) * t18 + mrSges(6,2) * t17) + m(6) * (t1 * t10 - t13 * t31 + t19 * t23 + t2 * t9 + t3 * t6 + t4 * t5) + t17 * t12 / 0.2e1 + t18 * t11 / 0.2e1 + t19 * t16 - (mrSges(5,1) * t22 - mrSges(5,2) * t30 - mrSges(4,3) * t33 + t169 * t86 + (Ifges(4,2) + Ifges(5,3)) * t85) * t135 - t154 * t152 + (Ifges(6,5) * t17 + Ifges(6,6) * t18) * t222 + (t1 * t60 + t10 * t15 - t14 * t9 - t17 * t5 + t18 * t6 - t2 * t61) * mrSges(6,3); -t58 * t171 - (Ifges(3,5) * t131 - Ifges(3,6) * t129) * t156 / 0.2e1 - t56 * t155 - Ifges(3,6) * t150 + (t11 - t211) * t193 + (-t14 * t76 + t15 * t77 - t213) * mrSges(6,3) + (-mrSges(4,3) * t153 + t219) * t86 - (t203 * t95 - t219 * t93) * qJD(2) / 0.2e1 + t220 * t32 + t208 * t37 + (t1 * t77 + t2 * t76 + t208 * t5 + t209 * t6 - t23 * t25) * m(6) + t209 * t36 + (-Ifges(4,2) * t95 + t202 - t90) * t190 + (-mrSges(5,2) * t119 - mrSges(4,3) * t176 + t203) * t85 - (-t205 * t93 - t179 + t46 + t89) * t95 / 0.2e1 + t166 * t62 - t167 * t63 + t198 * t159 + t218 * t195 - t216 - t112 * (mrSges(4,1) * t95 - mrSges(4,2) * t93) - t38 * (mrSges(5,1) * t95 + mrSges(5,3) * t93) + qJD(4) * t74 - t39 * t55 + t30 * mrSges(5,3) - t33 * mrSges(4,2) - t25 * t16 + (t119 * t30 + t121 * t32 + t221 * t54 - t38 * t39 - t45 * t62) * m(5) + ((t127 * t33 - t163 * t32) * pkin(2) - t112 * t155 + t58 * t62 - t59 * t63) * m(4) + t12 * t194 + (-mrSges(3,1) * t149 + mrSges(3,2) * t150) * pkin(6) - (-Ifges(3,2) * t159 + t103 + t123) * t158 / 0.2e1 + (-t129 * (Ifges(3,1) * t131 - t165) / 0.2e1 + pkin(1) * (mrSges(3,1) * t129 + mrSges(3,2) * t131)) * qJD(1) ^ 2 + t121 * t172 + t158 * t175 + t59 * t180 + t54 * t181 + t45 * t182 + t47 * t188 + (Ifges(5,3) * t95 - t177) * t191 + Ifges(3,5) * t149; t166 * t95 + t167 * t93 - t52 * t36 - t139 * t37 - t143 + t151 + t152 + (-t139 * t5 - t52 * t6 + t13) * m(6) + (-t45 * t95 + t54 * t93 + t22) * m(5) + (t58 * t95 + t59 * t93 + t144) * m(4); t172 + (-t16 + t55) * t95 + t140 * qJD(5) + (t128 * t15 - t130 * t14) * mrSges(6,3) + (-t140 - t74) * qJD(2) + (t1 * t128 + t130 * t2 - t23 * t95 + t124 * (-t128 * t5 + t130 * t6)) * m(6) + (-qJD(2) * t54 + t38 * t95 + t32) * m(5); t211 * t193 + t11 * t192 - t5 * t36 + t6 * t37 + t213 * mrSges(6,3) + (t12 - t218) * t195 + t216;];
tauc = t7(:);
