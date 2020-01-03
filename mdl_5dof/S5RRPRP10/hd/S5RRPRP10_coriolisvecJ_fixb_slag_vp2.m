% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:25
% EndTime: 2019-12-31 20:09:36
% DurationCPUTime: 5.00s
% Computational Cost: add. (2148->383), mult. (5189->505), div. (0->0), fcn. (2653->4), ass. (0->179)
t254 = Ifges(5,4) + Ifges(6,4);
t255 = Ifges(5,1) + Ifges(6,1);
t244 = Ifges(5,5) + Ifges(6,5);
t253 = Ifges(5,2) + Ifges(6,2);
t243 = Ifges(5,6) + Ifges(6,6);
t128 = sin(qJ(4));
t130 = cos(qJ(4));
t131 = cos(qJ(2));
t132 = -pkin(2) - pkin(7);
t129 = sin(qJ(2));
t165 = -qJ(3) * t129 - pkin(1);
t85 = t131 * t132 + t165;
t62 = t85 * qJD(1);
t188 = qJD(1) * t129;
t119 = pkin(6) * t188;
t95 = -pkin(3) * t188 - t119;
t64 = qJD(2) * t132 + qJD(3) - t95;
t19 = -t128 * t62 + t130 * t64;
t20 = t128 * t64 + t130 * t62;
t143 = t19 * t128 - t20 * t130;
t156 = mrSges(6,1) * t130 - mrSges(6,2) * t128;
t158 = mrSges(5,1) * t130 - mrSges(5,2) * t128;
t187 = qJD(1) * t131;
t89 = -qJD(2) * t128 - t130 * t187;
t11 = qJ(5) * t89 + t20;
t185 = qJD(2) * t130;
t90 = -t128 * t187 + t185;
t10 = -qJ(5) * t90 + t19;
t116 = qJD(4) + t188;
t9 = pkin(4) * t116 + t10;
t159 = t11 * t130 - t9 * t128;
t207 = -t130 / 0.2e1;
t209 = -t128 / 0.2e1;
t251 = t254 * t89;
t228 = t244 * t116 + t255 * t90 + t251;
t247 = t254 * t90;
t229 = t243 * t116 + t253 * t89 + t247;
t127 = qJD(2) * qJ(3);
t121 = pkin(6) * t187;
t96 = pkin(3) * t187 + t121;
t78 = t127 + t96;
t38 = -pkin(4) * t89 + qJD(5) + t78;
t252 = t143 * mrSges(5,3) - t159 * mrSges(6,3) + t38 * t156 + t78 * t158 + t207 * t229 + t209 * t228;
t105 = -t121 - t127;
t210 = t116 / 0.2e1;
t215 = t90 / 0.2e1;
t249 = t254 * t130;
t224 = t255 * t128 + t249;
t248 = t254 * t128;
t225 = t253 * t130 + t248;
t193 = Ifges(6,6) * t130;
t194 = Ifges(5,6) * t130;
t196 = Ifges(6,5) * t128;
t197 = Ifges(5,5) * t128;
t226 = t193 + t196 + t194 + t197;
t230 = qJD(2) / 0.2e1;
t231 = -qJD(2) / 0.2e1;
t232 = -qJD(1) / 0.2e1;
t233 = t89 / 0.2e1;
t103 = -pkin(2) * t131 + t165;
t79 = t103 * qJD(1);
t250 = t105 * mrSges(4,1) + Ifges(3,6) * t231 + (Ifges(3,4) * t129 + t131 * Ifges(3,2)) * t232 + Ifges(4,5) * t230 + (-Ifges(4,6) * t129 - t131 * Ifges(4,3)) * qJD(1) / 0.2e1 - t79 * mrSges(4,2) + t225 * t233 + t224 * t215 + t226 * t210 - t252;
t246 = qJD(3) + t119;
t214 = pkin(3) + pkin(6);
t245 = -mrSges(3,1) + mrSges(4,2);
t178 = qJD(1) * qJD(2);
t166 = t131 * t178;
t177 = qJD(2) * qJD(4);
t180 = qJD(4) * t131;
t186 = qJD(2) * t129;
t51 = -t128 * t177 + (t128 * t186 - t130 * t180) * qJD(1);
t168 = t128 * t180;
t139 = t129 * t185 + t168;
t52 = qJD(1) * t139 - t130 * t177;
t242 = t243 * t166 + t253 * t52 + t254 * t51;
t241 = t244 * t166 + t254 * t52 + t255 * t51;
t179 = t130 * qJD(5);
t182 = qJD(4) * t128;
t189 = qJ(5) - t132;
t120 = pkin(2) * t188;
t144 = pkin(7) * t129 - qJ(3) * t131;
t70 = qJD(1) * t144 + t120;
t31 = -t128 * t70 + t130 * t96;
t240 = t182 * t189 - t179 - (-qJ(5) * t128 * t129 + pkin(4) * t131) * qJD(1) - t31;
t101 = t189 * t130;
t32 = t128 * t96 + t130 * t70;
t239 = -qJ(5) * t130 * t188 - qJD(4) * t101 - t128 * qJD(5) - t32;
t169 = -pkin(4) * t130 - pkin(3);
t181 = qJD(4) * t130;
t238 = pkin(4) * t181 - t169 * t188 + t246;
t237 = -t253 * t128 + t249;
t236 = t255 * t130 - t248;
t167 = t129 * t178;
t115 = pkin(2) * t167;
t183 = qJD(3) * t129;
t136 = qJD(2) * t144 - t183;
t44 = qJD(1) * t136 + t115;
t184 = qJD(2) * t131;
t98 = t214 * t184;
t84 = qJD(1) * t98;
t4 = t128 * t84 + t130 * t44 + t64 * t181 - t182 * t62;
t5 = -qJD(4) * t20 - t128 * t44 + t130 * t84;
t223 = t4 * t128 + t5 * t130;
t102 = -qJD(2) * pkin(2) + t246;
t118 = Ifges(3,4) * t187;
t171 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t172 = -Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t173 = -Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1;
t195 = Ifges(4,6) * t131;
t221 = t171 * t116 - t172 * t89 - t173 * t90 + t102 * mrSges(4,1) + t19 * mrSges(5,1) + t9 * mrSges(6,1) + Ifges(3,1) * t188 / 0.2e1 + Ifges(3,5) * t230 + t118 / 0.2e1 + Ifges(4,4) * t231 + (-t129 * Ifges(4,2) - t195) * t232 - t11 * mrSges(6,2) - t20 * mrSges(5,2) - t79 * mrSges(4,3) + t243 * t233 + t244 * t215 + (Ifges(6,3) + Ifges(5,3)) * t210;
t218 = -t89 / 0.2e1;
t216 = -t90 / 0.2e1;
t213 = pkin(1) * mrSges(3,1);
t212 = pkin(1) * mrSges(3,2);
t211 = -t116 / 0.2e1;
t109 = t214 * t129;
t40 = t128 * t109 + t130 * t85;
t107 = -mrSges(4,1) * t187 - qJD(2) * mrSges(4,3);
t42 = -mrSges(5,1) * t89 + mrSges(5,2) * t90;
t192 = -t107 + t42;
t191 = qJD(2) * mrSges(3,2);
t190 = t130 * t131;
t110 = t214 * t131;
t176 = -0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * Ifges(4,6);
t175 = -Ifges(4,4) / 0.2e1 + Ifges(3,5) / 0.2e1;
t174 = Ifges(4,5) / 0.2e1 - Ifges(3,6) / 0.2e1;
t170 = m(4) * pkin(6) + mrSges(4,1);
t16 = -t52 * mrSges(6,1) + t51 * mrSges(6,2);
t164 = qJ(5) * t131 - t85;
t97 = t214 * t186;
t163 = m(4) * t102 + (mrSges(4,1) + mrSges(3,3)) * t188 + t245 * qJD(2);
t162 = m(4) * t105 - mrSges(3,3) * t187 + t107 + t191;
t1 = pkin(4) * t166 - qJ(5) * t51 - qJD(5) * t90 + t5;
t2 = qJ(5) * t52 + qJD(5) * t89 + t4;
t161 = -t1 * t130 - t2 * t128;
t157 = mrSges(5,1) * t128 + mrSges(5,2) * t130;
t155 = mrSges(6,1) * t128 + mrSges(6,2) * t130;
t123 = pkin(2) * t186;
t58 = t123 + t136;
t7 = t109 * t181 + t128 * t98 + t130 * t58 - t182 * t85;
t140 = -qJ(3) * t184 - t183;
t126 = qJD(2) * qJD(3);
t67 = -qJD(1) * t97 + t126;
t54 = -mrSges(6,2) * t116 + t89 * mrSges(6,3);
t55 = -mrSges(5,2) * t116 + mrSges(5,3) * t89;
t56 = mrSges(6,1) * t116 - t90 * mrSges(6,3);
t57 = mrSges(5,1) * t116 - mrSges(5,3) * t90;
t138 = (t54 + t55) * t130 + (-t56 - t57) * t128;
t137 = t5 * mrSges(5,1) + t1 * mrSges(6,1) - t4 * mrSges(5,2) - t2 * mrSges(6,2);
t117 = pkin(4) * t128 + qJ(3);
t114 = Ifges(5,3) * t166;
t113 = Ifges(6,3) * t166;
t100 = t189 * t128;
t99 = pkin(6) * t167 - t126;
t93 = (mrSges(4,2) * t131 - mrSges(4,3) * t129) * qJD(1);
t92 = t130 * t109;
t83 = t130 * t98;
t73 = pkin(4) * t190 + t110;
t72 = t123 + t140;
t60 = qJD(1) * t140 + t115;
t49 = Ifges(5,5) * t51;
t48 = Ifges(6,5) * t51;
t47 = Ifges(5,6) * t52;
t46 = Ifges(6,6) * t52;
t43 = -pkin(4) * t168 + (-pkin(6) + t169) * t186;
t41 = -mrSges(6,1) * t89 + mrSges(6,2) * t90;
t39 = -t128 * t85 + t92;
t36 = -mrSges(5,2) * t166 + mrSges(5,3) * t52;
t35 = -mrSges(6,2) * t166 + mrSges(6,3) * t52;
t34 = mrSges(5,1) * t166 - mrSges(5,3) * t51;
t33 = mrSges(6,1) * t166 - mrSges(6,3) * t51;
t30 = -qJ(5) * t190 + t40;
t23 = pkin(4) * t129 + t128 * t164 + t92;
t22 = -pkin(4) * t52 + t67;
t17 = -mrSges(5,1) * t52 + mrSges(5,2) * t51;
t8 = -qJD(4) * t40 - t128 * t58 + t83;
t6 = qJ(5) * t139 - t131 * t179 + t7;
t3 = pkin(4) * t184 + t83 + t164 * t181 + (-qJ(5) * t186 - qJD(4) * t109 + qJD(5) * t131 - t58) * t128;
t12 = [t110 * t17 + t73 * t16 + t23 * t33 + t3 * t56 + t30 * t35 + t39 * t34 + t40 * t36 + t43 * t41 - t97 * t42 + t6 * t54 + t7 * t55 + t8 * t57 + t72 * t93 + m(5) * (t110 * t67 + t19 * t8 + t20 * t7 + t39 * t5 + t4 * t40 - t78 * t97) + m(6) * (t1 * t23 + t11 * t6 + t2 * t30 + t22 * t73 + t3 * t9 + t38 * t43) + m(4) * (t103 * t60 + t72 * t79) + (-t60 * mrSges(4,3) + t113 / 0.2e1 + t114 / 0.2e1 + t49 / 0.2e1 + t47 / 0.2e1 + t48 / 0.2e1 + t46 / 0.2e1 - t172 * t52 - t173 * t51 + (t174 * qJD(2) + (-t103 * mrSges(4,2) + t129 * t176 - 0.2e1 * t213) * qJD(1) + t162 * pkin(6) + t250) * qJD(2) + t137) * t129 + (t67 * t158 + t22 * t156 + t60 * mrSges(4,2) - t170 * t99 + (t1 * t128 - t130 * t2) * mrSges(6,3) + (t128 * t5 - t130 * t4) * mrSges(5,3) + (-t78 * t157 - t38 * t155 + (t11 * t128 + t130 * t9) * mrSges(6,3) + (t128 * t20 + t130 * t19) * mrSges(5,3) + t237 * t218 + t236 * t216 + (t243 * t128 - t130 * t244) * t210 + t229 * t128 / 0.2e1) * qJD(4) + (t163 * pkin(6) + (-t103 * mrSges(4,3) - 0.2e1 * t212 + (-t196 / 0.2e1 - t193 / 0.2e1 - t197 / 0.2e1 - t194 / 0.2e1 - t176) * t131 + (0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + t170 * pkin(6) + t171) * t129) * qJD(1) + t175 * qJD(2) + t221) * qJD(2) - t224 * t51 / 0.2e1 - t225 * t52 / 0.2e1 + t241 * t209 + (qJD(4) * t228 + t242) * t207) * t131; t161 * mrSges(6,3) + t22 * t155 + t67 * t157 + ((-m(5) * t143 - t128 * t57 + t130 * t55) * t132 + t225 * t218 + t224 * t216 + t226 * t211 + t252) * qJD(4) + (-m(4) * t79 - t93) * (-qJ(3) * t187 + t120) + t241 * t130 / 0.2e1 + t242 * t209 + t239 * t54 + (-t1 * t101 - t100 * t2 + t11 * t239 + t117 * t22 + t238 * t38 + t240 * t9) * m(6) + t240 * t56 + t238 * t41 + t237 * t52 / 0.2e1 + t236 * t51 / 0.2e1 + m(5) * (t67 * qJ(3) + t78 * qJD(3) + t132 * t223) + t192 * qJD(3) + t117 * t16 - t99 * mrSges(4,3) - t100 * t35 - t101 * t33 - t95 * t42 - t32 * t55 - t31 * t57 + qJ(3) * t17 - m(5) * (t19 * t31 + t20 * t32 + t78 * t95) + (((t212 - t195 / 0.2e1) * qJD(1) - t118 / 0.2e1 + ((-m(4) * pkin(2) + t245) * qJD(2) - t163) * pkin(6) + (-pkin(2) * mrSges(4,1) + t172 * t128 - t173 * t130 + t175) * qJD(2) - t221) * t131 + ((Ifges(3,2) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1) * t187 + (-t162 + t191) * pkin(6) + (-qJ(3) * mrSges(4,1) + t174) * qJD(2) + (t213 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t129) * qJD(1) - t250) * t129) * qJD(1) + (t128 * t36 + t130 * t34) * t132 + m(4) * (-qJ(3) * t99 - qJD(3) * t105) - t223 * mrSges(5,3); (t33 + t34) * t130 + (t35 + t36) * t128 + (-t41 - t192) * qJD(2) + t138 * qJD(4) + (t170 * t184 + (t93 + t138) * t129) * qJD(1) - m(4) * (-qJD(2) * t105 - t188 * t79) + (-qJD(2) * t38 + t116 * t159 - t161) * m(6) + (-qJD(2) * t78 - t116 * t143 + t223) * m(5); t137 + (-t41 * t90 + t33) * pkin(4) + t113 + t114 + t49 + t47 + t48 + t46 - t78 * (mrSges(5,1) * t90 + mrSges(5,2) * t89) - t38 * (mrSges(6,1) * t90 + mrSges(6,2) * t89) - t10 * t54 - t19 * t55 + t11 * t56 + t20 * t57 + (t11 * t90 + t89 * t9) * mrSges(6,3) + (t19 * t89 + t20 * t90) * mrSges(5,3) + (-(t10 - t9) * t11 + (-t38 * t90 + t1) * pkin(4)) * m(6) + (t255 * t89 - t247) * t216 + t229 * t215 + (-t243 * t90 + t244 * t89) * t211 + (-t253 * t90 + t228 + t251) * t218; -t89 * t54 + t90 * t56 + 0.2e1 * (t22 / 0.2e1 + t11 * t218 + t9 * t215) * m(6) + t16;];
tauc = t12(:);
