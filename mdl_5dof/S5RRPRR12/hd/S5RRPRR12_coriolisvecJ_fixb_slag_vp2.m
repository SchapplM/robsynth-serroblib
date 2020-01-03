% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:03
% EndTime: 2019-12-31 20:29:15
% DurationCPUTime: 4.62s
% Computational Cost: add. (3822->396), mult. (8946->522), div. (0->0), fcn. (5392->6), ass. (0->177)
t145 = sin(qJ(2));
t194 = qJD(1) * t145;
t130 = pkin(6) * t194;
t97 = pkin(7) * t194 - t130;
t262 = qJD(3) - t97;
t247 = -qJD(2) / 0.2e1;
t261 = -mrSges(3,1) - mrSges(4,1);
t144 = sin(qJ(4));
t147 = cos(qJ(4));
t148 = cos(qJ(2));
t87 = -t144 * t148 + t145 * t147;
t260 = -Ifges(5,2) / 0.2e1;
t237 = qJD(2) - qJD(4);
t193 = qJD(1) * t148;
t80 = -t144 * t193 + t147 * t194;
t248 = t80 * mrSges(5,3);
t143 = sin(qJ(5));
t146 = cos(qJ(5));
t54 = -t143 * t80 - t146 * t237;
t55 = -t143 * t237 + t146 * t80;
t208 = mrSges(5,1) * t237 - mrSges(6,1) * t54 + mrSges(6,2) * t55 + t248;
t149 = -pkin(2) - pkin(3);
t101 = t147 * qJ(3) + t144 * t149;
t131 = pkin(6) * t193;
t174 = -pkin(7) * t193 + t131;
t259 = t101 * qJD(4) + t262 * t144 + t147 * t174;
t110 = Ifges(6,5) * t143 + Ifges(6,6) * t146;
t205 = Ifges(6,4) * t143;
t111 = Ifges(6,2) * t146 + t205;
t204 = Ifges(6,4) * t146;
t112 = Ifges(6,1) * t143 + t204;
t78 = -t144 * t194 - t147 * t193;
t197 = -qJD(5) + t78;
t249 = t78 * t260;
t82 = -qJD(1) * pkin(1) - pkin(2) * t193 - qJ(3) * t194;
t62 = pkin(3) * t193 - t82;
t32 = -pkin(4) * t78 - pkin(8) * t80 + t62;
t180 = t149 * qJD(2);
t63 = t180 + t262;
t141 = qJD(2) * qJ(3);
t81 = t141 + t174;
t40 = t144 * t63 + t147 * t81;
t37 = -pkin(8) * t237 + t40;
t8 = -t143 * t37 + t146 * t32;
t9 = t143 * t32 + t146 * t37;
t151 = t62 * mrSges(5,1) + t8 * mrSges(6,1) - t9 * mrSges(6,2) - Ifges(5,4) * t80 + Ifges(6,5) * t55 + Ifges(5,6) * t237 + Ifges(6,6) * t54 - Ifges(6,3) * t197 + t249;
t222 = pkin(6) - pkin(7);
t114 = t222 * t148;
t177 = qJD(2) * t114;
t162 = qJD(1) * t177;
t39 = -t144 * t81 + t147 * t63;
t140 = qJD(2) * qJD(3);
t192 = qJD(2) * t145;
t98 = t222 * t192;
t66 = -qJD(1) * t98 + t140;
t16 = qJD(4) * t39 + t144 * t162 + t147 * t66;
t86 = t144 * t145 + t147 * t148;
t51 = t237 * t86;
t44 = t51 * qJD(1);
t28 = -qJD(5) * t55 - t143 * t44;
t230 = t28 / 0.2e1;
t27 = qJD(5) * t54 + t146 * t44;
t231 = t27 / 0.2e1;
t169 = mrSges(6,1) * t143 + mrSges(6,2) * t146;
t36 = pkin(4) * t237 - t39;
t245 = t169 * t36;
t166 = Ifges(6,5) * t146 - Ifges(6,6) * t143;
t167 = -Ifges(6,2) * t143 + t204;
t168 = Ifges(6,1) * t146 - t205;
t226 = t55 / 0.2e1;
t257 = t167 * t54 / 0.2e1 + t168 * t226 - t166 * t197 / 0.2e1;
t178 = qJD(1) * t192;
t190 = qJD(2) * t148;
t240 = t87 * qJD(4) + t144 * t190;
t45 = t240 * qJD(1) - t147 * t178;
t258 = (t110 / 0.2e1 - Ifges(5,6)) * t45 - t16 * mrSges(5,2) + Ifges(5,5) * t44 + t111 * t230 + t112 * t231 + (t245 + t257) * qJD(5) - t151 * t80;
t102 = -qJD(2) * pkin(2) + qJD(3) + t130;
t203 = Ifges(4,5) * t148;
t250 = -Ifges(3,4) * t193 / 0.2e1;
t254 = -Ifges(3,1) / 0.2e1;
t256 = (m(4) * t102 + (mrSges(4,2) + mrSges(3,3)) * t194 + t261 * qJD(2)) * pkin(6) - t82 * mrSges(4,3) + (t145 * Ifges(4,1) - t203) * qJD(1) / 0.2e1 - t194 * t254 - t250 + t102 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t247;
t106 = t131 + t141;
t108 = mrSges(4,2) * t193 + qJD(2) * mrSges(4,3);
t206 = Ifges(3,4) * t145;
t251 = -Ifges(4,5) * t194 / 0.2e1;
t253 = Ifges(4,3) / 0.2e1;
t255 = -(m(4) * t106 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t193 + t108) * pkin(6) - t106 * mrSges(4,2) - t193 * t253 - t251 - (t148 * Ifges(3,2) + t206) * qJD(1) / 0.2e1 + t82 * mrSges(4,1) + (-Ifges(4,6) + Ifges(3,6)) * t247;
t69 = Ifges(5,4) * t78;
t252 = -t69 / 0.2e1;
t244 = (t260 + Ifges(5,1) / 0.2e1) * t80;
t172 = -t143 * t9 - t146 * t8;
t239 = t172 * mrSges(6,3) + t69 / 0.2e1;
t238 = -t143 * t8 + t146 * t9;
t176 = t145 * t180;
t125 = qJ(3) * t193;
t134 = t145 * qJD(3);
t196 = qJD(1) * t134 + qJD(2) * t125;
t52 = qJD(1) * t176 + t196;
t12 = pkin(4) * t45 - pkin(8) * t44 + t52;
t1 = qJD(5) * t8 + t12 * t143 + t146 * t16;
t2 = -qJD(5) * t9 + t12 * t146 - t143 * t16;
t236 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t27 + Ifges(6,6) * t28;
t218 = t146 / 0.2e1;
t219 = -t143 / 0.2e1;
t216 = Ifges(6,4) * t55;
t23 = Ifges(6,2) * t54 - Ifges(6,6) * t197 + t216;
t53 = Ifges(6,4) * t54;
t24 = Ifges(6,1) * t55 - Ifges(6,5) * t197 + t53;
t159 = t218 * t24 + t219 * t23;
t217 = Ifges(5,1) * t80;
t234 = Ifges(5,5) * t237 + t252 - t217 / 0.2e1 - t62 * mrSges(5,2) - t159 - t245;
t229 = -t54 / 0.2e1;
t227 = -t55 / 0.2e1;
t224 = t197 / 0.2e1;
t221 = pkin(1) * mrSges(3,1);
t220 = pkin(1) * mrSges(3,2);
t113 = t222 * t145;
t163 = t147 * t113 - t114 * t144;
t17 = qJD(4) * t40 + t144 * t66 - t147 * t162;
t209 = t17 * t163;
t170 = mrSges(6,1) * t146 - mrSges(6,2) * t143;
t200 = t170 + mrSges(5,1);
t195 = qJ(3) * t190 + t134;
t191 = qJD(2) * t147;
t189 = qJD(5) * t143;
t188 = qJD(5) * t146;
t103 = -t148 * pkin(2) - t145 * qJ(3) - pkin(1);
t185 = Ifges(3,5) / 0.2e1 + Ifges(4,4) / 0.2e1;
t184 = 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4);
t182 = -Ifges(3,6) / 0.2e1 + Ifges(4,6) / 0.2e1;
t181 = m(4) * pkin(6) + mrSges(4,2);
t85 = t148 * pkin(3) - t103;
t47 = pkin(4) * t80 - pkin(8) * t78;
t173 = t1 * t146 - t143 * t2;
t10 = mrSges(6,1) * t45 - mrSges(6,3) * t27;
t11 = -mrSges(6,2) * t45 + mrSges(6,3) * t28;
t165 = -t143 * t10 + t146 * t11;
t38 = pkin(4) * t86 - pkin(8) * t87 + t85;
t57 = t113 * t144 + t114 * t147;
t20 = -t143 * t57 + t146 * t38;
t21 = t143 * t38 + t146 * t57;
t100 = -t144 * qJ(3) + t147 * t149;
t68 = t149 * t194 + t125;
t34 = mrSges(6,2) * t197 + mrSges(6,3) * t54;
t35 = -mrSges(6,1) * t197 - mrSges(6,3) * t55;
t58 = mrSges(5,2) * t237 + t78 * mrSges(5,3);
t161 = -t143 * t35 + t146 * t34 + t58;
t60 = t176 + t195;
t155 = qJD(5) * t172 + t173;
t150 = t166 * t224 + t167 * t229 + t168 * t227 + t234;
t99 = -pkin(6) * t178 + t140;
t96 = -pkin(8) + t101;
t95 = pkin(4) - t100;
t89 = (-mrSges(4,1) * t148 - mrSges(4,3) * t145) * qJD(1);
t79 = t143 * t194 + t146 * t191;
t77 = -t143 * t191 + t146 * t194;
t71 = pkin(2) * t192 - t195;
t64 = t147 * qJD(3) + qJD(4) * t100;
t61 = pkin(2) * t178 - t196;
t50 = -t145 * t191 + t240;
t49 = t144 * t174 + t147 * t97;
t46 = -mrSges(5,1) * t78 + mrSges(5,2) * t80;
t43 = Ifges(6,3) * t45;
t33 = -t47 + t68;
t30 = qJD(4) * t57 - t144 * t98 - t147 * t177;
t29 = qJD(4) * t163 + t144 * t177 - t147 * t98;
t19 = t143 * t47 + t146 * t39;
t18 = -t143 * t39 + t146 * t47;
t15 = pkin(4) * t50 - pkin(8) * t51 + t60;
t14 = t143 * t33 + t146 * t49;
t13 = -t143 * t49 + t146 * t33;
t7 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t6 = Ifges(6,1) * t27 + Ifges(6,4) * t28 + Ifges(6,5) * t45;
t5 = t27 * Ifges(6,4) + t28 * Ifges(6,2) + t45 * Ifges(6,6);
t4 = -qJD(5) * t21 - t143 * t29 + t146 * t15;
t3 = qJD(5) * t20 + t143 * t15 + t146 * t29;
t22 = [t71 * t89 + t85 * (mrSges(5,1) * t45 + mrSges(5,2) * t44) - t163 * t7 + t29 * t58 + t60 * t46 + t3 * t34 + t4 * t35 + t20 * t10 + t21 * t11 + (t249 + t151) * t50 + (t217 / 0.2e1 - t150 + t239) * t51 + t208 * t30 + m(5) * (t16 * t57 + t29 * t40 - t30 * t39 + t52 * t85 + t60 * t62 - t209) + m(6) * (t1 * t21 + t2 * t20 + t3 * t9 + t30 * t36 + t4 * t8 - t209) + m(4) * (t103 * t61 + t71 * t82) + (-t163 * t44 - t39 * t51 - t40 * t50 - t45 * t57) * mrSges(5,3) + (-t61 * mrSges(4,3) + (t182 * qJD(2) + (t103 * mrSges(4,1) + t145 * t184 - 0.2e1 * t221) * qJD(1) + t255) * qJD(2)) * t145 + (-t16 * mrSges(5,3) + t43 / 0.2e1 - Ifges(5,4) * t44 + t52 * mrSges(5,1) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t45 + t236) * t86 + (t6 * t218 + t5 * t219 + t167 * t230 + t168 * t231 + Ifges(5,1) * t44 + t52 * mrSges(5,2) + (mrSges(5,3) + t169) * t17 + (-t1 * t143 - t146 * t2) * mrSges(6,3) + (t110 * t224 + t112 * t227 + t36 * t170 + t111 * t229 - t146 * t23 / 0.2e1 + t24 * t219 - t238 * mrSges(6,3)) * qJD(5) + (t166 / 0.2e1 - Ifges(5,4)) * t45) * t87 + (-t61 * mrSges(4,1) + t181 * t99 + (t185 * qJD(2) + (-t103 * mrSges(4,3) - 0.2e1 * t220 - t184 * t148 + (-0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) + t181 * pkin(6)) * t145) * qJD(1) + t256) * qJD(2)) * t148; -t258 + (-m(4) * t82 - t89) * (pkin(2) * t194 - t125) + t208 * t259 + (-t13 * t8 - t14 * t9 + t155 * t96 + t17 * t95 + t238 * t64 + t259 * t36) * m(6) + (-t100 * t17 + t101 * t16 - t62 * t68 + (t64 - t49) * t40 - t259 * t39) * m(5) + (-t234 + t239 + t244 + t257) * t78 + ((t250 + (t220 + t203 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t261) * pkin(6) + t185) * qJD(2) - t256) * t148 + (t251 + (t221 + t206 / 0.2e1) * qJD(1) + (Ifges(3,2) / 0.2e1 + t253 + t254 - Ifges(4,1) / 0.2e1) * t193 + (pkin(6) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t182) * qJD(2) - t255) * t145) * qJD(1) + qJD(3) * t108 + t95 * t7 + t99 * mrSges(4,3) - t68 * t46 + t161 * t64 - t49 * t58 - t14 * t34 - t13 * t35 + (-t100 * t44 - t101 * t45 - t39 * t78 - t40 * t80) * mrSges(5,3) + m(4) * (qJ(3) * t99 + qJD(3) * t106) + t200 * t17 + ((-t24 / 0.2e1 + t8 * mrSges(6,3) - t96 * t35) * t146 + (t23 / 0.2e1 + t9 * mrSges(6,3) - t96 * t34) * t143) * qJD(5) + (t2 * mrSges(6,3) - t96 * t10 - t6 / 0.2e1) * t143 + (-t1 * mrSges(6,3) + t96 * t11 - t5 / 0.2e1) * t146; -qJD(2) * t108 - t79 * t34 - t77 * t35 + ((-t46 + t89) * t145 + t181 * t190) * qJD(1) + (-t44 * mrSges(5,3) - qJD(2) * t58 + t161 * qJD(4) - t7) * t147 + (-t45 * mrSges(5,3) + (-t143 * t34 - t146 * t35) * qJD(5) + t165 - t237 * t208) * t144 - m(4) * (qJD(2) * t106 - t194 * t82) + ((qJD(4) * t238 - t17) * t147 - t77 * t8 - t79 * t9 + (-t237 * t36 + t155) * t144) * m(6) + (t144 * t16 - t147 * t17 - t62 * t194 - t237 * (-t144 * t39 + t147 * t40)) * m(5); t5 * t218 + t143 * t6 / 0.2e1 - t39 * t58 - t19 * t34 - t18 * t35 - pkin(4) * t7 + ((t197 * t8 + t1) * t146 + (t197 * t9 - t2) * t143) * mrSges(6,3) + t159 * qJD(5) + (t39 * mrSges(5,3) + t150 - t244 + t252) * t78 - m(6) * (t18 * t8 + t19 * t9) + (-m(6) * pkin(4) - t200) * t17 + (m(6) * (-t188 * t8 - t189 * t9 + t173) - t35 * t188 - t34 * t189 + t165) * pkin(8) + (-m(6) * t36 - t208 + t248) * t40 + t258; t43 - t36 * (mrSges(6,1) * t55 + mrSges(6,2) * t54) + (Ifges(6,1) * t54 - t216) * t227 + t23 * t226 + (Ifges(6,5) * t54 - Ifges(6,6) * t55) * t224 - t8 * t34 + t9 * t35 + (t54 * t8 + t55 * t9) * mrSges(6,3) + (-Ifges(6,2) * t55 + t24 + t53) * t229 + t236;];
tauc = t22(:);
