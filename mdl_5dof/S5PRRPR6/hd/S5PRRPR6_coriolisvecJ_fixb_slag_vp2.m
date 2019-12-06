% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:07
% EndTime: 2019-12-05 16:30:25
% DurationCPUTime: 4.79s
% Computational Cost: add. (2732->404), mult. (7352->605), div. (0->0), fcn. (5178->10), ass. (0->200)
t153 = sin(qJ(3));
t156 = cos(qJ(3));
t169 = pkin(3) * t153 - qJ(4) * t156;
t108 = qJD(3) * t169 - qJD(4) * t153;
t148 = sin(pkin(10));
t150 = cos(pkin(10));
t195 = qJD(3) * t153;
t189 = pkin(7) * t195;
t74 = t108 * t150 + t148 * t189;
t154 = sin(qJ(2));
t149 = sin(pkin(5));
t200 = qJD(1) * t149;
t157 = cos(qJ(2));
t203 = t156 * t157;
t90 = (-t148 * t203 + t150 * t154) * t200;
t255 = -t90 + t74;
t205 = t150 * t156;
t166 = pkin(4) * t153 - pkin(8) * t205;
t163 = t166 * qJD(3);
t254 = t163 + t255;
t102 = t148 * t108;
t209 = t148 * t156;
t191 = pkin(8) * t209;
t206 = t150 * t153;
t91 = (t148 * t154 + t150 * t203) * t200;
t253 = t102 + (-pkin(7) * t206 - t191) * qJD(3) - t91;
t252 = qJD(3) / 0.2e1;
t152 = sin(qJ(5));
t155 = cos(qJ(5));
t167 = t148 * t152 - t150 * t155;
t128 = t169 * qJD(2);
t187 = t154 * t200;
t131 = qJD(2) * pkin(7) + t187;
t120 = t153 * t131;
t151 = cos(pkin(5));
t204 = t151 * t156;
t184 = qJD(1) * t204;
t92 = -t120 + t184;
t50 = t128 * t150 - t148 * t92;
t35 = qJD(2) * t166 + t50;
t51 = t128 * t148 + t150 * t92;
t40 = -qJD(2) * t191 + t51;
t225 = pkin(8) + qJ(4);
t134 = t225 * t148;
t135 = t225 * t150;
t83 = -t134 * t155 - t135 * t152;
t251 = -qJD(4) * t167 + qJD(5) * t83 - t152 * t35 - t155 * t40;
t126 = t148 * t155 + t150 * t152;
t84 = -t134 * t152 + t135 * t155;
t250 = -qJD(4) * t126 - qJD(5) * t84 + t152 * t40 - t155 * t35;
t197 = qJD(2) * t153;
t249 = t197 / 0.2e1;
t181 = -Ifges(4,6) * qJD(3) / 0.2e1;
t248 = Ifges(4,5) * t252;
t133 = -pkin(3) * t156 - qJ(4) * t153 - pkin(2);
t122 = t150 * t133;
t69 = -pkin(8) * t206 + t122 + (-pkin(7) * t148 - pkin(4)) * t156;
t96 = pkin(7) * t205 + t133 * t148;
t80 = -pkin(8) * t148 * t153 + t96;
t27 = t152 * t69 + t155 * t80;
t247 = -qJD(5) * t27 - t152 * t253 + t155 * t254;
t26 = -t152 * t80 + t155 * t69;
t246 = qJD(5) * t26 + t152 * t254 + t155 * t253;
t193 = qJD(2) * qJD(3);
t179 = t156 * t193;
t175 = t148 * t179;
t224 = mrSges(5,2) * t150;
t101 = mrSges(5,1) * t175 + t179 * t224;
t164 = t167 * t156;
t161 = qJD(3) * t164;
t123 = qJD(3) * t150 - t148 * t197;
t124 = qJD(3) * t148 + t150 * t197;
t178 = t123 * t155 - t124 * t152;
t33 = -qJD(2) * t161 + qJD(5) * t178;
t165 = t126 * t156;
t162 = qJD(3) * t165;
t66 = t123 * t152 + t124 * t155;
t34 = -qJD(2) * t162 - qJD(5) * t66;
t11 = -t34 * mrSges(6,1) + mrSges(6,2) * t33;
t245 = t101 + t11;
t110 = t167 * qJD(5);
t198 = qJD(2) * t149;
t182 = t157 * t198;
t176 = qJD(1) * t182;
t201 = qJD(3) * t184 + t156 * t176;
t55 = (qJD(4) - t120) * qJD(3) + t201;
t72 = (t108 + t187) * qJD(2);
t19 = -t148 * t55 + t150 * t72;
t14 = qJD(2) * t163 + t19;
t20 = t148 * t72 + t150 * t55;
t15 = -pkin(8) * t175 + t20;
t196 = qJD(2) * t156;
t199 = qJD(1) * t153;
t185 = t151 * t199;
t93 = t131 * t156 + t185;
t89 = qJD(3) * qJ(4) + t93;
t186 = t157 * t200;
t94 = qJD(2) * t133 - t186;
t36 = -t148 * t89 + t150 * t94;
t21 = -pkin(4) * t196 - pkin(8) * t124 + t36;
t37 = t148 * t94 + t150 * t89;
t23 = pkin(8) * t123 + t37;
t7 = -t152 * t23 + t155 * t21;
t1 = qJD(5) * t7 + t14 * t152 + t15 * t155;
t8 = t152 * t21 + t155 * t23;
t2 = -qJD(5) * t8 + t14 * t155 - t15 * t152;
t244 = -t2 * mrSges(6,1) + t1 * mrSges(6,2) - Ifges(6,5) * t33 - Ifges(6,6) * t34;
t132 = -qJD(2) * pkin(2) - t186;
t147 = Ifges(4,4) * t196;
t220 = Ifges(5,2) * t148;
t222 = Ifges(5,4) * t150;
t170 = -t220 + t222;
t223 = Ifges(5,4) * t148;
t171 = Ifges(5,1) * t150 - t223;
t172 = mrSges(5,1) * t148 + t224;
t231 = t150 / 0.2e1;
t232 = -t148 / 0.2e1;
t85 = -qJD(3) * pkin(3) + qJD(4) - t92;
t243 = -(t148 * t37 + t150 * t36) * mrSges(5,3) + t132 * mrSges(4,2) + t85 * t172 + Ifges(4,1) * t249 + t147 / 0.2e1 + t248 + t123 * t170 / 0.2e1 + t124 * t171 / 0.2e1 + (Ifges(5,4) * t124 + Ifges(5,2) * t123 - Ifges(5,6) * t196) * t232 + (t124 * Ifges(5,1) + t123 * Ifges(5,4) - Ifges(5,5) * t196) * t231 - t92 * mrSges(4,3);
t242 = t33 / 0.2e1;
t241 = t34 / 0.2e1;
t240 = -t178 / 0.2e1;
t239 = t178 / 0.2e1;
t238 = -t66 / 0.2e1;
t237 = t66 / 0.2e1;
t104 = t126 * t153;
t236 = -t104 / 0.2e1;
t105 = t167 * t153;
t235 = -t105 / 0.2e1;
t144 = qJD(5) - t196;
t234 = -t144 / 0.2e1;
t233 = t144 / 0.2e1;
t230 = Ifges(6,4) * t66;
t229 = pkin(4) * t148;
t221 = Ifges(5,5) * t150;
t219 = Ifges(5,6) * t148;
t208 = t149 * t154;
t112 = t153 * t208 - t204;
t194 = qJD(3) * t156;
t58 = t131 * t194 + (qJD(3) * t151 + t182) * t199;
t218 = t112 * t58;
t111 = t126 * qJD(5);
t99 = qJD(2) * t165;
t214 = -t111 + t99;
t213 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t123 + mrSges(5,2) * t124 + mrSges(4,3) * t197;
t207 = t149 * t157;
t100 = qJD(2) * t164;
t202 = -t100 + t110;
t192 = pkin(7) * t153 * t58;
t22 = -mrSges(6,1) * t178 + mrSges(6,2) * t66;
t190 = t22 + t213;
t188 = pkin(7) + t229;
t183 = t154 * t198;
t180 = t153 * t193;
t177 = t153 * t186;
t173 = -mrSges(4,1) * t156 + mrSges(4,2) * t153;
t113 = t151 * t153 + t156 * t208;
t76 = -t113 * t148 - t150 * t207;
t77 = t113 * t150 - t148 * t207;
t28 = -t152 * t77 + t155 * t76;
t29 = t152 * t76 + t155 * t77;
t70 = t185 + (qJD(2) * t229 + t131) * t156;
t160 = t132 * mrSges(4,1) + t36 * mrSges(5,1) + t7 * mrSges(6,1) + t181 - (Ifges(4,4) * t153 + t156 * Ifges(4,2)) * qJD(2) / 0.2e1 + t144 * Ifges(6,3) + t66 * Ifges(6,5) + t178 * Ifges(6,6) - Ifges(5,3) * t196 / 0.2e1 + t124 * Ifges(5,5) + t123 * Ifges(5,6) - t37 * mrSges(5,2) - t8 * mrSges(6,2) - t93 * mrSges(4,3);
t158 = qJD(2) ^ 2;
t146 = -pkin(4) * t150 - pkin(3);
t142 = Ifges(6,3) * t180;
t139 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t196;
t129 = t188 * t153;
t127 = t173 * qJD(2);
t119 = t188 * t194;
t118 = (mrSges(4,1) * t153 + mrSges(4,2) * t156) * t193;
t107 = (mrSges(5,1) * t153 - mrSges(5,3) * t205) * t193;
t106 = (-mrSges(5,2) * t153 - mrSges(5,3) * t209) * t193;
t98 = -mrSges(5,1) * t196 - mrSges(5,3) * t124;
t97 = mrSges(5,2) * t196 + mrSges(5,3) * t123;
t95 = -pkin(7) * t209 + t122;
t82 = (Ifges(5,5) * t153 + t156 * t171) * t193;
t81 = (Ifges(5,6) * t153 + t156 * t170) * t193;
t79 = -qJD(3) * t112 + t156 * t182;
t78 = qJD(3) * t113 + t153 * t182;
t75 = -t150 * t189 + t102;
t62 = Ifges(6,4) * t178;
t57 = -t131 * t195 + t201;
t54 = t110 * t153 - t162;
t53 = -t111 * t153 - t161;
t52 = -pkin(4) * t123 + t85;
t49 = t148 * t183 + t150 * t79;
t48 = -t148 * t79 + t150 * t183;
t45 = mrSges(6,1) * t144 - mrSges(6,3) * t66;
t44 = -mrSges(6,2) * t144 + mrSges(6,3) * t178;
t41 = qJD(3) * t70 + t153 * t176;
t25 = -mrSges(6,2) * t180 + mrSges(6,3) * t34;
t24 = mrSges(6,1) * t180 - mrSges(6,3) * t33;
t18 = Ifges(6,1) * t66 + Ifges(6,5) * t144 + t62;
t17 = Ifges(6,2) * t178 + Ifges(6,6) * t144 + t230;
t10 = Ifges(6,1) * t33 + Ifges(6,4) * t34 + Ifges(6,5) * t180;
t9 = Ifges(6,4) * t33 + Ifges(6,2) * t34 + Ifges(6,6) * t180;
t4 = -qJD(5) * t29 - t152 * t49 + t155 * t48;
t3 = qJD(5) * t28 + t152 * t48 + t155 * t49;
t5 = [-t113 * mrSges(4,3) * t180 + t77 * t106 + t76 * t107 + t79 * t139 + t28 * t24 + t29 * t25 + t3 * t44 + t4 * t45 + t48 * t98 + t49 * t97 + t190 * t78 + (mrSges(4,3) * t179 + t245) * t112 + ((-mrSges(3,2) * t158 - t118) * t157 + (-mrSges(3,1) * t158 + qJD(2) * t127) * t154) * t149 + m(4) * (t218 + t113 * t57 - t78 * t92 + t79 * t93 + (t132 - t186) * t183) + m(5) * (t19 * t76 + t20 * t77 + t36 * t48 + t37 * t49 + t78 * t85 + t218) + m(6) * (t1 * t29 + t112 * t41 + t2 * t28 + t3 * t8 + t4 * t7 + t52 * t78); (-t91 + t75) * t97 + t129 * t11 - m(4) * (t132 * t154 + (-t153 * t92 + t156 * t93) * t157) * t200 - pkin(2) * t118 + t119 * t22 + t96 * t106 + t95 * t107 + t52 * (-mrSges(6,1) * t54 + mrSges(6,2) * t53) + t54 * t17 / 0.2e1 + t53 * t18 / 0.2e1 + t26 * t24 + t27 * t25 + t246 * t44 + t247 * t45 + (t1 * t27 + t129 * t41 + t2 * t26 + t246 * t8 + t247 * t7 + (t119 - t177) * t52) * m(6) + ((t160 + (-m(4) * t93 - t139) * pkin(7) + t181) * t153 + (t248 + (-m(4) * t92 + m(5) * t85 + t213) * pkin(7) + t243) * t156) * qJD(3) + (-t19 * mrSges(5,1) + t20 * mrSges(5,2) - t142 / 0.2e1 + t57 * mrSges(4,3) - t139 * t186 + t244) * t156 + (Ifges(6,5) * t53 + Ifges(6,6) * t54) * t233 + t10 * t235 + t9 * t236 + (Ifges(6,1) * t53 + Ifges(6,4) * t54) * t237 + (Ifges(6,4) * t53 + Ifges(6,2) * t54) * t239 + (pkin(7) * t101 + t82 * t231 + t81 * t232 + (mrSges(4,3) + t172) * t58 + (-t148 * t20 - t150 * t19) * mrSges(5,3) - t190 * t186) * t153 + ((-m(4) * pkin(2) + t173) * t187 + ((Ifges(6,5) * t235 + Ifges(6,6) * t236 + (-0.3e1 / 0.2e1 * Ifges(4,4) + t221 / 0.2e1 - t219 / 0.2e1) * t153) * t153 + ((-0.3e1 / 0.2e1 * t221 + 0.3e1 / 0.2e1 * t219 + 0.3e1 / 0.2e1 * Ifges(4,4)) * t156 + (-0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - Ifges(6,3) / 0.2e1 + Ifges(5,1) * t150 ^ 2 / 0.2e1 + (-t222 + t220 / 0.2e1) * t148) * t153) * t156) * qJD(3)) * qJD(2) + m(4) * (pkin(7) * t156 * t57 + t192) + m(5) * (t19 * t95 + t20 * t96 + t36 * t74 + t37 * t75 + t192) - t127 * t187 - m(5) * (t177 * t85 + t36 * t90 + t37 * t91) + (-t1 * t104 + t105 * t2 - t53 * t7 + t54 * t8) * mrSges(6,3) + t41 * (mrSges(6,1) * t104 - mrSges(6,2) * t105) + (-Ifges(6,4) * t105 - Ifges(6,2) * t104) * t241 + (-Ifges(6,1) * t105 - Ifges(6,4) * t104) * t242 + t255 * t98; t146 * t11 - t92 * t139 + t126 * t10 / 0.2e1 + t250 * t45 + (t1 * t84 + t146 * t41 + t2 * t83 + t250 * t7 + t251 * t8 - t52 * t70) * m(6) + t251 * t44 + ((-t147 / 0.2e1 + (-t219 + t221) * t196 / 0.2e1 + (Ifges(4,5) / 0.2e1 + (Ifges(5,1) * t148 + t222) * t231 + (Ifges(5,2) * t150 + t223) * t232) * qJD(3) - t243) * t156 + (-t160 + t181 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t196 + Ifges(4,4) * t249 + (Ifges(5,5) * t148 + Ifges(6,5) * t126 + Ifges(5,6) * t150 - Ifges(6,6) * t167) * t252) * t153) * qJD(2) - t167 * t9 / 0.2e1 + t41 * (mrSges(6,1) * t167 + mrSges(6,2) * t126) + (Ifges(6,4) * t126 - Ifges(6,2) * t167) * t241 + (Ifges(6,1) * t126 - Ifges(6,4) * t167) * t242 + (-t1 * t167 - t126 * t2 + t202 * t7 + t214 * t8) * mrSges(6,3) + (-t111 / 0.2e1 + t99 / 0.2e1) * t17 + (-Ifges(6,5) * t100 - Ifges(6,6) * t99) * t234 + (-Ifges(6,1) * t100 - Ifges(6,4) * t99) * t238 + (-Ifges(6,4) * t100 - Ifges(6,2) * t99) * t240 + (-t110 / 0.2e1 + t100 / 0.2e1) * t18 + (-t36 * t50 - t37 * t51 - t85 * t93 - pkin(3) * t58 + (-t148 * t36 + t150 * t37) * qJD(4) + (-t148 * t19 + t150 * t20) * qJ(4)) * m(5) - t51 * t97 - t50 * t98 - pkin(3) * t101 + t83 * t24 + t84 * t25 - t70 * t22 - t57 * mrSges(4,2) - t58 * mrSges(4,1) - t213 * t93 + (-mrSges(6,1) * t214 - mrSges(6,2) * t202) * t52 + (-t58 * mrSges(5,1) + t81 / 0.2e1 + qJD(4) * t97 + qJ(4) * t106 + t20 * mrSges(5,3)) * t150 + (t82 / 0.2e1 + t58 * mrSges(5,2) - qJD(4) * t98 - qJ(4) * t107 - t19 * mrSges(5,3)) * t148 + (-Ifges(6,5) * t110 - Ifges(6,6) * t111) * t233 + (-Ifges(6,1) * t110 - Ifges(6,4) * t111) * t237 + (-Ifges(6,4) * t110 - Ifges(6,2) * t111) * t239; -t123 * t97 + t124 * t98 - t178 * t44 + t66 * t45 + (-t178 * t8 + t66 * t7 + t41) * m(6) + (-t123 * t37 + t124 * t36 + t58) * m(5) + t245; t142 - t52 * (mrSges(6,1) * t66 + mrSges(6,2) * t178) + (Ifges(6,1) * t178 - t230) * t238 + t17 * t237 + (Ifges(6,5) * t178 - Ifges(6,6) * t66) * t234 + t8 * t45 - t7 * t44 + (t178 * t7 + t66 * t8) * mrSges(6,3) + (-Ifges(6,2) * t66 + t18 + t62) * t240 - t244;];
tauc = t5(:);
