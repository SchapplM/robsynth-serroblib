% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:00:12
% EndTime: 2019-03-09 02:00:20
% DurationCPUTime: 3.93s
% Computational Cost: add. (4712->382), mult. (11609->487), div. (0->0), fcn. (8059->8), ass. (0->177)
t246 = Ifges(6,1) + Ifges(7,1);
t242 = Ifges(7,4) + Ifges(6,5);
t248 = Ifges(6,6) - Ifges(7,6);
t247 = Ifges(5,2) / 0.2e1;
t122 = sin(pkin(10));
t124 = cos(pkin(10));
t127 = sin(qJ(4));
t129 = cos(qJ(4));
t111 = t122 * t127 - t129 * t124;
t104 = t111 * qJD(1);
t237 = qJD(5) + t104;
t112 = t122 * t129 + t124 * t127;
t105 = t112 * qJD(1);
t244 = t105 * Ifges(5,1) / 0.2e1;
t243 = t104 * t247;
t171 = qJD(1) * (t122 ^ 2 + t124 ^ 2);
t241 = mrSges(4,3) * t171;
t139 = -cos(pkin(9)) * pkin(1) - pkin(3) * t124 - pkin(2);
t102 = qJD(1) * t139 + qJD(3);
t126 = sin(qJ(5));
t128 = cos(qJ(5));
t115 = sin(pkin(9)) * pkin(1) + qJ(3);
t114 = t115 * qJD(1);
t118 = t124 * qJD(2);
t191 = pkin(7) * qJD(1);
t86 = t118 + (-t114 - t191) * t122;
t96 = t122 * qJD(2) + t124 * t114;
t87 = t124 * t191 + t96;
t50 = t127 * t86 + t129 * t87;
t47 = qJD(4) * pkin(8) + t50;
t57 = pkin(4) * t104 - pkin(8) * t105 + t102;
t17 = -t126 * t47 + t128 * t57;
t18 = t126 * t57 + t128 * t47;
t145 = t126 * t18 + t128 * t17;
t232 = qJD(6) - t17;
t11 = -pkin(5) * t237 + t232;
t12 = qJ(6) * t237 + t18;
t148 = t11 * t128 - t12 * t126;
t192 = Ifges(7,5) * t128;
t152 = Ifges(7,3) * t126 + t192;
t194 = Ifges(6,4) * t128;
t158 = -Ifges(6,2) * t126 + t194;
t163 = mrSges(7,1) * t126 - mrSges(7,3) * t128;
t165 = mrSges(6,1) * t126 + mrSges(6,2) * t128;
t140 = t128 * qJD(4) - t105 * t126;
t49 = -t127 * t87 + t129 * t86;
t46 = -qJD(4) * pkin(4) - t49;
t89 = qJD(4) * t126 + t105 * t128;
t21 = -pkin(5) * t140 - qJ(6) * t89 + t46;
t213 = t128 / 0.2e1;
t215 = t126 / 0.2e1;
t216 = -t126 / 0.2e1;
t220 = t89 / 0.2e1;
t222 = -t140 / 0.2e1;
t223 = t140 / 0.2e1;
t193 = Ifges(7,5) * t126;
t195 = Ifges(6,4) * t126;
t230 = t246 * t128 + t193 - t195;
t231 = -t126 * t248 + t242 * t128;
t209 = Ifges(7,5) * t140;
t85 = Ifges(6,4) * t140;
t235 = t242 * t237 + t246 * t89 - t209 + t85;
t84 = Ifges(7,5) * t89;
t35 = Ifges(7,6) * t237 - Ifges(7,3) * t140 + t84;
t210 = Ifges(6,4) * t89;
t38 = Ifges(6,2) * t140 + Ifges(6,6) * t237 + t210;
t131 = t148 * mrSges(7,2) - t145 * mrSges(6,3) + t152 * t222 + t158 * t223 + t163 * t21 + t165 * t46 + t215 * t35 + t216 * t38 + t230 * t220 + t231 * t237 / 0.2e1 + t235 * t213;
t240 = t102 * mrSges(5,2) - t49 * mrSges(5,3) - Ifges(5,4) * t104 + Ifges(5,5) * qJD(4) + t131 + t244;
t239 = t242 * t126 + t248 * t128;
t238 = t246 * t126 - t192 + t194;
t221 = -t89 / 0.2e1;
t218 = -t237 / 0.2e1;
t106 = t111 * qJD(4);
t97 = qJD(1) * t106;
t63 = qJD(5) * t140 - t128 * t97;
t64 = qJD(5) * t89 - t126 * t97;
t107 = t112 * qJD(4);
t98 = qJD(1) * t107;
t236 = t242 * t98 + (-Ifges(6,4) + Ifges(7,5)) * t64 + t246 * t63;
t149 = pkin(5) * t126 - qJ(6) * t128;
t234 = -qJD(6) * t126 + t237 * t149 - t50;
t202 = pkin(7) + t115;
t108 = t202 * t122;
t109 = t202 * t124;
t233 = -t129 * t108 - t109 * t127;
t177 = qJD(5) * t128;
t178 = qJD(5) * t126;
t137 = t111 * qJD(3);
t32 = -qJD(1) * t137 + qJD(4) * t49;
t72 = pkin(4) * t98 + pkin(8) * t97;
t3 = t126 * t72 + t128 * t32 + t57 * t177 - t178 * t47;
t4 = -qJD(5) * t18 - t126 * t32 + t128 * t72;
t168 = -t4 * t126 + t3 * t128;
t1 = qJ(6) * t98 + qJD(6) * t237 + t3;
t2 = -pkin(5) * t98 - t4;
t170 = t1 * t128 + t2 * t126;
t74 = pkin(4) * t111 - pkin(8) * t112 + t139;
t80 = -t108 * t127 + t109 * t129;
t196 = t126 * t74 + t128 * t80;
t52 = qJD(4) * t233 - t137;
t82 = pkin(4) * t107 + pkin(8) * t106;
t9 = -qJD(5) * t196 - t126 * t52 + t128 * t82;
t173 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t174 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t175 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t227 = t174 * t237 - t173 * t140 - t175 * t89 - t11 * mrSges(7,1) - t18 * mrSges(6,2) - t50 * mrSges(5,3) - Ifges(6,6) * t222 - Ifges(7,6) * t223 - Ifges(5,6) * qJD(4) - t105 * Ifges(5,4) + t243 + t102 * mrSges(5,1) + t12 * mrSges(7,3) + t17 * mrSges(6,1) - t242 * t221 - (Ifges(6,3) + Ifges(7,2)) * t218;
t226 = t63 / 0.2e1;
t225 = -t64 / 0.2e1;
t224 = t64 / 0.2e1;
t219 = t98 / 0.2e1;
t214 = -t128 / 0.2e1;
t212 = mrSges(6,3) * t140;
t211 = mrSges(6,3) * t89;
t138 = t112 * qJD(3);
t33 = qJD(1) * t138 + qJD(4) * t50;
t205 = t33 * t233;
t41 = mrSges(6,1) * t98 - mrSges(6,3) * t63;
t42 = -t98 * mrSges(7,1) + t63 * mrSges(7,2);
t201 = -t41 + t42;
t43 = -mrSges(6,2) * t98 - mrSges(6,3) * t64;
t44 = -mrSges(7,2) * t64 + mrSges(7,3) * t98;
t200 = t43 + t44;
t81 = pkin(4) * t105 + pkin(8) * t104;
t25 = t126 * t81 + t128 * t49;
t199 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t140 + mrSges(6,2) * t89 + t105 * mrSges(5,3);
t65 = mrSges(7,2) * t140 + mrSges(7,3) * t237;
t66 = -mrSges(6,2) * t237 + t212;
t198 = t65 + t66;
t67 = mrSges(6,1) * t237 - t211;
t68 = -mrSges(7,1) * t237 + mrSges(7,2) * t89;
t197 = t67 - t68;
t186 = t111 * t33;
t182 = t106 * t126;
t181 = t106 * t128;
t55 = -mrSges(7,1) * t140 - mrSges(7,3) * t89;
t176 = t55 + t199;
t172 = t98 * mrSges(5,1) - t97 * mrSges(5,2);
t169 = -t1 * t126 + t2 * t128;
t167 = -t126 * t3 - t128 * t4;
t166 = mrSges(6,1) * t128 - mrSges(6,2) * t126;
t164 = mrSges(7,1) * t128 + mrSges(7,3) * t126;
t157 = Ifges(6,2) * t128 + t195;
t151 = -Ifges(7,3) * t128 + t193;
t150 = pkin(5) * t128 + qJ(6) * t126;
t147 = t11 * t126 + t12 * t128;
t146 = -t122 * (-t114 * t122 + t118) + t124 * t96;
t144 = t126 * t17 - t128 * t18;
t24 = -t126 * t49 + t128 * t81;
t29 = -t126 * t80 + t128 * t74;
t8 = t126 * t82 + t128 * t52 + t74 * t177 - t178 * t80;
t136 = t126 * t201 + t128 * t200;
t135 = -t126 * t198 - t128 * t197;
t134 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t53 = qJD(4) * t80 + t138;
t113 = -pkin(4) - t150;
t94 = Ifges(7,2) * t98;
t93 = Ifges(6,3) * t98;
t90 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t104;
t62 = Ifges(7,4) * t63;
t61 = Ifges(6,5) * t63;
t60 = Ifges(6,6) * t64;
t59 = Ifges(7,6) * t64;
t54 = pkin(5) * t89 - qJ(6) * t140;
t34 = t112 * t149 - t233;
t27 = -pkin(5) * t111 - t29;
t26 = qJ(6) * t111 + t196;
t23 = mrSges(6,1) * t64 + mrSges(6,2) * t63;
t22 = mrSges(7,1) * t64 - mrSges(7,3) * t63;
t20 = -pkin(5) * t105 - t24;
t19 = qJ(6) * t105 + t25;
t14 = t63 * Ifges(6,4) - t64 * Ifges(6,2) + t98 * Ifges(6,6);
t13 = t63 * Ifges(7,5) + t98 * Ifges(7,6) + t64 * Ifges(7,3);
t10 = -t149 * t106 + (qJD(5) * t150 - qJD(6) * t128) * t112 + t53;
t7 = -pkin(5) * t107 - t9;
t6 = pkin(5) * t64 - qJ(6) * t63 - qJD(6) * t89 + t33;
t5 = qJ(6) * t107 + qJD(6) * t111 + t8;
t15 = [t34 * t22 + t29 * t41 + t27 * t42 + t196 * t43 + t26 * t44 + t10 * t55 + t5 * t65 + t8 * t66 + t9 * t67 + t7 * t68 - t233 * t23 + t52 * t90 + t139 * t172 + t199 * t53 + (t233 * t97 - t80 * t98) * mrSges(5,3) + m(5) * (t32 * t80 - t49 * t53 + t50 * t52 - t205) + m(6) * (t17 * t9 + t18 * t8 + t196 * t3 + t29 * t4 + t46 * t53 - t205) + m(7) * (t1 * t26 + t10 * t21 + t11 * t7 + t12 * t5 + t2 * t27 + t34 * t6) + (m(4) * (t115 * t171 + t146) + 0.2e1 * t241) * qJD(3) + (t243 + t227) * t107 - (t244 + t240) * t106 + (-t32 * mrSges(5,3) + t61 / 0.2e1 - t60 / 0.2e1 + t93 / 0.2e1 + t62 / 0.2e1 + t94 / 0.2e1 + t59 / 0.2e1 + Ifges(5,4) * t97 + t173 * t64 - t175 * t63 + (Ifges(5,2) + t174) * t98 + t134) * t111 + (-Ifges(5,1) * t97 - Ifges(5,4) * t98 + t13 * t215 + t6 * t163 + t152 * t224 + t158 * t225 + (mrSges(5,3) + t165) * t33 + t167 * mrSges(6,3) + t169 * mrSges(7,2) + (-t147 * mrSges(7,2) + mrSges(6,3) * t144 + t151 * t223 + t157 * t222 + t164 * t21 + t166 * t46 + t214 * t38 + t239 * t218 + t238 * t221) * qJD(5) + t230 * t226 + t231 * t219 + (qJD(5) * t235 + t14) * t216 + (qJD(5) * t35 + t236) * t213) * t112; (-mrSges(5,3) * t97 + t22 + t23) * t111 + t176 * t107 - (-t126 * t197 + t128 * t198 + t90) * t106 + m(5) * (-t106 * t50 - t107 * t49 + t186) + m(6) * (t107 * t46 + t17 * t182 - t18 * t181 + t186) + m(7) * (t107 * t21 - t11 * t182 + t111 * t6 - t12 * t181) + (-t98 * mrSges(5,3) + t135 * qJD(5) + m(5) * t32 + m(6) * (-t17 * t177 - t178 * t18 + t168) + m(7) * (t11 * t177 - t12 * t178 + t170) + t136) * t112; t104 * t90 - t176 * t105 + (t198 * t237 - t201) * t128 + (-t197 * t237 + t200) * t126 - m(5) * (-t104 * t50 - t105 * t49) + t172 + (-m(4) * t146 - t241) * qJD(1) + (-t21 * t105 + t147 * t237 - t169) * m(7) + (-t105 * t46 - t144 * t237 - t167) * m(6); (t136 + (-m(6) * t145 + m(7) * t148 + t135) * qJD(5) + m(6) * t168 + m(7) * t170) * pkin(8) + t151 * t224 + t157 * t225 + t14 * t213 + t13 * t214 - t227 * t105 + t236 * t215 + (-t11 * t20 + t113 * t6 - t12 * t19 + t21 * t234) * m(7) + t234 * t55 + (-pkin(4) * t33 - t17 * t24 - t18 * t25 - t46 * t50) * m(6) - ((t247 - Ifges(5,1) / 0.2e1) * t105 - t240) * t104 + t131 * qJD(5) - pkin(4) * t23 - t32 * mrSges(5,2) - t19 * t65 - t25 * t66 - t24 * t67 - t20 * t68 + t238 * t226 + t239 * t219 - t49 * t90 - Ifges(5,5) * t97 - Ifges(5,6) * t98 - t6 * t164 + (-mrSges(5,1) - t166) * t33 + t168 * mrSges(6,3) + t113 * t22 + t170 * mrSges(7,2) - t199 * t50; (-t11 * t140 + t12 * t89) * mrSges(7,2) + t38 * t220 + (Ifges(7,3) * t89 + t209) * t223 + t94 + t93 + t62 + t61 + t59 - t60 + t134 - pkin(5) * t42 + qJ(6) * t44 - t54 * t55 + qJD(6) * t65 - t21 * (mrSges(7,1) * t89 - mrSges(7,3) * t140) - t46 * (mrSges(6,1) * t89 + mrSges(6,2) * t140) + (t197 + t211) * t18 + (-t198 + t212) * t17 + (t242 * t140 - t248 * t89) * t218 + (-pkin(5) * t2 + qJ(6) * t1 - t11 * t18 + t12 * t232 - t21 * t54) * m(7) + (-Ifges(6,2) * t89 + t235 + t85) * t222 + (t246 * t140 - t210 + t35 + t84) * t221; -t237 * t65 + t89 * t55 + 0.2e1 * (t2 / 0.2e1 + t12 * t218 + t21 * t220) * m(7) + t42;];
tauc  = t15(:);
