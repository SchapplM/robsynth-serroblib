% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:43
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:42:55
% EndTime: 2018-11-23 15:43:00
% DurationCPUTime: 5.00s
% Computational Cost: add. (5872->435), mult. (13725->603), div. (0->0), fcn. (9854->8), ass. (0->196)
t172 = sin(pkin(9));
t174 = cos(pkin(9));
t245 = sin(qJ(4));
t246 = cos(qJ(4));
t144 = t245 * t172 - t246 * t174;
t280 = Ifges(5,1) / 0.2e1;
t137 = t144 * qJD(1);
t129 = t137 * qJD(4);
t251 = -t129 / 0.2e1;
t145 = t172 * t246 + t174 * t245;
t135 = t145 * qJD(1);
t279 = -t135 / 0.2e1;
t278 = Ifges(5,2) + Ifges(6,3);
t277 = t144 * qJD(3);
t171 = sin(pkin(10));
t173 = cos(pkin(10));
t118 = qJD(4) * t171 - t137 * t173;
t176 = sin(qJ(6));
t177 = cos(qJ(6));
t200 = t173 * qJD(4) + t137 * t171;
t276 = -t118 * t176 + t177 * t200;
t70 = t118 * t177 + t176 * t200;
t128 = qJD(4) * t135;
t185 = t171 * t176 - t173 * t177;
t32 = qJD(6) * t276 + t128 * t185;
t260 = t32 / 0.2e1;
t146 = t171 * t177 + t173 * t176;
t33 = -qJD(6) * t70 + t128 * t146;
t259 = t33 / 0.2e1;
t275 = Ifges(5,4) * t279;
t274 = t135 / 0.2e1;
t273 = t137 * t280;
t272 = Ifges(6,6) * t251;
t238 = pkin(8) + qJ(5);
t151 = t238 * t171;
t152 = t238 * t173;
t115 = -t151 * t176 + t152 * t177;
t242 = pkin(8) * t173;
t109 = -pkin(4) * t137 + qJ(5) * t135;
t175 = -pkin(1) - qJ(3);
t155 = qJD(1) * t175 + qJD(2);
t201 = -pkin(7) * qJD(1) + t155;
t130 = t201 * t172;
t123 = t245 * t130;
t131 = t201 * t174;
t96 = t131 * t246 - t123;
t42 = t173 * t109 - t171 * t96;
t27 = -pkin(5) * t137 + t135 * t242 + t42;
t215 = t135 * t171;
t43 = t171 * t109 + t173 * t96;
t36 = pkin(8) * t215 + t43;
t271 = -qJD(5) * t146 - qJD(6) * t115 + t176 * t36 - t177 * t27;
t114 = -t151 * t177 - t152 * t176;
t270 = -qJD(5) * t185 + qJD(6) * t114 - t176 * t27 - t177 * t36;
t11 = -t33 * mrSges(7,1) + t32 * mrSges(7,2);
t216 = t128 * t173;
t217 = t128 * t171;
t86 = -mrSges(6,1) * t217 - mrSges(6,2) * t216;
t269 = t11 + t86;
t139 = t146 * qJD(6);
t203 = qJD(4) * t245;
t204 = qJD(4) * t246;
t141 = -t172 * t203 + t174 * t204;
t268 = -t146 * qJD(1) - t139 * t145 - t141 * t185;
t138 = t185 * qJD(6);
t267 = t185 * qJD(1) + t138 * t145 - t141 * t146;
t92 = t146 * t135;
t221 = -t92 - t139;
t93 = t185 * t135;
t220 = -t93 - t138;
t239 = -pkin(7) + t175;
t148 = t239 * t172;
t149 = t239 * t174;
t112 = t245 * t148 - t246 * t149;
t211 = t172 ^ 2 + t174 ^ 2;
t199 = qJD(1) * t211;
t183 = t145 * qJD(3);
t181 = -qJD(1) * t183 + t131 * t204;
t58 = (-t123 + qJD(5)) * qJD(4) + t181;
t210 = qJD(1) * qJD(2);
t66 = -pkin(4) * t129 + qJ(5) * t128 + qJD(5) * t137 + t210;
t18 = -t171 * t58 + t173 * t66;
t14 = -pkin(5) * t129 + pkin(8) * t216 + t18;
t19 = t171 * t66 + t173 * t58;
t15 = pkin(8) * t217 + t19;
t97 = t130 * t246 + t131 * t245;
t85 = qJD(4) * qJ(5) + t97;
t170 = qJD(1) * qJ(2);
t165 = qJD(3) + t170;
t167 = t172 * pkin(3);
t150 = qJD(1) * t167 + t165;
t87 = pkin(4) * t135 + qJ(5) * t137 + t150;
t38 = -t171 * t85 + t173 * t87;
t17 = pkin(5) * t135 - pkin(8) * t118 + t38;
t39 = t171 * t87 + t173 * t85;
t26 = pkin(8) * t200 + t39;
t5 = t17 * t177 - t176 * t26;
t1 = qJD(6) * t5 + t14 * t176 + t15 * t177;
t6 = t17 * t176 + t177 * t26;
t2 = -qJD(6) * t6 + t14 * t177 - t15 * t176;
t264 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t32 + Ifges(7,6) * t33;
t187 = t171 * t39 + t173 * t38;
t223 = t171 * Ifges(6,2);
t234 = Ifges(6,4) * t173;
t191 = -t223 + t234;
t235 = Ifges(6,4) * t171;
t192 = Ifges(6,1) * t173 - t235;
t193 = mrSges(6,1) * t171 + mrSges(6,2) * t173;
t247 = t173 / 0.2e1;
t248 = -t171 / 0.2e1;
t80 = -qJD(4) * pkin(4) + qJD(5) - t96;
t263 = t191 * t200 / 0.2e1 + t192 * t118 / 0.2e1 + t80 * t193 + t150 * mrSges(5,2) - t187 * mrSges(6,3) + t275 + Ifges(5,5) * qJD(4) - t273 + (t118 * Ifges(6,4) + Ifges(6,2) * t200 + Ifges(6,6) * t135) * t248 + (t118 * Ifges(6,1) + Ifges(6,4) * t200 + Ifges(6,5) * t135) * t247;
t262 = Ifges(7,4) * t260 + Ifges(7,2) * t259 + Ifges(7,6) * t251;
t261 = Ifges(7,1) * t260 + Ifges(7,4) * t259 + Ifges(7,5) * t251;
t258 = -t128 * t191 / 0.2e1 + t272;
t257 = -t276 / 0.2e1;
t256 = t276 / 0.2e1;
t255 = -t70 / 0.2e1;
t254 = t70 / 0.2e1;
t133 = qJD(6) + t135;
t250 = -t133 / 0.2e1;
t249 = t133 / 0.2e1;
t244 = m(3) * qJ(2);
t243 = Ifges(7,4) * t70;
t140 = -t172 * t204 - t174 * t203;
t75 = pkin(4) * t141 - qJ(5) * t140 + qJD(5) * t144 + qJD(2);
t81 = -qJD(4) * t112 - t183;
t35 = t171 * t75 + t173 * t81;
t237 = m(4) * qJD(3);
t236 = Ifges(5,4) * t137;
t233 = Ifges(6,5) * t173;
t232 = Ifges(6,6) * t171;
t64 = -qJD(1) * t277 + qJD(4) * t97;
t231 = t112 * t64;
t228 = t129 * Ifges(6,5);
t224 = t144 * t64;
t161 = qJ(2) + t167;
t108 = pkin(4) * t145 + qJ(5) * t144 + t161;
t113 = t148 * t246 + t149 * t245;
t52 = t171 * t108 + t173 * t113;
t222 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t200 + mrSges(6,2) * t118 - mrSges(5,3) * t137;
t214 = t140 * t171;
t213 = t144 * t171;
t194 = mrSges(4,1) * t172 + mrSges(4,2) * t174;
t212 = mrSges(5,1) * t135 - mrSges(5,2) * t137 + qJD(1) * t194;
t28 = -mrSges(7,1) * t276 + mrSges(7,2) * t70;
t209 = -t28 - t222;
t208 = Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t34 = -t171 * t81 + t173 * t75;
t202 = -t129 * mrSges(5,1) - t128 * mrSges(5,2);
t51 = t173 * t108 - t113 * t171;
t190 = -t232 + t233;
t189 = t171 * t19 + t173 * t18;
t188 = -t171 * t18 + t173 * t19;
t186 = t171 * t38 - t173 * t39;
t37 = pkin(5) * t145 + t144 * t242 + t51;
t40 = pkin(8) * t213 + t52;
t12 = -t176 * t40 + t177 * t37;
t13 = t176 * t37 + t177 * t40;
t124 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t135;
t90 = -mrSges(6,2) * t135 + mrSges(6,3) * t200;
t91 = mrSges(6,1) * t135 - mrSges(6,3) * t118;
t184 = -t171 * t91 + t173 * t90 + t124;
t103 = t146 * t144;
t63 = -t130 * t203 + t181;
t182 = t140 * t96 + t141 * t97 + t145 * t63 + t224;
t82 = qJD(4) * t113 - t277;
t179 = t150 * mrSges(5,1) + t38 * mrSges(6,1) + t5 * mrSges(7,1) - Ifges(5,6) * qJD(4) + t236 / 0.2e1 + t133 * Ifges(7,3) + t70 * Ifges(7,5) + t276 * Ifges(7,6) + t118 * Ifges(6,5) + t200 * Ifges(6,6) - t39 * mrSges(6,2) - t6 * mrSges(7,2) + t278 * t274;
t178 = qJD(1) ^ 2;
t163 = -pkin(5) * t173 - pkin(4);
t127 = Ifges(7,3) * t129;
t105 = t185 * t144;
t104 = t185 * t145;
t102 = t146 * t145;
t95 = -mrSges(6,1) * t129 + mrSges(6,3) * t216;
t94 = mrSges(6,2) * t129 + mrSges(6,3) * t217;
t83 = -pkin(5) * t213 + t112;
t67 = Ifges(7,4) * t276;
t65 = -pkin(5) * t215 + t97;
t60 = -t128 * t192 - t228;
t57 = pkin(5) * t214 + t82;
t50 = -pkin(5) * t200 + t80;
t49 = mrSges(7,1) * t133 - mrSges(7,3) * t70;
t48 = -mrSges(7,2) * t133 + mrSges(7,3) * t276;
t47 = -t138 * t144 - t140 * t146;
t45 = qJD(6) * t103 - t140 * t185;
t41 = -pkin(5) * t217 + t64;
t25 = -pkin(8) * t214 + t35;
t24 = mrSges(7,2) * t129 + mrSges(7,3) * t33;
t23 = -mrSges(7,1) * t129 - mrSges(7,3) * t32;
t22 = Ifges(7,1) * t70 + Ifges(7,5) * t133 + t67;
t21 = Ifges(7,2) * t276 + Ifges(7,6) * t133 + t243;
t16 = pkin(5) * t141 - t140 * t242 + t34;
t4 = -qJD(6) * t13 + t16 * t177 - t176 * t25;
t3 = qJD(6) * t12 + t16 * t176 + t177 * t25;
t7 = [(-t127 / 0.2e1 - t19 * mrSges(6,2) + t18 * mrSges(6,1) - t190 * t128 - (Ifges(7,3) / 0.2e1 + t278) * t129 + t264) * t145 + 0.2e1 * qJD(3) * mrSges(4,3) * t199 + ((mrSges(5,1) * t145 - mrSges(5,2) * t144 + (2 * mrSges(3,3)) + t194 + 0.2e1 * t244) * qJD(1) + m(5) * (qJD(1) * t161 + t150) + m(4) * (t165 + t170) + t212) * qJD(2) + (t171 * t258 - t173 * t60 / 0.2e1 - t64 * t193 + t129 * t190 / 0.2e1 + t189 * mrSges(6,3) - (-Ifges(5,1) - t173 ^ 2 * Ifges(6,1) / 0.2e1 + (t234 - t223 / 0.2e1) * t171) * t128) * t144 + (-t112 * t128 + t113 * t129 - t182) * mrSges(5,3) + (t140 * t279 + t137 * t141 / 0.2e1 - t144 * t129 + t145 * t128) * Ifges(5,4) + (-t155 * t211 - t175 * t199) * t237 + m(7) * (t1 * t13 + t12 * t2 + t3 * t6 + t4 * t5 + t41 * t83 + t50 * t57) + (t1 * t103 - t105 * t2 - t45 * t5 + t47 * t6) * mrSges(7,3) + m(6) * (t18 * t51 + t19 * t52 + t34 * t38 + t35 * t39 + t80 * t82 + t231) + m(5) * (t113 * t63 + t81 * t97 - t82 * t96 + t231) + t222 * t82 + (t135 * t208 + t179) * t141 + t161 * t202 + (Ifges(7,4) * t45 + Ifges(7,2) * t47) * t256 + (Ifges(7,4) * t105 + Ifges(7,2) * t103) * t259 + (Ifges(7,1) * t105 + Ifges(7,4) * t103) * t260 + t105 * t261 + t103 * t262 + (t190 * t274 + t263 - t273) * t140 + (Ifges(7,5) * t45 + Ifges(7,6) * t47) * t249 + (Ifges(7,5) * t105 + Ifges(7,6) * t103) * t251 + (Ifges(7,1) * t45 + Ifges(7,4) * t47) * t254 + t12 * t23 + t13 * t24 + t45 * t22 / 0.2e1 + t47 * t21 / 0.2e1 + t3 * t48 + t4 * t49 + t50 * (-mrSges(7,1) * t47 + mrSges(7,2) * t45) + t57 * t28 + t83 * t11 + t35 * t90 + t34 * t91 + t52 * t94 + t51 * t95 + t41 * (-mrSges(7,1) * t103 + mrSges(7,2) * t105) + t112 * t86 + t81 * t124; -t102 * t23 - t104 * t24 + t267 * t49 + t268 * t48 + (-mrSges(3,3) - t244) * t178 + (t129 * mrSges(5,3) - t171 * t95 + t173 * t94) * t145 + (-t128 * mrSges(5,3) + t269) * t144 + t184 * t141 + t209 * t140 + m(5) * t182 + m(6) * (-t140 * t80 - t141 * t186 + t145 * t188 + t224) + (-m(4) * t165 - m(5) * t150 - m(6) * t187 - t171 * t90 - t173 * t91 - t211 * t237 - t212) * qJD(1) + (-t1 * t104 - t102 * t2 - t140 * t50 + t144 * t41 + t267 * t5 + t268 * t6) * m(7); -t185 * t23 + t146 * t24 + t171 * t94 + t173 * t95 + t221 * t49 + t220 * t48 + (m(4) + m(5)) * t210 - t211 * t178 * mrSges(4,3) - t209 * t137 + t184 * t135 - m(5) * (-t135 * t97 + t137 * t96) + m(4) * t155 * t199 + t202 + (t1 * t146 + t137 * t50 - t185 * t2 + t220 * t6 + t221 * t5) * m(7) + (-t135 * t186 + t137 * t80 + t189) * m(6); (t1 * t115 + t114 * t2 + t163 * t41 + t270 * t6 + t271 * t5 - t50 * t65) * m(7) + t271 * t49 + (-t64 * mrSges(6,1) + t19 * mrSges(6,3) + qJ(5) * t94 + qJD(5) * t90 + t258 + t272) * t173 + t270 * t48 + (-qJ(5) * t95 - qJD(5) * t91 - t18 * mrSges(6,3) + t60 / 0.2e1 + t64 * mrSges(6,2) - t228 / 0.2e1) * t171 + (-Ifges(7,4) * t138 - Ifges(7,2) * t139) * t256 + (-Ifges(7,5) * t138 - Ifges(7,6) * t139) * t249 + (-Ifges(7,1) * t138 - Ifges(7,4) * t139) * t254 + (-t92 / 0.2e1 - t139 / 0.2e1) * t21 + (-t93 / 0.2e1 - t138 / 0.2e1) * t22 + (-t1 * t185 - t146 * t2 - t220 * t5 + t221 * t6) * mrSges(7,3) + (Ifges(7,4) * t146 - Ifges(7,2) * t185) * t259 + (Ifges(7,1) * t146 - Ifges(7,4) * t185) * t260 + (Ifges(7,5) * t146 - Ifges(7,6) * t185) * t251 + t41 * (mrSges(7,1) * t185 + mrSges(7,2) * t146) - t185 * t262 + (-pkin(4) * t64 + qJ(5) * t188 - qJD(5) * t186 - t38 * t42 - t39 * t43 - t80 * t97) * m(6) - ((Ifges(6,1) * t171 + t234) * t247 + (Ifges(6,2) * t173 + t235) * t248 + Ifges(5,5)) * t128 + (-mrSges(7,1) * t221 + mrSges(7,2) * t220) * t50 - t222 * t97 + (Ifges(7,1) * t93 + Ifges(7,4) * t92) * t255 + (Ifges(7,4) * t93 + Ifges(7,2) * t92) * t257 + t146 * t261 + (-t96 * mrSges(5,3) + t275 + (t233 / 0.2e1 - t232 / 0.2e1) * t135 - (t280 - t208) * t137 + t263) * t135 - (-t236 / 0.2e1 - t179 + t97 * mrSges(5,3)) * t137 + (Ifges(7,5) * t93 + Ifges(7,6) * t92) * t250 - t63 * mrSges(5,2) - t64 * mrSges(5,1) - t65 * t28 - pkin(4) * t86 - t43 * t90 - t42 * t91 + t114 * t23 + t115 * t24 - t96 * t124 + Ifges(5,6) * t129 + t163 * t11; t118 * t91 - t200 * t90 - t276 * t48 + t70 * t49 + (-t276 * t6 + t5 * t70 + t41) * m(7) + (t118 * t38 - t200 * t39 + t64) * m(6) + t269; -t127 - t50 * (mrSges(7,1) * t70 + mrSges(7,2) * t276) + (Ifges(7,1) * t276 - t243) * t255 + t21 * t254 + (Ifges(7,5) * t276 - Ifges(7,6) * t70) * t250 - t5 * t48 + t6 * t49 + (t276 * t5 + t6 * t70) * mrSges(7,3) + (-Ifges(7,2) * t70 + t22 + t67) * t257 + t264;];
tauc  = t7(:);
