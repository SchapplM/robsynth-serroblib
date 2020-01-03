% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:36
% EndTime: 2019-12-31 19:31:51
% DurationCPUTime: 6.68s
% Computational Cost: add. (4512->443), mult. (12028->637), div. (0->0), fcn. (8664->8), ass. (0->209)
t171 = sin(pkin(8));
t176 = cos(qJ(2));
t213 = cos(pkin(8));
t190 = t213 * t176;
t174 = sin(qJ(2));
t203 = qJD(1) * t174;
t136 = -qJD(1) * t190 + t171 * t203;
t244 = -t136 / 0.2e1;
t276 = Ifges(4,6) * qJD(2) / 0.2e1;
t152 = t171 * t176 + t174 * t213;
t138 = t152 * qJD(1);
t170 = sin(pkin(9));
t172 = cos(pkin(9));
t116 = qJD(2) * t170 + t138 * t172;
t173 = sin(qJ(5));
t175 = cos(qJ(5));
t188 = t172 * qJD(2) - t138 * t170;
t275 = -t116 * t173 + t175 * t188;
t65 = t116 * t175 + t173 * t188;
t167 = -pkin(2) * t176 - pkin(1);
t204 = qJD(1) * t167;
t158 = qJD(3) + t204;
t222 = Ifges(5,6) * t188;
t224 = Ifges(5,5) * t116;
t227 = Ifges(4,4) * t138;
t231 = -qJ(3) - pkin(6);
t161 = t231 * t174;
t155 = qJD(1) * t161;
t149 = qJD(2) * pkin(2) + t155;
t162 = t231 * t176;
t156 = qJD(1) * t162;
t191 = t213 * t156;
t108 = t171 * t149 - t191;
t101 = qJD(2) * qJ(4) + t108;
t79 = t136 * pkin(3) - t138 * qJ(4) + t158;
t38 = -t101 * t170 + t172 * t79;
t39 = t172 * t101 + t170 * t79;
t20 = pkin(4) * t136 - pkin(7) * t116 + t38;
t28 = pkin(7) * t188 + t39;
t5 = -t173 * t28 + t175 * t20;
t6 = t173 * t20 + t175 * t28;
t274 = t39 * mrSges(5,2) + t6 * mrSges(6,2) + Ifges(4,2) * t244 + t227 / 0.2e1 + t276 - t158 * mrSges(4,1) - t38 * mrSges(5,1) - t5 * mrSges(6,1) - t222 / 0.2e1 - t224 / 0.2e1;
t178 = -t171 * t174 + t190;
t139 = t178 * qJD(2);
t130 = qJD(1) * t139;
t179 = t170 * t173 - t172 * t175;
t32 = qJD(5) * t275 - t130 * t179;
t253 = t32 / 0.2e1;
t153 = t170 * t175 + t172 * t173;
t33 = -qJD(5) * t65 - t130 * t153;
t252 = t33 / 0.2e1;
t137 = t152 * qJD(2);
t129 = qJD(1) * t137;
t247 = t129 / 0.2e1;
t236 = pkin(2) * t171;
t164 = qJ(4) + t236;
t232 = pkin(7) + t164;
t147 = t232 * t170;
t148 = t232 * t172;
t103 = -t147 * t173 + t148 * t175;
t207 = t136 * t172;
t142 = t171 * t156;
t110 = t155 * t213 + t142;
t199 = pkin(2) * t203;
t91 = pkin(3) * t138 + qJ(4) * t136 + t199;
t43 = -t110 * t170 + t172 * t91;
t27 = pkin(4) * t138 + pkin(7) * t207 + t43;
t208 = t136 * t170;
t44 = t172 * t110 + t170 * t91;
t37 = pkin(7) * t208 + t44;
t273 = -qJD(4) * t153 - qJD(5) * t103 + t173 * t37 - t175 * t27;
t102 = -t147 * t175 - t148 * t173;
t272 = -qJD(4) * t179 + qJD(5) * t102 - t173 * t27 - t175 * t37;
t140 = t179 * qJD(5);
t86 = t179 * t136;
t271 = t140 + t86;
t141 = t153 * qJD(5);
t85 = t153 * t136;
t270 = t141 + t85;
t230 = mrSges(4,3) * t138;
t214 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t188 + mrSges(5,2) * t116 + t230;
t269 = Ifges(4,5) * qJD(2);
t80 = -mrSges(5,2) * t136 + mrSges(5,3) * t188;
t81 = mrSges(5,1) * t136 - mrSges(5,3) * t116;
t265 = -t170 * t81 + t172 * t80;
t200 = qJD(1) * qJD(2);
t195 = t174 * t200;
t187 = pkin(2) * t195;
t57 = pkin(3) * t129 - qJ(4) * t130 - qJD(4) * t138 + t187;
t192 = qJD(2) * t231;
t134 = qJD(3) * t176 + t174 * t192;
t125 = t134 * qJD(1);
t135 = -t174 * qJD(3) + t176 * t192;
t177 = t135 * qJD(1);
t75 = t213 * t125 + t171 * t177;
t72 = qJD(2) * qJD(4) + t75;
t23 = -t170 * t72 + t172 * t57;
t24 = t170 * t57 + t172 * t72;
t264 = -t170 * t23 + t172 * t24;
t202 = qJD(1) * t176;
t211 = Ifges(3,6) * qJD(2);
t228 = Ifges(3,4) * t174;
t263 = pkin(6) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t202) + t211 / 0.2e1 + (t176 * Ifges(3,2) + t228) * qJD(1) / 0.2e1;
t209 = t130 * t172;
t14 = pkin(4) * t129 - pkin(7) * t209 + t23;
t210 = t130 * t170;
t15 = -pkin(7) * t210 + t24;
t1 = qJD(5) * t5 + t14 * t173 + t15 * t175;
t2 = -qJD(5) * t6 + t14 * t175 - t15 * t173;
t262 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t32 + Ifges(6,6) * t33;
t223 = Ifges(5,2) * t170;
t225 = Ifges(5,4) * t172;
t183 = -t223 + t225;
t226 = Ifges(5,4) * t170;
t184 = Ifges(5,1) * t172 - t226;
t185 = mrSges(5,1) * t170 + mrSges(5,2) * t172;
t240 = t172 / 0.2e1;
t241 = -t170 / 0.2e1;
t107 = t149 * t213 + t142;
t94 = -qJD(2) * pkin(3) + qJD(4) - t107;
t260 = (t116 * Ifges(5,1) + Ifges(5,4) * t188 + Ifges(5,5) * t136) * t240 + (t116 * Ifges(5,4) + Ifges(5,2) * t188 + t136 * Ifges(5,6)) * t241 + t158 * mrSges(4,2) + t94 * t185 + t116 * t184 / 0.2e1 + t188 * t183 / 0.2e1;
t259 = -0.2e1 * pkin(1);
t257 = Ifges(6,4) * t253 + Ifges(6,2) * t252 + Ifges(6,6) * t247;
t256 = Ifges(6,1) * t253 + Ifges(6,4) * t252 + Ifges(6,5) * t247;
t133 = qJD(5) + t136;
t239 = Ifges(6,4) * t65;
t18 = Ifges(6,2) * t275 + Ifges(6,6) * t133 + t239;
t255 = t18 / 0.2e1;
t62 = Ifges(6,4) * t275;
t19 = Ifges(6,1) * t65 + Ifges(6,5) * t133 + t62;
t254 = t19 / 0.2e1;
t251 = -t275 / 0.2e1;
t250 = t275 / 0.2e1;
t249 = -t65 / 0.2e1;
t248 = t65 / 0.2e1;
t246 = -t133 / 0.2e1;
t245 = t133 / 0.2e1;
t243 = t136 / 0.2e1;
t242 = -t138 / 0.2e1;
t238 = Ifges(6,5) * t65;
t237 = Ifges(6,6) * t275;
t235 = pkin(6) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t203);
t233 = pkin(7) * t172;
t201 = qJD(2) * t174;
t68 = pkin(2) * t201 + pkin(3) * t137 - qJ(4) * t139 - qJD(4) * t152;
t93 = t134 * t213 + t171 * t135;
t35 = t170 * t68 + t172 * t93;
t229 = Ifges(4,1) * t138;
t221 = Ifges(6,3) * t133;
t112 = -t213 * t161 - t162 * t171;
t74 = t125 * t171 - t213 * t177;
t220 = t112 * t74;
t219 = t136 * mrSges(4,3);
t106 = -pkin(3) * t178 - t152 * qJ(4) + t167;
t113 = t171 * t161 - t162 * t213;
t52 = t170 * t106 + t172 * t113;
t212 = Ifges(3,5) * qJD(2);
t206 = t139 * t170;
t205 = t152 * t170;
t78 = mrSges(5,1) * t210 + mrSges(5,2) * t209;
t196 = t213 * pkin(2);
t11 = -t33 * mrSges(6,1) + t32 * mrSges(6,2);
t194 = t176 * t200;
t34 = -t170 * t93 + t172 * t68;
t189 = t129 * mrSges(4,1) + t130 * mrSges(4,2);
t51 = t172 * t106 - t113 * t170;
t92 = t134 * t171 - t213 * t135;
t109 = t155 * t171 - t191;
t166 = -t196 - pkin(3);
t182 = Ifges(5,5) * t172 - Ifges(5,6) * t170;
t181 = t170 * t24 + t172 * t23;
t180 = t170 * t38 - t172 * t39;
t36 = -pkin(4) * t178 - t152 * t233 + t51;
t40 = -pkin(7) * t205 + t52;
t12 = -t173 * t40 + t175 * t36;
t13 = t173 * t36 + t175 * t40;
t168 = Ifges(3,4) * t202;
t157 = -t172 * pkin(4) + t166;
t146 = Ifges(3,1) * t203 + t168 + t212;
t132 = Ifges(4,4) * t136;
t128 = Ifges(6,3) * t129;
t123 = -qJD(2) * mrSges(4,2) - t219;
t104 = mrSges(4,1) * t136 + mrSges(4,2) * t138;
t98 = -t132 + t229 + t269;
t96 = t179 * t152;
t95 = t153 * t152;
t89 = pkin(4) * t205 + t112;
t88 = mrSges(5,1) * t129 - mrSges(5,3) * t209;
t87 = -mrSges(5,2) * t129 - mrSges(5,3) * t210;
t73 = -pkin(4) * t208 + t109;
t59 = pkin(4) * t206 + t92;
t58 = -pkin(4) * t188 + t94;
t56 = t129 * Ifges(5,5) + t130 * t184;
t55 = t129 * Ifges(5,6) + t130 * t183;
t48 = Ifges(5,3) * t136 + t222 + t224;
t47 = pkin(4) * t210 + t74;
t46 = mrSges(6,1) * t133 - mrSges(6,3) * t65;
t45 = -mrSges(6,2) * t133 + mrSges(6,3) * t275;
t42 = -t139 * t153 + t140 * t152;
t41 = -t139 * t179 - t141 * t152;
t26 = -mrSges(6,1) * t275 + mrSges(6,2) * t65;
t25 = -pkin(7) * t206 + t35;
t22 = -mrSges(6,2) * t129 + mrSges(6,3) * t33;
t21 = mrSges(6,1) * t129 - mrSges(6,3) * t32;
t17 = t221 + t237 + t238;
t16 = pkin(4) * t137 - t139 * t233 + t34;
t4 = -qJD(5) * t13 + t16 * t175 - t173 * t25;
t3 = qJD(5) * t12 + t16 * t173 + t175 * t25;
t7 = [(-t1 * t95 + t2 * t96 - t41 * t5 + t42 * t6) * mrSges(6,3) + (-Ifges(6,5) * t96 - Ifges(6,6) * t95) * t247 + (-Ifges(6,4) * t96 - Ifges(6,2) * t95) * t252 + (-Ifges(6,1) * t96 - Ifges(6,4) * t95) * t253 + t47 * (mrSges(6,1) * t95 - mrSges(6,2) * t96) + (-t107 * t139 - t108 * t137 + t112 * t130 - t113 * t129 + t152 * t74 + t178 * t75) * mrSges(4,3) + (-t152 * t129 + t130 * t178 + t137 * t242 + t139 * t244) * Ifges(4,4) - (t128 / 0.2e1 - t24 * mrSges(5,2) + t23 * mrSges(5,1) + t182 * t130 + (Ifges(4,2) + Ifges(6,3) / 0.2e1 + Ifges(5,3)) * t129 + t262) * t178 + (-t211 / 0.2e1 + (mrSges(3,1) * t259 - 0.3e1 / 0.2e1 * t228 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t176) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t178 + mrSges(4,2) * t152) + m(4) * (t158 + t204) + t104) * pkin(2) - t263) * t201 + m(6) * (t1 * t13 + t12 * t2 + t3 * t6 + t4 * t5 + t47 * t89 + t58 * t59) + (Ifges(6,1) * t41 + Ifges(6,4) * t42) * t248 + (Ifges(6,4) * t41 + Ifges(6,2) * t42) * t250 + t41 * t254 + t42 * t255 - t96 * t256 - t95 * t257 + (Ifges(6,5) * t41 + Ifges(6,6) * t42) * t245 + (t48 / 0.2e1 + t237 / 0.2e1 + t238 / 0.2e1 + t221 / 0.2e1 + t17 / 0.2e1 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t136 - t274) * t137 + (t229 / 0.2e1 + t98 / 0.2e1 + t182 * t243 + (-t39 * t170 - t38 * t172) * mrSges(5,3) + t260) * t139 + t167 * t189 + t3 * t45 + t4 * t46 + t58 * (-mrSges(6,1) * t42 + mrSges(6,2) * t41) + t59 * t26 + t35 * t80 + t34 * t81 + t52 * t87 + t51 * t88 + t89 * t11 + t112 * t78 + t93 * t123 + t214 * t92 + m(5) * (t23 * t51 + t24 * t52 + t34 * t38 + t35 * t39 + t92 * t94 + t220) + m(4) * (-t107 * t92 + t108 * t93 + t113 * t75 + t220) + (t55 * t241 + t56 * t240 + t74 * t185 + t182 * t247 - t181 * mrSges(5,3) + (Ifges(4,1) + Ifges(5,1) * t172 ^ 2 / 0.2e1 + (-t225 + t223 / 0.2e1) * t170) * t130) * t152 + (Ifges(4,5) * t139 / 0.2e1 - Ifges(4,6) * t137 / 0.2e1 + (t146 / 0.2e1 - t235 + t212 / 0.2e1 + (mrSges(3,2) * t259 + 0.3e1 / 0.2e1 * Ifges(3,4) * t176) * qJD(1)) * t176) * qJD(2) + t12 * t21 + t13 * t22; (Ifges(6,5) * t249 - Ifges(4,2) * t243 + Ifges(6,6) * t251 + Ifges(5,3) * t244 + Ifges(6,3) * t246 + t274 + t276) * t138 + (-Ifges(6,1) * t140 - Ifges(6,4) * t141) * t248 + (-Ifges(6,4) * t140 - Ifges(6,2) * t141) * t250 + (-Ifges(6,5) * t140 - Ifges(6,6) * t141) * t245 + t273 * t46 + (t1 * t103 + t102 * t2 + t157 * t47 + t272 * t6 + t273 * t5 - t58 * t73) * m(6) - (-Ifges(3,2) * t203 + t146 + t168) * t202 / 0.2e1 + (pkin(1) * (mrSges(3,1) * t174 + mrSges(3,2) * t176) - t174 * (Ifges(3,1) * t176 - t228) / 0.2e1) * qJD(1) ^ 2 + Ifges(3,5) * t194 + t263 * t203 + (Ifges(6,4) * t153 - Ifges(6,2) * t179) * t252 + (Ifges(6,1) * t153 - Ifges(6,4) * t179) * t253 + (Ifges(5,5) * t170 + Ifges(6,5) * t153 + Ifges(5,6) * t172 - Ifges(6,6) * t179) * t247 + t47 * (mrSges(6,1) * t179 + mrSges(6,2) * t153) + (-t1 * t179 - t153 * t2 - t270 * t6 + t271 * t5) * mrSges(6,3) - t179 * t257 + (Ifges(6,1) * t86 + Ifges(6,4) * t85) * t249 + ((t171 * t75 - t213 * t74) * pkin(2) + t107 * t109 - t108 * t110 - t158 * t199) * m(4) + (Ifges(6,5) * t86 + Ifges(6,6) * t85) * t246 + (-t132 + t98) * t243 + (-t129 * t236 - t130 * t196) * mrSges(4,3) - t214 * t109 + (mrSges(6,1) * t270 - mrSges(6,2) * t271) * t58 - t140 * t254 - t141 * t255 + t153 * t256 + t202 * t235 + t55 * t240 + t108 * t230 + t272 * t45 + (Ifges(6,4) * t86 + Ifges(6,2) * t85) * t251 + (-t182 * t244 + t269 / 0.2e1 + t260) * t136 + (-t207 * t38 - t208 * t39 + t264) * mrSges(5,3) + (-t180 * qJD(4) - t109 * t94 + t164 * t264 + t166 * t74 - t38 * t43 - t39 * t44) * m(5) + t265 * qJD(4) + (-t170 * t88 + t172 * t87) * t164 + (-mrSges(5,1) * t172 + mrSges(5,2) * t170 - mrSges(4,1)) * t74 + (-Ifges(4,1) * t136 + t17 - t227 + t48) * t242 - Ifges(3,6) * t195 - t73 * t26 - t75 * mrSges(4,2) - t44 * t80 - t43 * t81 - t85 * t18 / 0.2e1 - t86 * t19 / 0.2e1 - t104 * t199 + t102 * t21 + t103 * t22 - (Ifges(3,5) * t176 - Ifges(3,6) * t174) * t200 / 0.2e1 - t110 * t123 - Ifges(4,6) * t129 + Ifges(4,5) * t130 - t107 * t219 + (Ifges(5,1) * t170 + t225) * t209 / 0.2e1 - (Ifges(5,2) * t172 + t226) * t210 / 0.2e1 + t157 * t11 + t166 * t78 + (-mrSges(3,1) * t194 + mrSges(3,2) * t195) * pkin(6) + t170 * t56 / 0.2e1; -t179 * t21 + t153 * t22 + t170 * t87 + t172 * t88 - t270 * t46 - t271 * t45 + (-t26 - t214) * t138 + (t123 + t265) * t136 + t189 + (t1 * t153 - t58 * t138 - t179 * t2 - t270 * t5 - t271 * t6) * m(6) + (-t136 * t180 - t138 * t94 + t181) * m(5) + (t107 * t138 + t108 * t136 + t187) * m(4); t116 * t81 - t188 * t80 - t275 * t45 + t65 * t46 + t11 + t78 + (-t275 * t6 + t5 * t65 + t47) * m(6) + (t116 * t38 - t188 * t39 + t74) * m(5); t128 - t58 * (mrSges(6,1) * t65 + mrSges(6,2) * t275) + (Ifges(6,1) * t275 - t239) * t249 + t18 * t248 + (Ifges(6,5) * t275 - Ifges(6,6) * t65) * t246 - t5 * t45 + t6 * t46 + (t275 * t5 + t6 * t65) * mrSges(6,3) + (-Ifges(6,2) * t65 + t19 + t62) * t251 + t262;];
tauc = t7(:);
