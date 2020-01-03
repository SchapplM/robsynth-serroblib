% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:24
% EndTime: 2019-12-31 19:23:33
% DurationCPUTime: 4.26s
% Computational Cost: add. (4300->386), mult. (11240->544), div. (0->0), fcn. (10746->6), ass. (0->183)
t171 = sin(pkin(8));
t174 = cos(pkin(5));
t176 = cos(qJ(2));
t173 = cos(pkin(8));
t175 = sin(qJ(2));
t210 = t175 * t173;
t136 = t171 * t176 + t174 * t210;
t246 = -t136 / 0.2e1;
t212 = t174 * t175;
t137 = -t171 * t212 + t173 * t176;
t243 = t137 / 0.2e1;
t274 = Ifges(4,1) + Ifges(6,3);
t273 = Ifges(4,4) - Ifges(6,6);
t277 = Ifges(6,4) + Ifges(5,5);
t276 = Ifges(4,5) + Ifges(6,5);
t270 = Ifges(6,2) + Ifges(5,3);
t269 = -Ifges(5,6) + Ifges(6,6);
t172 = sin(pkin(5));
t215 = t172 * t175;
t275 = t276 * t215 / 0.2e1 + t274 * t243 + t273 * t246;
t266 = m(5) + m(6);
t267 = t270 * t136 + t269 * t137 + t277 * t215;
t211 = t174 * t176;
t204 = t171 * t211;
t135 = t204 + t210;
t232 = mrSges(5,3) * t135;
t161 = t173 * t211;
t217 = t171 * t175;
t134 = -t161 + t217;
t233 = mrSges(5,2) * t134;
t191 = -t232 - t233;
t71 = -mrSges(6,2) * t135 + mrSges(6,3) * t134;
t265 = t71 + t191;
t149 = -pkin(2) * t176 - qJ(3) * t215 - pkin(1);
t194 = qJ(3) * t174 + pkin(7);
t157 = t194 * t175;
t264 = t173 * (t149 * t172 - t157 * t174);
t214 = t172 * t176;
t103 = t134 * mrSges(5,1) + mrSges(5,3) * t214;
t104 = -t134 * mrSges(6,1) - mrSges(6,2) * t214;
t263 = t104 - t103;
t216 = t172 * t173;
t153 = -mrSges(5,1) * t216 - mrSges(5,3) * t174;
t154 = mrSges(6,1) * t216 + mrSges(6,2) * t174;
t262 = t154 - t153;
t261 = t172 ^ 2 + t174 ^ 2;
t259 = -Ifges(5,4) + t276;
t258 = -Ifges(4,6) + t277;
t255 = m(4) / 0.2e1;
t254 = -m(5) / 0.2e1;
t253 = m(5) / 0.2e1;
t252 = -m(6) / 0.2e1;
t251 = m(6) / 0.2e1;
t156 = t175 * pkin(2) - qJ(3) * t214;
t158 = t194 * t176;
t89 = t174 * t156 + t158 * t172;
t250 = m(4) * t89;
t249 = t134 / 0.2e1;
t248 = -t135 / 0.2e1;
t245 = t136 / 0.2e1;
t244 = -t137 / 0.2e1;
t140 = (-mrSges(6,2) * t171 - mrSges(6,3) * t173) * t172;
t242 = -t140 / 0.2e1;
t238 = t173 / 0.2e1;
t237 = t174 / 0.2e1;
t236 = t175 / 0.2e1;
t235 = pkin(3) + qJ(5);
t234 = m(6) * qJD(5);
t102 = t135 * mrSges(6,1) + mrSges(6,3) * t214;
t105 = t135 * mrSges(5,1) - mrSges(5,2) * t214;
t106 = mrSges(4,2) * t214 - t134 * mrSges(4,3);
t107 = -mrSges(4,1) * t214 - t135 * mrSges(4,3);
t143 = t171 * t158;
t209 = pkin(3) * t214 + t143;
t213 = t173 * t174;
t26 = t157 * t213 + t135 * pkin(4) + (qJ(5) * t176 - t149 * t173) * t172 + t209;
t218 = t171 * t174;
t219 = t171 * t172;
t45 = t149 * t219 - t157 * t218 + t173 * t158;
t39 = qJ(4) * t214 - t45;
t31 = -t134 * pkin(4) - t39;
t40 = t209 - t264;
t44 = -t143 + t264;
t96 = t210 * t261 + t204;
t97 = -t217 * t261 + t161;
t6 = (t106 + t263) * t97 + (t102 + t105 - t107) * t96 + m(6) * (t26 * t96 + t31 * t97) + m(4) * (-t44 * t96 + t45 * t97) + m(5) * (-t39 * t97 + t40 * t96);
t231 = qJD(1) * t6;
t228 = t136 * mrSges(6,1);
t100 = mrSges(6,2) * t215 - t228;
t101 = t137 * mrSges(5,1) + mrSges(5,2) * t215;
t108 = -mrSges(4,2) * t215 - mrSges(4,3) * t136;
t109 = mrSges(4,1) * t215 - mrSges(4,3) * t137;
t199 = -t214 / 0.2e1;
t84 = t174 * t149 + t157 * t172;
t184 = -qJ(4) * t135 + t84;
t27 = t134 * t235 + t184;
t142 = t171 * t157;
t192 = t158 * t213 - t142;
t220 = t156 * t173;
t28 = pkin(4) * t137 + (-t175 * t235 - t220) * t172 + t192;
t183 = -qJ(4) * t137 + t89;
t29 = t136 * t235 + t183;
t47 = t156 * t219 - t173 * t157 - t158 * t218;
t41 = -qJ(4) * t215 - t47;
t32 = -pkin(4) * t136 - t41;
t36 = pkin(3) * t134 + t184;
t37 = pkin(3) * t136 + t183;
t42 = (-pkin(3) * t175 - t220) * t172 + t192;
t46 = t142 + (t156 * t172 - t158 * t174) * t173;
t58 = Ifges(5,4) * t215 - t137 * Ifges(5,2) + t136 * Ifges(5,6);
t59 = Ifges(6,1) * t215 + Ifges(6,4) * t136 + Ifges(6,5) * t137;
t60 = Ifges(5,1) * t215 - Ifges(5,4) * t137 + Ifges(5,5) * t136;
t61 = Ifges(4,5) * t137 - Ifges(4,6) * t136 + Ifges(4,3) * t215;
t62 = t137 * Ifges(4,4) - t136 * Ifges(4,2) + Ifges(4,6) * t215;
t224 = t137 * mrSges(6,2);
t226 = t136 * mrSges(6,3);
t72 = -t224 + t226;
t225 = t137 * mrSges(4,2);
t229 = t136 * mrSges(4,1);
t73 = t225 + t229;
t223 = t137 * mrSges(5,3);
t227 = t136 * mrSges(5,2);
t74 = -t223 - t227;
t98 = t137 * mrSges(6,1) - mrSges(6,3) * t215;
t99 = mrSges(5,1) * t136 - mrSges(5,3) * t215;
t1 = (-t134 * t273 + t135 * t274 - t214 * t276) * t243 + (t134 * t270 + t135 * t269 - t214 * t277) * t245 + t37 * t191 + t135 * t275 + t267 * t249 + (-Ifges(5,4) * t214 - Ifges(5,2) * t135 + Ifges(5,6) * t134) * t244 + (Ifges(4,4) * t135 - Ifges(4,2) * t134 - Ifges(4,6) * t214) * t246 + t58 * t248 + (t61 + t60 + t59) * t199 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t175 + ((-Ifges(5,1) - Ifges(6,1) - Ifges(4,3)) * t214 + t259 * t135 + t258 * t134) * t172 / 0.2e1) * t175 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) - Ifges(3,2)) * t175 + Ifges(3,4) * t176) * t176 - t134 * t62 / 0.2e1 + t89 * (mrSges(4,1) * t134 + mrSges(4,2) * t135) + t32 * t104 + t42 * t105 + t47 * t106 + t46 * t107 + t45 * t108 + t44 * t109 + t26 * t98 + t39 * t99 + t31 * t100 + t40 * t101 + t28 * t102 + t41 * t103 + t84 * t73 + t29 * t71 + t27 * t72 + t36 * t74 + m(6) * (t26 * t28 + t27 * t29 + t31 * t32) + m(5) * (t36 * t37 + t39 * t41 + t40 * t42) + m(4) * (t44 * t46 + t45 * t47 + t84 * t89);
t230 = t1 * qJD(1);
t7 = (-m(5) * t39 + m(6) * t31 + t263) * t214 + (m(5) * t36 + m(6) * t27 + t265) * t135;
t222 = t7 * qJD(1);
t12 = t134 * t71 + m(6) * (t27 * t134 + t214 * t26) + t102 * t214;
t221 = t12 * qJD(1);
t141 = pkin(2) * t218 + qJ(3) * t216;
t208 = qJD(2) * t172;
t207 = qJD(2) * t174;
t206 = t172 * t234;
t205 = -mrSges(6,3) / 0.2e1 + mrSges(5,2) / 0.2e1;
t203 = mrSges(6,2) * t236;
t201 = -pkin(2) * t173 - pkin(3);
t198 = t214 / 0.2e1;
t197 = -t104 / 0.2e1 + t103 / 0.2e1;
t196 = -t153 / 0.2e1 + t154 / 0.2e1;
t195 = -qJ(4) * t171 - pkin(2);
t193 = t205 * t175;
t113 = -qJ(4) * t174 - t141;
t166 = qJ(3) * t219;
t117 = t174 * t201 + t166;
t138 = pkin(2) * t213 - t166;
t150 = mrSges(4,1) * t174 - mrSges(4,3) * t219;
t151 = -mrSges(4,2) * t174 + mrSges(4,3) * t216;
t152 = mrSges(6,1) * t219 - mrSges(6,3) * t174;
t155 = mrSges(5,1) * t219 + mrSges(5,2) * t174;
t88 = pkin(4) * t219 + t166 + (-qJ(5) + t201) * t174;
t94 = pkin(4) * t216 - t113;
t10 = ((t151 + t262) * t173 + (-t150 + t152 + t155) * t171 + m(6) * (t171 * t88 + t173 * t94) + m(4) * (-t138 * t171 + t141 * t173) + m(5) * (-t113 * t173 + t117 * t171)) * t172;
t177 = ((t106 / 0.2e1 - t197) * t173 + (-t107 / 0.2e1 + t105 / 0.2e1 + t102 / 0.2e1) * t171) * t172 + (t151 / 0.2e1 + t196) * t97 + (t152 / 0.2e1 - t150 / 0.2e1 + t155 / 0.2e1) * t96 + (-t138 * t96 + t141 * t97 + (-t171 * t44 + t173 * t45) * t172) * t255 + (-t113 * t97 + t117 * t96 + (t171 * t40 - t173 * t39) * t172) * t253 + (t88 * t96 + t94 * t97 + (t171 * t26 + t173 * t31) * t172) * t251;
t181 = -t250 / 0.2e1 + t37 * t254 + t29 * t252;
t3 = (mrSges(6,2) / 0.2e1 - mrSges(4,2) / 0.2e1 + mrSges(5,3) / 0.2e1) * t137 + (-mrSges(4,1) / 0.2e1 + t205) * t136 + t177 + t181;
t190 = qJD(1) * t3 + qJD(2) * t10;
t118 = (-pkin(3) * t173 + t195) * t172;
t139 = (mrSges(5,2) * t173 - mrSges(5,3) * t171) * t172;
t95 = (-t173 * t235 + t195) * t172;
t18 = (m(5) * t118 + m(6) * t95 + t139 + t140) * t219 + (m(5) * t113 - m(6) * t94 - t262) * t174;
t178 = (-t118 * t135 - t174 * t39 + (t113 * t176 - t171 * t36) * t172) * t254 + (-t95 * t135 + t174 * t31 + (-t171 * t27 - t176 * t94) * t172) * t252;
t180 = (mrSges(6,1) / 0.2e1 + mrSges(5,1) / 0.2e1) * t137 + t42 * t253 + t28 * t251;
t4 = t197 * t174 + (t140 / 0.2e1 + t139 / 0.2e1) * t135 + ((t71 / 0.2e1 - t233 / 0.2e1 - t232 / 0.2e1) * t171 + t196 * t176 + t193) * t172 + t178 + t180;
t189 = -qJD(1) * t4 - qJD(2) * t18;
t21 = t140 * t216 - m(6) * (-t174 * t88 - t216 * t95) + t174 * t152;
t179 = (t95 * t134 - t174 * t26 + (-t173 * t27 + t176 * t88) * t172) * t252 + t134 * t242 + t102 * t237;
t185 = t32 * t251 - t228 / 0.2e1;
t8 = (t203 + t71 * t238 - t176 * t152 / 0.2e1) * t172 + t179 + t185;
t188 = -qJD(1) * t8 - qJD(2) * t21;
t147 = t266 * t219;
t22 = 0.2e1 * t266 * (t96 / 0.4e1 + t135 / 0.4e1);
t186 = qJD(1) * t22 + qJD(2) * t147;
t48 = 0.2e1 * (t97 / 0.4e1 - t134 / 0.4e1) * m(6);
t182 = m(6) * t173 * t208 + qJD(1) * t48;
t148 = (qJD(1) * t214 - t207) * m(6);
t49 = m(6) * t249 + t251 * t97;
t23 = t266 * (t248 + t96 / 0.2e1);
t9 = -t71 * t216 / 0.2e1 + t152 * t198 + t172 * t203 - t179 + t185;
t5 = t104 * t237 + t135 * t242 + t153 * t198 + t139 * t248 + t154 * t199 - t174 * t103 / 0.2e1 + t172 * t193 - t178 + t180 - t265 * t219 / 0.2e1;
t2 = t226 / 0.2e1 - t224 / 0.2e1 + t225 / 0.2e1 + t229 / 0.2e1 - t223 / 0.2e1 - t227 / 0.2e1 + t177 - t181;
t11 = [qJD(2) * t1 + qJD(3) * t6 - qJD(4) * t7 + qJD(5) * t12, t230 + (-Ifges(3,6) * t175 + Ifges(3,5) * t176 + t28 * t152 + t41 * t153 + t32 * t154 + t42 * t155 + t37 * t139 + t29 * t140 + t141 * t108 + t46 * t150 + t47 * t151 + t138 * t109 + t118 * t74 + t113 * t99 + t117 * t101 + t95 * t72 + t88 * t98 + t94 * t100 + (-mrSges(3,1) * t176 + mrSges(3,2) * t175) * pkin(7)) * qJD(2) + t2 * qJD(3) + t5 * qJD(4) + t9 * qJD(5) + 0.2e1 * ((t138 * t46 + t141 * t47) * t255 + (t28 * t88 + t29 * t95 + t32 * t94) * t251 + (t113 * t41 + t117 * t42 + t118 * t37) * t253) * qJD(2) + (t59 / 0.2e1 + t60 / 0.2e1 + t61 / 0.2e1 + (Ifges(6,5) / 0.2e1 - Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t137 + (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1 - Ifges(4,6) / 0.2e1) * t136) * t207 + (t89 * (-mrSges(4,1) * t173 + mrSges(4,2) * t171) + (-Ifges(5,2) * t171 - Ifges(5,6) * t173) * t244 + (Ifges(4,4) * t171 + Ifges(4,2) * t173) * t246 - t171 * t58 / 0.2e1 + t62 * t238 + (t171 * t259 - t173 * t258) * t236 * t172 + (-t73 - t250) * pkin(2) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t212 + (t171 * t269 - t173 * t270) * t245 + (t171 * t274 + t173 * t273) * t243 + t171 * t275 - t267 * t173 / 0.2e1) * t208, qJD(2) * t2 + qJD(4) * t23 + qJD(5) * t49 + t231, qJD(2) * t5 + qJD(3) * t23 - t222, qJD(2) * t9 + qJD(3) * t49 + t221; qJD(3) * t3 - qJD(4) * t4 - qJD(5) * t8 - t230, qJD(3) * t10 - qJD(4) * t18 - qJD(5) * t21, t190, t189, t188; -qJD(2) * t3 - qJD(4) * t22 - qJD(5) * t48 - t231, -qJD(4) * t147 - t173 * t206 - t190, 0, -t186, -t182; t4 * qJD(2) + t22 * qJD(3) + t176 * t206 + t222, qJD(3) * t147 - t174 * t234 - t189, t186, 0, t148; -m(6) * qJD(4) * t214 + t8 * qJD(2) + t48 * qJD(3) - t221, (qJD(3) * t216 + qJD(4) * t174) * m(6) - t188, t182, -t148, 0;];
Cq = t11;
