% Calculate time derivative of joint inertia matrix for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:08
% EndTime: 2019-03-09 11:43:15
% DurationCPUTime: 3.87s
% Computational Cost: add. (7045->343), mult. (15277->489), div. (0->0), fcn. (15355->8), ass. (0->166)
t280 = Ifges(6,1) + Ifges(7,1);
t162 = sin(pkin(10));
t163 = cos(pkin(10));
t166 = sin(qJ(2));
t169 = cos(qJ(2));
t128 = -t162 * t166 + t163 * t169;
t129 = t162 * t169 + t163 * t166;
t165 = sin(qJ(4));
t168 = cos(qJ(4));
t179 = t128 * t168 - t129 * t165;
t241 = Ifges(7,4) + Ifges(6,5);
t279 = t179 * t241;
t164 = sin(qJ(5));
t167 = cos(qJ(5));
t278 = t164 ^ 2 + t167 ^ 2;
t229 = Ifges(7,5) * t164;
t231 = Ifges(6,4) * t164;
t277 = t167 * t280 + t229 - t231;
t228 = Ifges(7,5) * t167;
t230 = Ifges(6,4) * t167;
t276 = t164 * t280 - t228 + t230;
t275 = Ifges(6,6) * t167 + t164 * t241;
t274 = Ifges(7,2) + Ifges(6,3);
t210 = qJD(5) * t164;
t120 = t129 * qJD(2);
t121 = t128 * qJD(2);
t70 = qJD(4) * t179 - t120 * t165 + t121 * t168;
t221 = t167 * t70;
t98 = t128 * t165 + t129 * t168;
t177 = t210 * t98 - t221;
t209 = qJD(5) * t167;
t202 = t98 * t209;
t224 = t164 * t70;
t178 = t202 + t224;
t71 = qJD(4) * t98 + t120 * t168 + t121 * t165;
t273 = t241 * t71 + (-Ifges(6,4) + Ifges(7,5)) * t178 - t280 * t177;
t236 = t277 * t98 - t279;
t141 = -t167 * mrSges(6,1) + t164 * mrSges(6,2);
t272 = -mrSges(5,1) + t141;
t271 = t277 * qJD(5);
t270 = Ifges(7,6) * t210 + t209 * t241;
t151 = pkin(2) * t163 + pkin(3);
t248 = pkin(2) * t162;
t114 = t151 * t168 - t165 * t248;
t108 = t114 * qJD(4);
t268 = t278 * t108;
t240 = -qJ(3) - pkin(7);
t196 = qJD(2) * t240;
t117 = qJD(3) * t169 + t166 * t196;
t118 = -qJD(3) * t166 + t169 * t196;
t89 = -t117 * t162 + t118 * t163;
t176 = -pkin(8) * t121 + t89;
t142 = t240 * t166;
t143 = t240 * t169;
t100 = t142 * t163 + t143 * t162;
t91 = -pkin(8) * t129 + t100;
t101 = t142 * t162 - t143 * t163;
t92 = pkin(8) * t128 + t101;
t265 = -t165 * t92 + t168 * t91;
t90 = t117 * t163 + t118 * t162;
t77 = -pkin(8) * t120 + t90;
t26 = qJD(4) * t265 + t165 * t176 + t168 * t77;
t159 = qJD(2) * t166 * pkin(2);
t104 = pkin(3) * t120 + t159;
t36 = pkin(4) * t71 - pkin(9) * t70 + t104;
t152 = -pkin(2) * t169 - pkin(1);
t105 = -pkin(3) * t128 + t152;
t50 = -pkin(4) * t179 - pkin(9) * t98 + t105;
t52 = t165 * t91 + t168 * t92;
t6 = t164 * t36 + t167 * t26 + t209 * t50 - t210 * t52;
t235 = t164 * t50 + t167 * t52;
t7 = -qJD(5) * t235 - t164 * t26 + t167 * t36;
t267 = -t164 * t7 + t167 * t6;
t2 = qJ(6) * t71 - qJD(6) * t179 + t6;
t22 = -qJ(6) * t179 + t235;
t32 = -t164 * t52 + t167 * t50;
t23 = pkin(5) * t179 - t32;
t4 = -pkin(5) * t71 - t7;
t266 = t164 * t4 + t167 * t2 + t209 * t23 - t210 * t22;
t115 = t151 * t165 + t168 * t248;
t264 = (-mrSges(5,2) + (mrSges(7,2) + mrSges(6,3)) * t278) * t108;
t263 = 2 * m(5);
t262 = 2 * m(6);
t261 = 2 * m(7);
t260 = -2 * mrSges(5,3);
t259 = -2 * Ifges(5,4);
t27 = qJD(4) * t52 + t165 * t77 - t168 * t176;
t258 = 0.2e1 * t27;
t257 = -0.2e1 * t265;
t256 = 0.2e1 * t104;
t189 = t164 * mrSges(7,1) - t167 * mrSges(7,3);
t133 = t189 * qJD(5);
t255 = 0.2e1 * t133;
t254 = 0.2e1 * t152;
t253 = m(4) * pkin(2);
t252 = t71 / 0.2e1;
t145 = Ifges(6,2) * t167 + t231;
t250 = -t145 / 0.2e1;
t244 = t27 * t265;
t242 = t179 * Ifges(6,6);
t28 = mrSges(6,1) * t71 + mrSges(6,3) * t177;
t29 = -t71 * mrSges(7,1) - mrSges(7,2) * t177;
t239 = t28 - t29;
t30 = -mrSges(6,2) * t71 - mrSges(6,3) * t178;
t31 = -mrSges(7,2) * t178 + mrSges(7,3) * t71;
t238 = t30 + t31;
t185 = Ifges(7,3) * t164 + t228;
t42 = -Ifges(7,6) * t179 + t185 * t98;
t186 = -Ifges(6,2) * t164 + t230;
t43 = t186 * t98 - t242;
t237 = t42 - t43;
t223 = t164 * t98;
t72 = mrSges(6,2) * t179 - mrSges(6,3) * t223;
t75 = -mrSges(7,2) * t223 - mrSges(7,3) * t179;
t234 = t72 + t75;
t220 = t167 * t98;
t73 = -mrSges(6,1) * t179 - mrSges(6,3) * t220;
t74 = mrSges(7,1) * t179 + mrSges(7,2) * t220;
t233 = -t73 + t74;
t112 = pkin(9) + t115;
t232 = t268 * t112;
t227 = Ifges(6,6) * t164;
t109 = t115 * qJD(4);
t225 = t109 * t265;
t219 = qJD(5) * t98;
t216 = t108 * t164;
t215 = t108 * t167;
t213 = t268 * pkin(9);
t208 = qJD(6) * t167;
t207 = 0.2e1 * t169;
t206 = m(7) * t208;
t199 = t71 * mrSges(5,1) + mrSges(5,2) * t70;
t198 = -t210 / 0.2e1;
t195 = t272 * t109;
t194 = t120 * mrSges(4,1) + mrSges(4,2) * t121;
t193 = mrSges(7,2) * t208 + t270;
t190 = t164 * mrSges(6,1) + t167 * mrSges(6,2);
t140 = -t167 * mrSges(7,1) - t164 * mrSges(7,3);
t184 = -pkin(5) * t167 - qJ(6) * t164;
t183 = pkin(5) * t164 - qJ(6) * t167;
t180 = Ifges(7,6) * t178 + t221 * t241 + t274 * t71;
t139 = -pkin(4) + t184;
t116 = -pkin(5) * t210 + qJ(6) * t209 + qJD(6) * t164;
t175 = t164 * t233 + t167 * t234;
t174 = mrSges(7,2) * t184 - t227;
t135 = t185 * qJD(5);
t136 = t186 * qJD(5);
t144 = -Ifges(7,3) * t167 + t229;
t173 = (t144 - t145) * t210 + t276 * t209 + (-t135 + t136) * t167 + t271 * t164;
t172 = m(7) * t184 + t140 + t141;
t171 = t238 * t167 - t239 * t164 + (-t164 * t234 + t167 * t233) * qJD(5) + m(7) * t266 + m(6) * (-t209 * t32 - t210 * t235 + t267);
t134 = t190 * qJD(5);
t14 = -Ifges(7,5) * t177 + t71 * Ifges(7,6) + Ifges(7,3) * t178;
t15 = -Ifges(6,4) * t177 - Ifges(6,2) * t178 + t71 * Ifges(6,6);
t37 = t183 * t98 - t265;
t9 = (pkin(5) * t70 + qJ(6) * t219) * t164 + (-qJ(6) * t70 + (pkin(5) * qJD(5) - qJD(6)) * t98) * t167 + t27;
t170 = t202 * t250 + t43 * t198 - t26 * mrSges(5,2) + Ifges(5,5) * t70 - Ifges(5,6) * t71 + t37 * t133 - t265 * t134 + t9 * t140 + t42 * t210 / 0.2e1 + t272 * t27 + t275 * t252 - (-Ifges(6,6) * t210 + t270) * t179 / 0.2e1 + t273 * t164 / 0.2e1 + (t250 + t144 / 0.2e1) * t224 + (-t136 / 0.2e1 + t135 / 0.2e1) * t223 + t271 * t220 / 0.2e1 + (-Ifges(7,6) * t252 + t15 / 0.2e1 - t14 / 0.2e1) * t167 + ((-t164 * t235 - t167 * t32) * qJD(5) + t267) * mrSges(6,3) + (t144 * t98 + t236) * t209 / 0.2e1 + t266 * mrSges(7,2) + t276 * (t98 * t198 + t221 / 0.2e1);
t154 = mrSges(7,2) * t209;
t111 = -pkin(4) - t114;
t99 = -t114 + t139;
t96 = t109 - t116;
t64 = t190 * t98;
t63 = t189 * t98;
t19 = mrSges(6,1) * t178 - mrSges(6,2) * t177;
t18 = mrSges(7,1) * t178 + mrSges(7,3) * t177;
t1 = [t194 * t254 + t19 * t257 + t64 * t258 + (t2 * t22 + t23 * t4 + t37 * t9) * t261 + 0.2e1 * m(4) * (t100 * t89 + t101 * t90) - (mrSges(5,1) * t256 + t26 * t260 + (t259 - t227) * t70 + ((2 * Ifges(5,2)) + t274) * t71 + t180) * t179 - 0.2e1 * t128 * Ifges(4,2) * t120 + 0.2e1 * t129 * t121 * Ifges(4,1) + 0.2e1 * t235 * t30 + (mrSges(5,2) * t256 + mrSges(5,3) * t258 + 0.2e1 * Ifges(5,1) * t70 + t273 * t167 + (t14 - t15) * t164 + (t259 + t241 * t167 + (-Ifges(6,6) + Ifges(7,6)) * t164) * t71 + ((t237 + t242) * t167 + (-t236 + t279) * t164) * qJD(5)) * t98 + 0.2e1 * (-t120 * t129 + t121 * t128) * Ifges(4,4) + 0.2e1 * (-t100 * t121 - t101 * t120 + t128 * t90 - t129 * t89) * mrSges(4,3) + t52 * t71 * t260 + (t104 * t105 + t26 * t52 - t244) * t263 + (t235 * t6 + t32 * t7 - t244) * t262 + 0.2e1 * t105 * t199 + 0.2e1 * t23 * t29 + 0.2e1 * t22 * t31 + 0.2e1 * t32 * t28 + 0.2e1 * t37 * t18 + 0.2e1 * t9 * t63 + 0.2e1 * t6 * t72 + 0.2e1 * t7 * t73 + 0.2e1 * t4 * t74 + 0.2e1 * t2 * t75 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t169) * t207 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t128 + mrSges(4,2) * t129) + t253 * t254 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t166 + (Ifges(3,1) - Ifges(3,2)) * t207) * t166) * qJD(2) + (mrSges(5,3) * t257 + t237 * t164 + t236 * t167) * t70; t170 + (mrSges(5,3) * t179 + t175) * t108 + (Ifges(3,5) * t169 - Ifges(3,6) * t166 + (-mrSges(3,1) * t169 + mrSges(3,2) * t166) * pkin(7)) * qJD(2) + (t109 * t98 - t114 * t70 - t115 * t71) * mrSges(5,3) + m(7) * (t215 * t22 + t216 * t23 + t37 * t96 + t9 * t99) + m(6) * (t111 * t27 + t215 * t235 - t216 * t32 - t225) + t171 * t112 + m(5) * (t108 * t52 - t114 * t27 + t115 * t26 - t225) + (-t120 * t162 - t121 * t163) * pkin(2) * mrSges(4,3) + (t162 * t90 + t163 * t89) * t253 + t89 * mrSges(4,1) - t90 * mrSges(4,2) + t96 * t63 + t99 * t18 + t109 * t64 + t111 * t19 - Ifges(4,6) * t120 + Ifges(4,5) * t121; t173 + 0.2e1 * t195 + 0.2e1 * t264 + (t96 * t99 + t232) * t261 + (t109 * t111 + t232) * t262 + (t108 * t115 - t109 * t114) * t263 + t99 * t255 + 0.2e1 * t111 * t134 + 0.2e1 * t96 * t140; m(4) * t159 + t239 * t167 + t238 * t164 + t175 * qJD(5) + m(7) * (t164 * t2 - t167 * t4 + (t164 * t23 + t167 * t22) * qJD(5)) + m(6) * (t164 * t6 + t167 * t7 + (-t164 * t32 + t167 * t235) * qJD(5)) + m(5) * t104 + t194 + t199; 0; 0; t170 + m(7) * (-t116 * t37 + t139 * t9) + t171 * pkin(9) - t116 * t63 + t139 * t18 + (-m(6) * t27 - t19) * pkin(4); (t96 - t116) * t140 + (-pkin(4) + t111) * t134 + (t99 + t139) * t133 + t195 + m(7) * (-t116 * t99 + t139 * t96 + t213) + m(6) * (-pkin(4) * t109 + t213) + t264 + t173; 0; -0.2e1 * pkin(4) * t134 + t139 * t255 + 0.2e1 * (-m(7) * t139 - t140) * t116 + t173; -Ifges(6,6) * t224 - pkin(5) * t29 + m(7) * (-pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t22) + qJD(6) * t75 + qJ(6) * t31 + t2 * mrSges(7,3) - t6 * mrSges(6,2) + t7 * mrSges(6,1) - t4 * mrSges(7,1) - t275 * t219 + t180; t112 * t206 + (-m(7) * t183 - t189 - t190) * t108 + (t112 * t172 + t174) * qJD(5) + t193; m(7) * t116 + ((-mrSges(6,2) + mrSges(7,3)) * t167 + (-mrSges(6,1) - mrSges(7,1)) * t164) * qJD(5); t174 * qJD(5) + (qJD(5) * t172 + t206) * pkin(9) + t193; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t4 + t29; t154 + (t112 * t209 + t216) * m(7); m(7) * t210; m(7) * pkin(9) * t209 + t154; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
