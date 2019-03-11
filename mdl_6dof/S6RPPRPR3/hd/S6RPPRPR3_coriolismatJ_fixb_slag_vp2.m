% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:11
% EndTime: 2019-03-09 01:44:14
% DurationCPUTime: 2.55s
% Computational Cost: add. (6555->274), mult. (11748->385), div. (0->0), fcn. (12080->8), ass. (0->171)
t156 = sin(pkin(10));
t158 = sin(qJ(4));
t237 = cos(pkin(10));
t267 = cos(qJ(4));
t134 = t156 * t158 - t237 * t267;
t135 = t156 * t267 + t237 * t158;
t157 = sin(qJ(6));
t154 = t157 ^ 2;
t159 = cos(qJ(6));
t155 = t159 ^ 2;
t220 = t154 + t155;
t291 = 0.1e1 - t220;
t205 = m(7) * t291;
t295 = t205 * t134 * t135;
t296 = t295 * qJD(4);
t279 = m(7) / 0.2e1;
t294 = t135 * mrSges(6,1) - t134 * mrSges(6,2);
t143 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t129 = (-qJ(5) + t143) * t158;
t213 = t267 * t143;
t130 = -t267 * qJ(5) + t213;
t87 = t129 * t156 - t237 * t130;
t219 = t267 * pkin(4);
t96 = -t134 * pkin(5) + t135 * pkin(8) + t219;
t44 = t157 * t87 + t159 * t96;
t45 = t157 * t96 - t159 * t87;
t192 = -t44 * t157 + t45 * t159;
t293 = (-t192 - t87) * t279;
t144 = sin(pkin(9)) * pkin(1) + qJ(3);
t138 = t158 * pkin(4) + t144;
t292 = m(6) * t138;
t245 = t159 * mrSges(7,1);
t250 = t157 * mrSges(7,2);
t200 = t245 - t250;
t232 = t134 * t200;
t244 = t159 * mrSges(7,2);
t251 = t157 * mrSges(7,1);
t179 = t244 / 0.2e1 + t251 / 0.2e1;
t170 = t179 * t135;
t152 = Ifges(7,4) * t159;
t197 = Ifges(7,2) * t157 - t152;
t289 = Ifges(7,1) * t157 + t152;
t290 = t289 - t197;
t288 = t134 * t156 + t237 * t135;
t287 = t220 * mrSges(7,3);
t229 = t134 * t157;
t97 = -t135 * mrSges(7,2) + mrSges(7,3) * t229;
t248 = t157 * t97;
t255 = t134 * mrSges(7,3);
t286 = -t248 / 0.2e1 + t220 * t255 / 0.2e1;
t148 = -t237 * pkin(4) - pkin(5);
t285 = m(7) * t148 - t200;
t108 = t134 ^ 2;
t283 = t135 ^ 2;
t281 = m(6) / 0.2e1;
t280 = -m(7) / 0.2e1;
t278 = m(6) * pkin(4);
t277 = mrSges(7,1) / 0.2e1;
t276 = -mrSges(7,2) / 0.2e1;
t275 = -Ifges(7,3) / 0.2e1;
t273 = t135 / 0.2e1;
t140 = t244 + t251;
t272 = t140 / 0.2e1;
t147 = pkin(4) * t156 + pkin(8);
t271 = t147 / 0.2e1;
t270 = t157 / 0.2e1;
t269 = -t159 / 0.2e1;
t268 = t159 / 0.2e1;
t263 = mrSges(7,3) * t135;
t261 = Ifges(7,4) * t157;
t260 = Ifges(7,5) * t135;
t259 = Ifges(7,5) * t159;
t257 = Ifges(7,6) * t135;
t256 = Ifges(7,6) * t157;
t254 = t134 * t87;
t76 = t197 * t134 + t257;
t249 = t157 * t76;
t228 = t134 * t159;
t98 = t135 * mrSges(7,1) + mrSges(7,3) * t228;
t247 = t157 * t98;
t246 = t158 * mrSges(5,1);
t199 = Ifges(7,1) * t159 - t261;
t78 = -t199 * t134 + t260;
t243 = t159 * t78;
t242 = t159 * t97;
t241 = t159 * t98;
t238 = t87 * t140;
t165 = -t267 * mrSges(5,2) - t246 - t294;
t174 = t237 * t129 + t156 * t130;
t86 = pkin(5) * t135 + pkin(8) * t134 + t138;
t40 = -t157 * t174 + t159 * t86;
t41 = t157 * t86 + t159 * t174;
t194 = t41 * t157 + t40 * t159;
t17 = t248 + t241 + mrSges(4,3) + m(7) * t194 + t292 + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t144 - t165;
t236 = qJD(1) * t17;
t184 = t134 * mrSges(7,2) + t157 * t263;
t227 = t135 * t159;
t185 = -t134 * mrSges(7,1) + mrSges(7,3) * t227;
t161 = (t157 * t45 + t159 * t44) * t279 + t184 * t270 + t185 * t268 + t219 * t281;
t206 = t220 * t147;
t231 = t134 * t148;
t189 = (-t135 * t206 - t231) * t279 + (t237 * t134 - t135 * t156) * pkin(4) * t281;
t166 = t232 / 0.2e1 + t189;
t131 = t135 * mrSges(6,2);
t208 = t134 * mrSges(6,1) + t131;
t12 = -t161 + t166 + t208 - t220 * t263 / 0.2e1;
t235 = t12 * qJD(1);
t177 = t108 * t200;
t180 = t245 / 0.2e1 - t250 / 0.2e1;
t182 = t248 / 0.2e1 + t241 / 0.2e1;
t210 = t155 / 0.2e1 + t154 / 0.2e1;
t202 = mrSges(7,3) * t210;
t13 = t177 / 0.2e1 + (-t134 * t202 + t182) * t135 + t180;
t234 = t13 * qJD(1);
t226 = t147 * t159;
t141 = Ifges(7,2) * t159 + t261;
t225 = t157 * t141;
t176 = t135 * t200;
t16 = (t176 / 0.2e1 + t210 * t255 - t182) * t134;
t224 = t16 * qJD(1);
t181 = t247 / 0.2e1 - t242 / 0.2e1;
t18 = t170 + t181;
t223 = t18 * qJD(1);
t163 = (-t220 * t283 - t108) * t279 + (-t108 - t283) * t281;
t188 = t220 * t280 - m(6) / 0.2e1;
t30 = t163 + t188;
t222 = t30 * qJD(1);
t33 = -t295 / 0.2e1;
t221 = t33 * qJD(1);
t218 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t217 = t267 * mrSges(5,1);
t211 = -t227 / 0.2e1;
t196 = -t256 + t259;
t195 = Ifges(7,5) * t157 + Ifges(7,6) * t159;
t193 = t157 * t40 - t159 * t41;
t3 = t40 * t97 - t41 * t98 + (t195 * t273 - t87 * t200 + t78 * t270 + t76 * t268 + (t289 * t269 + t225 / 0.2e1) * t134 - t193 * mrSges(7,3)) * t134;
t191 = t3 * qJD(1) - t16 * qJD(2);
t9 = (mrSges(6,3) * t135 - t242 + t247) * t135 + (mrSges(6,3) + t140) * t108 + m(7) * (t193 * t135 - t254) + m(6) * (-t135 * t174 - t254);
t190 = -qJD(1) * t9 - qJD(2) * t33;
t187 = t45 * t276 + t44 * t277;
t183 = -t259 / 0.2e1 + t256 / 0.2e1;
t178 = t134 * t195;
t175 = t135 * t140;
t171 = t183 * t135;
t75 = -Ifges(7,6) * t134 + t197 * t135;
t77 = -Ifges(7,5) * t134 - t199 * t135;
t1 = t144 * t217 - t138 * t131 + t44 * t98 + t45 * t97 + m(7) * (t87 * t174 + t40 * t44 + t41 * t45) + (t77 * t269 + t75 * t270 - t138 * mrSges(6,1) + t41 * mrSges(7,2) - t40 * mrSges(7,1) - t174 * t140 + (-Ifges(6,4) - t183) * t134) * t134 + (-t243 / 0.2e1 + t249 / 0.2e1 + Ifges(6,4) * t135 - t238 + t171 + t194 * mrSges(7,3) + (-Ifges(6,2) + Ifges(6,1) - Ifges(7,3)) * t134) * t135 + (-t144 * mrSges(5,2) + Ifges(5,4) * t158) * t158 + ((-Ifges(5,1) + Ifges(5,2)) * t158 - Ifges(5,4) * t267) * t267 + (t292 + t294) * t219;
t162 = (t193 + t174) * t280 - t181;
t7 = t135 * t293 + (t175 / 0.2e1 + t162) * t134;
t8 = (-t162 - t170) * t135 + t134 * t293;
t173 = t1 * qJD(1) + t8 * qJD(2) - t7 * qJD(3);
t169 = -t140 / 0.2e1 + t179;
t21 = t283 * t205 / 0.2e1 - t291 * t279 * t108;
t168 = -t7 * qJD(1) + t21 * qJD(2) + qJD(3) * t295;
t167 = t8 * qJD(1) - qJD(2) * t295 + t21 * qJD(3);
t4 = -t238 / 0.2e1 + (0.3e1 / 0.4e1 * t257 + t76 / 0.4e1 + t97 * t271) * t157 + (t275 - t147 * t202 + (t148 * t276 - t289 / 0.4e1 - t152 / 0.4e1 - t218 * t157) * t157) * t134 + (-0.3e1 / 0.4e1 * t260 - t78 / 0.4e1 + t98 * t271 + (-0.3e1 / 0.4e1 * t261 + t148 * t277 - t141 / 0.4e1 + t218 * t159) * t134) * t159 + t187;
t47 = t199 * t270 - t225 / 0.2e1 + t148 * t140 + t290 * t268;
t49 = t169 * t135;
t50 = t169 * t134;
t164 = t4 * qJD(1) + t49 * qJD(2) + t50 * qJD(3) - t47 * qJD(4);
t52 = t135 * t272 + t170;
t51 = (t179 + t272) * t134;
t29 = t163 - t188;
t20 = t21 * qJD(4);
t19 = t170 - t181;
t15 = -t135 * t202 + t161 + t166;
t14 = -t177 / 0.2e1 + t98 * t211 + t180 + t286 * t135;
t6 = t7 * qJD(4);
t5 = t135 * t196 / 0.4e1 + t243 / 0.4e1 - t249 / 0.4e1 - t200 * t231 / 0.2e1 + t238 / 0.2e1 - t98 * t226 / 0.2e1 + t134 * t275 + t171 + t187 + (t141 / 0.2e1 - t199 / 0.4e1) * t228 + (t289 + t290) * t229 / 0.4e1 + t286 * t147;
t2 = qJD(4) * t8 + qJD(5) * t33 - qJD(6) * t16;
t10 = [qJD(3) * t17 + qJD(4) * t1 + qJD(5) * t9 + qJD(6) * t3, t2, qJD(5) * t29 + qJD(6) * t14 + t236 - t6 (-Ifges(6,5) * t135 + Ifges(6,6) * t134 - Ifges(5,5) * t158 - Ifges(5,6) * t267 - t178 / 0.2e1 + t77 * t270 + t75 * t268 - t148 * t175 + t289 * t211 + t225 * t273 - mrSges(5,2) * t213 - t143 * t246 + t184 * t226 - (t156 * t278 - mrSges(6,2)) * t87 + (m(7) * t192 - t157 * t185) * t147 + (-t237 * t278 - mrSges(6,1) + t285) * t174 + t288 * mrSges(6,3) * pkin(4) + t192 * mrSges(7,3)) * qJD(4) + t15 * qJD(5) + t5 * qJD(6) + t173, qJD(3) * t29 + qJD(4) * t15 + qJD(6) * t19 - t190, t14 * qJD(3) + t5 * qJD(4) + t19 * qJD(5) + (-mrSges(7,1) * t41 - mrSges(7,2) * t40 + t178) * qJD(6) + t191; t2, -t296, t20, t52 * qJD(6) + t167 + (t158 * mrSges(5,2) - t287 * t135 + 0.2e1 * t189 + t208 - t217 + t232) * qJD(4), t221, t52 * qJD(4) + qJD(6) * t232 - t224; qJD(5) * t30 - qJD(6) * t13 - t236 - t6, t20, t296 (-t288 * t278 + t165 + t285 * t135 + (-m(7) * t206 - t287) * t134) * qJD(4) + t51 * qJD(6) + t168, t222, t51 * qJD(4) - qJD(6) * t176 - t234; qJD(5) * t12 - qJD(6) * t4 - t173, -qJD(6) * t49 - t167, -qJD(6) * t50 - t168, t47 * qJD(6), t235 (-t147 * t200 + t196) * qJD(6) - t164; -qJD(3) * t30 - qJD(4) * t12 - qJD(6) * t18 + t190, -t221, -t222, -t235, 0, -qJD(6) * t140 - t223; qJD(3) * t13 + qJD(4) * t4 + qJD(5) * t18 - t191, qJD(4) * t49 + t224, qJD(4) * t50 + t234, t164, t223, 0;];
Cq  = t10;
