% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:52
% EndTime: 2019-03-09 01:29:55
% DurationCPUTime: 1.33s
% Computational Cost: add. (2534->249), mult. (5081->345), div. (0->0), fcn. (3713->6), ass. (0->161)
t152 = cos(qJ(6));
t229 = t152 * mrSges(7,2);
t150 = sin(qJ(6));
t233 = t150 * mrSges(7,1);
t118 = t229 + t233;
t151 = sin(qJ(5));
t211 = t151 * t118;
t263 = t211 / 0.2e1;
t153 = cos(qJ(5));
t230 = t152 * mrSges(7,1);
t232 = t150 * mrSges(7,2);
t181 = t230 - t232;
t262 = t181 * t153;
t146 = t150 ^ 2;
t148 = t152 ^ 2;
t183 = (-t148 / 0.2e1 - t146 / 0.2e1) * mrSges(7,3);
t200 = t146 + t148;
t261 = t200 * mrSges(7,3);
t145 = Ifges(7,4) * t152;
t120 = Ifges(7,1) * t150 + t145;
t247 = -t152 / 0.2e1;
t260 = t120 * t247;
t179 = -Ifges(7,2) * t150 + t145;
t241 = t151 * pkin(8);
t243 = pkin(5) * t153;
t123 = t241 + t243;
t182 = sin(pkin(9)) * pkin(1) + qJ(3);
t136 = -pkin(7) + t182;
t212 = t150 * t153;
t46 = t152 * t123 - t136 * t212;
t206 = t152 * t153;
t47 = t150 * t123 + t136 * t206;
t176 = -t46 * t150 + t47 * t152;
t259 = -mrSges(6,2) + t261;
t258 = t200 * t241 + t243;
t257 = m(7) / 0.2e1;
t256 = -pkin(8) / 0.2e1;
t255 = mrSges(7,2) / 0.2e1;
t147 = t151 ^ 2;
t149 = t153 ^ 2;
t38 = (-0.1e1 + t200) * (-t147 + t149);
t254 = -t38 / 0.2e1;
t253 = -t262 / 0.2e1;
t184 = t200 * t153;
t50 = (-t153 + t184) * t151;
t252 = m(7) * t50;
t251 = -t150 / 0.2e1;
t250 = t150 / 0.2e1;
t249 = t151 / 0.2e1;
t248 = t151 / 0.4e1;
t246 = t152 / 0.2e1;
t245 = -t153 / 0.2e1;
t244 = t153 / 0.2e1;
t242 = t151 * pkin(5);
t240 = m(7) * qJD(5);
t238 = Ifges(7,4) * t150;
t237 = Ifges(7,5) * t151;
t144 = Ifges(7,5) * t152;
t235 = Ifges(7,6) * t150;
t234 = Ifges(7,6) * t151;
t80 = t153 * t179 + t234;
t231 = t150 * t80;
t210 = t151 * t152;
t137 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t95 = -t153 * pkin(8) + t137 + t242;
t37 = t136 * t210 + t150 * t95;
t228 = t152 * t37;
t180 = Ifges(7,1) * t152 - t238;
t82 = t153 * t180 + t237;
t227 = t152 * t82;
t112 = t151 * mrSges(7,1) - mrSges(7,3) * t206;
t213 = t150 * t151;
t109 = -t153 * mrSges(7,2) + mrSges(7,3) * t213;
t209 = t152 * t109;
t111 = t153 * mrSges(7,1) + mrSges(7,3) * t210;
t214 = t150 * t111;
t98 = t153 * t118;
t161 = t98 / 0.2e1 - t214 / 0.2e1 + t209 / 0.2e1;
t168 = m(7) * t176;
t189 = t213 / 0.2e1;
t216 = t149 * t136;
t36 = -t136 * t213 + t152 * t95;
t193 = mrSges(7,3) * t212;
t110 = -t151 * mrSges(7,2) - t193;
t208 = t152 * t110;
t83 = t151 * t208;
t6 = t112 * t189 - t83 / 0.2e1 - t211 * t249 + (t147 * t136 - t37 * t210 + t36 * t213 - t216) * t257 + (t168 / 0.2e1 + t161) * t153;
t224 = t6 * qJD(1);
t167 = t112 * t251 + t208 / 0.2e1;
t177 = -t150 * t36 + t228;
t218 = t136 * t153;
t8 = (t177 * t257 + t167 + t263) * t153 + ((t176 - 0.2e1 * t218) * t257 + t161) * t151;
t223 = t8 * qJD(1);
t222 = -mrSges(6,1) - t181;
t207 = t152 * t112;
t215 = t150 * t110;
t165 = t207 / 0.2e1 + t215 / 0.2e1;
t15 = -t149 * t183 + t151 * t253 + t165 * t153;
t221 = qJD(1) * t15;
t17 = t151 * mrSges(6,1) + t153 * mrSges(6,2) + t215 + t207 + mrSges(5,3) + m(7) * (t150 * t37 + t152 * t36) + 0.4e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * t137;
t220 = qJD(1) * t17;
t156 = (t153 * t183 - t165) * t151 + t262 * t245;
t169 = -t232 / 0.2e1 + t230 / 0.2e1;
t10 = t156 - t169;
t219 = t10 * qJD(1);
t143 = t153 * mrSges(6,1);
t157 = (-t150 * t47 - t152 * t46) * t257 + t109 * t251 + t111 * t247;
t158 = t258 * t257 + t262 / 0.2e1;
t14 = -t143 + (mrSges(6,2) + t183) * t151 + t157 - t158;
t217 = t14 * qJD(1);
t188 = -t210 / 0.2e1;
t191 = -t233 / 0.2e1;
t201 = mrSges(7,2) * t188 + t151 * t191;
t19 = t167 + t201;
t204 = t19 * qJD(1);
t198 = t149 + t147;
t186 = m(6) * t198;
t159 = t186 / 0.2e1 + (t200 * t147 + t149) * t257;
t174 = -m(6) / 0.2e1 - m(7) * t200 / 0.2e1;
t28 = -m(5) - t159 + t174;
t203 = t28 * qJD(1);
t202 = mrSges(7,1) * t189 + t210 * t255;
t197 = qJD(5) * t151;
t196 = t50 * t240;
t195 = t252 / 0.2e1;
t194 = pkin(5) * t253;
t192 = Ifges(7,2) / 0.4e1 - Ifges(7,1) / 0.4e1;
t185 = t144 - t235;
t119 = Ifges(7,2) * t152 + t238;
t178 = Ifges(7,5) * t150 + Ifges(7,6) * t152;
t5 = t37 * t112 + (t136 * t262 + t80 * t246 + t82 * t250 + mrSges(7,3) * t228 + t178 * t249 + (t119 * t251 - t260) * t153) * t153 + (-t110 - t193) * t36;
t175 = -t5 * qJD(1) - t15 * qJD(2);
t173 = qJD(2) * t254 - t50 * qJD(4);
t172 = -t46 * mrSges(7,1) / 0.2e1 + t47 * t255;
t171 = t50 * qJD(2) + qJD(4) * t254;
t170 = t191 - t229 / 0.2e1;
t166 = t119 * t250 + t260;
t9 = t112 * t213 + t153 * t98 - mrSges(5,2) - mrSges(4,3) - t83 + t198 * mrSges(6,3) - m(7) * (t151 * t177 + t216) - t136 * t186 + (-m(5) - m(4)) * t182;
t164 = -t9 * qJD(1) + qJD(2) * t195;
t163 = t144 / 0.2e1 - t235 / 0.2e1 - Ifges(6,4);
t79 = Ifges(7,6) * t153 - t151 * t179;
t81 = Ifges(7,5) * t153 - t151 * t180;
t1 = m(7) * (t36 * t46 + t37 * t47) + t47 * t110 + t37 * t109 + t46 * t112 + t36 * t111 + t137 * t143 + (t136 * t211 + t163 * t153 + t81 * t246 + t79 * t251) * t153 + (t136 * t98 + t231 / 0.2e1 - t227 / 0.2e1 - t137 * mrSges(6,2) - t163 * t151 + (-m(7) * t136 ^ 2 - Ifges(6,1) + Ifges(6,2) + Ifges(7,3)) * t153) * t151;
t162 = t1 * qJD(1) + t6 * qJD(2) + t8 * qJD(4);
t26 = pkin(5) * t118 + t179 * t247 + t180 * t251 + t166;
t155 = -t136 * t118 / 0.2e1 + pkin(8) * t183 + (-t145 / 0.4e1 - t120 / 0.4e1 + t192 * t150) * t150 + (-t119 / 0.4e1 - 0.3e1 / 0.4e1 * t238 - t192 * t152) * t152;
t3 = t144 * t248 + t194 + (t112 * t256 + t82 / 0.4e1 + t237 / 0.2e1) * t152 + (t110 * t256 - t80 / 0.4e1 - 0.3e1 / 0.4e1 * t234) * t150 + (-Ifges(7,3) / 0.2e1 + t155) * t153 + t172;
t39 = -t211 / 0.2e1 + t202;
t40 = (t118 / 0.2e1 + t170) * t153;
t160 = t3 * qJD(1) - t39 * qJD(2) - t40 * qJD(4) - t26 * qJD(5);
t42 = t118 * t245 + t170 * t153;
t41 = t263 + t202;
t35 = t38 * t240 / 0.2e1;
t31 = t159 + t174;
t20 = -t167 + t201;
t16 = t249 * t261 + t157 + t158;
t11 = t156 + t169;
t7 = t8 * qJD(5);
t4 = qJD(3) * t195 + t6 * qJD(5) - t15 * qJD(6);
t2 = t227 / 0.4e1 - t231 / 0.4e1 + t185 * t248 + t194 + Ifges(7,6) * t189 + Ifges(7,5) * t188 + Ifges(7,3) * t244 - t165 * pkin(8) + t155 * t153 - t172;
t12 = [-qJD(3) * t9 + qJD(4) * t17 + qJD(5) * t1 - qJD(6) * t5, t4, t31 * qJD(4) + t16 * qJD(5) + t20 * qJD(6) + t164, qJD(3) * t31 + qJD(6) * t11 + t220 + t7, t16 * qJD(3) + t2 * qJD(6) + (-Ifges(6,5) + (-m(7) * pkin(5) + t222) * t136 + t166) * t197 + t162 + (-mrSges(6,2) * t218 - Ifges(6,6) * t153 + pkin(5) * t211 + t178 * t244 + t79 * t246 + t81 * t250 + (t168 + t209 - t214) * pkin(8) + t176 * mrSges(7,3)) * qJD(5), t20 * qJD(3) + t11 * qJD(4) + t2 * qJD(5) + (-t37 * mrSges(7,1) - t36 * mrSges(7,2) - t153 * t178) * qJD(6) + t175; t4, -t196, qJD(1) * t195, t35, t224 + (-t143 - t262) * qJD(5) + t41 * qJD(6) - t259 * t197 + (-t258 * qJD(5) - t171) * m(7), qJD(5) * t41 - qJD(6) * t262 - t221; t28 * qJD(4) + t14 * qJD(5) - t19 * qJD(6) - t164, -qJD(1) * t252 / 0.2e1, 0, t203, t217, qJD(6) * t118 - t204; -qJD(3) * t28 + qJD(6) * t10 - t220 + t7, t35, -t203, t196, t223 + t42 * qJD(6) + t222 * t197 + t259 * qJD(5) * t153 + ((pkin(8) * t184 - t242) * qJD(5) - t173) * m(7), -qJD(6) * t151 * t181 + t42 * qJD(5) + t219; -qJD(3) * t14 + qJD(6) * t3 - t162, m(7) * t171 - t39 * qJD(6) - t224, -t217, m(7) * t173 - t40 * qJD(6) - t223, -t26 * qJD(6) (-pkin(8) * t181 + t185) * qJD(6) + t160; qJD(3) * t19 - qJD(4) * t10 - qJD(5) * t3 - t175, qJD(5) * t39 + t221, t204, qJD(5) * t40 - t219, -t160, 0;];
Cq  = t12;
