% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:38
% EndTime: 2019-03-09 01:41:41
% DurationCPUTime: 1.98s
% Computational Cost: add. (6325->283), mult. (12080->383), div. (0->0), fcn. (12670->8), ass. (0->164)
t155 = sin(qJ(6));
t150 = t155 ^ 2;
t157 = cos(qJ(6));
t151 = t157 ^ 2;
t197 = t150 + t151;
t259 = t197 * mrSges(7,3);
t152 = sin(pkin(10));
t153 = cos(pkin(10));
t156 = sin(qJ(4));
t238 = cos(qJ(4));
t141 = t238 * t152 + t156 * t153;
t139 = t152 * t156 - t238 * t153;
t208 = qJ(5) * t139;
t99 = pkin(4) * t141 + t208;
t247 = m(6) * t99;
t204 = t139 * t157;
t97 = -t141 * mrSges(7,2) + mrSges(7,3) * t204;
t212 = t157 * t97;
t205 = t139 * t155;
t95 = t141 * mrSges(7,1) - mrSges(7,3) * t205;
t220 = t155 * t95;
t169 = t220 / 0.2e1 - t212 / 0.2e1;
t256 = mrSges(7,3) * (t150 / 0.2e1 + t151 / 0.2e1);
t258 = -t139 * t256 - t169;
t217 = t157 * mrSges(7,2);
t225 = t155 * mrSges(7,1);
t143 = t217 + t225;
t91 = t139 * t143;
t193 = -t91 / 0.2e1;
t235 = -mrSges(5,1) + mrSges(6,2);
t218 = t157 * mrSges(7,1);
t224 = t155 * mrSges(7,2);
t142 = -t218 + t224;
t90 = t141 * t142;
t245 = -t139 / 0.2e1;
t255 = t245 * t259;
t226 = Ifges(7,6) * t157;
t228 = Ifges(7,5) * t155;
t172 = t228 / 0.2e1 + t226 / 0.2e1;
t254 = Ifges(6,6) + Ifges(5,4) - t172;
t253 = 2 * m(7);
t252 = 0.2e1 * t141;
t251 = 2 * qJD(4);
t250 = m(6) / 0.2e1;
t249 = m(7) / 0.2e1;
t248 = pkin(4) + pkin(8);
t94 = t197 * t141;
t246 = m(7) * t94;
t244 = t141 / 0.2e1;
t243 = t142 / 0.2e1;
t242 = -t155 / 0.2e1;
t241 = t155 / 0.2e1;
t240 = -t157 / 0.2e1;
t239 = t157 / 0.2e1;
t237 = m(6) * t141;
t236 = m(7) * t141;
t234 = mrSges(5,3) + mrSges(6,1);
t146 = sin(pkin(9)) * pkin(1) + qJ(3);
t233 = pkin(7) + t146;
t232 = mrSges(7,3) * t141;
t231 = Ifges(7,4) * t155;
t230 = Ifges(7,4) * t157;
t229 = Ifges(7,5) * t141;
t227 = Ifges(7,6) * t141;
t134 = t233 * t153;
t191 = t233 * t152;
t80 = t238 * t134 - t156 * t191;
t56 = -t139 * pkin(5) + t80;
t70 = t248 * t141 + t208;
t37 = -t155 * t70 + t157 * t56;
t223 = t155 * t37;
t38 = t155 * t56 + t157 * t70;
t222 = t155 * t38;
t131 = Ifges(7,4) * t204;
t68 = Ifges(7,1) * t205 + t131 + t229;
t221 = t155 * t68;
t219 = t155 * t97;
t216 = t157 * t37;
t215 = t157 * t38;
t183 = Ifges(7,2) * t157 + t231;
t66 = t183 * t139 + t227;
t214 = t157 * t66;
t213 = t157 * t95;
t211 = t248 * t95;
t210 = t248 * t97;
t96 = -t139 * mrSges(7,1) - t155 * t232;
t98 = t139 * mrSges(7,2) + t157 * t232;
t168 = t98 * t239 + t96 * t242;
t135 = t139 * mrSges(6,3);
t136 = t139 * mrSges(5,2);
t199 = t136 - t135;
t53 = -t248 * t94 - t208;
t9 = t193 + (-t256 + t235) * t141 + (t53 / 0.4e1 + t223 / 0.4e1 - t215 / 0.4e1) * t253 - t247 - t168 + t199;
t209 = t9 * qJD(1);
t100 = -mrSges(6,2) * t139 - mrSges(6,3) * t141;
t175 = -cos(pkin(9)) * pkin(1) - t153 * pkin(3) - pkin(2);
t159 = -t141 * qJ(5) + t175;
t52 = t248 * t139 + t159;
t79 = t134 * t156 + t238 * t191;
t54 = pkin(5) * t141 + t79;
t31 = -t155 * t52 + t157 * t54;
t32 = t155 * t54 + t157 * t52;
t71 = pkin(4) * t139 + t159;
t11 = (m(7) * (t155 * t31 - t157 * t32) - t212 + t220 - m(6) * t71 - t100) * t141;
t207 = qJD(1) * t11;
t12 = t258 * t139 + t91 * t244;
t206 = qJD(1) * t12;
t162 = (t217 / 0.2e1 + t225 / 0.2e1) * t141;
t15 = t162 - t258;
t203 = t15 * qJD(1);
t170 = t224 / 0.2e1 - t218 / 0.2e1;
t163 = t170 * t141;
t167 = t219 / 0.2e1 + t213 / 0.2e1;
t17 = -t163 + t167;
t202 = t17 * qJD(1);
t164 = m(7) * (-t141 + t94) * t139;
t28 = t164 / 0.2e1;
t201 = t28 * qJD(1);
t188 = -m(7) * t197 / 0.4e1;
t34 = -t246 / 0.2e1 + (t188 - m(6) / 0.2e1) * t252;
t200 = t34 * qJD(1);
t196 = qJD(4) * t141;
t195 = -Ifges(7,2) / 0.4e1 + Ifges(7,1) / 0.4e1;
t145 = Ifges(7,1) * t157 - t231;
t184 = Ifges(7,1) * t155 + t230;
t67 = -Ifges(7,6) * t139 + t183 * t141;
t69 = -Ifges(7,5) * t139 + t184 * t141;
t89 = t142 * t139;
t1 = -t54 * t89 + t56 * t90 + t37 * t95 + t31 * t96 + t38 * t97 + t32 * t98 + t99 * t100 - t175 * t136 + m(7) * (t31 * t37 + t32 * t38 - t54 * t56) + (t214 / 0.2e1 + t221 / 0.2e1 + t175 * mrSges(5,1) - t254 * t141) * t141 + (t67 * t239 + t69 * t241 + (Ifges(6,3) - Ifges(5,1) - Ifges(6,2) + Ifges(5,2) - Ifges(7,3)) * t141 + t254 * t139) * t139 + (-t141 * mrSges(6,2) + t135 + t247) * t71;
t166 = t96 * t239 + t98 * t241;
t179 = t216 + t222;
t180 = t155 * t32 + t157 * t31;
t4 = ((t180 - t54) * t249 + t90 / 0.2e1 + t167) * t141 + ((t179 - t56) * t249 - t89 / 0.2e1 + t166) * t139;
t182 = t1 * qJD(1) + t4 * qJD(2);
t130 = Ifges(7,5) * t204;
t3 = t56 * t91 + t130 * t244 + t31 * t97 - t32 * t95 + ((t68 / 0.2e1 + t131 / 0.2e1 - t31 * mrSges(7,3)) * t157 + (-t66 / 0.2e1 - t227 / 0.2e1 - t32 * mrSges(7,3) + (-t231 / 0.2e1 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.2e1) * t157) * t139) * t155) * t139;
t178 = t3 * qJD(1) + t12 * qJD(2);
t177 = t4 * qJD(1) + qJD(2) * t164;
t7 = (t234 * t139 - t89) * t139 + (t234 * t141 + t213 + t219) * t141 + m(7) * (-t56 * t139 + t180 * t141) + (m(4) * t146 + mrSges(4,3)) * (t152 ^ 2 + t153 ^ 2) + (m(6) + m(5)) * (-t139 * t80 + t141 * t79);
t176 = qJD(1) * t7 + qJD(2) * t28;
t174 = t53 * t249 - t250 * t99;
t173 = t37 * mrSges(7,1) / 0.2e1 - t38 * mrSges(7,2) / 0.2e1;
t171 = qJ(5) * t193 + t56 * t243;
t144 = -Ifges(7,2) * t155 + t230;
t165 = t144 * t239 + t145 * t241;
t41 = -qJ(5) * t142 + t183 * t241 + t184 * t240 - t165;
t43 = (t243 - t170) * t141;
t5 = (-Ifges(7,3) / 0.2e1 - t248 * t256) * t139 + (0.3e1 / 0.4e1 * t229 + t68 / 0.4e1 + t131 / 0.4e1 - t211 / 0.2e1 + (t144 / 0.4e1 + t195 * t155) * t139) * t155 + (0.3e1 / 0.4e1 * t227 + t66 / 0.4e1 + t210 / 0.2e1 + (0.3e1 / 0.4e1 * t231 - t145 / 0.4e1 - t195 * t157) * t139) * t157 + t171 + t173;
t161 = t5 * qJD(1) + t43 * qJD(2) - t41 * qJD(4);
t133 = mrSges(6,3) + (m(6) + m(7)) * qJ(5) + t143;
t14 = t170 * t139 + (t56 / 0.4e1 - t222 / 0.4e1 - t216 / 0.4e1) * t253 - t166;
t77 = t246 / 0.2e1;
t49 = t77 - t236 / 0.2e1;
t160 = qJD(1) * t14 - qJD(2) * t49 + qJD(4) * t133;
t44 = -t90 / 0.2e1 - t163;
t40 = t237 + t236 / 0.2e1 + t77;
t33 = t237 / 0.2e1 + t77 + (t188 - m(6) / 0.4e1) * t252;
t18 = -t163 - t167;
t16 = t162 - t169 + t255;
t13 = m(6) * t80 + (-mrSges(6,1) + t170) * t139 + t166 + (t179 + t56) * t249;
t10 = (t215 - t223) * t249 + t247 / 0.2e1 + t193 - t141 * t256 + t168 + t174;
t6 = t141 * (-t226 - t228) / 0.4e1 - t221 / 0.4e1 - t214 / 0.4e1 - t155 * (-Ifges(7,2) * t205 + t131) / 0.4e1 - t210 * t239 - t211 * t242 + Ifges(7,3) * t245 + t172 * t141 - t171 + t173 + (t145 / 0.2e1 - t183 / 0.4e1) * t204 - t255 * t248 - (t184 + t144) * t205 / 0.4e1;
t2 = qJD(3) * t28 + qJD(4) * t4 + qJD(6) * t12;
t8 = [qJD(3) * t7 + qJD(4) * t1 + qJD(5) * t11 + qJD(6) * t3, t2, qJD(4) * t10 + qJD(5) * t33 + qJD(6) * t18 + t176, t10 * qJD(3) + t13 * qJD(5) + t6 * qJD(6) + ((-qJ(5) * t54 - t179 * t248) * t249 + (-pkin(4) * t80 - qJ(5) * t79) * t250) * t251 + (-qJ(5) * mrSges(6,1) + Ifges(6,5) - Ifges(5,6) + t165) * t196 + t182 + (qJ(5) * t90 - t54 * t143 + (t69 / 0.2e1 - t248 * t96 - t37 * mrSges(7,3)) * t157 + (-t67 / 0.2e1 - t248 * t98 - t38 * mrSges(7,3)) * t155 + (pkin(4) * mrSges(6,1) + Ifges(7,5) * t240 + Ifges(7,6) * t241 + Ifges(6,4) - Ifges(5,5)) * t139 + t235 * t80 + (mrSges(5,2) - mrSges(6,3)) * t79) * qJD(4), qJD(3) * t33 + qJD(4) * t13 + qJD(6) * t16 + t207, t18 * qJD(3) + t6 * qJD(4) + t16 * qJD(5) + (-mrSges(7,1) * t32 - mrSges(7,2) * t31 - Ifges(7,6) * t205 + t130) * qJD(6) + t178; t2, t164 * qJD(4), t201 (t199 - t91) * qJD(4) + t40 * qJD(5) + t44 * qJD(6) + t174 * t251 + (t235 - t259) * t196 + t177, t40 * qJD(4), qJD(4) * t44 - qJD(6) * t91 + t206; -qJD(4) * t9 + qJD(5) * t34 - qJD(6) * t17 - t176, -t201, 0, -t209, t200, qJD(6) * t142 - t202; qJD(3) * t9 + qJD(5) * t14 - qJD(6) * t5 - t182, -qJD(5) * t49 - qJD(6) * t43 - t177, t209, qJD(5) * t133 + qJD(6) * t41, t160 ((mrSges(7,2) * t248 - Ifges(7,6)) * t157 + (mrSges(7,1) * t248 - Ifges(7,5)) * t155) * qJD(6) - t161; -qJD(3) * t34 - qJD(4) * t14 - qJD(6) * t15 - t207, t49 * qJD(4), -t200, -t160, 0, -qJD(6) * t143 - t203; qJD(3) * t17 + qJD(4) * t5 + qJD(5) * t15 - t178, qJD(4) * t43 - t206, t202, t161, t203, 0;];
Cq  = t8;
