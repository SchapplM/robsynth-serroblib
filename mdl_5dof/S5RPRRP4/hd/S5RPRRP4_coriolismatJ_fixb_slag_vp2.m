% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:31:53
% EndTime: 2022-01-23 09:31:58
% DurationCPUTime: 2.18s
% Computational Cost: add. (5778->280), mult. (12661->361), div. (0->0), fcn. (12427->6), ass. (0->152)
t259 = Ifges(6,4) + Ifges(5,4);
t258 = Ifges(5,5) + Ifges(6,5);
t257 = Ifges(5,6) + Ifges(6,6);
t154 = cos(pkin(8));
t155 = sin(qJ(4));
t156 = sin(qJ(3));
t157 = cos(qJ(4));
t158 = cos(qJ(3));
t140 = -t155 * t158 - t157 * t156;
t153 = sin(pkin(8));
t126 = t140 * t153;
t211 = t126 * mrSges(6,3);
t100 = mrSges(6,2) * t154 + t211;
t212 = t126 * mrSges(5,3);
t101 = mrSges(5,2) * t154 + t212;
t256 = t100 + t101;
t193 = t158 * t153;
t196 = t155 * t156;
t125 = t153 * t196 - t157 * t193;
t213 = t125 * mrSges(6,3);
t102 = -mrSges(6,1) * t154 + t213;
t214 = t125 * mrSges(5,3);
t103 = -mrSges(5,1) * t154 + t214;
t255 = t102 + t103;
t245 = t259 * t125;
t253 = t154 / 0.2e1;
t252 = mrSges(5,1) + mrSges(6,1);
t228 = pkin(3) * t157;
t148 = pkin(4) + t228;
t251 = t148 - t228;
t181 = mrSges(6,3) / 0.2e1 + mrSges(5,3) / 0.2e1;
t165 = -t102 / 0.2e1 - t103 / 0.2e1 + t181 * t125;
t117 = t125 * qJ(5);
t141 = -pkin(2) * t154 - pkin(6) * t153 - pkin(1);
t135 = t158 * t141;
t185 = pkin(7) * t193;
t204 = qJ(2) * t156;
t86 = -t185 + t135 + (-pkin(3) - t204) * t154;
t112 = qJ(2) * t154 * t158 + t156 * t141;
t198 = t153 * t156;
t97 = -pkin(7) * t198 + t112;
t89 = t155 * t97;
t47 = t157 * t86 - t89;
t33 = t117 + t47;
t29 = -pkin(4) * t154 + t33;
t248 = m(6) * (t29 - t33);
t250 = -t248 / 0.2e1 + t165;
t249 = t47 * mrSges(5,3) + t29 * mrSges(6,3) - t259 * t126;
t139 = t157 * t158 - t196;
t219 = mrSges(5,2) + mrSges(6,2);
t167 = -t139 * t219 + t252 * t140;
t229 = pkin(3) * t155;
t247 = Ifges(5,1) + Ifges(6,1);
t246 = Ifges(5,2) + Ifges(6,2);
t171 = t257 * t125 + t258 * t126;
t244 = -t246 + t247;
t233 = t140 / 0.2e1;
t237 = -t125 / 0.2e1;
t243 = t256 * t139 / 0.2e1 + t255 * t233 + (mrSges(6,3) + mrSges(5,3)) * (-t139 * t126 / 0.2e1 + t140 * t237);
t242 = m(5) / 0.2e1;
t241 = m(6) / 0.2e1;
t240 = m(6) * pkin(4);
t203 = qJ(5) * t126;
t91 = t157 * t97;
t111 = -t154 * t204 + t135;
t96 = t111 - t185;
t56 = -t155 * t96 - t91;
t39 = t56 - t203;
t239 = -t39 / 0.2e1;
t48 = t155 * t86 + t91;
t34 = t48 + t203;
t238 = m(6) * t34;
t232 = -t148 / 0.2e1;
t230 = m(6) * t140;
t227 = t33 * mrSges(6,2);
t226 = t34 * mrSges(6,1);
t225 = t39 * mrSges(6,1);
t57 = t157 * t96 - t89;
t40 = t117 + t57;
t224 = t40 * mrSges(6,2);
t223 = t47 * mrSges(5,2);
t222 = t48 * mrSges(5,1);
t221 = t56 * mrSges(5,1);
t220 = t57 * mrSges(5,2);
t136 = t154 * mrSges(4,2) - mrSges(4,3) * t198;
t137 = -t154 * mrSges(4,1) - mrSges(4,3) * t193;
t138 = pkin(3) * t198 + t153 * qJ(2);
t118 = t126 * mrSges(6,2);
t172 = -t125 * mrSges(6,1) + t118;
t87 = -pkin(4) * t126 + t138;
t162 = -t34 * t213 - t48 * t214 - t87 * t172 - t138 * (-mrSges(5,1) * t125 + mrSges(5,2) * t126) + t171 * t253;
t74 = -mrSges(6,1) * t126 - mrSges(6,2) * t125;
t75 = -mrSges(5,1) * t126 - mrSges(5,2) * t125;
t98 = pkin(3) * t193 - t125 * pkin(4);
t1 = t111 * t136 - t112 * t137 + t40 * t100 + t57 * t101 + t39 * t102 + t56 * t103 + t98 * t74 + ((t111 * mrSges(4,3) + Ifges(4,5) * t154 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t156) * t153) * t156 + (-t112 * mrSges(4,3) + Ifges(4,6) * t154 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t158 + (-Ifges(4,1) + Ifges(4,2)) * t156) * t153 + (m(5) * t138 + t75) * pkin(3)) * t158) * t153 - t162 + m(5) * (t47 * t56 + t48 * t57) + m(6) * (t29 * t39 + t34 * t40 + t87 * t98) + t245 * t237 + (-t257 * t154 - t245) * t125 / 0.2e1 + (t247 * t237 - t258 * t253 + (t246 / 0.2e1 - t244 / 0.2e1) * t125 - t249) * t126;
t215 = t1 * qJD(1);
t210 = t156 * mrSges(4,1);
t209 = t158 * mrSges(4,2);
t2 = -t33 * t100 + t48 * t103 - t47 * t101 + (t245 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t154 + (m(6) * t87 + t74) * pkin(4)) * t125 + ((Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1) * t154 + t244 * t125 + t249) * t126 + t162 + (t102 + t248) * t34;
t208 = t2 * qJD(1);
t128 = t139 * t154;
t107 = t128 * t229;
t127 = t140 * t154;
t170 = (-mrSges(6,2) / 0.2e1 - mrSges(5,2) / 0.2e1) * t128;
t182 = mrSges(6,1) / 0.2e1 + mrSges(5,1) / 0.2e1;
t159 = t182 * t127 + t170 + (-t209 / 0.2e1 - t210 / 0.2e1) * t154 + (t127 * t228 + t107) * t242 + (t127 * t148 + t107) * t241;
t194 = t158 * t136;
t195 = t156 * t137;
t161 = -m(5) * ((t47 - t57) * t140 + (t48 + t56) * t139) / 0.2e1 - m(6) * ((t29 - t40) * t140 + (t34 + t39) * t139) / 0.2e1 + t195 / 0.2e1 - t194 / 0.2e1;
t173 = -t100 / 0.2e1 - t101 / 0.2e1;
t164 = (t126 * t181 + t173) * t139;
t5 = t140 * t165 + t159 + t161 + t164;
t207 = t5 * qJD(1);
t186 = t240 / 0.2e1;
t163 = t170 + (t186 + t182) * t127;
t7 = t250 * t140 + t163 + t164;
t206 = t7 * qJD(1);
t151 = t153 ^ 2;
t149 = t151 * qJ(2);
t152 = t154 ^ 2;
t168 = t209 + t210;
t9 = (t194 - t195) * t154 + t256 * t128 + t255 * t127 + (t152 + t151) * mrSges(3,3) + (t153 * t168 + t74 + t75) * t153 + m(6) * (t127 * t29 + t34 * t128 + t87 * t153) + m(5) * (t127 * t47 + t128 * t48 + t138 * t153) + m(4) * (t149 + (-t111 * t156 + t112 * t158) * t154) + m(3) * (qJ(2) * t152 + t149);
t205 = t9 * qJD(1);
t16 = t125 * t102 + t126 * t100 + m(6) * (t125 * t29 + t126 * t34);
t202 = qJD(1) * t16;
t201 = t125 * t139;
t200 = t126 * t140;
t199 = t148 * t125;
t197 = t155 * t125;
t184 = t126 * t229;
t20 = 0.2e1 * (t184 / 0.4e1 + t199 / 0.4e1 - t98 / 0.4e1) * m(6) - t172;
t192 = t20 * qJD(1);
t44 = (t201 / 0.2e1 - t200 / 0.2e1 - t153 / 0.2e1) * m(6);
t191 = t44 * qJD(1);
t179 = mrSges(6,1) + t240;
t59 = t125 * t179 - t118;
t190 = t59 * qJD(1);
t180 = t232 + pkin(4) / 0.2e1;
t176 = -t211 / 0.2e1;
t30 = (t228 / 0.2e1 + t180) * t230;
t160 = ((t238 / 0.2e1 - t212 / 0.2e1 - t173) * t157 + t250 * t155) * pkin(3);
t4 = t180 * t211 + (-t33 / 0.2e1 + t40 / 0.2e1) * mrSges(6,2) + (-t47 / 0.2e1 + t57 / 0.2e1) * mrSges(5,2) + (-t34 / 0.2e1 + t239) * mrSges(6,1) + (-t48 / 0.2e1 - t56 / 0.2e1) * mrSges(5,1) + (pkin(4) * t239 + t232 * t34) * m(6) + t160;
t72 = t219 * t228 + (t251 * m(6) + t252) * t229;
t166 = t4 * qJD(1) - t30 * qJD(2) - t72 * qJD(3);
t131 = t139 * t229;
t43 = (-t200 + t201 + t153) * t241;
t41 = (t184 + t199 + t98) * t241;
t19 = t167 + (pkin(4) + t251) * t230 / 0.2e1;
t8 = t233 * t248 + t163 + t243;
t6 = t159 - t161 + t243;
t3 = -t226 / 0.2e1 - t227 / 0.2e1 - t222 / 0.2e1 - t223 / 0.2e1 + t225 / 0.2e1 - t224 / 0.2e1 + t221 / 0.2e1 - t220 / 0.2e1 + pkin(4) * t176 + t39 * t186 + (-t238 / 0.2e1 + t176) * t148 + t160 + t171;
t10 = [qJD(2) * t9 + qJD(3) * t1 - qJD(4) * t2 + qJD(5) * t16, t205 + 0.2e1 * (t241 + t242) * (t127 * t139 - t128 * t140) * qJD(2) + t6 * qJD(3) + t8 * qJD(4) + t43 * qJD(5), t6 * qJD(2) + t3 * qJD(4) + t41 * qJD(5) + t215 + (-t112 * mrSges(4,1) - t111 * mrSges(4,2) - Ifges(4,5) * t198 - Ifges(4,6) * t193 + t171 - t220 + t221 - t224 + t225 + (mrSges(6,3) * t197 + (-t157 * t126 + t197) * mrSges(5,3) + m(6) * t155 * t40 + m(5) * (t155 * t57 + t157 * t56)) * pkin(3) + (m(6) * t39 - t211) * t148) * qJD(3), -t208 + t8 * qJD(2) + t3 * qJD(3) + (-t222 - t226 - t223 - t227 + (-t211 - t238) * pkin(4) + t171) * qJD(4), qJD(2) * t43 + qJD(3) * t41 + t202; -qJD(3) * t5 - qJD(4) * t7 + qJD(5) * t44 - t205, 0, -t207 + (t167 - t168) * qJD(3) + t19 * qJD(4) + 0.2e1 * ((t140 * t228 + t131) * t242 + (t140 * t148 + t131) * t241) * qJD(3), -t206 + t19 * qJD(3) + (pkin(4) * t230 + t167) * qJD(4), t191; qJD(2) * t5 + qJD(4) * t4 + qJD(5) * t20 - t215, -qJD(4) * t30 + t207, -t72 * qJD(4), (-t219 * t157 + (-mrSges(5,1) - t179) * t155) * qJD(4) * pkin(3) + t166, t192; qJD(2) * t7 - qJD(3) * t4 + qJD(5) * t59 + t208, qJD(3) * t30 + t206, -t166, 0, t190; -qJD(2) * t44 - qJD(3) * t20 - qJD(4) * t59 - t202, -t191, -t192, -t190, 0;];
Cq = t10;
