% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR14_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:33
% EndTime: 2019-12-31 18:34:37
% DurationCPUTime: 1.70s
% Computational Cost: add. (4393->233), mult. (8156->324), div. (0->0), fcn. (8377->6), ass. (0->142)
t135 = sin(pkin(8));
t137 = sin(qJ(3));
t188 = cos(pkin(8));
t216 = cos(qJ(3));
t111 = t135 * t137 - t188 * t216;
t112 = t135 * t216 + t188 * t137;
t235 = t112 * mrSges(5,1) - t111 * mrSges(5,2);
t127 = t137 * pkin(3) + qJ(2);
t234 = m(5) * t127;
t233 = m(6) * t112;
t138 = cos(qJ(5));
t130 = Ifges(6,5) * t138;
t136 = sin(qJ(5));
t209 = Ifges(6,6) * t136;
t232 = Ifges(5,4) - t130 / 0.2e1 + t209 / 0.2e1;
t131 = Ifges(6,4) * t138;
t212 = Ifges(6,1) * t136;
t121 = t131 + t212;
t133 = t136 ^ 2;
t134 = t138 ^ 2;
t179 = t133 + t134;
t231 = t111 * t135 + t188 * t112;
t178 = t216 * pkin(3);
t76 = -t111 * pkin(4) + t112 * pkin(7) + t178;
t139 = -pkin(1) - pkin(6);
t115 = (-qJ(4) + t139) * t137;
t175 = t216 * t139;
t116 = -t216 * qJ(4) + t175;
t81 = t115 * t135 - t188 * t116;
t35 = t136 * t81 + t138 * t76;
t36 = t136 * t76 - t138 * t81;
t163 = -t35 * t136 + t36 * t138;
t168 = t136 * Ifges(6,2) - t131;
t91 = t111 ^ 2;
t230 = t112 ^ 2;
t229 = m(5) / 0.2e1;
t228 = -m(6) / 0.2e1;
t227 = m(6) / 0.2e1;
t226 = m(5) * pkin(3);
t198 = t138 * mrSges(6,1);
t204 = t136 * mrSges(6,2);
t117 = -t198 + t204;
t74 = t111 * t117;
t225 = t74 / 0.2e1;
t223 = -t111 / 0.2e1;
t222 = -t112 / 0.2e1;
t197 = t138 * mrSges(6,2);
t205 = t136 * mrSges(6,1);
t118 = t197 + t205;
t221 = t118 / 0.2e1;
t220 = -t136 / 0.2e1;
t219 = t136 / 0.2e1;
t218 = -t138 / 0.2e1;
t217 = t138 / 0.2e1;
t211 = Ifges(6,4) * t136;
t210 = Ifges(6,2) * t138;
t150 = t188 * t115 + t135 * t116;
t75 = pkin(4) * t112 + pkin(7) * t111 + t127;
t34 = t136 * t75 + t138 * t150;
t196 = t138 * t34;
t33 = -t136 * t150 + t138 * t75;
t164 = t136 * t33 - t196;
t207 = t111 * mrSges(6,3);
t177 = t136 * t207;
t77 = -t112 * mrSges(6,2) + t177;
t194 = t138 * t77;
t78 = t112 * mrSges(6,1) + t138 * t207;
t200 = t136 * t78;
t206 = t111 * t81;
t7 = (mrSges(5,3) * t112 - t194 + t200) * t112 + (mrSges(5,3) + t118) * t91 + m(6) * (t164 * t112 - t206) + m(5) * (-t112 * t150 - t206);
t208 = qJD(1) * t7;
t60 = t112 * Ifges(6,6) + t168 * t111;
t202 = t136 * t60;
t201 = t136 * t77;
t199 = t137 * mrSges(4,1);
t122 = Ifges(6,1) * t138 - t211;
t62 = Ifges(6,5) * t112 - t111 * t122;
t195 = t138 * t62;
t193 = t138 * t78;
t119 = t210 + t211;
t167 = Ifges(6,5) * t136 + Ifges(6,6) * t138;
t4 = t34 * t78 - t81 * t74 + (t62 * t220 + t60 * t218 - mrSges(6,3) * t196 + t167 * t222 + (t119 * t220 + t121 * t217) * t111) * t111 + (-t77 + t177) * t33;
t190 = t4 * qJD(1);
t106 = t112 * mrSges(5,2);
t126 = -t188 * pkin(3) - pkin(4);
t125 = pkin(3) * t135 + pkin(7);
t170 = t179 * t125;
t141 = (-t126 * t111 - t112 * t170) * t227 + t117 * t223 + (t188 * t111 - t112 * t135) * t226 / 0.2e1 + t179 * mrSges(6,3) * t222;
t183 = t112 * t136;
t158 = mrSges(6,2) * t111 + mrSges(6,3) * t183;
t182 = t112 * t138;
t159 = -mrSges(6,1) * t111 + mrSges(6,3) * t182;
t142 = (t136 * t36 + t138 * t35) * t227 + t158 * t219 + t159 * t217 + t178 * t229;
t9 = t111 * mrSges(5,1) + t106 + t141 - t142;
t189 = t9 * qJD(1);
t147 = -t216 * mrSges(4,2) - t199 - t235;
t165 = t136 * t34 + t138 * t33;
t14 = t201 + t193 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(6) * t165 + t234 - t147;
t187 = qJD(1) * t14;
t155 = -t201 / 0.2e1 - t193 / 0.2e1;
t169 = mrSges(6,3) * (t134 / 0.2e1 + t133 / 0.2e1);
t145 = (t111 * t169 + t155) * t112 + t111 * t225;
t156 = -t204 / 0.2e1 + t198 / 0.2e1;
t10 = t145 - t156;
t186 = t10 * qJD(1);
t184 = t112 * t117;
t153 = t197 / 0.2e1 + t205 / 0.2e1;
t149 = t153 * t112;
t154 = t200 / 0.2e1 - t194 / 0.2e1;
t15 = t149 + t154;
t181 = t15 * qJD(1);
t146 = (-t179 * t230 - t91) * t227 + (-t91 - t230) * t229;
t161 = -m(5) / 0.2e1 + t179 * t228;
t21 = t146 + t161;
t180 = t21 * qJD(1);
t174 = t183 / 0.2e1;
t173 = -t182 / 0.2e1;
t171 = t130 - t209;
t59 = -Ifges(6,6) * t111 + t168 * t112;
t61 = -Ifges(6,5) * t111 - t122 * t112;
t1 = -t127 * t106 + t35 * t78 + t36 * t77 + m(6) * (t81 * t150 + t33 * t35 + t34 * t36) + (-t127 * mrSges(5,1) - t33 * mrSges(6,1) + t34 * mrSges(6,2) - t232 * t111 - t150 * t118 + t61 * t218 + t59 * t219) * t111 + (t202 / 0.2e1 - t195 / 0.2e1 - t81 * t118 + t165 * mrSges(6,3) + (-Ifges(5,2) + Ifges(5,1) - Ifges(6,3)) * t111 + t232 * t112) * t112 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t137) * t137 + (qJ(2) * mrSges(4,1) + (-Ifges(4,1) + Ifges(4,2)) * t137 - Ifges(4,4) * t216) * t216 + (t234 + t235) * t178;
t151 = t112 * t118;
t6 = -(t163 + t81) * t233 / 0.2e1 + (t151 / 0.2e1 + (t164 + t150) * t228 - t154) * t111;
t166 = t1 * qJD(1) - t6 * qJD(2);
t24 = (0.1e1 - t179) * t111 * t233;
t162 = -t6 * qJD(1) + t24 * qJD(2);
t160 = -t35 * mrSges(6,1) / 0.2e1 + t36 * mrSges(6,2) / 0.2e1;
t152 = t111 * t167;
t143 = t125 * t169 + (-t122 / 0.4e1 + t119 / 0.4e1 + t210 / 0.4e1) * t138 + (t121 / 0.4e1 - t168 / 0.4e1 + t131 / 0.2e1 + t212 / 0.4e1) * t136;
t144 = t155 * t125 + t126 * t225 - t202 / 0.4e1 + t195 / 0.4e1 + t81 * t221;
t3 = (-0.3e1 / 0.4e1 * t209 + 0.3e1 / 0.4e1 * t130) * t112 + (Ifges(6,3) / 0.2e1 + t143) * t111 + t144 + t160;
t38 = t126 * t118 + (t121 / 0.2e1 - t168 / 0.2e1) * t138 + (t122 / 0.2e1 - t119 / 0.2e1) * t136;
t42 = (-t118 / 0.2e1 + t153) * t111;
t148 = t3 * qJD(1) - t42 * qJD(2) + t38 * qJD(3);
t43 = (t153 + t221) * t111;
t20 = t146 - t161;
t16 = t149 - t154;
t12 = t141 + t142;
t11 = t145 + t156;
t5 = t6 * qJD(3);
t2 = t112 * t171 / 0.4e1 + Ifges(6,5) * t173 + Ifges(6,6) * t174 + Ifges(6,3) * t223 + t143 * t111 + t144 - t160;
t8 = [qJD(2) * t14 + qJD(3) * t1 + qJD(4) * t7 - qJD(5) * t4, qJD(4) * t20 + qJD(5) * t11 + t187 - t5, (-Ifges(5,5) * t112 + Ifges(5,6) * t111 - Ifges(4,5) * t137 - Ifges(4,6) * t216 + t121 * t173 + t119 * t174 - mrSges(4,2) * t175 - t139 * t199 + t61 * t219 + t59 * t217 - t126 * t151 - t152 / 0.2e1 - (t135 * t226 - mrSges(5,2)) * t81 + t231 * mrSges(5,3) * pkin(3) + (m(6) * t126 - t188 * t226 - mrSges(5,1) + t117) * t150 + (m(6) * t163 - t136 * t159 + t138 * t158) * t125 + t163 * mrSges(6,3)) * qJD(3) + t12 * qJD(4) + t2 * qJD(5) + t166, qJD(2) * t20 + qJD(3) * t12 + qJD(5) * t16 + t208, -t190 + t11 * qJD(2) + t2 * qJD(3) + t16 * qJD(4) + (-mrSges(6,1) * t34 - mrSges(6,2) * t33 + t152) * qJD(5); qJD(4) * t21 + qJD(5) * t10 - t187 - t5, t24 * qJD(3), (-t231 * t226 + m(6) * (-t111 * t170 + t126 * t112) + t184 + t147 - t179 * t207) * qJD(3) + t43 * qJD(5) + t162, t180, t43 * qJD(3) + qJD(5) * t184 + t186; qJD(4) * t9 + qJD(5) * t3 - t166, -qJD(5) * t42 - t162, t38 * qJD(5), t189, (t117 * t125 + t171) * qJD(5) + t148; -qJD(2) * t21 - qJD(3) * t9 - qJD(5) * t15 - t208, -t180, -t189, 0, -qJD(5) * t118 - t181; -qJD(2) * t10 - qJD(3) * t3 + qJD(4) * t15 + t190, qJD(3) * t42 - t186, -t148, t181, 0;];
Cq = t8;
