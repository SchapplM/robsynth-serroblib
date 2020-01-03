% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:57
% EndTime: 2019-12-31 18:19:01
% DurationCPUTime: 1.71s
% Computational Cost: add. (4299->223), mult. (8397->328), div. (0->0), fcn. (8627->8), ass. (0->133)
t181 = sin(pkin(9));
t182 = cos(pkin(9));
t210 = sin(qJ(3));
t211 = cos(qJ(3));
t111 = t181 * t210 - t182 * t211;
t108 = t111 * mrSges(5,2);
t112 = -t181 * t211 - t182 * t210;
t224 = -t112 * mrSges(5,1) - t108;
t135 = sin(qJ(5));
t136 = cos(qJ(5));
t158 = mrSges(6,1) * t136 - t135 * mrSges(6,2);
t177 = t112 * t158;
t130 = Ifges(6,4) * t136;
t118 = Ifges(6,1) * t135 + t130;
t132 = t135 ^ 2;
t133 = t136 ^ 2;
t173 = -t132 - t133;
t126 = sin(pkin(8)) * pkin(1) + pkin(6);
t166 = t210 * t126;
t110 = -t210 * qJ(4) - t166;
t120 = t211 * t126;
t174 = t211 * qJ(4) + t120;
t78 = -t182 * t110 + t181 * t174;
t172 = t210 * pkin(3);
t83 = -t112 * pkin(4) + t111 * pkin(7) + t172;
t34 = t135 * t78 + t136 * t83;
t35 = t135 * t83 - t136 * t78;
t153 = -t34 * t135 + t35 * t136;
t157 = Ifges(6,2) * t135 - t130;
t223 = t181 * t110 + t182 * t174;
t222 = t112 ^ 2;
t221 = m(5) / 0.2e1;
t220 = m(6) / 0.2e1;
t219 = -mrSges(6,1) / 0.2e1;
t218 = mrSges(6,2) / 0.2e1;
t204 = mrSges(6,3) * t112;
t170 = t136 * t204;
t196 = t111 * mrSges(6,1);
t86 = t170 + t196;
t217 = -t86 / 0.2e1;
t190 = t136 * mrSges(6,2);
t194 = t135 * mrSges(6,1);
t115 = t190 + t194;
t216 = t115 / 0.2e1;
t215 = -t135 / 0.2e1;
t214 = t135 / 0.2e1;
t213 = -t136 / 0.2e1;
t212 = t136 / 0.2e1;
t209 = t35 * mrSges(6,2);
t205 = mrSges(6,3) * t111;
t202 = Ifges(6,4) * t135;
t201 = Ifges(6,5) * t112;
t129 = Ifges(6,5) * t136;
t199 = Ifges(6,2) * t136;
t198 = Ifges(6,6) * t112;
t197 = Ifges(6,6) * t135;
t165 = t181 * pkin(3);
t125 = t165 + pkin(7);
t195 = t125 * mrSges(6,3);
t60 = Ifges(6,6) * t111 + t157 * t112;
t193 = t135 * t60;
t178 = t111 * t136;
t85 = -t112 * mrSges(6,1) + mrSges(6,3) * t178;
t192 = t135 * t85;
t191 = t135 * t86;
t168 = -cos(pkin(8)) * pkin(1) - pkin(2);
t113 = -t211 * pkin(3) + t168;
t72 = t111 * pkin(4) + t112 * pkin(7) + t113;
t31 = t135 * t72 + t136 * t223;
t189 = t136 * t31;
t119 = Ifges(6,1) * t136 - t202;
t62 = Ifges(6,5) * t111 - t112 * t119;
t188 = t136 * t62;
t179 = t111 * t135;
t84 = mrSges(6,2) * t112 + mrSges(6,3) * t179;
t187 = t136 * t84;
t184 = t78 * t112;
t138 = (t135 * t35 + t136 * t34) * t220 + t84 * t214 + t85 * t212 + t172 * t221;
t167 = t182 * pkin(3);
t127 = -t167 - pkin(4);
t149 = (t173 * t125 * t111 - t112 * t127) * t220 + (-t181 * t111 + t182 * t112) * pkin(3) * t221;
t162 = t132 / 0.2e1 + t133 / 0.2e1;
t139 = -t162 * t205 + t149;
t9 = t108 + (t158 / 0.2e1 + mrSges(5,1)) * t112 - t138 + t139;
t183 = t9 * qJD(1);
t169 = t196 / 0.2e1;
t11 = ((-t132 / 0.2e1 + t162) * t204 + (t217 + t169) * t136) * t112;
t180 = t11 * qJD(1);
t13 = mrSges(6,2) * t178 + (t169 + t86 / 0.2e1 - t170 / 0.2e1) * t135;
t176 = t13 * qJD(1);
t21 = m(6) * (-0.1e1 - t173) * t112 * t111;
t22 = t21 / 0.2e1;
t175 = t22 * qJD(1);
t171 = t135 * t204;
t164 = t179 / 0.2e1;
t163 = -t178 / 0.2e1;
t160 = t129 - t197;
t116 = t199 + t202;
t156 = Ifges(6,5) * t135 + Ifges(6,6) * t136;
t143 = t210 * mrSges(4,1) + t211 * mrSges(4,2);
t30 = -t135 * t223 + t136 * t72;
t59 = t157 * t111 - t198;
t61 = -t119 * t111 - t201;
t82 = t111 * t115;
t1 = -t78 * t82 + t31 * t84 + t30 * t85 + t34 * t86 + t168 * t143 + m(6) * (t78 * t223 + t30 * t34 + t31 * t35) + (mrSges(5,1) * t172 + t193 / 0.2e1 - t188 / 0.2e1 - t209 + (Ifges(5,4) - t129 / 0.2e1 + t197 / 0.2e1) * t111) * t111 + (-mrSges(5,2) * t172 - Ifges(5,4) * t112 + (-t61 / 0.2e1 - t223 * mrSges(6,2) + t201 / 0.2e1) * t136 + (t59 / 0.2e1 + t35 * mrSges(6,3) - t223 * mrSges(6,1) - t198 / 0.2e1) * t135 + (-Ifges(5,2) + Ifges(5,1) - Ifges(6,3)) * t111) * t112 + (-Ifges(4,2) + Ifges(4,1)) * t211 * t210 + (-t210 ^ 2 + t211 ^ 2) * Ifges(4,4) + (m(5) * t172 + t224) * t113;
t154 = t135 * t30 - t189;
t147 = -t111 * mrSges(6,2) + t171;
t142 = t136 * t147;
t63 = t111 * t142;
t64 = t222 * t115;
t6 = -t63 / 0.2e1 + t64 / 0.2e1 + (t191 / 0.2e1 - t82 / 0.2e1 + (t154 + t223) * t220) * t111 + (-t187 / 0.2e1 + t192 / 0.2e1 + (-t153 - t78) * t220) * t112;
t155 = t1 * qJD(1) + t6 * qJD(2);
t3 = t158 * t184 + t31 * t86 + (t62 * t215 + t60 * t213 - t111 * t156 / 0.2e1 - mrSges(6,3) * t189 + (t116 * t215 + t118 * t212) * t112) * t112 + (t171 - t147) * t30;
t152 = -t3 * qJD(1) - t11 * qJD(2);
t151 = t6 * qJD(1) + t21 * qJD(2);
t7 = -t86 * t179 + t63 - t64 - m(6) * (t154 * t111 - t184) - m(5) * (-t111 * t223 - t184) + (-t111 ^ 2 - t222) * mrSges(5,3);
t150 = -qJD(1) * t7 + qJD(2) * t22;
t148 = t34 * t219 + t209 / 0.2e1;
t146 = -t193 / 0.4e1 + t78 * t216;
t145 = t190 / 0.2e1 + t194 / 0.2e1;
t144 = t112 * t156;
t141 = t145 * t111;
t38 = t127 * t115 + (t118 / 0.2e1 - t157 / 0.2e1) * t136 + (t119 / 0.2e1 - t116 / 0.2e1) * t135;
t42 = (-t115 / 0.2e1 + t145) * t111;
t137 = t162 * t195 + (-t119 / 0.4e1 + t116 / 0.4e1 + t127 * t219 + t199 / 0.4e1) * t136 + (t118 / 0.4e1 - t157 / 0.4e1 + t130 / 0.2e1 + t127 * t218 + (-t195 / 0.2e1 + Ifges(6,1) / 0.4e1) * t135) * t135;
t5 = (t125 * t217 + t62 / 0.4e1) * t136 + (0.3e1 / 0.4e1 * t129 + (t125 * t218 - 0.3e1 / 0.4e1 * Ifges(6,6)) * t135) * t111 + (Ifges(6,3) / 0.2e1 + t137) * t112 + t146 + t148;
t140 = t5 * qJD(1) - t42 * qJD(2) + t38 * qJD(3);
t43 = t111 * t216 + t141;
t14 = -t191 / 0.2e1 + t142 / 0.2e1 + t141;
t10 = t177 / 0.2e1 + t138 + t139;
t4 = t188 / 0.4e1 + t111 * t160 / 0.4e1 + Ifges(6,5) * t163 + Ifges(6,6) * t164 + (mrSges(6,2) * t164 + t86 * t213) * t125 + t146 - t148 + (-Ifges(6,3) / 0.2e1 + t137) * t112;
t2 = qJD(3) * t6 + qJD(4) * t22 - qJD(5) * t11;
t8 = [qJD(3) * t1 - qJD(4) * t7 - qJD(5) * t3, t2, (-Ifges(5,5) * t111 + Ifges(5,6) * t112 + Ifges(4,5) * t211 - Ifges(4,6) * t210 + t118 * t163 + t116 * t164 + mrSges(4,2) * t166 - mrSges(4,1) * t120 + t61 * t214 + t59 * t212 - t223 * t158 - t144 / 0.2e1 + m(5) * (-t181 * t78 - t182 * t223) * pkin(3) + t78 * mrSges(5,2) - t223 * mrSges(5,1) + (m(6) * t223 - t82) * t127 + (m(6) * t153 + t187 - t192) * t125 + t153 * mrSges(6,3) + (t111 * t167 + t112 * t165) * mrSges(5,3)) * qJD(3) + t10 * qJD(4) + t4 * qJD(5) + t155, qJD(3) * t10 + qJD(5) * t14 + t150, t4 * qJD(3) + t14 * qJD(4) + (-t31 * mrSges(6,1) - t30 * mrSges(6,2) + t144) * qJD(5) + t152; t2, t21 * qJD(3), t43 * qJD(5) + t151 + (t173 * t205 - t143 + 0.2e1 * t149 + t177 - t224) * qJD(3), t175, t43 * qJD(3) + qJD(5) * t177 - t180; qJD(4) * t9 + qJD(5) * t5 - t155, -qJD(5) * t42 - t151, t38 * qJD(5), t183, (-t158 * t125 + t160) * qJD(5) + t140; -qJD(3) * t9 - qJD(5) * t13 - t150, -t175, -t183, 0, -qJD(5) * t115 - t176; -qJD(3) * t5 + qJD(4) * t13 - t152, qJD(3) * t42 + t180, -t140, t176, 0;];
Cq = t8;
