% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:39
% EndTime: 2019-12-31 18:48:44
% DurationCPUTime: 1.90s
% Computational Cost: add. (5310->198), mult. (10430->250), div. (0->0), fcn. (11608->6), ass. (0->115)
t139 = sin(pkin(8));
t140 = cos(pkin(8));
t218 = sin(qJ(3));
t220 = cos(qJ(3));
t129 = -t139 * t218 + t140 * t220;
t141 = sin(qJ(4));
t155 = t139 * t220 + t140 * t218;
t219 = cos(qJ(4));
t117 = -t129 * t219 + t141 * t155;
t233 = t129 * t141 + t155 * t219;
t186 = t233 ^ 2;
t188 = t117 ^ 2;
t261 = -(Ifges(6,5) - Ifges(5,4)) * (t186 - t188) + (Ifges(5,1) - Ifges(5,2) + Ifges(6,1) - Ifges(6,3)) * t117 * t233;
t230 = m(6) / 0.2e1;
t246 = t233 * mrSges(6,1);
t171 = t117 * mrSges(6,3) + t246;
t255 = t117 * mrSges(5,2);
t172 = t233 * mrSges(5,1) - t255;
t259 = t246 / 0.2e1 - t255 / 0.2e1;
t223 = t117 / 0.2e1;
t256 = pkin(4) * t117;
t254 = t117 * mrSges(6,2);
t212 = pkin(6) + qJ(2);
t169 = t218 * t212;
t170 = t220 * t212;
t84 = (-pkin(7) * t218 - t169) * t140 + (-pkin(7) * t220 - t170) * t139;
t121 = t129 * t212;
t85 = t129 * pkin(7) + t121;
t64 = -t141 * t85 + t219 * t84;
t252 = t64 * t117;
t222 = -t233 / 0.2e1;
t250 = t233 / 0.2e1;
t251 = t250 + t222;
t165 = pkin(4) * t233 + qJ(5) * t117;
t153 = t155 * pkin(3);
t249 = m(5) * t153;
t190 = qJ(5) * t233;
t157 = t141 * t84 + t219 * t85;
t244 = t157 * t233;
t136 = m(6) * qJ(5) + mrSges(6,3);
t243 = qJD(4) * t136;
t242 = t136 * qJD(5);
t241 = mrSges(6,1) + mrSges(5,1);
t239 = Ifges(5,5) + Ifges(6,4);
t236 = -Ifges(5,6) + Ifges(6,6);
t182 = t219 * pkin(3);
t134 = -t182 - pkin(4);
t235 = t134 * t157;
t234 = -t244 + t252;
t231 = t129 ^ 2;
t229 = m(5) * pkin(3);
t228 = m(6) * pkin(3);
t227 = t157 / 0.2e1;
t226 = m(6) * t157;
t224 = -t117 / 0.2e1;
t217 = pkin(3) * t141;
t216 = t157 * mrSges(5,1);
t215 = t157 * mrSges(6,1);
t214 = t64 * mrSges(5,2);
t213 = t64 * mrSges(6,3);
t211 = mrSges(6,1) * t117;
t210 = mrSges(5,3) * t233;
t209 = Ifges(6,4) * t117;
t208 = Ifges(5,5) * t117;
t207 = Ifges(5,6) * t233;
t206 = Ifges(6,6) * t233;
t120 = -t139 * t170 - t140 * t169;
t145 = t155 ^ 2;
t5 = (t145 + t231) * mrSges(4,3) + m(4) * (-t120 * t155 + t121 * t129) + (m(3) * qJ(2) + mrSges(3,3)) * (t139 ^ 2 + t140 ^ 2) + (m(6) + m(5)) * (-t117 * t157 - t233 * t64) + (mrSges(5,3) + mrSges(6,2)) * (t186 + t188);
t205 = qJD(1) * t5;
t147 = t153 + t165;
t151 = mrSges(4,1) * t155 + mrSges(4,2) * t129;
t179 = -pkin(2) * t140 - pkin(1);
t160 = pkin(3) * t129 - t179;
t154 = -t160 + t256;
t152 = -t154 + t190;
t167 = -mrSges(6,3) * t233 + t211;
t1 = -t179 * t151 + t157 * t210 - t64 * t254 + t152 * t171 - (mrSges(5,1) * t117 + mrSges(5,2) * t233) * t153 + (t145 - t231) * Ifges(4,4) + (-Ifges(4,1) + Ifges(4,2)) * t129 * t155 + (t234 - t252) * mrSges(5,3) + (t234 + t244) * mrSges(6,2) + (t172 + t249) * t160 + (m(6) * t152 - t167) * t147 + t261;
t204 = t1 * qJD(1);
t156 = m(6) * t165;
t2 = t160 * t172 - t165 * t167 + (t156 + t171) * t152 + t261;
t195 = t2 * qJD(1);
t142 = t147 * t230 + t249 / 0.2e1;
t133 = qJ(5) + t217;
t144 = (-t117 * t133 + t134 * t233) * t230 + (-t117 * t141 - t219 * t233) * t229 / 0.2e1;
t8 = -t142 + t144 - t151 - t171 - t172;
t192 = t8 * qJD(1);
t148 = t156 / 0.2e1 + t259;
t150 = t165 * t230 + t259;
t9 = 0.2e1 * mrSges(5,1) * t250 + 0.2e1 * mrSges(6,3) * t223 + t148 + t150;
t191 = t9 * qJD(1);
t17 = (-m(6) * t154 + t136 * t233 - t211) * t233;
t189 = qJD(1) * t17;
t69 = t233 * m(6);
t184 = t69 * qJD(1);
t176 = t224 + t223;
t143 = (t235 + t157 * t182 - (-t133 + t217) * t64) * t230;
t146 = t133 * t222 + t134 * t224 + (t141 * t250 + t219 * t224) * pkin(3);
t158 = m(6) * (-pkin(4) * t157 + qJ(5) * t64);
t3 = -t158 / 0.2e1 + (t190 / 0.2e1 - t256 / 0.2e1 + t146) * mrSges(6,2) + t143 + t241 * (-t157 / 0.2e1 + t227) + t239 * t176 + t236 * t251;
t149 = -t241 * t217 + (-mrSges(5,2) + mrSges(6,3)) * t182;
t74 = -(t133 * t219 + t134 * t141) * t228 - t149;
t164 = t3 * qJD(1) - t74 * qJD(3);
t130 = m(6) * t133 + mrSges(6,3);
t16 = t176 * mrSges(6,2);
t161 = -qJD(1) * t16 - qJD(3) * t130;
t159 = -qJD(3) * t136 - t243;
t127 = 0.2e1 * t217 * t230 + t136;
t70 = t251 * m(6);
t20 = -t254 + t226;
t15 = -t254 + t226 / 0.2e1 + m(6) * t227;
t14 = t142 + t144;
t10 = t148 - t150;
t4 = -t209 / 0.2e1 - t208 / 0.2e1 - t207 / 0.2e1 + t206 / 0.2e1 + pkin(4) * t254 / 0.2e1 - t214 / 0.2e1 + t213 / 0.2e1 - t216 / 0.2e1 - t215 / 0.2e1 + t158 / 0.2e1 + (-mrSges(6,1) / 0.2e1 - mrSges(5,1) / 0.2e1) * t157 - (mrSges(5,2) / 0.2e1 - mrSges(6,3) / 0.2e1) * t64 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t233 + (-Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1) * t117 + t143 + (-t190 / 0.2e1 + t146) * mrSges(6,2);
t6 = [qJD(2) * t5 - qJD(3) * t1 - qJD(4) * t2 + qJD(5) * t17, qJD(3) * t14 + qJD(4) * t10 + qJD(5) * t70 + t205, -t204 + t14 * qJD(2) + (t117 * mrSges(5,3) * t182 - t134 * t254 - t133 * t233 * mrSges(6,2) - t210 * t217 - t214 - t215 + t213 - t208 - t207 + m(6) * (t133 * t64 + t235) - t209 + t206 - t216 + (t141 * t64 - t157 * t219) * t229 - t120 * mrSges(4,2) + Ifges(4,5) * t129 - Ifges(4,6) * t155 - t121 * mrSges(4,1)) * qJD(3) + t4 * qJD(4) + t15 * qJD(5), t10 * qJD(2) + t4 * qJD(3) + t20 * qJD(5) - t195 + ((-m(6) * pkin(4) - t241) * t157 - (mrSges(5,2) - t136) * t64 + (-qJ(5) * mrSges(6,2) + t236) * t233 + (pkin(4) * mrSges(6,2) - t239) * t117) * qJD(4), qJD(2) * t70 + qJD(3) * t15 + qJD(4) * t20 + t189; -qJD(3) * t8 + qJD(4) * t9 - qJD(5) * t69 - t205, 0, -t192, t191, -t184; qJD(2) * t8 + qJD(4) * t3 + qJD(5) * t16 + t204, t192, -qJD(4) * t74 + qJD(5) * t130, ((-pkin(4) * t141 + qJ(5) * t219) * t228 + t149) * qJD(4) + t127 * qJD(5) + t164, qJD(4) * t127 - t161; -qJD(2) * t9 - qJD(3) * t3 + t195, -t191, -t164 + t242, t242, -t159; qJD(2) * t69 - qJD(3) * t16 - t189, t184, t161 - t243, t159, 0;];
Cq = t6;
