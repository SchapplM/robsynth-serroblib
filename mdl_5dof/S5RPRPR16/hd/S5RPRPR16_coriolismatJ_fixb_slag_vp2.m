% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR16_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:46
% EndTime: 2019-12-31 18:38:49
% DurationCPUTime: 1.19s
% Computational Cost: add. (2056->229), mult. (3931->305), div. (0->0), fcn. (2990->4), ass. (0->133)
t222 = -mrSges(6,1) / 0.2e1;
t221 = mrSges(6,2) / 0.2e1;
t130 = cos(qJ(3));
t117 = t130 * qJ(4);
t128 = sin(qJ(3));
t166 = -t128 * pkin(3) + t117;
t98 = -t128 * mrSges(5,2) - t130 * mrSges(5,3);
t220 = m(5) * t166 - t98;
t94 = qJ(2) - t166;
t219 = m(5) * t94 + t98;
t203 = m(6) * t128;
t132 = -pkin(1) - pkin(6);
t217 = -pkin(4) + t132;
t129 = cos(qJ(5));
t196 = Ifges(6,6) * t129;
t127 = sin(qJ(5));
t200 = Ifges(6,5) * t127;
t216 = Ifges(5,6) + Ifges(4,4) - t200 / 0.2e1 - t196 / 0.2e1;
t125 = t127 ^ 2;
t126 = t129 ^ 2;
t165 = t125 + t126;
t167 = pkin(3) * t130 + qJ(4) * t128;
t84 = pkin(7) * t130 + t167;
t95 = t217 * t128;
t30 = -t127 * t84 + t95 * t129;
t188 = t129 * t30;
t31 = t95 * t127 + t129 * t84;
t194 = t127 * t31;
t215 = -t194 - t188;
t173 = t129 * t128;
t174 = t127 * t128;
t169 = t173 * t222 + t174 * t221;
t182 = t130 * mrSges(6,2);
t87 = mrSges(6,3) * t173 - t182;
t185 = t129 * t87;
t184 = t130 * mrSges(6,1);
t85 = -mrSges(6,3) * t174 + t184;
t192 = t127 * t85;
t75 = pkin(7) * t128 + t94;
t96 = t217 * t130;
t28 = -t127 * t75 - t129 * t96;
t29 = -t127 * t96 + t129 * t75;
t214 = m(6) * (t127 * t28 - t129 * t29) - t185 + t192 - t219;
t213 = m(6) / 0.2e1;
t189 = t129 * mrSges(6,2);
t195 = t127 * mrSges(6,1);
t97 = t189 + t195;
t77 = t97 * t128;
t212 = t77 / 0.2e1;
t211 = t95 / 0.2e1;
t209 = (0.1e1 - t165) * t130 * t203;
t112 = Ifges(6,4) * t173;
t208 = -t112 / 0.4e1;
t207 = -t127 / 0.2e1;
t206 = t127 / 0.2e1;
t204 = t129 / 0.2e1;
t202 = Ifges(6,4) * t127;
t201 = Ifges(6,4) * t129;
t199 = Ifges(6,5) * t129;
t198 = Ifges(6,2) * t129;
t197 = Ifges(6,6) * t127;
t181 = t130 * Ifges(6,5);
t61 = Ifges(6,1) * t174 + t112 + t181;
t193 = t127 * t61;
t170 = t130 * t129;
t88 = mrSges(6,2) * t128 + mrSges(6,3) * t170;
t191 = t127 * t88;
t190 = t128 * mrSges(6,3);
t154 = t198 + t202;
t180 = t130 * Ifges(6,6);
t59 = t154 * t128 + t180;
t187 = t129 * t59;
t171 = t130 * t127;
t86 = -mrSges(6,1) * t128 - mrSges(6,3) * t171;
t186 = t129 * t86;
t183 = t130 * mrSges(4,2);
t101 = Ifges(6,1) * t129 - t202;
t151 = t127 * t29 + t129 * t28;
t153 = -t197 + t199;
t4 = -t95 * t77 - t28 * t87 + t29 * t85 + (t59 * t206 - t130 * t153 / 0.2e1 + (t101 * t207 + t198 * t206) * t128 + t151 * mrSges(6,3) - (t112 + t61) * t129 / 0.2e1) * t128;
t179 = t4 * qJD(1);
t142 = t185 / 0.2e1 - t192 / 0.2e1;
t158 = t125 / 0.2e1 + t126 / 0.2e1;
t147 = t158 * t190;
t133 = (t147 - t142) * t130 + t128 * t212;
t143 = t189 / 0.2e1 + t195 / 0.2e1;
t8 = t133 + t143;
t178 = t8 * qJD(1);
t157 = -t128 * mrSges(4,1) - t183;
t11 = mrSges(3,3) + (m(4) + m(3)) * qJ(2) - t157 - t214;
t177 = qJD(1) * t11;
t12 = t214 * t130;
t176 = qJD(1) * t12;
t175 = qJD(5) * t97;
t13 = (t182 / 0.2e1 - t87 / 0.2e1) * t129 + (t184 / 0.2e1 + t85 / 0.2e1) * t127 + t147;
t172 = t13 * qJD(1);
t164 = qJD(3) * t128;
t161 = -Ifges(6,2) / 0.4e1 + Ifges(6,1) / 0.4e1;
t160 = qJ(4) * t212;
t156 = mrSges(6,1) * t129 - mrSges(6,2) * t127;
t155 = Ifges(6,1) * t127 + t201;
t60 = -Ifges(6,6) * t128 + t154 * t130;
t62 = -Ifges(6,5) * t128 + t155 * t130;
t76 = t130 * t156;
t1 = t31 * t87 + t29 * t88 - t95 * t76 + t30 * t85 + t28 * t86 + m(6) * (t28 * t30 + t29 * t31 + t95 * t96) + t219 * t167 + (t187 / 0.2e1 + t193 / 0.2e1 - t94 * mrSges(5,2) + qJ(2) * mrSges(4,1) - t216 * t130) * t130 + (t60 * t204 + t62 * t206 - t96 * t156 + t94 * mrSges(5,3) - qJ(2) * mrSges(4,2) + (-Ifges(4,1) + Ifges(4,2) - Ifges(5,2) + Ifges(5,3) - Ifges(6,3)) * t130 + t216 * t128) * t128;
t144 = -t191 / 0.2e1 - t186 / 0.2e1;
t6 = ((t215 + t95) * t213 + t144) * t130 + (-t76 + (t151 + t96) * t213 + t87 * t206 + t85 * t204) * t128;
t152 = t1 * qJD(1) + t6 * qJD(2);
t131 = -pkin(3) - pkin(7);
t148 = t158 * t131 * mrSges(6,3);
t146 = t31 * t221 + t30 * t222;
t145 = t6 * qJD(1) + qJD(2) * t209;
t100 = -Ifges(6,2) * t127 + t201;
t141 = t100 * t204 + t101 * t206;
t140 = t128 * t156;
t139 = -t100 / 0.4e1 - t161 * t127;
t138 = t101 / 0.4e1 + t161 * t129;
t21 = -qJ(4) * t156 + t154 * t207 + t155 * t204 + t141;
t3 = t160 + (Ifges(6,3) / 0.2e1 - t148) * t128 + (-t131 * t85 / 0.2e1 + t208 - t95 * mrSges(6,2) / 0.2e1 - t61 / 0.4e1 - 0.3e1 / 0.4e1 * t181 + t139 * t128) * t127 + (t131 * t87 / 0.2e1 - t59 / 0.4e1 + mrSges(6,1) * t211 - 0.3e1 / 0.4e1 * t180 + (-0.3e1 / 0.4e1 * t202 + t138) * t128) * t129 + t146;
t32 = -t140 / 0.2e1 - t169;
t135 = t3 * qJD(1) - t32 * qJD(2) - t21 * qJD(3);
t16 = 0.2e1 * (t95 / 0.4e1 - t194 / 0.4e1 - t188 / 0.4e1) * m(6) + t144 + t169;
t41 = (-0.1e1 / 0.2e1 + t158) * t203;
t73 = mrSges(5,3) + (m(5) + m(6)) * qJ(4) + t97;
t134 = qJD(1) * t16 - qJD(2) * t41 + qJD(3) * t73;
t34 = t203 / 0.2e1 + (t165 * t213 + m(5)) * t128;
t33 = t140 / 0.2e1 - t169;
t14 = t143 * t130 + t142 - t165 * t190 / 0.2e1;
t10 = m(6) * t211 - t215 * t213 + (m(5) * t132 - mrSges(5,1)) * t128 - t144 + t169;
t9 = t133 - t143;
t5 = t6 * qJD(3);
t2 = -t187 / 0.4e1 + t127 * t208 + t156 * t211 - t193 / 0.4e1 + t160 + t130 * (-t196 - t200) / 0.4e1 + Ifges(6,6) * t170 / 0.2e1 + Ifges(6,5) * t171 / 0.2e1 + t142 * t131 - t146 + (-Ifges(6,3) / 0.2e1 - t148 + t138 * t129 + (-0.3e1 / 0.4e1 * t201 + t139) * t127) * t128;
t7 = [qJD(2) * t11 + qJD(3) * t1 + qJD(4) * t12 - qJD(5) * t4, qJD(5) * t9 + t177 + t5, t10 * qJD(4) + t2 * qJD(5) + (Ifges(5,4) - Ifges(4,5) + pkin(3) * mrSges(5,1) - t199 / 0.2e1 + t197 / 0.2e1) * t164 + t152 + (t60 * t207 + t96 * t97 + t62 * t204 + (Ifges(5,5) - Ifges(4,6) + t141) * t130 + (t157 + t220) * t132 + (m(6) * t96 - t130 * mrSges(5,1) - t76) * qJ(4) + (-m(6) * t215 + t186 + t191) * t131 + t215 * mrSges(6,3)) * qJD(3), qJD(3) * t10 + qJD(5) * t14 + t176, -t179 + t9 * qJD(2) + t2 * qJD(3) + t14 * qJD(4) + (-mrSges(6,1) * t29 - mrSges(6,2) * t28 + t128 * t153) * qJD(5); qJD(5) * t8 - t177 + t5, qJD(3) * t209, t34 * qJD(4) + t33 * qJD(5) + (-t165 * mrSges(6,3) - mrSges(4,1)) * t164 + t145 + (t130 * t97 - t183 + 0.2e1 * (t165 * t131 * t128 + t117) * t213 + t220) * qJD(3), t34 * qJD(3), t33 * qJD(3) + t130 * t175 + t178; qJD(4) * t16 + qJD(5) * t3 - t152, -t41 * qJD(4) - t32 * qJD(5) - t145, qJD(4) * t73 - qJD(5) * t21, t134, ((-mrSges(6,2) * t131 - Ifges(6,6)) * t129 + (-mrSges(6,1) * t131 - Ifges(6,5)) * t127) * qJD(5) + t135; -qJD(3) * t16 - qJD(5) * t13 - t176, t41 * qJD(3), -t134, 0, -t172 - t175; -qJD(2) * t8 - qJD(3) * t3 + qJD(4) * t13 + t179, qJD(3) * t32 - t178, -t135, t172, 0;];
Cq = t7;
