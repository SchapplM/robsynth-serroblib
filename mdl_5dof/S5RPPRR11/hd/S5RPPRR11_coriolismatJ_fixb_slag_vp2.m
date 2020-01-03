% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:33
% EndTime: 2019-12-31 18:05:36
% DurationCPUTime: 1.14s
% Computational Cost: add. (1715->205), mult. (3567->280), div. (0->0), fcn. (2567->4), ass. (0->118)
t137 = cos(qJ(5));
t127 = Ifges(6,5) * t137;
t135 = sin(qJ(5));
t216 = -t135 * Ifges(6,6) + t127;
t136 = sin(qJ(4));
t190 = t137 * mrSges(6,1);
t195 = t135 * mrSges(6,2);
t162 = t190 - t195;
t220 = t136 * t162;
t219 = Ifges(5,4) - t216;
t131 = qJ(2) - pkin(6);
t218 = t131 * t136;
t189 = t137 * mrSges(6,2);
t196 = t135 * mrSges(6,1);
t109 = t189 + t196;
t217 = t136 * t109;
t129 = t135 ^ 2;
t130 = t137 ^ 2;
t171 = t129 + t130;
t202 = pkin(7) * t136;
t138 = cos(qJ(4));
t203 = pkin(4) * t138;
t114 = t202 + t203;
t179 = t135 * t138;
t55 = t114 * t137 - t131 * t179;
t174 = t137 * t138;
t56 = t114 * t135 + t131 * t174;
t215 = -t135 * t55 + t137 * t56;
t128 = Ifges(6,4) * t137;
t199 = Ifges(6,2) * t135;
t214 = t199 - t128;
t201 = Ifges(6,4) * t135;
t110 = Ifges(6,2) * t137 + t201;
t111 = Ifges(6,1) * t135 + t128;
t206 = -t137 / 0.2e1;
t207 = -t135 / 0.2e1;
t213 = (Ifges(6,1) * t137 - t201) * t207 + t135 * t110 / 0.2e1 + t111 * t206;
t212 = t136 ^ 2;
t211 = t138 ^ 2;
t210 = m(6) / 0.2e1;
t209 = m(6) * (-0.1e1 + t171) * t138 * t136;
t192 = t136 * mrSges(6,1);
t104 = -mrSges(6,3) * t174 + t192;
t208 = -t104 / 0.2e1;
t205 = -t138 / 0.2e1;
t204 = pkin(4) * t136;
t132 = pkin(1) + qJ(3);
t200 = Ifges(6,5) * t136;
t198 = Ifges(6,6) * t136;
t105 = -pkin(7) * t138 + t132 + t204;
t180 = t135 * t136;
t45 = t105 * t137 - t131 * t180;
t178 = t136 * t137;
t46 = t105 * t135 + t131 * t178;
t158 = -t135 * t45 + t137 * t46;
t165 = m(5) * (t211 + t212);
t169 = mrSges(6,3) * t179;
t191 = t136 * mrSges(6,2);
t102 = -t169 - t191;
t176 = t137 * t102;
t181 = t135 * t104;
t7 = -mrSges(4,2) - mrSges(3,3) + (mrSges(5,3) + t109) * t211 + (-m(6) * t211 - t165) * t131 + (-m(6) * t158 + mrSges(5,3) * t136 - t176 + t181) * t136 + (-m(4) - m(3)) * qJ(2);
t197 = qJD(1) * t7;
t85 = t162 * t138;
t4 = t46 * t104 + (t131 * t85 + (-Ifges(6,4) * t179 + t200) * t135 + (Ifges(6,4) * t174 + t46 * mrSges(6,3) + t198 + (Ifges(6,1) - Ifges(6,2)) * t179) * t137) * t138 + (-t102 - t169) * t45;
t187 = t4 * qJD(1);
t175 = t137 * t104;
t183 = t135 * t102;
t151 = -t183 / 0.2e1 - t175 / 0.2e1;
t164 = (-t130 / 0.2e1 - t129 / 0.2e1) * mrSges(6,3);
t140 = (t138 * t164 + t151) * t136 + t85 * t205;
t153 = -t195 / 0.2e1 + t190 / 0.2e1;
t8 = t140 - t153;
t186 = t8 * qJD(1);
t163 = -t136 * mrSges(5,1) - t138 * mrSges(5,2);
t13 = t183 + t175 + mrSges(4,3) + m(6) * (t135 * t46 + t137 * t45) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t132 - t163;
t185 = qJD(1) * t13;
t126 = t138 * mrSges(5,1);
t101 = -mrSges(6,2) * t138 + mrSges(6,3) * t180;
t103 = mrSges(6,1) * t138 + mrSges(6,3) * t178;
t142 = (-t135 * t56 - t137 * t55) * t210 + t101 * t207 + t103 * t206;
t143 = (t171 * t202 + t203) * t210 - t162 * t205;
t10 = -t126 + (mrSges(5,2) + t164) * t136 + t142 - t143;
t184 = t10 * qJD(1);
t182 = t135 * t103;
t177 = t137 * t101;
t14 = (-t191 / 0.2e1 + t102 / 0.2e1) * t137 + (-t192 / 0.2e1 + t208) * t135;
t173 = t14 * qJD(1);
t141 = t165 / 0.2e1 + (t171 * t212 + t211) * t210;
t157 = -m(5) / 0.2e1 - m(6) * t171 / 0.2e1;
t24 = -m(4) - t141 + t157;
t172 = t24 * qJD(1);
t170 = -pkin(4) * t85 / 0.2e1;
t166 = t171 * t138;
t160 = -Ifges(6,5) * t135 - Ifges(6,6) * t137;
t1 = -m(6) * (t45 * t55 + t46 * t56) - t56 * t102 - t46 * t101 - t55 * t104 - t45 * t103 - t132 * (-t136 * mrSges(5,2) + t126) - t219 * t212 + (t219 * t138 + (t130 * Ifges(6,1) + Ifges(5,1) - Ifges(5,2) - Ifges(6,3) + (t199 - 0.2e1 * t128) * t135) * t136 + (m(6) * t218 - 0.2e1 * t217) * t131) * t138;
t150 = t181 / 0.2e1 - t176 / 0.2e1;
t152 = m(6) * t215;
t6 = (t152 / 0.2e1 + t177 / 0.2e1 - t182 / 0.2e1) * t136 + ((t158 - 0.2e1 * t218) * t210 + t217 - t150) * t138;
t159 = -t1 * qJD(1) + t6 * qJD(3);
t156 = -t55 * mrSges(6,1) / 0.2e1 + t56 * mrSges(6,2) / 0.2e1;
t155 = t6 * qJD(1) + qJD(3) * t209;
t154 = -t196 / 0.2e1 - t189 / 0.2e1;
t147 = -t110 / 0.4e1 + (Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.4e1) * t137;
t19 = pkin(4) * t109 - t214 * t206 + t213;
t144 = -t131 * t109 / 0.2e1 + pkin(7) * t164;
t145 = -t111 / 0.4e1 - t128 / 0.4e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.4e1) * t135;
t2 = t170 + t151 * pkin(7) + t216 * t136 + (-Ifges(6,3) / 0.2e1 + t147 * t137 + (-0.5e1 / 0.4e1 * t128 + t145) * t135 + t144) * t138 + t156;
t35 = (t109 / 0.2e1 + t154) * t138;
t146 = t2 * qJD(1) - t35 * qJD(3) - t19 * qJD(4);
t36 = t109 * t205 + t154 * t138;
t25 = t141 + t157;
t15 = t136 * t154 + t150;
t11 = t142 + t143 + t171 * t136 * mrSges(6,3) / 0.2e1;
t9 = t140 + t153;
t5 = t6 * qJD(4);
t3 = t170 + t136 * t127 / 0.4e1 - Ifges(6,5) * t178 / 0.2e1 + Ifges(6,6) * t180 / 0.2e1 + (-pkin(7) * t102 / 0.2e1 - t198 / 0.2e1) * t135 + (pkin(7) * t208 + t200 / 0.4e1) * t137 - t156 + (Ifges(6,3) / 0.2e1 + t144 + t145 * t135 + (-0.5e1 / 0.4e1 * t201 + t147) * t137) * t138;
t12 = [-qJD(2) * t7 + qJD(3) * t13 - qJD(4) * t1 - qJD(5) * t4, qJD(3) * t25 + qJD(4) * t11 + qJD(5) * t15 - t197, qJD(2) * t25 + qJD(5) * t9 + t185 + t5, t11 * qJD(2) + t3 * qJD(5) + t159 + (pkin(4) * t217 + (t152 + t177 - t182) * pkin(7) + (t137 * t214 / 0.2e1 - Ifges(5,5) + (-m(6) * pkin(4) - mrSges(5,1) - t162) * t131 + t213) * t136 + (-t131 * mrSges(5,2) - Ifges(5,6) - t160) * t138 + t215 * mrSges(6,3)) * qJD(4), -t187 + t15 * qJD(2) + t9 * qJD(3) + t3 * qJD(4) + (-mrSges(6,1) * t46 - mrSges(6,2) * t45 + t138 * t160) * qJD(5); qJD(3) * t24 + qJD(4) * t10 - qJD(5) * t14 + t197, 0, t172, t184, qJD(5) * t109 - t173; -qJD(2) * t24 + qJD(5) * t8 - t185 + t5, -t172, qJD(4) * t209, (-t220 + m(6) * (pkin(7) * t166 - t204) + mrSges(6,3) * t166 + t163) * qJD(4) + t36 * qJD(5) + t155, t36 * qJD(4) - qJD(5) * t220 + t186; -qJD(2) * t10 + qJD(5) * t2 - t159, -t184, -t35 * qJD(5) - t155, -t19 * qJD(5), (-pkin(7) * t162 + t216) * qJD(5) + t146; qJD(2) * t14 - qJD(3) * t8 - qJD(4) * t2 + t187, t173, qJD(4) * t35 - t186, -t146, 0;];
Cq = t12;
