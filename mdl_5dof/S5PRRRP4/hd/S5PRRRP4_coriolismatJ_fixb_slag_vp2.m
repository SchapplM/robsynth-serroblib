% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:39
% EndTime: 2019-12-05 16:45:42
% DurationCPUTime: 1.20s
% Computational Cost: add. (1548->145), mult. (3647->193), div. (0->0), fcn. (2938->6), ass. (0->101)
t120 = cos(qJ(4));
t118 = t120 ^ 2;
t218 = Ifges(6,5) - Ifges(5,4);
t119 = sin(qJ(4));
t217 = t119 ^ 2 + t118;
t161 = qJ(5) * t120;
t186 = t119 * pkin(4);
t203 = t161 - t186;
t82 = -t120 * mrSges(6,1) - t119 * mrSges(6,3);
t131 = -t203 * t82 - t218 * t118 + ((-Ifges(5,2) - Ifges(6,3) + Ifges(5,1) + Ifges(6,1)) * t120 + t218 * t119) * t119;
t191 = cos(qJ(3));
t132 = t217 * t191;
t212 = t203 / 0.2e1;
t211 = mrSges(6,2) + mrSges(5,3);
t194 = m(6) * t203;
t152 = -t194 / 0.2e1;
t151 = t191 * pkin(2);
t109 = -t151 - pkin(3);
t169 = t120 * mrSges(5,2);
t175 = t119 * mrSges(5,1);
t85 = t169 + t175;
t178 = t109 * t85;
t159 = t119 * qJ(5);
t140 = t120 * pkin(4) + t159;
t80 = -pkin(3) - t140;
t65 = -t151 + t80;
t182 = t65 + t80;
t193 = pkin(3) * t85;
t168 = t120 * mrSges(6,3);
t174 = t119 * mrSges(6,1);
t84 = -t168 + t174;
t195 = t84 / 0.2e1;
t210 = t193 / 0.2e1 - t178 / 0.2e1 - (t152 + t195) * t182 - t131;
t209 = -t120 * mrSges(5,1) + t119 * mrSges(5,2) + t82;
t208 = t132 * pkin(2);
t189 = sin(qJ(3));
t190 = sin(qJ(2));
t192 = cos(qJ(2));
t77 = t189 * t190 - t191 * t192;
t207 = t217 * t77;
t122 = (m(6) * t140 - t209) * qJD(4);
t206 = mrSges(4,1) - t209;
t149 = t189 * t77;
t78 = -t189 * t192 - t191 * t190;
t205 = (-t132 * t78 + t149) * pkin(2);
t202 = m(5) / 0.4e1 + m(6) / 0.4e1;
t199 = 2 * qJD(3);
t198 = m(5) / 0.2e1;
t197 = m(6) / 0.2e1;
t196 = m(6) * pkin(2);
t188 = m(6) * t119;
t187 = m(6) * t120;
t185 = t65 * t78;
t150 = t189 * pkin(2);
t108 = t150 + pkin(7);
t184 = t207 * t108;
t183 = t207 * pkin(7);
t179 = t109 * t78;
t111 = t120 * mrSges(6,2);
t7 = 0.4e1 * t202 * (-0.1e1 + t217) * t78 * t77;
t163 = t7 * qJD(1);
t148 = t84 - t194;
t8 = t148 * t65 + t131 + t178;
t162 = t8 * qJD(2);
t38 = (m(6) * t65 + t82) * t119;
t160 = qJD(2) * t38;
t13 = (-t85 / 0.2e1 - t84 / 0.2e1 + (mrSges(5,2) / 0.2e1 - mrSges(6,3) / 0.2e1) * t120 + (mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1) * t119 + (t186 / 0.2e1 - t161 / 0.2e1 + t212) * m(6)) * t77;
t158 = t13 * qJD(1);
t157 = t208 * pkin(7);
t110 = m(6) * qJ(5) + mrSges(6,3);
t156 = t110 * qJD(4);
t155 = qJD(2) + qJD(3);
t14 = (t85 / 0.2e1 + t195 + 0.2e1 * t152 + t169 / 0.2e1 + t175 / 0.2e1 + t174 / 0.2e1 - t168 / 0.2e1) * t77;
t36 = t77 * t188;
t154 = t14 * qJD(4) - t36 * qJD(5) + t163;
t153 = -t13 * qJD(4) - t163;
t142 = t151 / 0.2e1;
t126 = -mrSges(4,2) * t151 - t206 * t150 + t211 * t208;
t129 = t132 * t108;
t16 = -(t189 * t65 + t129) * t196 - m(5) * (t189 * t109 + t129) * pkin(2) - t126;
t133 = (-t78 * t80 - t183) * t197 + (pkin(3) * t78 - t183) * t198;
t134 = m(5) * (-t179 - t184);
t135 = m(6) * (-t184 - t185);
t2 = -t134 / 0.2e1 - t135 / 0.2e1 - 0.2e1 * t202 * t205 + t133;
t138 = -t2 * qJD(1) - t16 * qJD(2);
t15 = t148 * t80 + t131 - t193;
t121 = t142 * t168 + t196 * t191 * t212 - (t169 + (mrSges(5,1) + mrSges(6,1)) * t119) * t151 / 0.2e1;
t3 = t121 + t210;
t137 = -t3 * qJD(2) + t15 * qJD(3);
t21 = (t82 + (t142 + t65 / 0.2e1 + t80 / 0.2e1) * m(6)) * t119;
t39 = (m(6) * t80 + t82) * t119;
t136 = qJD(2) * t21 + qJD(3) * t39;
t130 = t77 * mrSges(4,2) + t206 * t78 - t211 * t207;
t125 = (-mrSges(6,2) * t159 - pkin(4) * t111 + (Ifges(6,4) + Ifges(5,5)) * t120 + (-Ifges(5,6) + Ifges(6,6)) * t119) * qJD(4) - t158;
t124 = t205 - t184;
t81 = pkin(7) * t187 + t111;
t64 = t108 * t187 + t111;
t22 = -t182 * t188 / 0.2e1 + (m(6) * t142 - t82) * t119;
t4 = t121 - t210;
t1 = (t124 - t185) * t197 + (t124 - t179) * t198 + t130 + t133;
t5 = [t155 * t7, (-t192 * mrSges(3,2) - t190 * mrSges(3,1) + m(4) * (t191 * t78 - t149) * pkin(2) + t134 + t135 + t130) * qJD(2) + t1 * qJD(3) + t154, t1 * qJD(2) + qJD(3) * t130 + t133 * t199 + t154, t155 * t14 + (-qJD(5) * t187 + t122) * t78, -qJD(4) * t78 * t187 - t155 * t36; -qJD(3) * t2 + t153, -qJD(3) * t16 + qJD(4) * t8 - qJD(5) * t38, t126 * qJD(3) + t4 * qJD(4) + t22 * qJD(5) + ((t150 * t80 + t157) * t197 + (-pkin(3) * t150 + t157) * t198) * t199 + t138, t4 * qJD(3) + t64 * qJD(5) - t108 * t122 + t125 + t162, qJD(3) * t22 + qJD(4) * t64 - t160; qJD(2) * t2 + t153, -qJD(4) * t3 - qJD(5) * t21 - t138, qJD(4) * t15 - qJD(5) * t39, -pkin(7) * t122 + t81 * qJD(5) + t125 + t137, qJD(4) * t81 - t136; t155 * t13, qJD(3) * t3 + t158 - t162, -t137 + t158, t110 * qJD(5), t156; 0, qJD(3) * t21 + t160, t136, -t156, 0;];
Cq = t5;
