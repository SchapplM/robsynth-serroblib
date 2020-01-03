% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:44
% EndTime: 2019-12-31 17:57:46
% DurationCPUTime: 1.22s
% Computational Cost: add. (4037->182), mult. (7996->260), div. (0->0), fcn. (8254->8), ass. (0->112)
t129 = sin(pkin(9));
t130 = cos(pkin(9));
t132 = sin(qJ(4));
t194 = cos(qJ(4));
t112 = t129 * t132 - t130 * t194;
t131 = sin(qJ(5));
t127 = t131 ^ 2;
t133 = cos(qJ(5));
t128 = t133 ^ 2;
t208 = t128 + t127;
t158 = t208 * t112;
t210 = mrSges(6,3) * t158;
t156 = t133 * mrSges(6,1) - t131 * mrSges(6,2);
t209 = -t156 - mrSges(5,1);
t121 = sin(pkin(8)) * pkin(1) + qJ(3);
t190 = pkin(6) + t121;
t110 = t190 * t130;
t161 = t190 * t129;
t84 = t110 * t132 + t161 * t194;
t114 = t129 * t194 + t132 * t130;
t192 = pkin(4) * t114;
t90 = t112 * pkin(7) + t192;
t38 = t131 * t84 + t133 * t90;
t39 = t131 * t90 - t133 * t84;
t152 = -t131 * t38 + t133 * t39;
t124 = Ifges(6,4) * t133;
t188 = Ifges(6,2) * t131;
t207 = t188 - t124;
t206 = t112 ^ 2;
t204 = m(6) / 0.2e1;
t203 = -mrSges(6,1) / 0.2e1;
t202 = mrSges(6,2) / 0.2e1;
t71 = -pkin(7) * t158 - t192;
t201 = m(6) * t71;
t200 = t114 / 0.2e1;
t182 = t133 * mrSges(6,2);
t185 = t131 * mrSges(6,1);
t116 = t182 + t185;
t199 = t116 / 0.2e1;
t198 = -t131 / 0.2e1;
t197 = t131 / 0.2e1;
t196 = -t133 / 0.2e1;
t195 = t133 / 0.2e1;
t193 = m(6) * t114;
t191 = pkin(4) * t116;
t189 = Ifges(6,4) * t131;
t123 = Ifges(6,5) * t133;
t187 = Ifges(6,6) * t131;
t85 = t110 * t194 - t132 * t161;
t186 = t114 * t85;
t183 = t131 * t39;
t181 = t133 * t38;
t179 = t84 * t114;
t176 = t112 * t131;
t142 = -t114 * mrSges(6,2) + mrSges(6,3) * t176;
t175 = t112 * t133;
t144 = t114 * mrSges(6,1) + mrSges(6,3) * t175;
t135 = t142 * t197 + t144 * t195;
t111 = t112 * mrSges(5,2);
t159 = t114 * mrSges(5,1) - t111;
t163 = t156 * t200;
t11 = -t163 + 0.2e1 * (t71 / 0.4e1 - t183 / 0.4e1 - t181 / 0.4e1) * m(6) - t135 - t159 - t210 / 0.2e1;
t177 = t11 * qJD(1);
t174 = t114 * t131;
t14 = t112 * t116;
t173 = t14 * qJD(1);
t22 = (0.1e1 - t208) * t112 * t193;
t23 = t22 / 0.2e1;
t172 = t23 * qJD(1);
t169 = pkin(7) * t203;
t168 = pkin(7) * t202;
t167 = pkin(7) * mrSges(6,3) / 0.2e1;
t166 = mrSges(6,3) * t174;
t165 = -Ifges(6,1) / 0.4e1 + Ifges(6,2) / 0.4e1;
t160 = Ifges(5,4) - t123;
t157 = mrSges(6,3) * (-t128 / 0.2e1 - t127 / 0.2e1);
t120 = Ifges(6,1) * t133 - t189;
t155 = -Ifges(6,5) * t131 - Ifges(6,6) * t133;
t139 = -cos(pkin(8)) * pkin(1) - pkin(2) - t130 * pkin(3);
t143 = t112 * mrSges(6,2) + t166;
t145 = -t133 * t114 * mrSges(6,3) + t112 * mrSges(6,1);
t80 = t112 * pkin(4) - t114 * pkin(7) + t139;
t34 = -t131 * t85 + t133 * t80;
t35 = t131 * t80 + t133 * t85;
t1 = mrSges(5,3) * t186 - m(6) * (t34 * t38 + t35 * t39 + t84 * t85) + t39 * t143 - t35 * t142 - t38 * t145 - t34 * t144 - t139 * t159 + (-t160 - t187) * t206 + (-t85 * mrSges(5,3) + t160 * t114 + (t128 * Ifges(6,1) + Ifges(5,1) - Ifges(5,2) - Ifges(6,3)) * t112 + (Ifges(6,6) * t114 + (t188 - 0.2e1 * t124) * t112) * t131) * t114 + (t84 * t112 - t186) * t116;
t153 = t131 * t34 - t133 * t35;
t6 = (t152 + t84) * t193 / 0.2e1 + (t153 + t85) * t112 * t204;
t154 = -t1 * qJD(1) + t6 * qJD(2);
t3 = -t156 * t179 + t35 * t145 + ((-Ifges(6,4) * t174 + Ifges(6,5) * t112) * t131 + (t35 * mrSges(6,3) + Ifges(6,6) * t112 + (t124 + (Ifges(6,1) - Ifges(6,2)) * t131) * t114) * t133) * t114 + (t143 - t166) * t34;
t151 = t3 * qJD(1);
t150 = t6 * qJD(1) + t22 * qJD(2);
t7 = m(6) * (t112 * t153 + t179) + m(5) * (-t85 * t112 + t179) + (m(4) * t121 + mrSges(4,3)) * (t129 ^ 2 + t130 ^ 2) + (t114 ^ 2 + t206) * (mrSges(5,3) + t116);
t149 = -qJD(1) * t7 - qJD(2) * t23;
t148 = pkin(7) * t157;
t147 = t202 * t39 + t203 * t38;
t141 = t182 / 0.2e1 + t185 / 0.2e1;
t140 = t112 * t123 / 0.4e1 + t84 * t199;
t138 = t141 * t112;
t119 = Ifges(6,1) * t131 + t124;
t134 = (t207 / 0.4e1 - t119 / 0.4e1 - t124 + pkin(4) * t202 + (t167 + t165) * t131) * t114;
t117 = Ifges(6,2) * t133 + t189;
t136 = (-t117 / 0.4e1 + t120 / 0.4e1 + pkin(4) * t203 + (t167 - t165) * t133) * t114;
t4 = (-Ifges(6,3) / 0.2e1 + t148) * t114 + ((t169 + 0.3e1 / 0.4e1 * Ifges(6,5)) * t112 + t136) * t133 + ((t168 - Ifges(6,6)) * t112 + t134) * t131 + t140 + t147;
t51 = -t191 + (t119 / 0.2e1 - t207 / 0.2e1) * t133 + (t120 / 0.2e1 - t117 / 0.2e1) * t131;
t52 = (-t116 / 0.2e1 + t141) * t112;
t137 = t4 * qJD(1) - t52 * qJD(2) + t51 * qJD(4);
t53 = t112 * t199 + t138;
t15 = t143 * t196 + t145 * t198 + t138;
t12 = (t181 + t183) * t204 + t201 / 0.2e1 - t163 + t112 * t157 + t135;
t5 = Ifges(6,6) * t176 / 0.2e1 - Ifges(6,5) * t175 / 0.2e1 + Ifges(6,3) * t200 + t114 * t148 + ((t169 + Ifges(6,5) / 0.4e1) * t112 + t136) * t133 + ((t168 - Ifges(6,6) / 0.2e1) * t112 + t134) * t131 + t140 - t147;
t2 = qJD(3) * t23 + qJD(4) * t6;
t8 = [qJD(3) * t7 - qJD(4) * t1 - qJD(5) * t3, t2, qJD(4) * t12 + qJD(5) * t15 - t149, t12 * qJD(3) + t5 * qJD(5) + t154 + (t84 * mrSges(5,2) + (t117 * t197 + t119 * t196 + t120 * t198 + t207 * t195 - Ifges(5,5) + t191) * t112 + (-pkin(7) * t116 - Ifges(5,6) - t155) * t114 + (-m(6) * pkin(4) + t209) * t85 + (m(6) * pkin(7) + mrSges(6,3)) * t152) * qJD(4), t15 * qJD(3) + t5 * qJD(4) + (-t35 * mrSges(6,1) - t34 * mrSges(6,2) + t114 * t155) * qJD(5) - t151; t2, t22 * qJD(4), t172, (t209 * t114 + t111 + t201 - t210) * qJD(4) + t53 * qJD(5) + t150, -qJD(5) * t114 * t156 + t53 * qJD(4); -qJD(4) * t11 - qJD(5) * t14 + t149, -t172, 0, -t177, -qJD(5) * t116 - t173; qJD(3) * t11 + qJD(5) * t4 - t154, -qJD(5) * t52 - t150, t177, t51 * qJD(5), (-pkin(7) * t156 + t123 - t187) * qJD(5) + t137; qJD(3) * t14 - qJD(4) * t4 + t151, qJD(4) * t52, t173, -t137, 0;];
Cq = t8;
