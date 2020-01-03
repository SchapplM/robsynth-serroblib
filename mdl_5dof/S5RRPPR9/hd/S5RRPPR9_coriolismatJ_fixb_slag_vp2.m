% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:43
% EndTime: 2019-12-31 19:40:46
% DurationCPUTime: 1.19s
% Computational Cost: add. (2253->243), mult. (4210->333), div. (0->0), fcn. (3175->4), ass. (0->131)
t208 = pkin(6) - qJ(4);
t128 = cos(qJ(5));
t113 = t128 * mrSges(6,1);
t126 = sin(qJ(5));
t184 = t126 * mrSges(6,2);
t207 = -t113 + t184;
t121 = t126 ^ 2;
t122 = t128 ^ 2;
t160 = t121 + t122;
t115 = Ifges(6,6) * t126;
t190 = Ifges(6,5) * t128;
t206 = t190 / 0.2e1 - t115 / 0.2e1;
t205 = Ifges(4,5) - Ifges(3,4) - Ifges(5,4) + t206;
t129 = cos(qJ(2));
t111 = t129 * qJ(3);
t127 = sin(qJ(2));
t130 = -pkin(2) - pkin(3);
t72 = t127 * t130 + t111;
t204 = m(5) / 0.2e1;
t203 = m(6) / 0.2e1;
t202 = m(5) * t72;
t120 = -pkin(7) + t130;
t201 = -t120 / 0.2e1;
t200 = -t126 / 0.2e1;
t199 = t126 / 0.2e1;
t198 = t127 / 0.2e1;
t197 = -t128 / 0.2e1;
t196 = t128 / 0.2e1;
t195 = t129 / 0.2e1;
t194 = mrSges(5,3) - mrSges(4,2);
t193 = mrSges(6,3) * t129;
t192 = Ifges(6,4) * t126;
t191 = Ifges(6,4) * t128;
t189 = Ifges(6,6) * t127;
t85 = t208 * t129;
t170 = t129 * t85;
t166 = t126 * t129;
t179 = t127 * mrSges(6,2);
t74 = mrSges(6,3) * t166 - t179;
t174 = t128 * t74;
t165 = t128 * t129;
t180 = t127 * mrSges(6,1);
t76 = mrSges(6,3) * t165 + t180;
t181 = t126 * t76;
t79 = -t129 * pkin(2) - t127 * qJ(3) - pkin(1);
t68 = t129 * pkin(3) - t79;
t35 = pkin(4) * t127 + pkin(7) * t129 + t68;
t81 = t208 * t127;
t25 = -t126 * t81 + t128 * t35;
t26 = t126 * t35 + t128 * t81;
t176 = t128 * mrSges(6,2);
t185 = t126 * mrSges(6,1);
t146 = t176 + t185;
t65 = t146 * t129;
t7 = (t129 * mrSges(5,3) + t65) * t129 + (t127 * mrSges(5,3) - t174 + t181) * t127 + m(6) * (-t170 + (t126 * t25 - t128 * t26) * t127) + m(5) * (-t127 * t81 - t170);
t188 = qJD(1) * t7;
t154 = m(4) * t79 - t129 * mrSges(4,1) - t127 * mrSges(4,3);
t173 = t128 * t76;
t182 = t126 * t74;
t86 = mrSges(5,1) * t127 - mrSges(5,2) * t129;
t8 = (m(6) * (t126 * t26 + t128 * t25) + t182 + t173 + m(5) * t68 + t86 - t154) * t127;
t187 = qJD(1) * t8;
t161 = t129 * mrSges(5,1) + t127 * mrSges(5,2);
t36 = pkin(4) * t129 + t120 * t127 + t111;
t27 = -t126 * t85 + t128 * t36;
t28 = t126 * t36 + t128 * t85;
t143 = -Ifges(6,2) * t126 + t191;
t43 = Ifges(6,6) * t129 + t143 * t127;
t44 = -t143 * t129 + t189;
t145 = Ifges(6,1) * t128 - t192;
t45 = Ifges(6,5) * t129 + t145 * t127;
t106 = Ifges(6,4) * t166;
t177 = t127 * Ifges(6,5);
t46 = -Ifges(6,1) * t165 + t106 + t177;
t64 = t146 * t127;
t171 = t129 * mrSges(6,2);
t178 = t127 * mrSges(6,3);
t73 = -t126 * t178 - t171;
t172 = t129 * mrSges(6,1);
t75 = -t128 * t178 + t172;
t1 = t72 * t86 + t28 * t74 + t25 * t75 + t27 * t76 + t81 * t65 + t85 * t64 + t26 * t73 + t154 * (pkin(2) * t127 - t111) + m(6) * (t25 * t27 + t26 * t28 - t81 * t85) + (-pkin(1) * mrSges(3,2) - t79 * mrSges(4,3) - t205 * t129 + t45 * t197 + t43 * t199) * t129 + (t79 * mrSges(4,1) - pkin(1) * mrSges(3,1) + t46 * t196 + t44 * t200 + (Ifges(4,1) - Ifges(4,3) + Ifges(3,1) - Ifges(3,2) + Ifges(6,3) - Ifges(5,1) + Ifges(5,2)) * t129 + t205 * t127) * t127 + (t161 + t202) * t68;
t186 = t1 * qJD(1);
t183 = t126 * t27;
t175 = t128 * t28;
t162 = Ifges(6,5) * t166 + Ifges(6,6) * t165;
t63 = t207 * t129;
t66 = Ifges(6,2) * t165 + t106;
t144 = Ifges(6,1) * t126 + t191;
t67 = t144 * t129;
t4 = t85 * t63 + t162 * t198 + t25 * t74 - t26 * t76 + ((-t67 / 0.2e1 + t44 / 0.2e1 + t26 * mrSges(6,3)) * t128 + (t46 / 0.2e1 + t66 / 0.2e1 - t25 * mrSges(6,3)) * t126) * t129;
t169 = t4 * qJD(1);
t123 = qJ(3) + pkin(4);
t153 = -t122 / 0.2e1 - t121 / 0.2e1;
t148 = mrSges(6,3) * t153;
t151 = t160 * t127;
t131 = -t127 * t148 - t72 * t204 + (-t120 * t151 - t123 * t129) * t203 + t207 * t195;
t132 = t202 / 0.2e1 + (t126 * t28 + t128 * t27) * t203 + t73 * t199 + t75 * t196;
t5 = t131 - t132 - t161;
t168 = t5 * qJD(1);
t149 = t179 / 0.2e1 - t74 / 0.2e1;
t150 = t180 / 0.2e1 + t76 / 0.2e1;
t10 = -t149 * t126 + t150 * t128 + t153 * t193;
t167 = t10 * qJD(1);
t13 = t150 * t126 + t149 * t128;
t164 = t13 * qJD(1);
t29 = (-0.2e1 * t153 * m(6) + m(5)) * t127;
t163 = t29 * qJD(1);
t159 = -t144 / 0.4e1 - t143 / 0.4e1;
t87 = -Ifges(6,2) * t128 - t192;
t158 = t145 / 0.4e1 + t87 / 0.4e1;
t152 = -m(4) * qJ(3) - mrSges(4,3);
t142 = t175 - t183;
t20 = -t123 * t146 + (t144 / 0.2e1 + t143 / 0.2e1) * t128 + (t145 / 0.2e1 + t87 / 0.2e1) * t126;
t133 = t123 * t63 / 0.2e1 + t127 * t115 / 0.4e1 - t85 * t146 / 0.2e1;
t136 = t74 * t201 - t67 / 0.4e1 + t44 / 0.4e1;
t137 = -t46 / 0.4e1 - t66 / 0.4e1 + t76 * t201;
t138 = -t27 * mrSges(6,1) / 0.2e1 + t28 * mrSges(6,2) / 0.2e1;
t139 = t120 * t148;
t3 = (-Ifges(6,3) / 0.2e1 - t139) * t129 + (t189 / 0.2e1 + t159 * t129 + t136) * t126 + (-0.3e1 / 0.4e1 * t177 + t158 * t129 + t137) * t128 + t133 + t138;
t141 = -t3 * qJD(1) - t20 * qJD(2);
t12 = (-t171 / 0.2e1 - t73 / 0.2e1) * t128 + (-t172 / 0.2e1 + t75 / 0.2e1) * t126 + 0.2e1 * (t85 / 0.4e1 + t183 / 0.4e1 - t175 / 0.4e1) * m(6);
t32 = -m(5) * qJ(3) - m(6) * t123 - mrSges(5,1) + t152 + t207;
t140 = qJD(1) * t12 - qJD(2) * t32;
t135 = t176 / 0.2e1 + t185 / 0.2e1;
t30 = m(6) * t160 * t198 - t151 * t203;
t14 = t174 / 0.2e1 - t181 / 0.2e1 + t135 * t127;
t11 = -t182 / 0.2e1 - t173 / 0.2e1 + (-t184 / 0.2e1 + t113 / 0.2e1) * t127 + t160 * t193 / 0.2e1;
t9 = t142 * t203 + t73 * t196 + t75 * t200 + 0.2e1 * (m(6) / 0.4e1 + t204) * t85 + (m(4) * pkin(6) - t135 - t194) * t129;
t6 = t131 + t132;
t2 = Ifges(6,3) * t195 + (-t177 / 0.4e1 + t137) * t128 + t136 * t126 + (t159 * t126 + t158 * t128 - t139) * t129 + t133 - t138 + t206 * t127;
t15 = [qJD(2) * t1 + qJD(3) * t8 + qJD(4) * t7 + qJD(5) * t4, t9 * qJD(3) + t6 * qJD(4) + t2 * qJD(5) + t186 + (-t81 * mrSges(5,1) + t85 * mrSges(5,2) + t123 * t64 + t81 * t207 + (-t43 / 0.2e1 + t120 * t73 - t28 * mrSges(6,3)) * t128 + (-t45 / 0.2e1 - t120 * t75 + t27 * mrSges(6,3)) * t126 + 0.2e1 * (t120 * t142 - t123 * t81) * t203 + 0.2e1 * (-qJ(3) * t81 + t130 * t85) * t204 + (Ifges(6,5) * t200 + Ifges(6,6) * t197 + Ifges(3,5) + Ifges(5,6) + Ifges(4,4) - pkin(2) * mrSges(4,2) - t130 * mrSges(5,3) + (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1)) * pkin(6)) * t129 + (Ifges(4,6) - Ifges(3,6) - Ifges(5,5) - t144 * t196 + t87 * t200 + t194 * qJ(3) + (mrSges(3,2) + t152) * pkin(6)) * t127) * qJD(2), qJD(2) * t9 + qJD(4) * t30 + qJD(5) * t11 + t187, qJD(2) * t6 + qJD(3) * t30 + qJD(5) * t14 + t188, t169 + t2 * qJD(2) + t11 * qJD(3) + t14 * qJD(4) + (-mrSges(6,1) * t26 - mrSges(6,2) * t25 + t162) * qJD(5); qJD(3) * t12 + qJD(4) * t5 + qJD(5) * t3 - t186, -qJD(3) * t32 + qJD(5) * t20, t140, t168, (t120 * t207 + t115 - t190) * qJD(5) - t141; -qJD(2) * t12 - qJD(4) * t29 - qJD(5) * t10 - t187, -t140, 0, -t163, qJD(5) * t207 - t167; -qJD(2) * t5 + qJD(3) * t29 - qJD(5) * t13 - t188, -t168, t163, 0, -qJD(5) * t146 - t164; -qJD(2) * t3 + qJD(3) * t10 + qJD(4) * t13 - t169, t141, t167, t164, 0;];
Cq = t15;
