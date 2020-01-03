% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:24
% EndTime: 2020-01-03 11:22:26
% DurationCPUTime: 1.05s
% Computational Cost: add. (3435->199), mult. (8207->311), div. (0->0), fcn. (8565->8), ass. (0->123)
t114 = sin(pkin(9));
t117 = cos(pkin(9));
t119 = cos(pkin(7));
t140 = t119 * t117;
t116 = sin(pkin(7));
t118 = cos(pkin(8));
t143 = t118 * t116;
t93 = t114 * t143 + t140;
t168 = t93 * mrSges(6,2);
t120 = sin(qJ(5));
t115 = sin(pkin(8));
t121 = cos(qJ(5));
t145 = t115 * t121;
t141 = t119 * t114;
t94 = t117 * t143 - t141;
t75 = t116 * t145 - t94 * t120;
t171 = t75 * mrSges(6,3);
t50 = -t168 + t171;
t187 = -t50 / 0.2e1;
t186 = t75 / 0.2e1;
t148 = t115 * t116;
t101 = -t119 * pkin(2) - t116 * qJ(3) - pkin(1);
t142 = t118 * t119;
t85 = qJ(2) * t142 + t115 * t101;
t77 = -t119 * qJ(4) + t85;
t90 = (pkin(3) * t115 - qJ(4) * t118 + qJ(2)) * t116;
t43 = -t114 * t77 + t117 * t90;
t37 = -pkin(4) * t148 - t43;
t185 = m(6) * t37;
t146 = t115 * t120;
t76 = t116 * t146 + t94 * t121;
t167 = -mrSges(5,1) * t148 - t75 * mrSges(6,1) + t76 * mrSges(6,2) + t94 * mrSges(5,3);
t184 = t114 ^ 2;
t183 = t115 ^ 2;
t112 = t116 ^ 2;
t182 = t118 ^ 2;
t181 = m(4) / 0.4e1;
t180 = m(5) / 0.2e1;
t179 = m(6) / 0.2e1;
t178 = mrSges(6,1) / 0.2e1;
t177 = -mrSges(6,2) / 0.2e1;
t176 = t93 / 0.2e1;
t97 = t117 * t145 - t120 * t118;
t175 = -t97 / 0.2e1;
t174 = t116 / 0.2e1;
t173 = -t120 / 0.2e1;
t172 = -t121 / 0.2e1;
t170 = t76 * mrSges(6,3);
t169 = t93 * mrSges(6,1);
t166 = Ifges(6,5) * t75 - Ifges(6,6) * t76;
t44 = t114 * t90 + t117 * t77;
t147 = t115 * t119;
t84 = -qJ(2) * t147 + t118 * t101;
t78 = t119 * pkin(3) - t84;
t36 = t93 * pkin(4) - t94 * pkin(6) + t78;
t38 = pkin(6) * t148 + t44;
t22 = -t120 * t38 + t121 * t36;
t23 = t120 * t36 + t121 * t38;
t39 = t76 * mrSges(6,1) + t75 * mrSges(6,2);
t51 = t169 - t170;
t1 = t166 * t176 - t23 * t51 + t22 * t50 + t37 * t39 + (-t22 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,5) * t176) * t75 + (-t23 * mrSges(6,3) - Ifges(6,6) * t93 / 0.2e1 - Ifges(6,4) * t76 + (Ifges(6,1) - Ifges(6,2)) * t75) * t76;
t165 = t1 * qJD(1);
t144 = t117 * t116;
t92 = t118 * t141 - t144;
t164 = t114 * t92;
t163 = t114 * t93;
t162 = t117 * t92;
t161 = t117 * t94;
t160 = t120 * mrSges(6,1);
t159 = t120 * t51;
t158 = t121 * mrSges(6,2);
t157 = t121 * t50;
t110 = t112 * qJ(2);
t113 = t119 ^ 2;
t151 = t119 * mrSges(4,1) + t93 * mrSges(5,1) + t94 * mrSges(5,2) + mrSges(4,3) * t143;
t95 = t114 * t116 + t118 * t140;
t73 = t119 * t145 - t120 * t95;
t74 = t119 * t146 + t121 * t95;
t79 = -mrSges(5,2) * t148 - t93 * mrSges(5,3);
t99 = t119 * mrSges(4,2) - mrSges(4,3) * t148;
t2 = t74 * t50 + t73 * t51 + t95 * t79 + t167 * t92 + (mrSges(4,1) * t115 + mrSges(4,2) * t118) * t112 + (t113 + t112) * mrSges(3,3) + (t151 * t115 + t118 * t99) * t119 + m(6) * (t22 * t73 + t23 * t74 + t37 * t92) + m(5) * (t78 * t147 - t43 * t92 + t44 * t95) + m(4) * (t110 + (-t115 * t84 + t118 * t85) * t119) + m(3) * (t113 * qJ(2) + t110);
t156 = t2 * qJD(1);
t149 = t114 * t115;
t129 = t117 * t146 + t121 * t118;
t88 = t129 * t116;
t89 = t97 * t116;
t3 = t88 * t51 + m(6) * (t22 * t88 - t23 * t89) - t89 * t50 + (t151 * t118 + (-t167 * t114 - t117 * t79 - t99) * t115 - t149 * t185 + m(5) * (-t115 * t117 * t44 + t118 * t78 + t43 * t149) + m(4) * (-t115 * t85 - t118 * t84)) * t116;
t155 = t3 * qJD(1);
t126 = (t129 * t186 + t76 * t175) * mrSges(6,3) + t129 * t187 + t51 * t175 + t39 * t149 / 0.2e1;
t131 = t74 * t177 + t73 * t178;
t4 = t126 - t131;
t154 = t4 * qJD(1);
t6 = (m(5) * t43 - t167 - t185) * t94 + (-t159 + t157 + t79 - m(6) * (t120 * t22 - t121 * t23) + m(5) * t44) * t93;
t153 = t6 * qJD(1);
t124 = (t51 * t172 + t50 * t173 + (t120 * t186 + t76 * t172) * mrSges(6,3)) * t114 - t117 * t39 / 0.2e1;
t130 = -t89 * t177 + t88 * t178;
t8 = t124 - t130;
t152 = t8 * qJD(1);
t10 = (t168 / 0.2e1 + t171 / 0.2e1 + t187) * t121 + (t169 / 0.2e1 + t51 / 0.2e1 + t170 / 0.2e1) * t120;
t150 = t10 * qJD(1);
t136 = m(4) * t174;
t122 = t136 + (t114 * t95 - t162) * t180 + (-t162 + (-t120 * t73 + t121 * t74) * t114) * t179;
t132 = -t129 * t88 - t97 * t89;
t133 = -t117 ^ 2 - t184;
t13 = -m(6) * t132 / 0.2e1 + 0.2e1 * ((m(5) / 0.4e1 + t181) * t182 + (m(6) * t184 / 0.4e1 - m(5) * t133 / 0.4e1 + t181) * t183) * t116 + t122;
t139 = t13 * qJD(1);
t135 = t115 * t180;
t123 = (t94 * t149 + (-t120 * t129 - t121 * t97) * t93) * t179 + (t114 * t94 - t117 * t93) * t135;
t128 = (t120 * t74 + t121 * t73) * t179 + t119 * t135;
t18 = t123 - t128;
t138 = t18 * qJD(1);
t125 = (-t161 + (-t120 ^ 2 - t121 ^ 2) * t163) * t179 + (-t161 - t163) * t180;
t134 = m(5) * t174;
t127 = (-t120 * t89 + t121 * t88) * t179 + t118 * t134;
t20 = t125 - t127;
t137 = t20 * qJD(1);
t19 = t125 + t127;
t17 = t123 + t128;
t14 = (-t184 * t183 * t116 + t132) * t179 + (t133 * t183 - t182) * t134 + (-t182 - t183) * t136 + t122;
t11 = -t159 / 0.2e1 + t171 * t172 + t157 / 0.2e1 + t170 * t173 + (t158 / 0.2e1 + t160 / 0.2e1) * t93;
t9 = t124 + t130;
t5 = t126 + t131;
t7 = [t2 * qJD(2) + t3 * qJD(3) - t6 * qJD(4) + t1 * qJD(5), t156 + 0.2e1 * ((-t129 * t73 + t97 * t74) * t179 + (t164 * t179 + (t117 * t95 - t142 + t164) * t180) * t115) * qJD(2) + t14 * qJD(3) + t17 * qJD(4) + t5 * qJD(5), t155 + t14 * qJD(2) + t19 * qJD(4) + t9 * qJD(5) + m(6) * (t115 * t144 - t120 * t88 - t121 * t89) * qJD(3) * t114, t17 * qJD(2) + t19 * qJD(3) + t11 * qJD(5) - t153, t165 + t5 * qJD(2) + t9 * qJD(3) + t11 * qJD(4) + (-t23 * mrSges(6,1) - t22 * mrSges(6,2) + t166) * qJD(5); -t13 * qJD(3) + t18 * qJD(4) + t4 * qJD(5) - t156, 0, -t139, t138, t154 + (-t97 * mrSges(6,1) + mrSges(6,2) * t129) * qJD(5); t13 * qJD(2) + t20 * qJD(4) + t8 * qJD(5) - t155, t139, 0, t137, t152 + (-mrSges(6,1) * t121 + mrSges(6,2) * t120) * qJD(5) * t114; -t18 * qJD(2) - t20 * qJD(3) - t10 * qJD(5) + t153, -t138, -t137, 0, -t150 + (-t158 - t160) * qJD(5); -t4 * qJD(2) - t8 * qJD(3) + t10 * qJD(4) - t165, -t154, -t152, t150, 0;];
Cq = t7;
