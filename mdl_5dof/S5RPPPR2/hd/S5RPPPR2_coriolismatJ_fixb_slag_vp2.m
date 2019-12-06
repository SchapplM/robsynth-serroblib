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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:31:02
% EndTime: 2019-12-05 17:31:06
% DurationCPUTime: 1.16s
% Computational Cost: add. (3435->193), mult. (8207->303), div. (0->0), fcn. (8565->8), ass. (0->121)
t109 = sin(pkin(9));
t112 = cos(pkin(9));
t114 = cos(pkin(7));
t140 = t114 * t112;
t111 = sin(pkin(7));
t113 = cos(pkin(8));
t144 = t111 * t113;
t88 = t109 * t144 + t140;
t167 = t88 * mrSges(6,2);
t115 = sin(qJ(5));
t110 = sin(pkin(8));
t116 = cos(qJ(5));
t145 = t110 * t116;
t141 = t114 * t109;
t143 = t112 * t111;
t89 = t113 * t143 - t141;
t75 = t111 * t145 - t89 * t115;
t170 = t75 * mrSges(6,3);
t51 = -t167 + t170;
t188 = -t51 / 0.2e1;
t187 = t75 / 0.2e1;
t148 = t110 * t111;
t146 = t110 * t115;
t76 = t111 * t146 + t89 * t116;
t186 = -mrSges(5,1) * t148 - t75 * mrSges(6,1) + t76 * mrSges(6,2) + t89 * mrSges(5,3);
t185 = -t114 * mrSges(4,1) - t88 * mrSges(5,1) - t89 * mrSges(5,2) - mrSges(4,3) * t144;
t184 = t109 ^ 2;
t183 = t110 ^ 2;
t107 = t111 ^ 2;
t182 = t113 ^ 2;
t181 = m(5) / 0.2e1;
t180 = m(6) / 0.2e1;
t179 = mrSges(6,1) / 0.2e1;
t178 = -mrSges(6,2) / 0.2e1;
t177 = t88 / 0.2e1;
t92 = t112 * t145 - t115 * t113;
t176 = -t92 / 0.2e1;
t175 = -t115 / 0.2e1;
t174 = -t116 / 0.2e1;
t173 = m(4) * t111;
t172 = m(5) * t111;
t171 = m(6) * t111;
t169 = t76 * mrSges(6,3);
t168 = t88 * mrSges(6,1);
t166 = Ifges(6,5) * t75 - Ifges(6,6) * t76;
t142 = t113 * t114;
t96 = -t114 * pkin(2) - t111 * qJ(3) - pkin(1);
t151 = qJ(2) * t142 + t110 * t96;
t77 = -t114 * qJ(4) + t151;
t85 = (pkin(3) * t110 - qJ(4) * t113 + qJ(2)) * t111;
t45 = t109 * t85 + t112 * t77;
t147 = t110 * t114;
t131 = -qJ(2) * t147 + t113 * t96;
t130 = t114 * pkin(3) - t131;
t37 = t88 * pkin(4) - t89 * pkin(6) + t130;
t39 = pkin(6) * t148 + t45;
t22 = -t115 * t39 + t116 * t37;
t23 = t115 * t37 + t116 * t39;
t44 = -t109 * t77 + t112 * t85;
t38 = -pkin(4) * t148 - t44;
t40 = t76 * mrSges(6,1) + t75 * mrSges(6,2);
t52 = t168 - t169;
t1 = t38 * t40 + t22 * t51 - t23 * t52 + t166 * t177 + (-t22 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,5) * t177) * t75 + (-t23 * mrSges(6,3) - Ifges(6,6) * t88 / 0.2e1 - Ifges(6,4) * t76 + (Ifges(6,1) - Ifges(6,2)) * t75) * t76;
t165 = t1 * qJD(1);
t87 = t113 * t141 - t143;
t164 = t109 * t87;
t163 = t109 * t88;
t162 = t112 * t87;
t161 = t112 * t89;
t160 = t115 * mrSges(6,1);
t159 = t115 * t52;
t158 = t116 * mrSges(6,2);
t157 = t116 * t51;
t105 = t107 * qJ(2);
t108 = t114 ^ 2;
t90 = t109 * t111 + t113 * t140;
t73 = t114 * t145 - t115 * t90;
t74 = t114 * t146 + t116 * t90;
t78 = -mrSges(5,2) * t148 - t88 * mrSges(5,3);
t94 = t114 * mrSges(4,2) - mrSges(4,3) * t148;
t2 = (mrSges(4,1) * t110 + mrSges(4,2) * t113) * t107 + m(6) * (t22 * t73 + t23 * t74) + t74 * t51 + t73 * t52 + m(4) * (t105 + (-t131 * t110 + t151 * t113) * t114) + t94 * t142 + m(3) * (t108 * qJ(2) + t105) + (m(5) * t45 + t78) * t90 + (-m(5) * t44 + m(6) * t38 + t186) * t87 + (m(5) * t130 - t185) * t147 + (t108 + t107) * mrSges(3,3);
t156 = t2 * qJD(1);
t127 = t112 * t146 + t116 * t113;
t125 = t127 * t111;
t126 = t92 * t111;
t134 = t110 * t143;
t149 = t109 * t110;
t3 = -(t22 * t127 - t38 * t149 - t23 * t92) * t171 + t51 * t126 - t52 * t125 + t78 * t134 - (t130 * t113 + (t44 * t109 - t45 * t112) * t110) * t172 - (-t151 * t110 - t131 * t113) * t173 + t185 * t144 + (t186 * t109 + t94) * t148;
t155 = t3 * qJD(1);
t122 = (t127 * t187 + t76 * t176) * mrSges(6,3) + t127 * t188 + t52 * t176 + t40 * t149 / 0.2e1;
t128 = t74 * t178 + t73 * t179;
t4 = t122 - t128;
t154 = t4 * qJD(1);
t6 = t186 * t89 + (-t157 - t78 + t159) * t88 + m(6) * (t38 * t89 + (t115 * t22 - t116 * t23) * t88) + m(5) * (-t44 * t89 - t45 * t88);
t153 = t6 * qJD(1);
t120 = (t51 * t175 + t52 * t174 + (t115 * t187 + t76 * t174) * mrSges(6,3)) * t109 - t112 * t40 / 0.2e1;
t8 = (mrSges(6,2) * t176 - t127 * mrSges(6,1) / 0.2e1) * t111 + t120;
t152 = t8 * qJD(1);
t10 = (t167 / 0.2e1 + t188 + t170 / 0.2e1) * t116 + (t168 / 0.2e1 + t169 / 0.2e1 + t52 / 0.2e1) * t115;
t150 = t10 * qJD(1);
t117 = -(-t182 - t183) * t173 / 0.2e1 - (-t182 + (-t112 ^ 2 - t184) * t183) * t172 / 0.2e1 - (-t127 ^ 2 - t184 * t183 - t92 ^ 2) * t171 / 0.2e1;
t118 = t173 / 0.2e1 + (t109 * t90 - t162) * t181 + (-t162 + (-t115 * t73 + t116 * t74) * t109) * t180;
t13 = t117 + t118;
t139 = t13 * qJD(1);
t136 = t110 * t181;
t119 = (t89 * t149 + (-t115 * t127 - t116 * t92) * t88) * t180 + (t109 * t89 - t112 * t88) * t136;
t124 = (t115 * t74 + t116 * t73) * t180 + t114 * t136;
t18 = t119 - t124;
t138 = t18 * qJD(1);
t132 = t115 ^ 2 + t116 ^ 2;
t121 = (-t132 * t163 - t161) * t180 + (-t161 - t163) * t181;
t123 = -(m(6) * t132 + m(5)) * t144 / 0.2e1;
t20 = t121 + t123;
t137 = t20 * qJD(1);
t19 = t121 - t123;
t17 = t119 + t124;
t14 = -t117 + t118;
t11 = t157 / 0.2e1 + t169 * t175 - t159 / 0.2e1 + t170 * t174 + (t158 / 0.2e1 + t160 / 0.2e1) * t88;
t9 = t125 * t179 - t126 * t178 + t120;
t5 = t122 + t128;
t7 = [t2 * qJD(2) - t3 * qJD(3) + t6 * qJD(4) + t1 * qJD(5), t156 + 0.2e1 * ((-t127 * t73 + t92 * t74) * t180 + (t164 * t180 + (t112 * t90 - t142 + t164) * t181) * t110) * qJD(2) + t14 * qJD(3) + t17 * qJD(4) + t5 * qJD(5), -t155 + t14 * qJD(2) + m(6) * (-t132 + 0.1e1) * t109 * qJD(3) * t134 + t19 * qJD(4) + t9 * qJD(5), t17 * qJD(2) + t19 * qJD(3) + t11 * qJD(5) + t153, t165 + t5 * qJD(2) + t9 * qJD(3) + t11 * qJD(4) + (-t23 * mrSges(6,1) - t22 * mrSges(6,2) + t166) * qJD(5); -t13 * qJD(3) + t18 * qJD(4) + t4 * qJD(5) - t156, 0, -t139, t138, t154 + (-t92 * mrSges(6,1) + mrSges(6,2) * t127) * qJD(5); t13 * qJD(2) + t20 * qJD(4) + t8 * qJD(5) + t155, t139, 0, t137, t152 + (-mrSges(6,1) * t116 + mrSges(6,2) * t115) * qJD(5) * t109; -t18 * qJD(2) - t20 * qJD(3) - t10 * qJD(5) - t153, -t138, -t137, 0, -t150 + (-t158 - t160) * qJD(5); -t4 * qJD(2) - t8 * qJD(3) + t10 * qJD(4) - t165, -t154, -t152, t150, 0;];
Cq = t7;
