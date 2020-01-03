% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:45
% EndTime: 2019-12-31 18:23:47
% DurationCPUTime: 1.19s
% Computational Cost: add. (2076->214), mult. (4102->305), div. (0->0), fcn. (3197->6), ass. (0->128)
t104 = sin(qJ(5));
t101 = t104 ^ 2;
t106 = cos(qJ(5));
t102 = t106 ^ 2;
t144 = t101 + t102;
t107 = cos(qJ(3));
t147 = t106 * t107;
t105 = sin(qJ(3));
t160 = t105 * mrSges(6,2);
t75 = -mrSges(6,3) * t147 - t160;
t153 = t106 * t75;
t149 = t104 * t107;
t161 = t105 * mrSges(6,1);
t73 = mrSges(6,3) * t149 + t161;
t163 = t104 * t73;
t116 = -t163 / 0.2e1 + t153 / 0.2e1;
t172 = mrSges(6,3) * t107;
t196 = t116 + t144 * t172 / 0.2e1;
t178 = m(6) * t105;
t174 = mrSges(6,2) * t106;
t176 = mrSges(6,1) * t104;
t79 = t174 + t176;
t195 = mrSges(5,3) + t79;
t168 = Ifges(6,6) * t106;
t169 = Ifges(6,5) * t104;
t117 = t169 / 0.2e1 + t168 / 0.2e1;
t194 = Ifges(4,4) + Ifges(5,6) - t117;
t180 = t106 / 0.2e1;
t183 = t104 / 0.2e1;
t170 = Ifges(6,4) * t106;
t81 = -Ifges(6,2) * t104 + t170;
t171 = Ifges(6,4) * t104;
t82 = Ifges(6,1) * t106 - t171;
t114 = t81 * t180 + t82 * t183;
t96 = sin(pkin(8)) * pkin(1) + pkin(6);
t179 = pkin(4) + t96;
t68 = t179 * t107;
t146 = t107 * qJ(4);
t80 = -t105 * pkin(3) + t146;
t71 = pkin(7) * t105 - t80;
t25 = -t104 * t71 + t68 * t106;
t156 = t106 * t25;
t26 = t68 * t104 + t106 * t71;
t165 = t104 * t26;
t192 = -t165 - t156;
t191 = m(6) / 0.2e1;
t190 = Ifges(6,3) / 0.2e1;
t189 = t68 / 0.2e1;
t188 = t73 / 0.2e1;
t187 = -t75 / 0.2e1;
t108 = -pkin(3) - pkin(7);
t186 = (0.1e1 - t144) * t107 * t178;
t185 = qJ(4) / 0.2e1;
t184 = -t104 / 0.2e1;
t182 = t105 / 0.2e1;
t181 = -t106 / 0.2e1;
t177 = pkin(3) * t107;
t175 = mrSges(6,1) * t107;
t173 = mrSges(6,2) * t107;
t134 = -cos(pkin(8)) * pkin(1) - pkin(2);
t151 = qJ(4) * t105;
t120 = t134 - t151;
t65 = t120 - t177;
t78 = t107 * mrSges(5,2) - t105 * mrSges(5,3);
t135 = -m(5) * t65 - t78;
t49 = t108 * t107 + t120;
t67 = t179 * t105;
t21 = -t104 * t49 + t106 * t67;
t22 = t104 * t67 + t106 * t49;
t9 = (m(6) * (t104 * t21 - t106 * t22) - t153 + t163 + t135) * t105;
t167 = qJD(1) * t9;
t166 = t104 * mrSges(6,2);
t127 = Ifges(6,1) * t104 + t170;
t112 = t127 * t107;
t159 = t105 * Ifges(6,5);
t53 = -t112 + t159;
t164 = t104 * t53;
t148 = t105 * t106;
t74 = mrSges(6,3) * t148 - t173;
t162 = t104 * t74;
t158 = t105 * Ifges(6,6);
t157 = t106 * mrSges(6,1);
t126 = Ifges(6,2) * t106 + t171;
t111 = t126 * t107;
t51 = -t111 + t158;
t155 = t106 * t51;
t72 = -mrSges(6,3) * t104 * t105 + t175;
t154 = t106 * t72;
t133 = t101 / 0.2e1 + t102 / 0.2e1;
t8 = (t133 * t172 + t79 * t182 + t116) * t107;
t152 = t8 * qJD(1);
t150 = qJD(5) * t79;
t129 = mrSges(6,3) * t133;
t11 = (t160 / 0.2e1 + t187) * t106 + (t161 / 0.2e1 + t188) * t104 - t107 * t129;
t145 = t11 * qJD(1);
t143 = qJD(3) * t105;
t142 = qJD(3) * t107;
t141 = Ifges(6,2) / 0.4e1 - Ifges(6,1) / 0.4e1;
t132 = m(6) * t144;
t128 = t157 - t166;
t50 = Ifges(6,6) * t107 + t126 * t105;
t52 = Ifges(6,5) * t107 + t127 * t105;
t66 = t105 * t128;
t1 = t21 * t72 + t25 * t73 + t22 * t74 + t26 * t75 - t68 * t66 + m(6) * (t21 * t25 + t22 * t26 - t67 * t68) + t135 * t80 + (t134 * mrSges(4,1) - t65 * mrSges(5,2) + t164 / 0.2e1 + t155 / 0.2e1 - t194 * t105) * t105 + (-t67 * t128 + t134 * mrSges(4,2) - t65 * mrSges(5,3) + t50 * t181 + t52 * t184 + (Ifges(4,1) - Ifges(4,2) + Ifges(6,3) + Ifges(5,2) - Ifges(5,3)) * t105 + t194 * t107) * t107;
t113 = t157 / 0.2e1 - t166 / 0.2e1;
t115 = -t162 / 0.2e1 - t154 / 0.2e1;
t123 = t104 * t22 + t106 * t21;
t6 = ((t192 + t68) * t191 + t113 * t107 + t115) * t107 + (t75 * t183 + t73 * t180 - t66 / 0.2e1 + (t123 - t67) * t191) * t105;
t125 = t1 * qJD(1) + t6 * qJD(2);
t95 = Ifges(6,6) * t149;
t4 = t95 * t182 + t21 * t75 - t22 * t73 + (-t68 * t79 - Ifges(6,5) * t148 / 0.2e1 + t53 * t181 + t51 * t183 + t114 * t107 + t123 * mrSges(6,3)) * t107;
t124 = t4 * qJD(1) - t8 * qJD(2);
t18 = qJ(4) * t128 + t126 * t183 + t127 * t181 - t114;
t119 = -t26 * mrSges(6,2) / 0.2e1 + t25 * mrSges(6,1) / 0.2e1;
t2 = (-t108 * t129 + t190) * t107 + (0.3e1 / 0.4e1 * t158 - t68 * mrSges(6,1) / 0.2e1 + t51 / 0.4e1 + t108 * t187 + (mrSges(6,2) * t185 + t82 / 0.4e1 - t141 * t106) * t107) * t106 + (0.3e1 / 0.4e1 * t159 + mrSges(6,2) * t189 + t53 / 0.4e1 + t108 * t188 + (mrSges(6,1) * t185 - t170 - t81 / 0.4e1 + t141 * t104) * t107) * t104 + t119;
t121 = t2 * qJD(1) - t18 * qJD(3);
t118 = t6 * qJD(1) + qJD(2) * t186;
t13 = (t175 / 0.2e1 - t72 / 0.2e1) * t106 + (-t173 / 0.2e1 - t74 / 0.2e1) * t104 + 0.2e1 * (t68 / 0.4e1 - t165 / 0.4e1 - t156 / 0.4e1) * m(6);
t34 = (-0.1e1 / 0.2e1 + t133) * t178;
t64 = (m(5) + m(6)) * qJ(4) + t195;
t109 = qJD(1) * t13 - qJD(2) * t34 + qJD(3) * t64;
t29 = t178 / 0.2e1 + (m(5) + t132 / 0.2e1) * t105;
t28 = t113 * t105 + t128 * t182;
t12 = (t174 / 0.2e1 + t176 / 0.2e1) * t105 + t196;
t10 = m(6) * t189 - t192 * t191 + (m(5) * t96 + mrSges(5,1) + t113) * t107 - t115;
t5 = qJD(3) * t6 - qJD(5) * t8;
t3 = -t79 * t146 / 0.2e1 + t105 * (-t168 - t169) / 0.4e1 + t128 * t189 - t164 / 0.4e1 - t82 * t147 / 0.2e1 - t155 / 0.4e1 + t81 * t149 / 0.2e1 + t104 * t112 / 0.4e1 + t106 * t111 / 0.4e1 + t107 * t190 + t117 * t105 + t119 + t196 * t108;
t7 = [qJD(3) * t1 + qJD(4) * t9 + qJD(5) * t4, t5, t10 * qJD(4) + t3 * qJD(5) + (-pkin(3) * mrSges(5,1) + Ifges(6,5) * t180 + Ifges(6,6) * t184 - Ifges(5,4) + Ifges(4,5)) * t142 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + t114) * t143 + t125 + (-t67 * t79 + t52 * t180 + t50 * t184 + (m(5) * (-t151 - t177) + t105 * mrSges(4,2) - t107 * mrSges(4,1) + t78) * t96 + (-m(6) * t67 - t66) * qJ(4) + (-m(6) * t192 + t154 + t162) * t108 + t192 * mrSges(6,3)) * qJD(3), qJD(3) * t10 + qJD(5) * t12 + t167, t3 * qJD(3) + t12 * qJD(4) + (-mrSges(6,1) * t22 - mrSges(6,2) * t21 - Ifges(6,5) * t147 + t95) * qJD(5) + t124; t5, qJD(3) * t186, m(5) * t80 * qJD(3) + t29 * qJD(4) + t28 * qJD(5) + (m(6) * qJ(4) - mrSges(4,2) + t195) * t142 + (-t144 * mrSges(6,3) + t108 * t132 - mrSges(4,1) + mrSges(5,2)) * t143 + t118, t29 * qJD(3), t28 * qJD(3) + t107 * t150 - t152; qJD(4) * t13 - qJD(5) * t2 - t125, -t34 * qJD(4) - t118, qJD(4) * t64 + qJD(5) * t18, t109, ((-mrSges(6,2) * t108 - Ifges(6,6)) * t106 + (-mrSges(6,1) * t108 - Ifges(6,5)) * t104) * qJD(5) - t121; -qJD(3) * t13 - qJD(5) * t11 - t167, t34 * qJD(3), -t109, 0, -t145 - t150; qJD(3) * t2 + qJD(4) * t11 - t124, t152, t121, t145, 0;];
Cq = t7;
