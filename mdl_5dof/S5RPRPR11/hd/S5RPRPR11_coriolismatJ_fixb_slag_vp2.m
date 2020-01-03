% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR11
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:13
% EndTime: 2019-12-31 18:27:15
% DurationCPUTime: 1.24s
% Computational Cost: add. (3704->159), mult. (7157->208), div. (0->0), fcn. (7678->6), ass. (0->98)
t101 = cos(qJ(5));
t100 = sin(qJ(3));
t152 = cos(qJ(3));
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t87 = t100 * t97 - t152 * t98;
t89 = t100 * t98 + t152 * t97;
t99 = sin(qJ(5));
t114 = t101 * t89 + t99 * t87;
t137 = pkin(6) + qJ(2);
t120 = t137 * t97;
t93 = t137 * t98;
t74 = -t100 * t120 + t152 * t93;
t104 = t87 * pkin(7) + t74;
t73 = t100 * t93 + t152 * t120;
t43 = -t89 * pkin(7) + t73;
t165 = t101 * t104 + t99 * t43;
t28 = -t101 * t43 + t104 * t99;
t67 = t101 * t87 - t99 * t89;
t188 = -t165 * mrSges(6,1) + t28 * mrSges(6,2) + Ifges(6,5) * t67 - Ifges(6,6) * t114;
t189 = t188 * qJD(5);
t102 = -pkin(3) - pkin(4);
t132 = t87 * qJ(4);
t41 = t102 * t89 - t132;
t187 = m(6) * t41;
t178 = t67 * mrSges(6,2);
t124 = -t178 / 0.2e1;
t119 = -t114 * mrSges(6,1) - t178;
t121 = -t98 * pkin(2) - pkin(1);
t109 = t89 * qJ(4) - t121;
t36 = t102 * t87 + t109;
t185 = t36 * t119;
t64 = Ifges(6,4) * t67;
t184 = Ifges(6,2) * t114 - t64;
t182 = t101 * t165 + t99 * t28;
t181 = -t67 / 0.2e1;
t180 = t67 / 0.2e1;
t177 = Ifges(4,4) - Ifges(5,5);
t60 = t87 * pkin(3) - t109;
t175 = m(5) * t60 + t87 * mrSges(5,1) - t89 * mrSges(5,3);
t115 = t99 * mrSges(6,1) + t101 * mrSges(6,2);
t173 = t115 * qJD(5);
t151 = Ifges(6,4) * t114;
t172 = Ifges(6,1) * t67 - t151;
t170 = -t114 / 0.2e1;
t169 = t114 / 0.2e1;
t164 = t87 ^ 2;
t163 = 0.2e1 * t89;
t162 = -m(5) / 0.2e1;
t161 = m(5) / 0.2e1;
t160 = -m(6) / 0.4e1;
t159 = m(6) / 0.2e1;
t118 = t89 * mrSges(4,1) - t87 * mrSges(4,2);
t32 = -t67 * mrSges(6,1) + mrSges(6,2) * t114;
t33 = Ifges(6,2) * t67 + t151;
t34 = Ifges(6,1) * t114 + t64;
t70 = pkin(3) * t89 + t132;
t80 = t87 * mrSges(5,3);
t1 = t41 * t32 + t185 + t34 * t181 + t184 * t180 + t121 * t118 + t60 * t80 + t36 * t187 + (t60 * mrSges(5,1) + (-Ifges(4,1) - Ifges(5,1) + Ifges(4,2) + Ifges(5,3)) * t87 - t177 * t89) * t89 + t175 * t70 + t177 * t164 + (t33 - t172) * t169;
t135 = t1 * qJD(1);
t5 = (-t114 ^ 2 - t67 ^ 2) * mrSges(6,3) + m(6) * (-t114 * t28 - t165 * t67) + (t89 ^ 2 + t164) * (mrSges(4,3) + mrSges(5,2)) + (m(3) * qJ(2) + mrSges(3,3)) * (t97 ^ 2 + t98 ^ 2) + (m(5) + m(4)) * (t73 * t89 - t74 * t87);
t133 = t5 * qJD(1);
t91 = -t99 * qJ(4) + t101 * t102;
t92 = t101 * qJ(4) + t99 * t102;
t103 = t70 * t162 + (t114 * t91 - t67 * t92) * t159;
t108 = -t187 / 0.2e1 + t70 * t161;
t9 = t89 * mrSges(5,1) - t103 + t108 + t118 - t119 + t80;
t131 = t9 * qJD(1);
t12 = (m(6) * t36 - t175 + t32) * t89;
t130 = t12 * qJD(1);
t17 = 0.2e1 * t170 * mrSges(6,1) + 0.2e1 * t124;
t129 = t17 * qJD(1);
t106 = m(6) * (t101 * t114 - t67 * t99);
t20 = -t106 / 0.2e1 + (t160 + t162) * t163;
t128 = t20 * qJD(1);
t31 = t92 * mrSges(6,1) + t91 * mrSges(6,2);
t127 = t31 * qJD(5);
t126 = t99 * t114 * mrSges(6,3);
t123 = t180 + t181;
t122 = t169 + t170;
t3 = t123 * Ifges(6,5) + t122 * Ifges(6,6);
t113 = t3 * qJD(1) + t31 * qJD(3);
t40 = mrSges(5,3) + m(5) * qJ(4) + m(6) * (t92 * t101 - t91 * t99) + t115;
t105 = t182 * t159;
t107 = m(6) * t182;
t7 = t105 - t107 / 0.2e1 + (t123 * t101 - t122 * t99) * mrSges(6,3);
t112 = -t7 * qJD(1) + t40 * qJD(3);
t48 = t126 / 0.2e1;
t11 = -t126 / 0.2e1 + t48;
t2 = -t185 + t172 * t169 + t33 * t170 + (t34 - t184) * t180;
t111 = t2 * qJD(1) + t11 * qJD(4);
t110 = -t11 * qJD(1) - qJD(3) * t115;
t19 = t89 * t161 + t106 / 0.2e1 + (t160 - m(5) / 0.4e1) * t163;
t18 = t124 + t178 / 0.2e1;
t13 = t103 + t108;
t10 = t11 * qJD(5);
t6 = t48 - t87 * mrSges(5,2) + t107 / 0.2e1 + m(5) * t74 + t105 + (0.2e1 * t180 * t101 + t99 * t169) * mrSges(6,3);
t4 = [t5 * qJD(2) + t1 * qJD(3) + t12 * qJD(4) + t2 * qJD(5), t13 * qJD(3) + t19 * qJD(4) + t18 * qJD(5) + t133, t13 * qJD(2) + t6 * qJD(4) - t189 + t135 + (-t74 * mrSges(4,1) - t74 * mrSges(5,1) + t73 * mrSges(4,2) - t73 * mrSges(5,3) + 0.2e1 * (t165 * t91 + t92 * t28) * t159 + 0.2e1 * (-pkin(3) * t74 - qJ(4) * t73) * t161 + (-qJ(4) * mrSges(5,2) - Ifges(4,6) + Ifges(5,6)) * t89 + (pkin(3) * mrSges(5,2) - Ifges(5,4) - Ifges(4,5)) * t87 + (t114 * t92 + t67 * t91) * mrSges(6,3) + t188) * qJD(3), t19 * qJD(2) + t6 * qJD(3) + t10 + t130, t18 * qJD(2) - qJD(3) * t188 + t111 + t189; t9 * qJD(3) + t20 * qJD(4) + t17 * qJD(5) - t133, 0, t131, t128, t129; -t9 * qJD(2) - t7 * qJD(4) + t3 * qJD(5) - t135, -t131, t40 * qJD(4) + t127, t112, t113 - t127; -t20 * qJD(2) + t7 * qJD(3) + t10 - t130, -t128, -t112 + t173, 0, -t110 - t173; -t17 * qJD(2) - t3 * qJD(3) - t111, -t129, -qJD(4) * t115 - t113, t110, 0;];
Cq = t4;
