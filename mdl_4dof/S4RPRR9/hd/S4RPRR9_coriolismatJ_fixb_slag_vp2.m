% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR9_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:09
% EndTime: 2019-12-31 16:56:11
% DurationCPUTime: 0.83s
% Computational Cost: add. (1248->154), mult. (2725->220), div. (0->0), fcn. (2065->4), ass. (0->87)
t103 = sin(qJ(4));
t105 = cos(qJ(4));
t99 = Ifges(5,5) * t105;
t163 = -Ifges(5,6) * t103 + t99;
t104 = sin(qJ(3));
t142 = t105 * mrSges(5,1);
t145 = t103 * mrSges(5,2);
t122 = t142 - t145;
t165 = t104 * t122;
t164 = Ifges(4,4) - t163;
t106 = cos(qJ(3));
t107 = -pkin(1) - pkin(5);
t130 = t106 * t107;
t91 = pkin(3) * t106 + pkin(6) * t104;
t44 = -t103 * t130 + t105 * t91;
t45 = t103 * t91 + t105 * t130;
t162 = -t103 * t44 + t105 * t45;
t100 = Ifges(5,4) * t105;
t149 = Ifges(5,2) * t103;
t161 = t149 - t100;
t151 = Ifges(5,4) * t103;
t155 = -t105 / 0.2e1;
t156 = -t103 / 0.2e1;
t87 = Ifges(5,2) * t105 + t151;
t88 = Ifges(5,1) * t103 + t100;
t160 = (Ifges(5,1) * t105 - t151) * t156 + t103 * t87 / 0.2e1 + t88 * t155;
t102 = t105 ^ 2;
t159 = m(5) / 0.2e1;
t158 = -pkin(6) / 0.2e1;
t101 = t103 ^ 2;
t129 = t101 + t102;
t157 = m(5) * (-0.1e1 + t129) * t106 * t104;
t154 = t105 / 0.2e1;
t153 = -t106 / 0.2e1;
t152 = t104 * pkin(3);
t150 = Ifges(5,5) * t104;
t147 = Ifges(5,6) * t104;
t146 = t103 * mrSges(5,1);
t134 = t103 * t106;
t127 = mrSges(5,3) * t134;
t79 = -mrSges(5,2) * t104 - t127;
t143 = t103 * t79;
t141 = t105 * mrSges(5,2);
t131 = t105 * t106;
t80 = mrSges(5,1) * t104 - mrSges(5,3) * t131;
t139 = t105 * t80;
t132 = t104 * t107;
t84 = -pkin(6) * t106 + qJ(2) + t152;
t40 = -t103 * t132 + t105 * t84;
t41 = t103 * t84 + t105 * t132;
t67 = t122 * t106;
t4 = t41 * t80 + (t107 * t67 + (-Ifges(5,4) * t134 + t150) * t103 + (Ifges(5,4) * t131 + t41 * mrSges(5,3) + t147 + (Ifges(5,1) - Ifges(5,2)) * t134) * t105) * t106 + (-t79 - t127) * t40;
t138 = t4 * qJD(1);
t114 = -t143 / 0.2e1 - t139 / 0.2e1;
t124 = mrSges(5,3) * (-t102 / 0.2e1 - t101 / 0.2e1);
t108 = (t106 * t124 + t114) * t104 + t67 * t153;
t115 = -t145 / 0.2e1 + t142 / 0.2e1;
t7 = t108 - t115;
t137 = t7 * qJD(1);
t123 = -t104 * mrSges(4,1) - t106 * mrSges(4,2);
t10 = m(5) * (t103 * t41 + t105 * t40) + t143 + t139 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) - t123;
t136 = qJD(1) * t10;
t135 = t103 * t104;
t133 = t104 * t105;
t128 = -pkin(3) * t67 / 0.2e1;
t125 = t129 * t106;
t86 = t141 + t146;
t120 = -Ifges(5,5) * t103 - Ifges(5,6) * t105;
t68 = t86 * t104;
t1 = -m(5) * (t40 * t44 + t41 * t45) - t45 * t79 - t44 * t80 + (-t133 * t40 - t135 * t41) * mrSges(5,3) + (-qJ(2) * mrSges(4,1) - t40 * mrSges(5,1) + t41 * mrSges(5,2) + t164 * t106 - t107 * t68) * t106 + (-t86 * t130 + qJ(2) * mrSges(4,2) + (m(5) * t107 ^ 2 + t102 * Ifges(5,1) + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) + (t149 - 0.2e1 * t100) * t103) * t106 - t164 * t104) * t104;
t6 = ((-t103 * t40 + t105 * t41) * t159 + t79 * t154 + t68 / 0.2e1 + t80 * t156) * t106 + (-0.2e1 * t130 + t162) * t104 * t159;
t119 = -t1 * qJD(1) + t6 * qJD(2);
t118 = -t44 * mrSges(5,1) / 0.2e1 + t45 * mrSges(5,2) / 0.2e1;
t117 = t6 * qJD(1) + qJD(2) * t157;
t116 = -t146 / 0.2e1 - t141 / 0.2e1;
t112 = -t87 / 0.4e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.4e1) * t105;
t13 = pkin(3) * t86 - t155 * t161 + t160;
t109 = -t107 * t86 / 0.2e1 + pkin(6) * t124;
t110 = -t88 / 0.4e1 - t100 / 0.4e1 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.4e1) * t103;
t2 = t128 + t114 * pkin(6) + t163 * t104 + (-Ifges(5,3) / 0.2e1 + t112 * t105 + (-0.5e1 / 0.4e1 * t100 + t110) * t103 + t109) * t106 + t118;
t23 = (t86 / 0.2e1 + t116) * t106;
t111 = t2 * qJD(1) - t23 * qJD(2) - t13 * qJD(3);
t24 = t106 * t116 + t153 * t86;
t8 = t108 + t115;
t5 = t6 * qJD(3);
t3 = t128 + t104 * t99 / 0.4e1 - Ifges(5,5) * t133 / 0.2e1 + Ifges(5,6) * t135 / 0.2e1 + (t79 * t158 - t147 / 0.2e1) * t103 + (t150 / 0.4e1 + t80 * t158) * t105 - t118 + (Ifges(5,3) / 0.2e1 + t109 + t110 * t103 + (-0.5e1 / 0.4e1 * t151 + t112) * t105) * t106;
t9 = [qJD(2) * t10 - qJD(3) * t1 - qJD(4) * t4, qJD(4) * t8 + t136 + t5, t3 * qJD(4) + t119 + (pkin(3) * t68 + (t161 * t154 - Ifges(4,5) + (-m(5) * pkin(3) - mrSges(4,1) - t122) * t107 + t160) * t104 + (-t107 * mrSges(4,2) - pkin(6) * t86 - Ifges(4,6) - t120) * t106 + (m(5) * pkin(6) + mrSges(5,3)) * t162) * qJD(3), -t138 + t8 * qJD(2) + t3 * qJD(3) + (-t41 * mrSges(5,1) - t40 * mrSges(5,2) + t106 * t120) * qJD(4); qJD(4) * t7 - t136 + t5, qJD(3) * t157, (-t165 + m(5) * (pkin(6) * t125 - t152) + mrSges(5,3) * t125 + t123) * qJD(3) + t24 * qJD(4) + t117, t24 * qJD(3) - qJD(4) * t165 + t137; qJD(4) * t2 - t119, -t23 * qJD(4) - t117, -t13 * qJD(4), (-pkin(6) * t122 + t163) * qJD(4) + t111; -qJD(2) * t7 - qJD(3) * t2 + t138, qJD(3) * t23 - t137, -t111, 0;];
Cq = t9;
