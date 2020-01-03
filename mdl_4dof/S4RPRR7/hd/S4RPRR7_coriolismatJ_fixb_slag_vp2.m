% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:44
% EndTime: 2019-12-31 16:53:45
% DurationCPUTime: 0.73s
% Computational Cost: add. (2401->157), mult. (5045->226), div. (0->0), fcn. (5175->6), ass. (0->97)
t97 = cos(qJ(4));
t145 = -t97 / 0.2e1;
t90 = Ifges(5,4) * t97;
t95 = sin(qJ(4));
t87 = t95 * Ifges(5,1) + t90;
t157 = t145 * t87;
t138 = Ifges(5,6) * t95;
t89 = Ifges(5,5) * t97;
t106 = -t89 / 0.2e1 + t138 / 0.2e1;
t156 = Ifges(4,4) + t106;
t109 = -Ifges(5,2) * t95 + t90;
t153 = t97 ^ 2;
t154 = t95 ^ 2;
t155 = t153 + t154;
t152 = m(5) / 0.2e1;
t151 = -mrSges(5,1) / 0.2e1;
t150 = mrSges(5,2) / 0.2e1;
t149 = -mrSges(5,3) / 0.2e1;
t141 = cos(qJ(3));
t93 = sin(pkin(7));
t94 = cos(pkin(7));
t96 = sin(qJ(3));
t83 = t141 * t93 + t96 * t94;
t148 = t83 / 0.2e1;
t147 = -t95 / 0.2e1;
t146 = t95 / 0.2e1;
t144 = t97 / 0.2e1;
t143 = pkin(3) * t83;
t81 = -t141 * t94 + t93 * t96;
t142 = pkin(6) * t81;
t140 = Ifges(5,4) * t95;
t127 = pkin(5) + qJ(2);
t113 = t127 * t93;
t84 = t127 * t94;
t57 = t113 * t141 + t84 * t96;
t137 = t57 * t83;
t136 = t95 * mrSges(5,1);
t135 = t95 * mrSges(5,3);
t36 = Ifges(5,6) * t81 + t109 * t83;
t133 = t95 * t36;
t130 = t97 * mrSges(5,3);
t55 = t81 * mrSges(5,1) - t130 * t83;
t132 = t95 * t55;
t131 = t97 * mrSges(5,2);
t110 = Ifges(5,1) * t97 - t140;
t38 = Ifges(5,5) * t81 + t110 * t83;
t129 = t97 * t38;
t120 = t83 * t135;
t53 = -t81 * mrSges(5,2) - t120;
t128 = t97 * t53;
t114 = -pkin(2) * t94 - pkin(1);
t51 = pkin(3) * t81 - pkin(6) * t83 + t114;
t58 = -t113 * t96 + t141 * t84;
t19 = t51 * t97 - t95 * t58;
t20 = t95 * t51 + t58 * t97;
t85 = t131 + t136;
t50 = t85 * t83;
t5 = (mrSges(4,3) * t83 + t50) * t83 + (mrSges(4,3) * t81 - t128 + t132) * t81 + m(5) * (t137 + (t19 * t95 - t20 * t97) * t81) + m(4) * (-t58 * t81 + t137) + (m(3) * qJ(2) + mrSges(3,3)) * (t93 ^ 2 + t94 ^ 2);
t125 = qJD(1) * t5;
t56 = t142 + t143;
t26 = t56 * t97 + t95 * t57;
t27 = t95 * t56 - t57 * t97;
t35 = Ifges(5,6) * t83 - t109 * t81;
t37 = Ifges(5,5) * t83 - t110 * t81;
t49 = t85 * t81;
t52 = -t83 * mrSges(5,2) + t135 * t81;
t54 = t83 * mrSges(5,1) + t130 * t81;
t74 = t83 * mrSges(4,1);
t1 = t27 * t53 + t20 * t52 + t26 * t55 + t19 * t54 + m(5) * (t19 * t26 + t20 * t27 + t57 * t58) + t58 * t50 - t57 * t49 + t114 * t74 + (-t129 / 0.2e1 + t133 / 0.2e1 - t114 * mrSges(4,2) + t156 * t81) * t81 + (t37 * t144 + t35 * t147 - t156 * t83 + (Ifges(5,3) - Ifges(4,1) + Ifges(4,2)) * t81) * t83;
t124 = t1 * qJD(1);
t108 = Ifges(5,5) * t95 + Ifges(5,6) * t97;
t111 = t97 * mrSges(5,1) - t95 * mrSges(5,2);
t86 = t97 * Ifges(5,2) + t140;
t2 = t20 * t55 - t111 * t137 + (t38 * t146 + t36 * t144 + t20 * t130 + t81 * t108 / 0.2e1 + (t86 * t147 - t157) * t83) * t83 + (-t53 - t120) * t19;
t123 = t2 * qJD(1);
t101 = (t26 * t97 + t95 * t27) * t152 + t52 * t146 + t54 * t144;
t99 = (-t155 * t142 - t143) * t152 - t83 * t111 / 0.2e1;
t6 = t74 + (-mrSges(4,2) + (t153 / 0.2e1 + t154 / 0.2e1) * mrSges(5,3)) * t81 - t99 + t101;
t122 = t6 * qJD(1);
t102 = (t131 / 0.2e1 + t136 / 0.2e1) * t81;
t103 = -t128 / 0.2e1 + t132 / 0.2e1;
t9 = t102 + t103;
t121 = t9 * qJD(1);
t119 = pkin(6) * t149;
t115 = Ifges(5,1) / 0.4e1 - Ifges(5,2) / 0.4e1;
t112 = t89 - t138;
t104 = t146 * t86 + t157;
t25 = pkin(3) * t85 + t109 * t145 + t110 * t147 + t104;
t100 = (t145 * t55 + t147 * t53) * pkin(6) + t57 * t85 / 0.2e1 - t133 / 0.4e1 + t129 / 0.4e1;
t105 = t150 * t27 + t151 * t26;
t98 = (-t87 / 0.4e1 - t90 / 0.4e1 + pkin(3) * t150 + (t119 - t115) * t95) * t95 + (-0.3e1 / 0.4e1 * t140 - t86 / 0.4e1 + pkin(3) * t151 + (t119 + t115) * t97) * t97;
t4 = (-0.3e1 / 0.4e1 * t138 + 0.3e1 / 0.4e1 * t89) * t81 + (-Ifges(5,3) / 0.2e1 + t98) * t83 + t100 + t105;
t107 = -t4 * qJD(1) + t25 * qJD(3);
t10 = t102 - t103;
t7 = t155 * t81 * t149 + t101 + t99;
t3 = Ifges(5,3) * t148 + t98 * t83 + t100 - t105 + (t112 / 0.4e1 + t106) * t81;
t8 = [qJD(2) * t5 + qJD(3) * t1 - qJD(4) * t2, qJD(3) * t7 + qJD(4) * t10 + t125, t7 * qJD(2) + t3 * qJD(4) + t124 + (t27 * t130 - t26 * t135 + t108 * t148 + t37 * t146 + t35 * t144 + pkin(3) * t49 - Ifges(4,6) * t83 + t57 * mrSges(4,2) + (-Ifges(4,5) + t104) * t81 + (-m(5) * pkin(3) - mrSges(4,1) - t111) * t58 + (m(5) * (-t26 * t95 + t27 * t97) + t97 * t52 - t95 * t54) * pkin(6)) * qJD(3), -t123 + t10 * qJD(2) + t3 * qJD(3) + (-t20 * mrSges(5,1) - t19 * mrSges(5,2) - t108 * t83) * qJD(4); qJD(3) * t6 - qJD(4) * t9 - t125, 0, t122, -qJD(4) * t85 - t121; -qJD(2) * t6 + qJD(4) * t4 - t124, -t122, -t25 * qJD(4), (-pkin(6) * t111 + t112) * qJD(4) - t107; qJD(2) * t9 - qJD(3) * t4 + t123, t121, t107, 0;];
Cq = t8;
