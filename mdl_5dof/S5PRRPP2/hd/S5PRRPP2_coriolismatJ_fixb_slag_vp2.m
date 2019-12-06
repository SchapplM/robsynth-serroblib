% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:47
% EndTime: 2019-12-05 16:08:51
% DurationCPUTime: 1.06s
% Computational Cost: add. (1664->140), mult. (4078->192), div. (0->0), fcn. (3896->6), ass. (0->85)
t150 = m(6) / 0.2e1;
t125 = cos(pkin(8));
t92 = sin(qJ(3));
t111 = t125 * t92;
t124 = sin(pkin(8));
t94 = cos(qJ(3));
t100 = t124 * t94 + t111;
t144 = t92 * pkin(3);
t109 = t124 * t92;
t68 = -t125 * t94 + t109;
t35 = pkin(4) * t100 + qJ(5) * t68 + t144;
t166 = t35 * t150;
t84 = -pkin(3) * t94 - pkin(2);
t30 = pkin(4) * t68 - qJ(5) * t100 + t84;
t43 = mrSges(6,1) * t68 - mrSges(6,3) * t100;
t165 = m(6) * t30 + t43;
t164 = t100 * mrSges(5,1) - t68 * mrSges(5,2);
t162 = 0.2e1 * t150;
t161 = qJ(4) + pkin(6);
t160 = m(5) + m(6);
t159 = mrSges(5,1) + mrSges(6,1);
t157 = -Ifges(6,5) + Ifges(5,4);
t63 = t68 * mrSges(6,3);
t156 = mrSges(6,1) * t100 + t164 + t63;
t155 = -t94 * mrSges(4,1) + t92 * mrSges(4,2);
t119 = t124 * pkin(3);
t78 = t119 + qJ(5);
t74 = m(6) * t78 + mrSges(6,3);
t130 = t161 * t94;
t46 = t111 * t161 + t124 * t130;
t154 = -t109 * t161 + t125 * t130;
t149 = m(5) * pkin(3);
t120 = t125 * pkin(3);
t83 = -t120 - pkin(4);
t153 = m(6) * t83 - t125 * t149 - t159;
t152 = t124 * t149 - mrSges(5,2) + t74;
t151 = m(5) / 0.2e1;
t95 = cos(qJ(2));
t57 = t100 * t95;
t148 = -t57 / 0.2e1;
t146 = -t95 / 0.2e1;
t93 = sin(qJ(2));
t55 = t68 * t93;
t145 = m(6) * t55;
t143 = mrSges(4,2) * t94;
t142 = t68 * mrSges(6,2);
t141 = t68 * mrSges(5,3);
t140 = t100 * mrSges(6,2);
t139 = t100 * mrSges(5,3);
t137 = t92 * t95;
t136 = t93 * t95;
t129 = t92 ^ 2 + t94 ^ 2;
t56 = t100 * t93;
t58 = t95 * t68;
t6 = m(4) * (-0.1e1 + t129) * t136 + t160 * (t55 * t58 + t56 * t57 - t136);
t128 = t6 * qJD(1);
t122 = t149 / 0.2e1;
t102 = t92 * t122 + t166;
t98 = (t100 * t83 - t68 * t78) * t150 + (-t100 * t125 - t124 * t68) * t122;
t8 = t102 - t98 + t156;
t127 = t8 * qJD(2);
t32 = m(6) * t100;
t123 = t32 * qJD(2);
t121 = t100 * t146;
t115 = -t154 * t58 + t46 * t57;
t44 = mrSges(5,1) * t68 + mrSges(5,2) * t100;
t1 = t30 * t63 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t92 + pkin(3) * t44) * t92 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) * t94 + (Ifges(4,1) - Ifges(4,2)) * t92) * t94 - (-t30 * mrSges(6,1) + t100 * t157) * t100 + (-(Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3)) * t100 + t157 * t68) * t68 + (m(5) * t144 + t164) * t84 + t165 * t35;
t96 = -pkin(3) * t137 * t151 - t95 * t166 + (mrSges(4,1) * t92 + t143 + t156) * t146;
t97 = -mrSges(4,1) * t137 / 0.2e1 + t143 * t146 + (-t125 * t122 + t83 * t150) * t57 - (t78 * t150 - mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1 + t124 * t122) * t58 + t159 * t148;
t2 = -t96 + t97;
t108 = -t2 * qJD(1) + t1 * qJD(2);
t99 = (t151 + t150) * (t100 * t56 + t68 * t55);
t11 = (-m(5) / 0.2e1 - m(6) / 0.2e1) * t93 + t99;
t4 = (t100 ^ 2 + t68 ^ 2) * (mrSges(5,3) + mrSges(6,2)) + t160 * (t100 * t46 - t154 * t68);
t107 = qJD(1) * t11 + qJD(2) * t4;
t13 = t165 * t100;
t24 = (t148 - t121) * m(6);
t106 = qJD(1) * t24 - qJD(2) * t13;
t105 = qJD(3) * t74;
t23 = -m(6) * t121 + t57 * t150;
t14 = t154 * t162 - t142;
t12 = t98 + t102;
t10 = t99 + t160 * t93 / 0.2e1;
t3 = t96 + t97;
t5 = [qJD(2) * t6, t3 * qJD(3) + t10 * qJD(4) + t23 * qJD(5) + t128 + ((t129 * mrSges(4,3) - mrSges(3,2)) * t95 + (-mrSges(3,1) + t43 + t44 + t155) * t93 + m(4) * (t129 * t95 * pkin(6) - t93 * pkin(2)) + 0.2e1 * (t84 * t93 + t115) * t151 + (t30 * t93 + t115) * t162 - (-t141 - t142) * t58 + (t139 + t140) * t57) * qJD(2), t3 * qJD(2) + (-t152 * t56 - t153 * t55 + t155 * t93) * qJD(3) - qJD(5) * t145, t10 * qJD(2), t23 * qJD(2) - qJD(3) * t145; -qJD(3) * t2 + qJD(4) * t11 + qJD(5) * t24 - t128, qJD(3) * t1 + qJD(4) * t4 - qJD(5) * t13, (Ifges(4,5) * t94 - Ifges(4,6) * t92 - t119 * t139 + t120 * t141 - t78 * t140 - t83 * t142 - (-Ifges(6,6) + Ifges(5,6)) * t100 + (-Ifges(6,4) - Ifges(5,5)) * t68 - t152 * t46 + t153 * t154 + t155 * pkin(6)) * qJD(3) + t12 * qJD(4) + t14 * qJD(5) + t108, qJD(3) * t12 + t107, qJD(3) * t14 + t106; qJD(2) * t2, -qJD(4) * t8 - t108, t74 * qJD(5), -t127, t105; -t11 * qJD(2), qJD(3) * t8 - qJD(5) * t32 - t107, t127, 0, -t123; -t24 * qJD(2), qJD(4) * t32 - t106, -t105, t123, 0;];
Cq = t5;
