% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:36
% EndTime: 2019-12-05 17:39:39
% DurationCPUTime: 0.99s
% Computational Cost: add. (4322->121), mult. (7343->163), div. (0->0), fcn. (8242->6), ass. (0->79)
t146 = sin(qJ(4));
t94 = sin(pkin(8));
t95 = cos(pkin(8));
t99 = cos(qJ(4));
t80 = -t146 * t95 - t99 * t94;
t81 = -t146 * t94 + t99 * t95;
t97 = sin(qJ(5));
t98 = cos(qJ(5));
t111 = t98 * t80 - t97 * t81;
t181 = t111 * mrSges(6,1);
t117 = t181 / 0.2e1;
t109 = t97 * t80 + t98 * t81;
t96 = -pkin(1) - qJ(3);
t147 = -pkin(6) + t96;
t84 = t147 * t94;
t85 = t147 * t95;
t69 = -t146 * t84 + t99 * t85;
t105 = -t81 * pkin(7) + t69;
t70 = t146 * t85 + t99 * t84;
t45 = t80 * pkin(7) + t70;
t160 = t105 * t98 - t97 * t45;
t32 = t105 * t97 + t98 * t45;
t5 = -t32 * mrSges(6,1) - t160 * mrSges(6,2) + Ifges(6,5) * t111 - Ifges(6,6) * t109;
t185 = t5 * qJD(5);
t173 = t109 * mrSges(6,1);
t113 = t111 * mrSges(6,2) + t173;
t145 = Ifges(6,4) * t109;
t151 = -t109 / 0.2e1;
t167 = t111 / 0.2e1;
t175 = t109 / 0.2e1;
t87 = t94 * pkin(3) + qJ(2);
t73 = -t80 * pkin(4) + t87;
t2 = (Ifges(6,2) * t111 + t145) * t151 + (Ifges(6,1) * t111 - t145) * t175 + (0.2e1 * Ifges(6,4) * t111 + (Ifges(6,1) - Ifges(6,2)) * t109) * t167 + t73 * t113;
t164 = t109 * mrSges(6,2);
t184 = t181 - t164;
t156 = -m(6) / 0.2e1;
t183 = m(6) * t73;
t116 = t173 / 0.2e1;
t179 = t151 + t175;
t178 = t109 ^ 2 + t111 ^ 2;
t168 = t109 * t97 + t111 * t98;
t158 = t81 ^ 2;
t157 = -m(5) / 0.2e1;
t155 = m(6) * pkin(4);
t4 = (t109 * t167 - t111 * t175) * mrSges(6,3);
t149 = t4 * qJD(4);
t136 = t109 * t98;
t132 = t111 * t97;
t125 = t94 ^ 2 + t95 ^ 2;
t124 = t2 * qJD(1);
t114 = m(4) * t125;
t115 = t80 ^ 2 + t158;
t9 = t178 * mrSges(6,3) + t115 * mrSges(5,3) + t125 * mrSges(4,3) + m(6) * (-t109 * t160 + t111 * t32) + m(5) * (-t69 * t81 + t70 * t80) - t96 * t114;
t123 = t9 * qJD(1);
t103 = t178 * t156 + t115 * t157 - t114 / 0.2e1;
t110 = -m(4) / 0.2e1 + t157 + t156;
t15 = t103 + t110;
t122 = t15 * qJD(1);
t78 = t81 * mrSges(5,1);
t18 = t80 * mrSges(5,2) + t78 + (t81 / 0.2e1 + t136 / 0.2e1 - t132 / 0.2e1) * t155 + t113;
t121 = t18 * qJD(1);
t106 = t80 * mrSges(5,1) - t81 * mrSges(5,2) + t184;
t20 = t94 * mrSges(4,1) + t95 * mrSges(4,2) + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(5) * t87 + t183 - t106;
t120 = t20 * qJD(1);
t21 = 0.2e1 * t167 * mrSges(6,2) + 0.2e1 * t116;
t119 = t21 * qJD(1);
t1 = t87 * t78 - t158 * Ifges(5,4) + (t87 * mrSges(5,2) + Ifges(5,4) * t80 + (Ifges(5,1) - Ifges(5,2)) * t81) * t80 + (t183 - t184) * t81 * pkin(4) + t2;
t108 = t1 * qJD(1) - t4 * qJD(2);
t107 = t4 * qJD(1);
t12 = t117 - t181 / 0.2e1 + t179 * mrSges(6,2);
t6 = t179 * Ifges(6,6) + (t167 - t111 / 0.2e1) * Ifges(6,5);
t86 = (mrSges(6,1) * t97 + mrSges(6,2) * t98) * pkin(4);
t104 = -t6 * qJD(1) - t12 * qJD(2) + t86 * qJD(4);
t82 = t86 * qJD(5);
t33 = (t81 + t132 - t136) * t155 / 0.2e1;
t22 = t116 - t173 / 0.2e1;
t14 = t103 - t110;
t13 = -t164 + 0.2e1 * t117;
t3 = [t20 * qJD(2) + t9 * qJD(3) + t1 * qJD(4) + t2 * qJD(5), t14 * qJD(3) + t120 - t149, t14 * qJD(2) + t33 * qJD(4) + t22 * qJD(5) + t123, t33 * qJD(3) + t185 + t108 + (-t70 * mrSges(5,1) - t69 * mrSges(5,2) + Ifges(5,5) * t80 - Ifges(5,6) * t81 + (m(6) * (t160 * t97 - t32 * t98) - t168 * mrSges(6,3)) * pkin(4) + t5) * qJD(4), t22 * qJD(3) + t5 * qJD(4) + t124 + t185; t15 * qJD(3) - t120 - t149, 0, t122, (t168 * t155 + t106) * qJD(4) + t13 * qJD(5) - t107, t13 * qJD(4) + t184 * qJD(5); -t15 * qJD(2) + t18 * qJD(4) + t21 * qJD(5) - t123, -t122, 0, t121, t119; -t18 * qJD(3) + t6 * qJD(5) - t108, t12 * qJD(5) + t107, -t121, -t82, -t104 - t82; -t21 * qJD(3) - t6 * qJD(4) - t124, -t12 * qJD(4), -t119, t104, 0;];
Cq = t3;
