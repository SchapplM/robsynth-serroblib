% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRR4
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:15
% EndTime: 2019-12-31 16:50:17
% DurationCPUTime: 0.59s
% Computational Cost: add. (1280->148), mult. (2881->225), div. (0->0), fcn. (2224->6), ass. (0->92)
t84 = sin(qJ(4));
t86 = cos(qJ(4));
t105 = t86 * mrSges(5,1) - t84 * mrSges(5,2);
t145 = -mrSges(4,1) - t105;
t125 = Ifges(5,6) * t84;
t78 = Ifges(5,5) * t86;
t144 = Ifges(4,4) - t78 / 0.2e1 + t125 / 0.2e1;
t81 = t84 ^ 2;
t82 = t86 ^ 2;
t143 = mrSges(5,3) * (t82 / 0.2e1 + t81 / 0.2e1);
t79 = Ifges(5,4) * t86;
t142 = -Ifges(5,2) * t84 + t79;
t66 = Ifges(5,1) * t84 + t79;
t141 = -pkin(6) / 0.2e1;
t118 = t86 * mrSges(5,2);
t121 = t84 * mrSges(5,1);
t64 = t118 + t121;
t140 = t64 / 0.2e1;
t139 = -t84 / 0.2e1;
t137 = t84 / 0.2e1;
t85 = sin(qJ(3));
t136 = t85 / 0.2e1;
t135 = -t86 / 0.2e1;
t134 = t86 / 0.2e1;
t87 = cos(qJ(3));
t133 = -t87 / 0.2e1;
t132 = t87 / 0.2e1;
t113 = t81 + t82;
t131 = m(5) * (-0.1e1 + t113) * t87 * t85;
t130 = t85 * pkin(3);
t128 = Ifges(5,4) * t84;
t127 = Ifges(5,5) * t87;
t124 = Ifges(5,6) * t87;
t120 = t84 * t85;
t69 = -t87 * pkin(6) + t130;
t76 = sin(pkin(7)) * pkin(1) + pkin(5);
t20 = t76 * t120 + t86 * t69;
t123 = t20 * t84;
t122 = t76 * t85;
t119 = t84 * t87;
t117 = t86 * mrSges(5,3);
t114 = t87 * mrSges(5,3);
t60 = -t85 * mrSges(5,2) - t84 * t114;
t116 = t86 * t60;
t115 = t86 * t87;
t48 = t105 * t85;
t59 = t87 * mrSges(5,2) - mrSges(5,3) * t120;
t61 = -t87 * mrSges(5,1) - t85 * t117;
t8 = t48 * t132 + (t61 * t134 + t59 * t137 + t85 * t143) * t85;
t112 = t8 * qJD(1);
t111 = -t119 / 0.2e1;
t110 = t115 / 0.2e1;
t108 = -cos(pkin(7)) * pkin(1) - pkin(2);
t107 = t113 * t87;
t104 = Ifges(5,1) * t86 - t128;
t65 = Ifges(5,2) * t86 + t128;
t102 = Ifges(5,5) * t84 + Ifges(5,6) * t86;
t52 = -t87 * pkin(3) - t85 * pkin(6) + t108;
t14 = -t76 * t119 + t86 * t52;
t15 = t76 * t115 + t84 * t52;
t101 = t14 * t84 - t15 * t86;
t21 = -t86 * t122 + t84 * t69;
t100 = t21 * t86 - t123;
t37 = t142 * t85 - t124;
t38 = Ifges(5,6) * t85 + t142 * t87;
t39 = t104 * t85 - t127;
t40 = Ifges(5,5) * t85 + t104 * t87;
t49 = t87 * t64;
t62 = t85 * mrSges(5,1) - t86 * t114;
t1 = m(5) * (t14 * t20 + t15 * t21) + t21 * t59 + t15 * t60 + t20 * t61 + t14 * t62 + (t108 * mrSges(4,1) + t40 * t134 + t38 * t139 - t144 * t85 + t76 * t49) * t85 + (t108 * mrSges(4,2) + t39 * t134 + t37 * t139 + (-Ifges(5,3) + Ifges(4,1) - Ifges(4,2) + (m(5) * t76 + t64) * t76) * t85 + t144 * t87) * t87;
t6 = t85 ^ 2 * t140 + t49 * t133 + m(5) * ((-t76 * t87 - t101) * t87 + (t100 + t122) * t85) / 0.2e1 + t59 * t110 + t116 * t136 + t61 * t111 - t62 * t120 / 0.2e1;
t99 = t1 * qJD(1) + t6 * qJD(2);
t50 = t85 * t65;
t51 = t85 * t66;
t4 = t14 * t59 - t15 * t61 + (t101 * mrSges(5,3) + t102 * t132 - t51 * t134 + t37 * t135 + t76 * t48 + (t39 - t50) * t139) * t85;
t98 = t4 * qJD(1) - t8 * qJD(2);
t97 = -pkin(3) * t48 / 0.2e1 - t87 * t78 / 0.4e1;
t96 = -t37 / 0.4e1 - t51 / 0.4e1 + t59 * t141;
t95 = t61 * t141 + t39 / 0.4e1 - t50 / 0.4e1;
t93 = -t20 * mrSges(5,1) / 0.2e1 + t21 * mrSges(5,2) / 0.2e1;
t92 = -t121 / 0.2e1 - t118 / 0.2e1;
t91 = t66 * t134 + t65 * t139;
t90 = t6 * qJD(1) + qJD(2) * t131;
t11 = pkin(3) * t64 + t104 * t139 + t135 * t142 - t91;
t16 = (t140 + t92) * t87;
t88 = -pkin(6) * t143 + t76 * t140 + (t104 / 0.4e1 - t65 / 0.4e1) * t86 - (t66 + t142) * t84 / 0.4e1;
t3 = (-t127 / 0.2e1 + t95) * t86 + (0.3e1 / 0.4e1 * t124 + t96) * t84 + (-Ifges(5,3) / 0.2e1 + t88) * t85 + t93 + t97;
t89 = t3 * qJD(1) - t16 * qJD(2) - t11 * qJD(3);
t17 = t64 * t133 + t92 * t87;
t5 = t6 * qJD(3) - t8 * qJD(4);
t2 = Ifges(5,5) * t110 + Ifges(5,6) * t111 + Ifges(5,3) * t136 + (t124 / 0.4e1 + t96) * t84 + t95 * t86 + t88 * t85 - t93 + t97;
t7 = [t1 * qJD(3) + t4 * qJD(4), t5, t2 * qJD(4) + t99 + (mrSges(4,2) * t122 - mrSges(5,3) * t123 - Ifges(4,6) * t85 - pkin(3) * t49 + t102 * t136 + t21 * t117 + t38 * t134 + t40 * t137 + (m(5) * t100 - t84 * t62 + t116) * pkin(6) + (Ifges(4,5) + (-m(5) * pkin(3) + t145) * t76 + t91) * t87) * qJD(3), t2 * qJD(3) + (-t15 * mrSges(5,1) - t14 * mrSges(5,2) - t85 * t102) * qJD(4) + t98; t5, qJD(3) * t131, (-t87 * mrSges(4,2) + m(5) * (pkin(6) * t107 - t130) + mrSges(5,3) * t107 + t145 * t85) * qJD(3) + t17 * qJD(4) + t90, t17 * qJD(3) - t48 * qJD(4) - t112; t3 * qJD(4) - t99, -t16 * qJD(4) - t90, -t11 * qJD(4), (-t105 * pkin(6) - t125 + t78) * qJD(4) + t89; -t3 * qJD(3) - t98, t16 * qJD(3) + t112, -t89, 0;];
Cq = t7;
