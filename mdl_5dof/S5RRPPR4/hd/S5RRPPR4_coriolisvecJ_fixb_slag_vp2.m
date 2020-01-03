% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:42
% DurationCPUTime: 0.78s
% Computational Cost: add. (1147->144), mult. (1819->207), div. (0->0), fcn. (807->6), ass. (0->81)
t61 = sin(qJ(5));
t63 = cos(qJ(5));
t58 = qJD(1) + qJD(2);
t65 = -pkin(2) - pkin(3);
t103 = pkin(1) * qJD(1);
t64 = cos(qJ(2));
t93 = t64 * t103;
t77 = qJD(3) - t93;
t33 = t58 * t65 + t77;
t62 = sin(qJ(2));
t94 = t62 * t103;
t45 = qJ(3) * t58 + t94;
t59 = sin(pkin(8));
t60 = cos(pkin(8));
t13 = t33 * t59 + t45 * t60;
t10 = -pkin(7) * t58 + t13;
t7 = qJD(4) * t63 - t10 * t61;
t8 = qJD(4) * t61 + t10 * t63;
t78 = -t61 * t7 + t63 * t8;
t127 = m(6) * t78;
t109 = t58 * t61;
t44 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t109;
t107 = t63 * mrSges(6,3);
t46 = -qJD(5) * mrSges(6,2) - t107 * t58;
t74 = -t61 * t44 + t63 * t46;
t130 = t127 + t74;
t102 = pkin(1) * qJD(2);
t87 = qJD(1) * t102;
t80 = t64 * t87;
t99 = qJD(3) * t58;
t40 = t80 + t99;
t81 = t62 * t87;
t15 = t40 * t59 - t60 * t81;
t70 = (-mrSges(6,1) * t61 - mrSges(6,2) * t63) * qJD(5);
t30 = t58 * t70;
t129 = m(6) * t15 + t30;
t12 = t33 * t60 - t45 * t59;
t128 = m(4) * t45 + m(5) * (-t12 * t59 + t13 * t60);
t126 = mrSges(3,1) + mrSges(4,1);
t76 = mrSges(6,1) * t63 - mrSges(6,2) * t61;
t36 = t76 * t58;
t125 = mrSges(5,1) * t58 + t36;
t113 = mrSges(5,2) * t58;
t124 = t113 + t74;
t16 = t40 * t60 + t59 * t81;
t3 = -qJD(5) * t8 - t16 * t61;
t120 = t3 * t61;
t2 = qJD(5) * t7 + t16 * t63;
t67 = t2 * t63 - t120 + (-t61 * t8 - t63 * t7) * qJD(5);
t97 = qJD(5) * t63;
t98 = qJD(5) * t61;
t123 = m(6) * t67 - t44 * t97 - t46 * t98;
t9 = pkin(4) * t58 - t12;
t122 = m(5) * t12 - m(6) * t9 - t125;
t121 = m(5) * t13 + t124 + t127;
t119 = t59 * t9;
t116 = t7 * mrSges(6,3);
t115 = t8 * mrSges(6,3);
t112 = Ifges(6,4) * t61;
t111 = Ifges(6,4) * t63;
t110 = t15 * t60;
t88 = -pkin(1) * t64 - pkin(2);
t54 = -pkin(3) + t88;
t55 = pkin(1) * t62 + qJ(3);
t105 = t54 * t59 + t55 * t60;
t104 = qJ(3) * t60 + t59 * t65;
t101 = Ifges(6,5) * qJD(5);
t100 = Ifges(6,6) * qJD(5);
t92 = t62 * t102;
t91 = t64 * t102;
t83 = t58 * t93;
t73 = t54 * t60 - t55 * t59;
t72 = -qJ(3) * t59 + t60 * t65;
t71 = (-Ifges(6,2) * t63 - t112) * t58;
t28 = t71 + t100;
t52 = t58 * t111;
t29 = -Ifges(6,1) * t109 + t101 - t52;
t66 = -mrSges(3,2) * t80 - t107 * t2 + mrSges(6,3) * t120 + t16 * mrSges(5,2) + qJD(5) ^ 2 * (-Ifges(6,5) * t63 + Ifges(6,6) * t61) / 0.2e1 + t40 * mrSges(4,3) + t97 * t116 + t9 * t70 - t126 * t81 + (-t58 * (-Ifges(6,1) * t63 + t112) + t115) * t98 + (t28 + t71) * t98 / 0.2e1 + (mrSges(5,1) + t76) * t15 - (t29 + (-0.3e1 * t111 + (-Ifges(6,1) + 0.2e1 * Ifges(6,2)) * t61) * t58) * t97 / 0.2e1;
t53 = qJD(3) + t91;
t41 = -pkin(2) * t58 + t77;
t1 = [t66 + m(5) * (t105 * t16 - t15 * t73) + m(4) * (t40 * t55 + t45 * t53 + (qJD(1) * t88 + t41) * t92) + t129 * (pkin(4) - t73) - t122 * (t53 * t59 - t60 * t92) + t123 * (-pkin(7) + t105) + t121 * (t53 * t60 + t59 * t92) + (-mrSges(3,2) * t91 + t53 * mrSges(4,3) - t126 * t92) * t58; m(5) * (t104 * t16 - t15 * t72) + mrSges(3,2) * t83 + t66 + t126 * t58 * t94 + t123 * (-pkin(7) + t104) + (-t83 + t99) * mrSges(4,3) + (-(t41 * t62 + t45 * t64) * t103 - pkin(2) * t81 + t40 * qJ(3)) * m(4) + t122 * (t59 * t93 - t60 * t94) - t121 * (t59 * t62 + t60 * t64) * t103 + t129 * (pkin(4) - t72) + (m(6) * (t60 * t78 + t119) + t125 * t59 + t124 * t60 + t128) * qJD(3); m(4) * t81 - t60 * t30 + (-t44 * t63 - t46 * t61) * t59 * qJD(5) + m(5) * (t16 * t59 - t110) + m(6) * (t67 * t59 - t110) + (-t59 * t36 + (-mrSges(5,1) * t59 - mrSges(4,3)) * t58 - m(6) * t119 + (-t113 - t130) * t60 - t128) * t58; m(6) * (t2 * t61 + t3 * t63) + ((t61 ^ 2 + t63 ^ 2) * t58 * mrSges(6,3) + t130) * qJD(5); t3 * mrSges(6,1) - t2 * mrSges(6,2) + t8 * t44 - t7 * t46 + ((-t101 / 0.2e1 + t9 * mrSges(6,2) + t29 / 0.2e1 - t52 / 0.2e1 - t116) * t63 + (t100 / 0.2e1 + t9 * mrSges(6,1) - t28 / 0.2e1 - t115 + (t112 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t63) * t58) * t61) * t58;];
tauc = t1(:);
