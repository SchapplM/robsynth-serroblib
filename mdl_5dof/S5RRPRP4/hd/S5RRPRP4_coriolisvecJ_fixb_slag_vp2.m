% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:40
% DurationCPUTime: 1.08s
% Computational Cost: add. (1007->172), mult. (1622->234), div. (0->0), fcn. (588->4), ass. (0->96)
t163 = Ifges(6,4) + Ifges(5,5);
t162 = mrSges(5,1) + mrSges(6,1);
t66 = sin(qJ(4));
t138 = Ifges(6,5) * t66;
t140 = Ifges(5,4) * t66;
t68 = cos(qJ(4));
t141 = Ifges(6,1) * t68;
t63 = qJD(1) + qJD(2);
t161 = (Ifges(5,1) * t68 + t138 - t140 + t141) * t63;
t88 = mrSges(6,1) * t66 - mrSges(6,3) * t68;
t40 = t88 * t63;
t127 = pkin(1) * qJD(1);
t67 = sin(qJ(2));
t115 = t67 * t127;
t52 = t66 * pkin(4) - qJ(5) * t68 + qJ(3);
t8 = t52 * t63 + t115;
t158 = m(6) * t8 + t40;
t160 = t163 * qJD(4) + t161;
t90 = mrSges(5,1) * t66 + mrSges(5,2) * t68;
t41 = t90 * t63;
t152 = mrSges(4,3) * t63 + t41;
t159 = -t152 - t158;
t118 = qJD(4) * qJ(5);
t70 = -pkin(2) - pkin(7);
t69 = cos(qJ(2));
t93 = -t69 * t127 + qJD(3);
t38 = t70 * t63 + t93;
t14 = t38 * t66 + t118;
t128 = t68 * mrSges(5,3);
t129 = t68 * mrSges(6,2);
t153 = (-t128 - t129) * t63 + t162 * qJD(4);
t130 = t66 * mrSges(5,3);
t131 = t66 * mrSges(6,2);
t116 = t63 * t131;
t51 = qJD(4) * mrSges(6,3) - t116;
t154 = -qJD(4) * mrSges(5,2) - t63 * t130 + t51;
t9 = -qJD(4) * pkin(4) - t38 * t68 + qJD(5);
t157 = m(6) * (t14 * t66 - t68 * t9) + t153 * t68 + t154 * t66;
t92 = t14 * t68 + t9 * t66;
t155 = t92 * mrSges(6,2);
t122 = qJD(4) * t66;
t126 = pkin(1) * qJD(2);
t109 = qJD(1) * t126;
t99 = t67 * t109;
t6 = t38 * t122 - t68 * t99;
t144 = t6 * t68;
t121 = qJD(4) * t68;
t7 = t38 * t121 + t66 * t99;
t5 = qJD(4) * qJD(5) + t7;
t95 = t5 * t66 - t144;
t150 = (t66 ^ 2 + t68 ^ 2) * t38;
t146 = pkin(1) * t67;
t139 = Ifges(5,4) * t68;
t136 = Ifges(6,3) * t66;
t48 = qJ(3) * t63 + t115;
t133 = t48 * t69;
t132 = t63 * t68;
t125 = qJD(2) * t67;
t124 = qJD(3) * t63;
t120 = t14 * qJD(4);
t119 = -qJD(2) + t63;
t114 = pkin(1) * t125;
t113 = t69 * t126;
t110 = -pkin(1) * t69 - pkin(2);
t103 = t63 * t114;
t100 = -qJD(5) * t68 + qJD(3);
t98 = t69 * t109;
t94 = -t66 * t7 + t144;
t91 = mrSges(5,1) * t68 - mrSges(5,2) * t66;
t89 = mrSges(6,1) * t68 + mrSges(6,3) * t66;
t87 = pkin(4) * t68 + qJ(5) * t66;
t45 = t98 + t124;
t57 = qJD(3) + t113;
t59 = qJ(3) + t146;
t86 = t45 * t59 + t48 * t57;
t85 = t45 * qJ(3) + t48 * qJD(3);
t37 = pkin(4) * t121 + t66 * t118 + t100;
t84 = t66 * (-Ifges(5,2) * t68 - t140);
t83 = t66 * (Ifges(6,3) * t68 - t138);
t82 = t68 * (-Ifges(5,1) * t66 - t139);
t79 = (-Ifges(5,2) * t66 + t139) * t63;
t78 = t91 * qJD(4);
t77 = t89 * qJD(4);
t72 = m(6) * t92 - t153 * t66 + t154 * t68;
t56 = Ifges(6,5) * t132;
t27 = Ifges(6,6) * qJD(4) + t63 * t136 + t56;
t28 = Ifges(5,6) * qJD(4) + t79;
t3 = t98 + (t87 * qJD(4) + t100) * t63;
t71 = mrSges(4,2) * t99 + t6 * t129 - t5 * t131 + t3 * t88 + t45 * t90 + t48 * t78 + t8 * t77 + ((-Ifges(5,6) + Ifges(6,6)) * t68 - t163 * t66) * qJD(4) ^ 2 / 0.2e1 - (t28 + t79) * t121 / 0.2e1 - (t160 + t161) * t122 / 0.2e1 + (-t84 + t82 + t83) * qJD(4) * t63 + (t27 + (-0.2e1 * Ifges(6,1) * t66 + 0.3e1 * Ifges(6,5) * t68 + t136) * t63) * t121 / 0.2e1;
t58 = -pkin(7) + t110;
t46 = -pkin(2) * t63 + t93;
t44 = t52 + t146;
t32 = t63 * t78;
t31 = t63 * t77;
t16 = t37 + t113;
t1 = [t6 * t128 + t71 + m(6) * (t16 * t8 + t3 * t44 + (t92 * qJD(4) + t95) * t58) - t7 * t130 - t120 * t129 + t59 * t32 + t44 * t31 + t45 * mrSges(4,3) + t16 * t40 + m(4) * t86 + m(5) * (-t94 * t58 + t86) + mrSges(4,2) * t103 + t152 * t57 + t154 * t58 * t121 + (-t63 * t113 - t98) * mrSges(3,2) + (-t103 - t99) * mrSges(3,1) + (-t9 * mrSges(6,2) - t153 * t58) * t122 + (m(4) * (t110 * qJD(1) + t46) + m(5) * t150 + t157) * t114; (t45 + t124) * mrSges(4,3) + t71 + (t72 * t70 - t155) * qJD(4) + m(5) * (-t94 * t70 + t85) + t94 * mrSges(5,3) + qJD(3) * t41 + t52 * t31 + qJ(3) * t32 + t37 * t40 + (-m(5) * (t67 * t150 + t133) + (t119 * mrSges(3,2) + t159) * t69 + (t119 * mrSges(3,1) - t63 * mrSges(4,2) - t157) * t67 + (-pkin(2) * t125 - t46 * t67 - t133) * m(4)) * t127 + m(4) * t85 + m(6) * (t3 * t52 + t8 * t37 + t95 * t70); m(4) * t99 - m(5) * t94 + m(6) * t95 + t72 * qJD(4) + ((-m(4) - m(5)) * t48 + t159) * t63; -t7 * mrSges(5,2) + t5 * mrSges(6,3) + qJD(5) * t51 - t162 * t6 + m(6) * (-t6 * pkin(4) + t5 * qJ(5) + t14 * qJD(5)) - t72 * t38 + (-t48 * t91 - t8 * t89 + t68 * t28 / 0.2e1 + (-t83 / 0.2e1 - t82 / 0.2e1 + t84 / 0.2e1) * t63 + t155 + ((Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1 - qJ(5) * mrSges(6,2)) * t68 + (-Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1 + pkin(4) * mrSges(6,2)) * t66) * qJD(4) + (t141 * t63 + t160) * t66 / 0.2e1 - (t56 + t27) * t68 / 0.2e1) * t63 - t158 * t87 * t63; (-t51 - t116) * qJD(4) + 0.2e1 * (t6 / 0.2e1 - t120 / 0.2e1) * m(6) + t158 * t132;];
tauc = t1(:);
