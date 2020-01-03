% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:14
% DurationCPUTime: 1.78s
% Computational Cost: add. (1278->210), mult. (3542->293), div. (0->0), fcn. (2170->4), ass. (0->104)
t124 = Ifges(4,4) + Ifges(5,4);
t144 = Ifges(4,1) + Ifges(5,1);
t148 = Ifges(4,5) + Ifges(5,5);
t143 = Ifges(4,2) + Ifges(5,2);
t147 = Ifges(4,6) + Ifges(5,6);
t102 = sin(qJ(3));
t103 = sin(qJ(2));
t104 = cos(qJ(3));
t105 = cos(qJ(2));
t84 = -t102 * t103 + t104 * t105;
t72 = t84 * qJD(1);
t150 = t124 * t72;
t85 = t102 * t105 + t104 * t103;
t73 = t85 * qJD(1);
t149 = t124 * t73;
t101 = qJD(2) + qJD(3);
t146 = t147 * t101 + t143 * t72 + t149;
t145 = t148 * t101 + t144 * t73 + t150;
t136 = -pkin(6) - pkin(5);
t97 = t136 * t105;
t90 = qJD(1) * t97;
t74 = t102 * t90;
t96 = t136 * t103;
t89 = qJD(1) * t96;
t80 = qJD(2) * pkin(2) + t89;
t47 = t104 * t80 + t74;
t67 = t73 * qJ(4);
t21 = t47 - t67;
t140 = -t72 / 0.2e1;
t137 = t73 / 0.2e1;
t99 = -pkin(2) * t105 - pkin(1);
t95 = qJD(1) * t99;
t135 = m(4) * t95;
t134 = pkin(1) * mrSges(3,1);
t133 = pkin(1) * mrSges(3,2);
t132 = pkin(3) * t73;
t129 = mrSges(4,3) * t72;
t128 = mrSges(5,3) * t72;
t127 = t73 * mrSges(4,3);
t50 = t104 * t89 + t74;
t55 = t102 * t96 - t104 * t97;
t123 = Ifges(3,4) * t103;
t122 = qJ(4) * t72;
t52 = t101 * t85;
t42 = t52 * qJD(1);
t121 = t102 * t42;
t77 = t104 * t90;
t120 = Ifges(3,5) * qJD(2);
t119 = Ifges(3,6) * qJD(2);
t118 = qJD(2) * mrSges(3,1);
t117 = qJD(2) * mrSges(3,2);
t116 = qJD(1) * t103;
t115 = qJD(1) * t105;
t114 = qJD(2) * t103;
t113 = qJD(3) * t102;
t112 = qJD(3) * t104;
t111 = pkin(2) * t114;
t110 = qJD(2) * t136;
t109 = t120 / 0.2e1;
t108 = -t119 / 0.2e1;
t49 = -t102 * t89 + t77;
t54 = t102 * t97 + t104 * t96;
t107 = qJD(1) * t110;
t48 = t102 * t80 - t77;
t81 = t103 * t107;
t82 = t105 * t107;
t8 = t102 * t82 + t104 * t81 + t112 * t80 + t113 * t90;
t91 = t103 * t110;
t92 = t105 * t110;
t10 = t102 * t92 + t104 * t91 + t112 * t96 + t113 * t97;
t9 = -qJD(3) * t48 - t102 * t81 + t104 * t82;
t11 = -qJD(3) * t55 - t102 * t91 + t104 * t92;
t51 = t101 * t84;
t16 = pkin(3) * t101 + t21;
t2 = -qJ(4) * t42 + qJD(4) * t72 + t8;
t41 = t51 * qJD(1);
t3 = -qJ(4) * t41 - qJD(4) * t73 + t9;
t53 = -t72 * pkin(3) + qJD(4) + t95;
t106 = t9 * mrSges(4,1) + t3 * mrSges(5,1) - t8 * mrSges(4,2) - t2 * mrSges(5,2) + t16 * t128 - t53 * (mrSges(5,1) * t73 + mrSges(5,2) * t72) - t95 * (mrSges(4,1) * t73 + mrSges(4,2) * t72) + t47 * t129 - t147 * t42 + t148 * t41 - (t144 * t72 - t149) * t73 / 0.2e1 + t146 * t137 - (-t147 * t73 + t148 * t72) * t101 / 0.2e1 + (-t143 * t73 + t145 + t150) * t140;
t100 = Ifges(3,4) * t115;
t98 = pkin(2) * t104 + pkin(3);
t94 = mrSges(3,3) * t115 - t117;
t93 = -mrSges(3,3) * t116 + t118;
t71 = Ifges(3,1) * t116 + t100 + t120;
t70 = t119 + (Ifges(3,2) * t105 + t123) * qJD(1);
t61 = -pkin(3) * t84 + t99;
t60 = mrSges(4,1) * t101 - t127;
t59 = mrSges(5,1) * t101 - t73 * mrSges(5,3);
t58 = -mrSges(4,2) * t101 + t129;
t57 = -mrSges(5,2) * t101 + t128;
t56 = pkin(2) * t116 + t132;
t46 = -mrSges(4,1) * t72 + mrSges(4,2) * t73;
t45 = -mrSges(5,1) * t72 + mrSges(5,2) * t73;
t36 = t41 * mrSges(5,2);
t35 = pkin(3) * t52 + t111;
t34 = qJ(4) * t84 + t55;
t33 = -qJ(4) * t85 + t54;
t26 = pkin(3) * t42 + qJD(1) * t111;
t25 = -t67 + t50;
t24 = t49 - t122;
t22 = t48 + t122;
t5 = -qJ(4) * t51 - qJD(4) * t85 + t11;
t4 = -qJ(4) * t52 + qJD(4) * t84 + t10;
t1 = [(-t16 * t51 + t2 * t84 - t22 * t52 - t3 * t85) * mrSges(5,3) + (-t47 * t51 - t48 * t52 + t8 * t84 - t85 * t9) * mrSges(4,3) + t61 * t36 + m(5) * (t16 * t5 + t2 * t34 + t22 * t4 + t26 * t61 + t3 * t33 + t35 * t53) + m(4) * (t10 * t48 + t11 * t47 + t54 * t9 + t55 * t8) + (-t70 / 0.2e1 - pkin(5) * t94 + t108 + (-0.2e1 * t134 - 0.3e1 / 0.2e1 * t123 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t105) * qJD(1) + (0.2e1 * t135 + t46 + qJD(1) * (-mrSges(4,1) * t84 + mrSges(4,2) * t85)) * pkin(2)) * t114 - (-mrSges(4,1) * t99 - mrSges(5,1) * t61 + mrSges(4,3) * t55 + mrSges(5,3) * t34 + t124 * t85 + t143 * t84) * t42 + (mrSges(4,2) * t99 - mrSges(4,3) * t54 - mrSges(5,3) * t33 + t124 * t84 + t144 * t85) * t41 + (t71 / 0.2e1 - pkin(5) * t93 + t109 + (-0.2e1 * t133 + 0.3e1 / 0.2e1 * Ifges(3,4) * t105) * qJD(1)) * t105 * qJD(2) + t95 * (mrSges(4,1) * t52 + mrSges(4,2) * t51) + t26 * (-t84 * mrSges(5,1) + mrSges(5,2) * t85) + t4 * t57 + t10 * t58 + t5 * t59 + t11 * t60 + t53 * (mrSges(5,1) * t52 + mrSges(5,2) * t51) + t35 * t45 + t145 * t51 / 0.2e1 - t146 * t52 / 0.2e1 + (t124 * t51 - t143 * t52) * t72 / 0.2e1 + (-t124 * t52 + t144 * t51) * t137 + (-t147 * t52 + t148 * t51) * t101 / 0.2e1; ((t109 - t100 / 0.2e1 - t71 / 0.2e1 + qJD(1) * t133 + (t93 - t118) * pkin(5)) * t105 + (t108 + t70 / 0.2e1 + (t134 + t123 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t105) * qJD(1) + (t94 + t117) * pkin(5) + (-t46 - t135) * pkin(2)) * t103) * qJD(1) + (-mrSges(5,3) * t121 + (-t104 * t41 - t121) * mrSges(4,3) + ((t57 + t58) * t104 + (-t59 - t60) * t102) * qJD(3) + m(4) * (t102 * t8 + t104 * t9 + t112 * t48 - t113 * t47)) * pkin(2) + (t22 * t73 - t41 * t98) * mrSges(5,3) - m(4) * (t47 * t49 + t48 * t50) + t106 + t48 * t127 - t56 * t45 - t25 * t57 - t50 * t58 - t24 * t59 - t49 * t60 + ((t102 * t2 + t112 * t22 - t113 * t16) * pkin(2) + t3 * t98 - t16 * t24 - t22 * t25 - t53 * t56) * m(5); (-t53 * t132 - (-t16 + t21) * t22 + t3 * pkin(3)) * m(5) + t106 + (t48 * mrSges(4,3) + t22 * mrSges(5,3) - pkin(3) * t45) * t73 - t21 * t57 - t47 * t58 + t22 * t59 + t48 * t60 - pkin(3) * t41 * mrSges(5,3); t42 * mrSges(5,1) - t72 * t57 + t73 * t59 + t36 + 0.2e1 * (t26 / 0.2e1 + t16 * t137 + t22 * t140) * m(5);];
tauc = t1(:);
