% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP3
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
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:03
% DurationCPUTime: 0.74s
% Computational Cost: add. (673->147), mult. (1328->206), div. (0->0), fcn. (518->4), ass. (0->77)
t125 = mrSges(4,1) + mrSges(5,1);
t101 = pkin(1) * qJD(1);
t59 = qJD(1) + qJD(2);
t63 = sin(qJ(2));
t50 = t59 * pkin(6) + t63 * t101;
t62 = sin(qJ(3));
t65 = cos(qJ(2));
t100 = pkin(1) * qJD(2);
t87 = qJD(1) * t100;
t82 = t65 * t87;
t64 = cos(qJ(3));
t97 = qJD(3) * t64;
t10 = t50 * t97 + t62 * t82;
t111 = t10 * t62;
t53 = t64 * t82;
t98 = qJD(3) * t62;
t9 = -t50 * t98 + t53;
t124 = t64 * t9 + t111;
t123 = (t62 ^ 2 + t64 ^ 2) * t50;
t104 = t62 * t50;
t7 = t53 + (qJD(4) - t104) * qJD(3);
t122 = t64 * t7 + t111;
t106 = t59 * t62;
t103 = (-mrSges(5,2) - mrSges(4,3)) * t106 + t125 * qJD(3);
t105 = t59 * t64;
t90 = mrSges(5,2) * t105;
t49 = qJD(3) * mrSges(5,3) + t90;
t102 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t105 + t49;
t26 = -qJD(3) * pkin(3) + qJD(4) + t104;
t92 = qJD(3) * qJ(4);
t37 = t64 * t50 + t92;
t119 = t103 * t62 - m(5) * (t26 * t62 + t37 * t64) - t102 * t64;
t116 = t65 * pkin(1);
t115 = Ifges(4,4) * t62;
t113 = Ifges(5,5) * t62;
t112 = Ifges(5,5) * t64;
t110 = t26 * mrSges(5,2);
t99 = qJD(2) * t63;
t96 = t37 * qJD(3);
t95 = t62 * qJD(4);
t94 = -qJD(1) - t59;
t93 = -qJD(2) + t59;
t89 = t65 * t101;
t88 = t106 / 0.2e1;
t85 = t98 / 0.2e1;
t83 = t63 * t87;
t81 = t59 * t85;
t79 = -t64 * mrSges(4,1) + t62 * mrSges(4,2);
t78 = -t64 * mrSges(5,1) - t62 * mrSges(5,3);
t77 = pkin(3) * t62 - qJ(4) * t64;
t52 = -t64 * pkin(3) - t62 * qJ(4) - pkin(2);
t76 = (Ifges(5,1) * t62 - t112) * t59;
t75 = (Ifges(4,2) * t64 + t115) * t59;
t74 = (mrSges(4,1) * t62 + mrSges(4,2) * t64) * qJD(3);
t73 = (mrSges(5,1) * t62 - mrSges(5,3) * t64) * qJD(3);
t38 = pkin(3) * t98 - t64 * t92 - t95;
t68 = (-m(5) * t26 + t103) * t64 + (m(5) * t37 + t102) * t62;
t54 = Ifges(5,5) * t106;
t27 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t105 + t54;
t28 = Ifges(4,6) * qJD(3) + t75;
t29 = Ifges(5,4) * qJD(3) + t76;
t3 = t83 + (t77 * qJD(3) - t95) * t59;
t55 = Ifges(4,4) * t105;
t30 = Ifges(4,1) * t106 + Ifges(4,5) * qJD(3) + t55;
t51 = -t59 * pkin(2) - t89;
t8 = t52 * t59 - t89;
t67 = t79 * t83 + (-Ifges(5,3) * t64 + t113) * t81 + t27 * t85 + t8 * t73 + t51 * t74 + t3 * t78 + (-t59 * (Ifges(5,3) * t62 + t112) + t110) * t97 + ((Ifges(5,4) + Ifges(4,5)) * t64 + (-Ifges(4,6) + Ifges(5,6)) * t62) * qJD(3) ^ 2 / 0.2e1 - (t75 + t28) * t98 / 0.2e1 + 0.2e1 * (t113 - t115 + (Ifges(4,1) + Ifges(5,1)) * t64) * t81 + t124 * mrSges(4,3) + (-t96 * t62 + t122) * mrSges(5,2) + ((0.3e1 * Ifges(4,4) * t64 + (Ifges(4,1) - 0.2e1 * Ifges(4,2)) * t62) * t59 + t76 + t30 + t29) * t97 / 0.2e1;
t66 = m(4) * t124 + m(5) * t122 - t68 * qJD(3);
t57 = -pkin(2) - t116;
t45 = t52 - t116;
t41 = t77 * t59;
t40 = t79 * t59;
t39 = t78 * t59;
t32 = t59 * t74;
t31 = t59 * t73;
t21 = pkin(1) * t99 + t38;
t1 = [((m(4) * (qJD(1) * t57 + t51) + t40 + t94 * mrSges(3,1)) * t63 + (m(4) * t123 + t94 * mrSges(3,2) - t119) * t65) * t100 + m(5) * (t8 * t21 + t3 * t45) + t67 + t66 * (t63 * pkin(1) + pkin(6)) + t21 * t39 + t45 * t31 + t57 * t32; ((-m(5) * t8 + t93 * mrSges(3,1) - t39 - t40) * t63 + (t93 * mrSges(3,2) + t119) * t65 + (-pkin(2) * t99 - t65 * t123 - t51 * t63) * m(4)) * t101 + m(5) * (t3 * t52 + t8 * t38) + t67 + t66 * pkin(6) - pkin(2) * t32 + t38 * t39 + t52 * t31; -t9 * mrSges(4,2) + t7 * mrSges(5,3) + qJD(4) * t49 - t41 * t39 - t125 * t10 + t68 * t50 + ((-t29 / 0.2e1 - t30 / 0.2e1 - t55 / 0.2e1 - t110 + t8 * mrSges(5,3) - t51 * mrSges(4,2) + Ifges(5,5) * t105 / 0.2e1) * t64 + (-t27 / 0.2e1 + t28 / 0.2e1 - t54 / 0.2e1 + t37 * mrSges(5,2) - t8 * mrSges(5,1) - t51 * mrSges(4,1) + Ifges(4,4) * t88 + (-Ifges(5,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1) * t105) * t62 + ((-pkin(3) * mrSges(5,2) + Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * t64 + (-qJ(4) * mrSges(5,2) + Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * t62) * qJD(3)) * t59 + (-t10 * pkin(3) + t7 * qJ(4) + t37 * qJD(4) - t8 * t41) * m(5); t39 * t106 + (-t49 + t90) * qJD(3) + 0.2e1 * (t10 / 0.2e1 + t8 * t88 - t96 / 0.2e1) * m(5);];
tauc = t1(:);
