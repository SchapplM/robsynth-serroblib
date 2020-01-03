% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:03
% EndTime: 2019-12-31 17:52:06
% DurationCPUTime: 1.08s
% Computational Cost: add. (743->181), mult. (1591->251), div. (0->0), fcn. (695->4), ass. (0->91)
t55 = -pkin(1) - pkin(2);
t38 = qJD(1) * t55 + qJD(2);
t50 = sin(pkin(7));
t51 = cos(pkin(7));
t89 = qJ(2) * qJD(1);
t18 = t50 * t38 + t51 * t89;
t14 = -qJD(1) * pkin(6) + t18;
t54 = cos(qJ(4));
t47 = t54 * qJD(3);
t53 = sin(qJ(4));
t10 = -t14 * t53 + t47;
t11 = qJD(3) * t53 + t14 * t54;
t86 = qJD(1) * qJD(2);
t77 = t51 * t86;
t91 = qJD(4) * t53;
t3 = qJD(4) * t47 - t14 * t91 + t54 * t77;
t57 = t11 * qJD(4);
t4 = -t53 * t77 - t57;
t69 = t3 * t54 - t4 * t53;
t90 = qJD(4) * t54;
t119 = t10 * t90 + t11 * t91 - t69;
t88 = qJ(5) * qJD(1);
t8 = t53 * t88 + t10;
t5 = qJD(4) * pkin(4) + t8;
t9 = -t54 * t88 + t11;
t118 = m(6) * (t5 * t53 - t54 * t9);
t116 = -m(3) * qJ(2) - mrSges(3,3);
t17 = t38 * t51 - t50 * t89;
t13 = qJD(1) * pkin(3) - t17;
t94 = qJD(1) * t54;
t12 = pkin(4) * t94 + qJD(5) + t13;
t67 = t54 * mrSges(6,1) - t53 * mrSges(6,2);
t29 = t67 * qJD(1);
t115 = m(6) * t12 + t29;
t87 = qJ(5) * qJD(4);
t1 = (-qJD(5) * t54 + t53 * t87) * qJD(1) + t3;
t92 = qJD(2) * t51;
t70 = -qJD(5) + t92;
t64 = t70 * t53;
t2 = -t57 + (t54 * t87 - t64) * qJD(1);
t114 = t1 * t54 - t2 * t53 - t5 * t90 - t9 * t91;
t95 = qJD(1) * t53;
t33 = qJD(4) * mrSges(6,1) + mrSges(6,3) * t95;
t34 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t95;
t100 = t33 + t34;
t35 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t94;
t36 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t94;
t99 = t35 + t36;
t113 = t100 * t53 - t99 * t54;
t93 = qJD(2) * t50;
t32 = -pkin(4) * t91 + t93;
t26 = t32 * qJD(1);
t112 = m(6) * t26;
t106 = Ifges(5,4) * t53;
t105 = Ifges(5,4) * t54;
t104 = Ifges(6,4) * t53;
t103 = Ifges(6,4) * t54;
t102 = t51 * t53;
t101 = t51 * t54;
t98 = t51 * qJ(2) + t50 * t55;
t28 = -pkin(6) + t98;
t96 = qJ(5) - t28;
t85 = qJD(1) * qJD(4);
t78 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t76 = t53 * t85;
t75 = t54 * t85;
t72 = -t50 * qJ(2) + t51 * t55;
t71 = qJD(4) * t96;
t27 = pkin(3) - t72;
t66 = -t10 * t53 + t11 * t54;
t65 = -t17 * t50 + t18 * t51;
t30 = (mrSges(5,1) * t54 - mrSges(5,2) * t53) * qJD(1);
t63 = (-mrSges(5,1) * t53 - mrSges(5,2) * t54) * qJD(4);
t62 = (-t53 * mrSges(6,1) - t54 * mrSges(6,2)) * qJD(4);
t61 = (-Ifges(5,1) * t53 - t105) * qJD(1);
t60 = (-Ifges(6,1) * t53 - t103) * qJD(1);
t59 = (-Ifges(5,2) * t54 - t106) * qJD(1);
t58 = (-Ifges(6,2) * t54 - t104) * qJD(1);
t56 = qJD(1) ^ 2;
t25 = qJD(1) * t63;
t24 = qJD(1) * t62;
t23 = Ifges(5,5) * qJD(4) + t61;
t22 = Ifges(6,5) * qJD(4) + t60;
t21 = Ifges(5,6) * qJD(4) + t59;
t20 = Ifges(6,6) * qJD(4) + t58;
t19 = pkin(4) * t54 + t27;
t16 = t96 * t54;
t15 = t96 * t53;
t7 = t54 * t71 - t64;
t6 = t53 * t71 + t54 * t70;
t31 = [0.2e1 * t30 * t93 + 0.2e1 * mrSges(4,2) * t77 + m(6) * (-t1 * t16 + t12 * t32 + t15 * t2 + t19 * t26 + t5 * t7 + t6 * t9) + t26 * t67 + t32 * t29 + t7 * t33 + t6 * t35 + t19 * t24 + t27 * t25 + t12 * t62 + t13 * t63 + (-t53 * t34 + t54 * t36) * t92 + ((-Ifges(5,5) - Ifges(6,5)) * t54 + (Ifges(5,6) + Ifges(6,6)) * t53) * qJD(4) ^ 2 / 0.2e1 + (-t104 - t106 + (Ifges(5,1) + Ifges(6,1)) * t54) * t76 + (t103 + t105 + (-Ifges(5,2) - Ifges(6,2)) * t53) * t75 + (m(5) * (t66 * t51 + (qJD(1) * t27 + t13) * t50) + m(4) * ((-t50 * t72 + t51 * t98) * qJD(1) + t65)) * qJD(2) + (m(5) * ((-t10 * t54 - t11 * t53) * qJD(4) + t69) - t34 * t90 - t36 * t91) * t28 + 0.2e1 * (t50 * mrSges(4,1) - t116) * t86 + (t59 + t58 + t21 + t20) * t91 / 0.2e1 - (t61 + t60 + t23 + t22) * t90 / 0.2e1 + t119 * mrSges(5,3) + (t15 * t75 - t16 * t76 - t114) * mrSges(6,3); t116 * t56 + (-t56 * mrSges(4,2) - t112 - t24 - t25) * t51 + (t113 * t51 - m(4) * t65 - m(5) * (-t10 * t102 + t101 * t11) - m(6) * (t101 * t9 - t102 * t5)) * qJD(1) + (-t56 * mrSges(4,1) + (-t100 * t54 - t53 * t99) * qJD(4) + m(6) * t114 + (-t115 - t30) * qJD(1) + ((-t13 - t92) * qJD(1) - t119) * m(5)) * t50; m(5) * (t3 * t53 + t4 * t54) + m(6) * (t1 * t53 + t2 * t54) + (m(5) * t66 - t118 + (mrSges(5,3) + mrSges(6,3)) * qJD(1) * (t53 ^ 2 + t54 ^ 2) - t113) * qJD(4); t4 * mrSges(5,1) - t3 * mrSges(5,2) - t1 * mrSges(6,2) - t10 * t36 + t11 * t34 - t8 * t35 + (m(6) * pkin(4) + mrSges(6,1)) * t2 + (-m(6) * (-t5 + t8) + t33) * t9 + ((t22 / 0.2e1 + t23 / 0.2e1 + t12 * mrSges(6,2) + t13 * mrSges(5,2) - t10 * mrSges(5,3) - t5 * mrSges(6,3) - t78 * t94) * t54 + (-t20 / 0.2e1 - t21 / 0.2e1 + t12 * mrSges(6,1) + t13 * mrSges(5,1) - t9 * mrSges(6,3) - t11 * mrSges(5,3) + t78 * t95 + t115 * pkin(4) + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t94) * t53 + ((Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t53 + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1 + pkin(4) * mrSges(6,3)) * t54) * qJD(4)) * qJD(1); t112 + (-t53 * t33 + t54 * t35 - t118 + t62) * qJD(1);];
tauc = t31(:);
