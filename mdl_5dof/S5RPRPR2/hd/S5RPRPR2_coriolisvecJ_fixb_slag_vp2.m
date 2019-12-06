% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:20
% EndTime: 2019-12-05 17:49:23
% DurationCPUTime: 0.81s
% Computational Cost: add. (1501->151), mult. (3033->218), div. (0->0), fcn. (1825->8), ass. (0->82)
t78 = cos(pkin(8)) * pkin(1) + pkin(2);
t73 = t78 * qJD(1);
t112 = pkin(1) * sin(pkin(8));
t102 = qJD(1) * t112;
t92 = sin(qJ(3));
t75 = t92 * t102;
t94 = cos(qJ(3));
t55 = t94 * t73 - t75;
t98 = qJD(4) - t55;
t87 = sin(pkin(9));
t89 = cos(pkin(9));
t105 = t87 ^ 2 + t89 ^ 2;
t104 = qJD(3) * t94;
t49 = -qJD(3) * t75 + t73 * t104;
t86 = qJD(1) + qJD(3);
t39 = t86 * qJD(4) + t49;
t101 = t105 * t39;
t91 = sin(qJ(5));
t93 = cos(qJ(5));
t68 = -t91 * t87 + t93 * t89;
t65 = t68 * qJD(5);
t123 = t65 / 0.2e1;
t69 = t93 * t87 + t91 * t89;
t66 = t69 * qJD(5);
t122 = -t66 / 0.2e1;
t52 = t86 * t65;
t121 = t68 * t52;
t53 = t86 * t66;
t120 = t69 * t53;
t71 = (-pkin(7) - qJ(4)) * t87;
t82 = t89 * pkin(7);
t72 = t89 * qJ(4) + t82;
t43 = t91 * t71 + t93 * t72;
t119 = -t43 * qJD(5) - t98 * t69;
t42 = t93 * t71 - t91 * t72;
t118 = t42 * qJD(5) + t98 * t68;
t116 = t105 * mrSges(5,3);
t59 = t68 * t86;
t60 = t69 * t86;
t97 = -t89 * mrSges(5,1) + t87 * mrSges(5,2);
t117 = t59 * mrSges(6,1) - t60 * mrSges(6,2) - t97 * t86;
t106 = t94 * t112 + t92 * t78;
t56 = t94 * t102 + t92 * t73;
t41 = t86 * qJ(4) + t56;
t81 = t89 * qJD(2);
t28 = t81 + (-pkin(7) * t86 - t41) * t87;
t33 = t87 * qJD(2) + t89 * t41;
t29 = t86 * t82 + t33;
t8 = t91 * t28 + t93 * t29;
t3 = -t8 * qJD(5) - t69 * t39;
t7 = t93 * t28 - t91 * t29;
t115 = -t3 * t69 - t7 * t65;
t113 = t60 / 0.2e1;
t109 = t89 * pkin(4);
t107 = Ifges(6,4) * t60;
t76 = t92 * t112;
t79 = -pkin(3) - t109;
t23 = t53 * mrSges(6,1) + t52 * mrSges(6,2);
t100 = t94 * t78 - t76;
t99 = -pkin(3) - t100;
t96 = -(-t87 * t41 + t81) * t87 + t33 * t89;
t63 = qJ(4) + t106;
t44 = (-pkin(7) - t63) * t87;
t45 = t89 * t63 + t82;
t14 = t93 * t44 - t91 * t45;
t15 = t91 * t44 + t93 * t45;
t61 = -qJD(3) * t76 + t78 * t104;
t2 = t7 * qJD(5) + t68 * t39;
t24 = Ifges(6,2) * t59 + Ifges(6,6) * qJD(5) + t107;
t57 = Ifges(6,4) * t59;
t25 = Ifges(6,1) * t60 + Ifges(6,5) * qJD(5) + t57;
t34 = t79 * t86 + t98;
t50 = t56 * qJD(3);
t95 = -t49 * mrSges(4,2) + t25 * t123 + t24 * t122 + qJD(5) * (Ifges(6,5) * t65 - Ifges(6,6) * t66) / 0.2e1 + t34 * (t66 * mrSges(6,1) + t65 * mrSges(6,2)) + mrSges(5,3) * t101 + (t59 * t122 - t68 * t53) * Ifges(6,2) + (t65 * t113 + t69 * t52) * Ifges(6,1) + (t2 * t68 - t8 * t66) * mrSges(6,3) + (-t68 * mrSges(6,1) + t69 * mrSges(6,2) - mrSges(4,1) + t97) * t50 + (-t66 * t113 + t59 * t123 - t120 + t121) * Ifges(6,4);
t58 = qJD(4) + t61;
t54 = t99 - t109;
t48 = qJD(5) * mrSges(6,1) - t60 * mrSges(6,3);
t47 = -qJD(5) * mrSges(6,2) + t59 * mrSges(6,3);
t40 = -t86 * pkin(3) + t98;
t5 = -t15 * qJD(5) - t69 * t58;
t4 = t14 * qJD(5) + t68 * t58;
t1 = [m(5) * (t63 * t101 + t50 * t99 + t96 * t58) + t95 + m(4) * (-t50 * t100 + t49 * t106 + t56 * t61) + m(6) * (t3 * t14 + t2 * t15 + t8 * t4 + t7 * t5 + t50 * t54) - t61 * t86 * mrSges(4,2) + t4 * t47 + t5 * t48 + t54 * t23 + t58 * t86 * t116 + (-m(4) * t55 + m(5) * t40 + m(6) * t34 - t86 * mrSges(4,1) - t117) * t106 * qJD(3) + (-t14 * t52 - t15 * t53 + t115) * mrSges(6,3); t65 * t47 - t66 * t48 + m(6) * (t2 * t69 + t3 * t68 + t8 * t65 - t7 * t66) + (-t120 - t121) * mrSges(6,3); t95 + (t56 * mrSges(4,1) + t55 * mrSges(4,2) + t98 * t116) * t86 + t117 * t56 + t119 * t48 + t118 * t47 + (-t42 * t52 - t43 * t53 + t115) * mrSges(6,3) + t79 * t23 + (t118 * t8 + t119 * t7 + t2 * t43 + t3 * t42 - t34 * t56 + t50 * t79) * m(6) + (-t50 * pkin(3) + qJ(4) * t101 - t40 * t56 + t98 * t96) * m(5); -t86 ^ 2 * t116 - t59 * t47 + t60 * t48 + t23 + (-t8 * t59 + t7 * t60 + t50) * m(6) + (-t96 * t86 + t50) * m(5); Ifges(6,5) * t52 - Ifges(6,6) * t53 - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t34 * (t60 * mrSges(6,1) + t59 * mrSges(6,2)) - t60 * (Ifges(6,1) * t59 - t107) / 0.2e1 + t24 * t113 - qJD(5) * (Ifges(6,5) * t59 - Ifges(6,6) * t60) / 0.2e1 - t7 * t47 + t8 * t48 + (t7 * t59 + t8 * t60) * mrSges(6,3) - (-Ifges(6,2) * t60 + t25 + t57) * t59 / 0.2e1;];
tauc = t1(:);
