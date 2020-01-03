% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:40
% EndTime: 2019-12-31 16:34:43
% DurationCPUTime: 1.05s
% Computational Cost: add. (838->167), mult. (2248->263), div. (0->0), fcn. (1375->6), ass. (0->98)
t119 = -Ifges(4,1) / 0.2e1;
t73 = cos(qJ(3));
t98 = qJD(2) * t73;
t118 = -Ifges(4,4) * t98 / 0.2e1;
t74 = cos(qJ(2));
t100 = qJD(1) * t74;
t111 = -pkin(6) - pkin(5);
t70 = sin(qJ(3));
t59 = t111 * t70;
t60 = t111 * t73;
t69 = sin(qJ(4));
t72 = cos(qJ(4));
t31 = t59 * t72 + t60 * t69;
t89 = qJD(3) * t111;
t55 = t70 * t89;
t56 = t73 * t89;
t78 = t69 * t70 - t72 * t73;
t117 = t31 * qJD(4) + t78 * t100 + t55 * t72 + t56 * t69;
t32 = t59 * t69 - t60 * t72;
t53 = t69 * t73 + t70 * t72;
t77 = t53 * t74;
t116 = qJD(1) * t77 - t32 * qJD(4) - t55 * t69 + t56 * t72;
t71 = sin(qJ(2));
t43 = t78 * t71;
t47 = t78 * qJD(2);
t48 = t53 * qJD(2);
t26 = mrSges(5,1) * t47 + mrSges(5,2) * t48;
t64 = -pkin(3) * t73 - pkin(2);
t50 = t64 * qJD(2) - t100;
t90 = m(5) * t50 + t26;
t66 = qJD(3) + qJD(4);
t102 = Ifges(4,6) * qJD(3);
t108 = Ifges(4,4) * t70;
t63 = -qJD(2) * pkin(2) - t100;
t115 = -t90 * pkin(3) + t102 / 0.2e1 + (Ifges(4,2) * t73 + t108) * qJD(2) / 0.2e1 - t63 * mrSges(4,1);
t67 = t70 ^ 2;
t68 = t73 ^ 2;
t113 = -t47 / 0.2e1;
t112 = t48 / 0.2e1;
t109 = mrSges(5,3) * t47;
t107 = Ifges(5,4) * t48;
t101 = qJD(1) * t71;
t62 = qJD(2) * pkin(5) + t101;
t84 = pkin(6) * qJD(2) + t62;
t41 = t84 * t73;
t106 = t41 * t69;
t105 = t41 * t72;
t104 = t48 * mrSges(5,3);
t103 = Ifges(4,5) * qJD(3);
t99 = qJD(2) * t70;
t97 = qJD(3) * t70;
t96 = qJD(3) * t73;
t95 = qJD(4) * t69;
t94 = qJD(4) * t72;
t93 = qJD(1) * qJD(2);
t92 = qJD(2) * qJD(3);
t91 = t62 * t96;
t88 = (t67 + t68) * t62;
t87 = t74 * t93;
t86 = t103 / 0.2e1;
t85 = -t102 / 0.2e1;
t82 = t99 * t119 - t103 / 0.2e1 + t118 - t63 * mrSges(4,2);
t61 = t73 * t87;
t35 = -t62 * t97 + t61;
t36 = -t70 * t87 - t91;
t81 = t35 * t73 - t36 * t70;
t40 = t84 * t70;
t37 = qJD(3) * pkin(3) - t40;
t14 = t37 * t72 - t106;
t15 = t37 * t69 + t105;
t57 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t99;
t58 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t98;
t80 = t73 * t57 + t70 * t58;
t79 = t70 * t57 - t73 * t58;
t19 = -Ifges(5,2) * t47 + Ifges(5,6) * t66 + t107;
t29 = -qJD(3) * t40 + t61;
t30 = -t91 + (-pkin(6) * t96 - t70 * t100) * qJD(2);
t2 = t14 * qJD(4) + t29 * t72 + t30 * t69;
t44 = Ifges(5,4) * t47;
t20 = Ifges(5,1) * t48 + Ifges(5,5) * t66 - t44;
t27 = t66 * t78;
t23 = t27 * qJD(2);
t28 = t66 * t53;
t24 = t28 * qJD(2);
t3 = -t15 * qJD(4) - t29 * t69 + t30 * t72;
t76 = t3 * mrSges(5,1) - t2 * mrSges(5,2) - t14 * t109 + t19 * t112 - t50 * (mrSges(5,1) * t48 - mrSges(5,2) * t47) - t48 * (-Ifges(5,1) * t47 - t107) / 0.2e1 - t66 * (-Ifges(5,5) * t47 - Ifges(5,6) * t48) / 0.2e1 - Ifges(5,6) * t24 - Ifges(5,5) * t23 + (-Ifges(5,2) * t48 + t20 - t44) * t47 / 0.2e1;
t75 = qJD(2) ^ 2;
t51 = (pkin(3) * t97 + t101) * qJD(2);
t49 = (mrSges(4,1) * t70 + mrSges(4,2) * t73) * t92;
t42 = t53 * t71;
t34 = mrSges(5,1) * t66 - t104;
t33 = -mrSges(5,2) * t66 - t109;
t17 = -t40 * t72 - t106;
t16 = t40 * t69 - t105;
t7 = -qJD(2) * t77 + t66 * t43;
t6 = -t28 * t71 - t74 * t47;
t4 = mrSges(5,1) * t24 - mrSges(5,2) * t23;
t1 = [t6 * t33 + t7 * t34 + m(5) * (t14 * t7 + t15 * t6 - t2 * t43 - t3 * t42) + (-t23 * t42 + t24 * t43) * mrSges(5,3) + (-t49 - t4 - m(5) * t51 - t75 * mrSges(3,2) + (m(4) * t88 - t79) * qJD(2)) * t74 + (-t75 * mrSges(3,1) - t80 * qJD(3) + m(4) * (t81 - t87) + (m(4) * t63 + qJD(2) * (-mrSges(4,1) * t73 + mrSges(4,2) * t70) + t90) * qJD(2)) * t71; t66 * (-Ifges(5,5) * t27 - Ifges(5,6) * t28) / 0.2e1 + t64 * t4 - pkin(2) * t49 + t50 * (t28 * mrSges(5,1) - t27 * mrSges(5,2)) + t51 * (mrSges(5,1) * t78 + t53 * mrSges(5,2)) - t27 * t20 / 0.2e1 - t28 * t19 / 0.2e1 + t116 * t34 + t117 * t33 + (-t28 * t113 + t24 * t78) * Ifges(5,2) + (-t27 * t112 - t23 * t53) * Ifges(5,1) + t81 * mrSges(4,3) + (-t71 * t26 + t79 * t74) * qJD(1) + (t14 * t27 - t15 * t28 - t2 * t78 + t31 * t23 - t32 * t24 - t3 * t53) * mrSges(5,3) + ((-pkin(5) * t57 - t82 + t86) * t73 + (-pkin(5) * t58 + t85 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t98 - t115) * t70) * qJD(3) + (-t28 * t112 - t27 * t113 + t23 * t78 - t24 * t53) * Ifges(5,4) + (-0.3e1 / 0.2e1 * t67 + 0.3e1 / 0.2e1 * t68) * Ifges(4,4) * t92 + (-t50 * t101 + t116 * t14 + t117 * t15 + t2 * t32 + t3 * t31 + t51 * t64) * m(5) + (-pkin(2) * t71 * t93 + t81 * pkin(5) - (t63 * t71 + t74 * t88) * qJD(1)) * m(4); (m(5) * (-t14 * t95 + t15 * t94 + t2 * t69 + t3 * t72) + t33 * t94 - t34 * t95 + (t23 * t72 - t24 * t69) * mrSges(5,3)) * pkin(3) + ((t86 + t118 + t82) * t73 + (t85 + (t108 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t119) * t73) * qJD(2) + t115) * t70) * qJD(2) - m(5) * (t14 * t16 + t15 * t17) + t80 * t62 + t76 + t15 * t104 - t17 * t33 - t16 * t34 - t35 * mrSges(4,2) + t36 * mrSges(4,1); -t14 * t33 + (t34 + t104) * t15 + t76;];
tauc = t1(:);
