% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:38
% EndTime: 2019-12-05 15:42:44
% DurationCPUTime: 1.55s
% Computational Cost: add. (2002->215), mult. (5221->316), div. (0->0), fcn. (3872->6), ass. (0->98)
t87 = cos(pkin(9));
t91 = cos(qJ(4));
t105 = t87 * t91;
t86 = sin(pkin(9));
t89 = sin(qJ(4));
t71 = -t86 * t89 + t105;
t63 = t71 * qJD(2);
t72 = t86 * t91 + t87 * t89;
t64 = t72 * qJD(2);
t88 = sin(qJ(5));
t90 = cos(qJ(5));
t45 = t63 * t88 + t64 * t90;
t113 = Ifges(6,4) * t45;
t97 = t90 * t63 - t64 * t88;
t39 = Ifges(6,4) * t97;
t85 = qJD(4) + qJD(5);
t13 = Ifges(6,1) * t45 + Ifges(6,5) * t85 + t39;
t65 = t71 * qJD(4);
t58 = qJD(2) * t65;
t66 = t72 * qJD(4);
t59 = qJD(2) * t66;
t17 = qJD(5) * t97 + t58 * t90 - t59 * t88;
t18 = -qJD(5) * t45 - t58 * t88 - t59 * t90;
t101 = qJD(4) * t91;
t102 = qJD(3) * t86;
t104 = pkin(6) + qJ(3);
t76 = t104 * t86;
t82 = t87 * qJD(1);
t60 = -qJD(2) * t76 + t82;
t100 = qJ(3) * qJD(2);
t74 = t86 * qJD(1) + t87 * t100;
t61 = pkin(6) * qJD(2) * t87 + t74;
t79 = qJD(3) * t105;
t24 = t60 * t101 + qJD(2) * t79 + (-qJD(2) * t102 - qJD(4) * t61) * t89;
t19 = -pkin(7) * t59 + t24;
t36 = t60 * t89 + t61 * t91;
t93 = t72 * qJD(3);
t25 = -qJD(2) * t93 - qJD(4) * t36;
t20 = -pkin(7) * t58 + t25;
t30 = pkin(7) * t63 + t36;
t112 = t30 * t88;
t35 = t91 * t60 - t61 * t89;
t29 = -pkin(7) * t64 + t35;
t26 = qJD(4) * pkin(4) + t29;
t6 = t26 * t90 - t112;
t2 = qJD(5) * t6 + t19 * t90 + t20 * t88;
t111 = t30 * t90;
t7 = t26 * t88 + t111;
t3 = -qJD(5) * t7 - t19 * t88 + t20 * t90;
t99 = -pkin(3) * t87 - pkin(2);
t75 = qJD(2) * t99 + qJD(3);
t48 = -pkin(4) * t63 + t75;
t130 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t18 - (Ifges(6,5) * t97 - Ifges(6,6) * t45) * t85 / 0.2e1 + (t45 * t7 + t6 * t97) * mrSges(6,3) - (-Ifges(6,2) * t45 + t13 + t39) * t97 / 0.2e1 - t48 * (mrSges(6,1) * t45 + mrSges(6,2) * t97) - (Ifges(6,1) * t97 - t113) * t45 / 0.2e1;
t96 = (t86 ^ 2 + t87 ^ 2) * qJD(2);
t129 = mrSges(4,3) * t96;
t12 = Ifges(6,2) * t97 + Ifges(6,6) * t85 + t113;
t127 = t12 / 0.2e1;
t77 = t104 * t87;
t50 = -t89 * t76 + t91 * t77;
t121 = t97 / 0.2e1;
t119 = t45 / 0.2e1;
t117 = t65 / 0.2e1;
t116 = -t66 / 0.2e1;
t114 = Ifges(5,4) * t64;
t46 = t71 * t90 - t72 * t88;
t110 = t46 * t17;
t47 = t71 * t88 + t72 * t90;
t109 = t47 * t18;
t107 = t71 * t58;
t106 = t72 * t59;
t98 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t49 = -t91 * t76 - t77 * t89;
t37 = -pkin(7) * t72 + t49;
t38 = pkin(7) * t71 + t50;
t10 = t37 * t90 - t38 * t88;
t11 = t37 * t88 + t38 * t90;
t95 = -(-t100 * t86 + t82) * t86 + t74 * t87;
t33 = -t76 * t101 + t79 + (-qJD(4) * t77 - t102) * t89;
t34 = -qJD(4) * t50 - t93;
t62 = Ifges(5,4) * t63;
t56 = t58 * mrSges(5,2);
t55 = -pkin(4) * t71 + t99;
t53 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t64;
t52 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t63;
t41 = Ifges(5,1) * t64 + Ifges(5,5) * qJD(4) + t62;
t40 = Ifges(5,2) * t63 + Ifges(5,6) * qJD(4) + t114;
t32 = mrSges(6,1) * t85 - mrSges(6,3) * t45;
t31 = -mrSges(6,2) * t85 + mrSges(6,3) * t97;
t28 = -pkin(7) * t65 + t34;
t27 = -pkin(7) * t66 + t33;
t23 = -qJD(5) * t47 - t65 * t88 - t66 * t90;
t22 = qJD(5) * t46 + t65 * t90 - t66 * t88;
t21 = -mrSges(6,1) * t97 + mrSges(6,2) * t45;
t9 = t29 * t90 - t112;
t8 = -t29 * t88 - t111;
t5 = -qJD(5) * t11 - t27 * t88 + t28 * t90;
t4 = qJD(5) * t10 + t27 * t90 + t28 * t88;
t1 = [t22 * t31 + t23 * t32 + t65 * t52 - t66 * t53 + (t109 - t110) * mrSges(6,3) + (-t106 - t107) * mrSges(5,3) + m(5) * (t24 * t72 + t25 * t71 - t35 * t66 + t36 * t65) + m(6) * (t2 * t47 + t22 * t7 + t23 * t6 + t3 * t46); (-t10 * t17 + t11 * t18 + t2 * t46 - t22 * t6 + t23 * t7 - t3 * t47) * mrSges(6,3) + m(5) * (t24 * t50 + t25 * t49 + t33 * t36 + t34 * t35) + t75 * (mrSges(5,1) * t66 + mrSges(5,2) * t65) + qJD(4) * (Ifges(5,5) * t65 - Ifges(5,6) * t66) / 0.2e1 + t85 * (Ifges(6,5) * t22 + Ifges(6,6) * t23) / 0.2e1 + t48 * (-mrSges(6,1) * t23 + mrSges(6,2) * t22) + t33 * t52 + t34 * t53 + t4 * t31 + t5 * t32 + t22 * t13 / 0.2e1 + (m(4) * (qJ(3) * t96 + t95) + 0.2e1 * t129) * qJD(3) + (t119 * t22 + t17 * t47) * Ifges(6,1) + (t121 * t23 + t46 * t18) * Ifges(6,2) + (t119 * t23 + t121 * t22 + t109 + t110) * Ifges(6,4) + (t117 * t64 + t58 * t72) * Ifges(5,1) + t55 * t98 + t40 * t116 + t41 * t117 + (t24 * t71 - t25 * t72 - t35 * t65 - t36 * t66 - t49 * t58 - t50 * t59) * mrSges(5,3) + m(6) * (t10 * t3 + t11 * t2 + t4 * t7 + t5 * t6 + (t48 * t66 + t55 * t59) * pkin(4)) + (t59 * (-mrSges(6,1) * t46 + mrSges(6,2) * t47) + t66 * t21) * pkin(4) + (t116 * t63 - t71 * t59) * Ifges(5,2) + t99 * (t59 * mrSges(5,1) + t56) + (t116 * t64 + t117 * t63 - t106 + t107) * Ifges(5,4) + t23 * t127; -t97 * t31 + t45 * t32 - t63 * t52 + t64 * t53 + t56 - (-m(6) * pkin(4) - mrSges(5,1)) * t59 - m(5) * (-t35 * t64 + t36 * t63) - m(6) * (-t45 * t6 + t7 * t97) + t98 + (-m(4) * t95 - t129) * qJD(2); (-t64 * t21 + (t90 * t31 - t88 * t32) * qJD(5) + (-t17 * t90 + t18 * t88) * mrSges(6,3) + (t2 * t88 + t3 * t90 - t48 * t64 + (-t6 * t88 + t7 * t90) * qJD(5)) * m(6)) * pkin(4) + (t35 * t63 + t36 * t64) * mrSges(5,3) - t75 * (mrSges(5,1) * t64 + mrSges(5,2) * t63) - qJD(4) * (Ifges(5,5) * t63 - Ifges(5,6) * t64) / 0.2e1 + t64 * t40 / 0.2e1 + Ifges(5,5) * t58 - Ifges(5,6) * t59 - t35 * t52 + t36 * t53 - t9 * t31 - t8 * t32 - t24 * mrSges(5,2) + t25 * mrSges(5,1) - t64 * (Ifges(5,1) * t63 - t114) / 0.2e1 - m(6) * (t6 * t8 + t7 * t9) - (-Ifges(5,2) * t64 + t41 + t62) * t63 / 0.2e1 + t45 * t127 + t130; t12 * t119 - t6 * t31 + t7 * t32 + t130;];
tauc = t1(:);
