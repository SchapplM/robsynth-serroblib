% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR2
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:47
% EndTime: 2019-12-05 15:44:49
% DurationCPUTime: 0.71s
% Computational Cost: add. (1310->143), mult. (2923->219), div. (0->0), fcn. (2034->8), ass. (0->75)
t61 = sin(pkin(9));
t62 = cos(pkin(9));
t65 = sin(qJ(2));
t68 = cos(qJ(2));
t51 = t61 * t68 + t62 * t65;
t45 = t51 * qJD(1);
t50 = -t61 * t65 + t62 * t68;
t47 = t50 * qJD(1);
t64 = sin(qJ(4));
t67 = cos(qJ(4));
t106 = pkin(2) * t61;
t59 = pkin(2) * t62 + pkin(3);
t78 = -t64 * t106 + t59 * t67;
t95 = -t78 * qJD(4) - t45 * t64 + t47 * t67;
t93 = t67 * t106 + t64 * t59;
t94 = t93 * qJD(4) - t67 * t45 - t47 * t64;
t56 = qJD(2) * pkin(2) + qJD(1) * t68;
t90 = qJD(1) * t65;
t27 = t62 * t56 - t61 * t90;
t25 = qJD(2) * pkin(3) + t27;
t28 = t56 * t61 + t62 * t90;
t16 = t25 * t64 + t28 * t67;
t60 = qJD(2) + qJD(4);
t14 = pkin(7) * t60 + t16;
t63 = sin(qJ(5));
t66 = cos(qJ(5));
t11 = qJD(3) * t66 - t14 * t63;
t12 = qJD(3) * t63 + t14 * t66;
t81 = -t11 * t63 + t12 * t66;
t46 = t51 * qJD(2);
t40 = qJD(1) * t46;
t48 = t50 * qJD(2);
t41 = qJD(1) * t48;
t7 = t16 * qJD(4) + t67 * t40 + t64 * t41;
t80 = t50 * t67 - t51 * t64;
t105 = t80 * t7;
t15 = t25 * t67 - t28 * t64;
t6 = t15 * qJD(4) - t64 * t40 + t67 * t41;
t2 = qJD(5) * t11 + t66 * t6;
t104 = t2 * t66;
t89 = qJD(5) * t12;
t3 = -t63 * t6 - t89;
t103 = t3 * t63;
t102 = Ifges(6,4) * t63;
t100 = t11 * mrSges(6,3);
t97 = t60 * t63;
t96 = t60 * t66;
t92 = Ifges(6,5) * qJD(5);
t91 = Ifges(6,6) * qJD(5);
t88 = qJD(5) * t63;
t87 = qJD(5) * t66;
t85 = t87 / 0.2e1;
t82 = -mrSges(6,1) * t66 + mrSges(6,2) * t63;
t44 = t82 * t60;
t84 = t60 * mrSges(5,1) - t44;
t20 = t50 * t64 + t51 * t67;
t53 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t97;
t54 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t96;
t79 = -t63 * t53 + t66 * t54;
t77 = (Ifges(6,2) * t66 + t102) * t60;
t76 = (mrSges(6,1) * t63 + mrSges(6,2) * t66) * qJD(5);
t75 = t60 * mrSges(5,2) - t79;
t74 = -t103 + (-t11 * t66 - t12 * t63) * qJD(5);
t71 = t74 + t104;
t13 = -pkin(4) * t60 - t15;
t33 = t77 + t91;
t57 = Ifges(6,4) * t96;
t34 = Ifges(6,1) * t97 + t57 + t92;
t70 = -t6 * mrSges(5,2) + mrSges(6,3) * t104 + t13 * t76 + t34 * t85 + qJD(5) ^ 2 * (Ifges(6,5) * t66 - Ifges(6,6) * t63) / 0.2e1 - (t33 + t77) * t88 / 0.2e1 + (-mrSges(5,1) + t82) * t7 + ((Ifges(6,1) * t66 - t102) * t88 + (0.3e1 * Ifges(6,4) * t66 + (Ifges(6,1) - 0.2e1 * Ifges(6,2)) * t63) * t85) * t60;
t43 = pkin(7) + t93;
t42 = -pkin(4) - t78;
t37 = t60 * t76;
t9 = t20 * qJD(4) + t67 * t46 + t64 * t48;
t8 = t80 * qJD(4) - t46 * t64 + t48 * t67;
t1 = [-t80 * t37 - t84 * t9 + (-t53 * t66 - t54 * t63) * t20 * qJD(5) - t75 * t8 + m(5) * (-t15 * t9 + t16 * t8 + t20 * t6 - t105) + m(4) * (-t27 * t46 + t28 * t48 - t40 * t50 + t41 * t51) + m(6) * (t13 * t9 + t71 * t20 + t81 * t8 - t105) + (-mrSges(4,1) * t46 - mrSges(4,2) * t48 + (-mrSges(3,1) * t65 - mrSges(3,2) * t68) * qJD(2)) * qJD(2); (-qJD(5) * t43 * t54 + t95 * t53 + (-t3 - t89) * mrSges(6,3)) * t63 + (-t94 * mrSges(5,1) + t95 * mrSges(5,2)) * t60 + (-t95 * t54 + (-t43 * t53 - t100) * qJD(5)) * t66 + t94 * t44 + (qJD(2) * t47 - t41) * mrSges(4,2) + (qJD(2) * t45 - t40) * mrSges(4,1) + t42 * t37 + t70 + (-t94 * t15 - t95 * t16 + t6 * t93 - t7 * t78) * m(5) + ((-t40 * t62 + t41 * t61) * pkin(2) + t27 * t45 - t28 * t47) * m(4) + (t94 * t13 + t7 * t42 + t71 * t43 - t95 * t81) * m(6); m(6) * (t2 * t63 + t3 * t66) + (m(6) * t81 + (-t63 ^ 2 - t66 ^ 2) * t60 * mrSges(6,3) + t79) * qJD(5); t74 * mrSges(6,3) + t84 * t16 + t75 * t15 - pkin(4) * t37 + t70 + (-t53 * t87 - t54 * t88) * pkin(7) + (-t7 * pkin(4) + (-t11 * t87 - t12 * t88 - t103 + t104) * pkin(7) - t13 * t16 - t81 * t15) * m(6); t3 * mrSges(6,1) - t2 * mrSges(6,2) - t11 * t54 + t12 * t53 + ((t92 / 0.2e1 - t13 * mrSges(6,2) - t34 / 0.2e1 - t57 / 0.2e1 + t100) * t66 + (-t91 / 0.2e1 - t13 * mrSges(6,1) + t33 / 0.2e1 + t12 * mrSges(6,3) + (t102 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t66) * t60) * t63) * t60;];
tauc = t1(:);
