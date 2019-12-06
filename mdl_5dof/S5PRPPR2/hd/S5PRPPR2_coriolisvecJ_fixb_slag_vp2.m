% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:58
% EndTime: 2019-12-05 15:24:03
% DurationCPUTime: 1.07s
% Computational Cost: add. (871->151), mult. (2233->217), div. (0->0), fcn. (1583->8), ass. (0->82)
t64 = sin(pkin(9));
t66 = cos(pkin(9));
t86 = t64 ^ 2 + t66 ^ 2;
t97 = mrSges(5,3) * t86;
t65 = sin(pkin(8));
t69 = sin(qJ(2));
t84 = qJD(1) * t69;
t58 = t65 * t84;
t67 = cos(pkin(8));
t71 = cos(qJ(2));
t83 = qJD(1) * t71;
t81 = t67 * t83;
t42 = -t58 + t81;
t101 = qJD(4) - t42;
t68 = sin(qJ(5));
t70 = cos(qJ(5));
t73 = t64 * t68 - t66 * t70;
t41 = t73 * qJD(2);
t100 = t41 / 0.2e1;
t59 = pkin(2) * t65 + qJ(4);
t89 = pkin(6) + t59;
t47 = t89 * t64;
t48 = t89 * t66;
t19 = -t47 * t68 + t48 * t70;
t52 = t64 * t70 + t66 * t68;
t99 = -qJD(5) * t19 - t101 * t52;
t18 = -t47 * t70 - t48 * t68;
t98 = qJD(5) * t18 - t101 * t73;
t45 = t73 * qJD(5);
t36 = qJD(2) * t45;
t96 = t73 * t36;
t46 = t52 * qJD(5);
t37 = qJD(2) * t46;
t95 = t52 * t37;
t44 = t52 * qJD(2);
t75 = -mrSges(5,1) * t66 + mrSges(5,2) * t64;
t94 = mrSges(6,1) * t41 + mrSges(6,2) * t44 + t75 * qJD(2);
t91 = t44 / 0.2e1;
t90 = pkin(2) * t67;
t88 = Ifges(6,4) * t44;
t51 = t65 * t71 + t67 * t69;
t40 = t51 * qJD(2);
t34 = qJD(1) * t40;
t49 = t65 * t69 - t67 * t71;
t23 = t34 * t49;
t57 = qJD(2) * pkin(2) + t83;
t33 = t57 * t65 + t67 * t84;
t27 = qJD(2) * qJ(4) + t33;
t22 = qJD(3) * t64 + t27 * t66;
t85 = pkin(6) * qJD(2);
t82 = t42 * qJD(2);
t79 = -pkin(4) * t66 - pkin(3);
t56 = qJD(2) * t81;
t28 = t56 + (qJD(4) - t58) * qJD(2);
t78 = t86 * t28;
t11 = t37 * mrSges(6,1) - mrSges(6,2) * t36;
t32 = t57 * t67 - t58;
t76 = qJD(4) - t32;
t61 = t66 * qJD(3);
t16 = t61 + (-t27 - t85) * t64;
t17 = t66 * t85 + t22;
t5 = t16 * t70 - t17 * t68;
t6 = t16 * t68 + t17 * t70;
t74 = -(-t27 * t64 + t61) * t64 + t22 * t66;
t72 = qJD(2) ^ 2;
t55 = t79 - t90;
t43 = t49 * qJD(2);
t38 = Ifges(6,4) * t41;
t35 = -qJD(2) * t58 + t56;
t30 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t44;
t29 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t41;
t26 = -qJD(2) * pkin(3) + t76;
t24 = qJD(2) * t79 + t76;
t15 = Ifges(6,1) * t44 + Ifges(6,5) * qJD(5) - t38;
t14 = -Ifges(6,2) * t41 + Ifges(6,6) * qJD(5) + t88;
t13 = t73 * t51;
t12 = t52 * t51;
t4 = t43 * t52 + t45 * t51;
t3 = t43 * t73 - t46 * t51;
t2 = -qJD(5) * t6 - t28 * t52;
t1 = qJD(5) * t5 - t28 * t73;
t7 = [t49 * t11 + t3 * t29 + t4 * t30 + (-mrSges(3,1) * t69 - mrSges(3,2) * t71) * t72 + t94 * t40 + (-t12 * t36 + t13 * t37) * mrSges(6,3) + m(4) * (-t32 * t40 - t33 * t43 + t35 * t51 + t23) + m(6) * (-t1 * t13 - t12 * t2 + t24 * t40 + t3 * t6 + t4 * t5 + t23) + m(5) * (t26 * t40 - t43 * t74 + t51 * t78 + t23) + (-t40 * mrSges(4,1) - (-mrSges(4,2) + t97) * t43) * qJD(2); t55 * t11 - t45 * t15 / 0.2e1 + t24 * (t46 * mrSges(6,1) - t45 * mrSges(6,2)) + qJD(5) * (-Ifges(6,5) * t45 - Ifges(6,6) * t46) / 0.2e1 - t46 * t14 / 0.2e1 + t99 * t30 + t98 * t29 + (t100 * t46 + t37 * t73) * Ifges(6,2) + (-t36 * t52 - t45 * t91) * Ifges(6,1) + (t82 - t35) * mrSges(4,2) + (t1 * t19 + t18 * t2 + t5 * t99 + t6 * t98) * m(6) + (t101 * t74 + t59 * t78) * m(5) + (m(6) * t55 + m(5) * (-pkin(3) - t90) + t75 + t73 * mrSges(6,1) + t52 * mrSges(6,2) - mrSges(4,1)) * t34 + ((-t34 * t67 + t35 * t65) * pkin(2) - t33 * t42) * m(4) + (m(4) * t32 - m(5) * t26 - m(6) * t24 + qJD(2) * mrSges(4,1) - t94) * t51 * qJD(1) + (-t1 * t73 + t18 * t36 - t19 * t37 - t2 * t52 + t45 * t5 - t46 * t6) * mrSges(6,3) + (t100 * t45 - t46 * t91 - t95 + t96) * Ifges(6,4) + (qJD(2) * qJD(4) + t28 - t82) * t97; m(6) * (t1 * t52 - t2 * t73 - t45 * t6 - t46 * t5) - t45 * t29 - t46 * t30 + (-t95 - t96) * mrSges(6,3); -t72 * t97 + t41 * t29 + t44 * t30 + t11 + (t41 * t6 + t44 * t5 + t34) * m(6) + (-qJD(2) * t74 + t34) * m(5); -Ifges(6,5) * t36 - Ifges(6,6) * t37 - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t24 * (mrSges(6,1) * t44 - mrSges(6,2) * t41) - t44 * (-Ifges(6,1) * t41 - t88) / 0.2e1 + t14 * t91 - qJD(5) * (-Ifges(6,5) * t41 - Ifges(6,6) * t44) / 0.2e1 - t5 * t29 + t6 * t30 + (-t41 * t5 + t44 * t6) * mrSges(6,3) + (-Ifges(6,2) * t44 + t15 - t38) * t100;];
tauc = t7(:);
