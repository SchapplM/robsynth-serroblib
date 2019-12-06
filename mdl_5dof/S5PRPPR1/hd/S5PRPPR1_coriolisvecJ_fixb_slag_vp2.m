% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR1
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:47
% EndTime: 2019-12-05 15:21:51
% DurationCPUTime: 1.01s
% Computational Cost: add. (913->162), mult. (2489->260), div. (0->0), fcn. (1716->6), ass. (0->87)
t65 = sin(pkin(9));
t69 = sin(qJ(5));
t67 = cos(pkin(9));
t70 = cos(qJ(5));
t92 = t67 * t70;
t74 = t65 * t69 - t92;
t44 = t74 * qJD(5);
t68 = cos(pkin(8));
t87 = qJD(2) * t68;
t102 = t74 * t87 - t44;
t53 = t65 * t70 + t67 * t69;
t45 = t53 * qJD(5);
t72 = qJD(2) * t53;
t101 = t68 * t72 - t45;
t66 = sin(pkin(8));
t63 = t66 ^ 2;
t64 = t68 ^ 2;
t100 = (t63 + t64) * mrSges(4,3) * qJD(2);
t99 = -m(6) / 0.2e1;
t84 = t66 * qJD(2);
t79 = t65 * t84;
t30 = -t69 * t79 + t84 * t92;
t97 = t30 / 0.2e1;
t96 = -t68 / 0.2e1;
t95 = Ifges(6,4) * t30;
t94 = t65 * t66;
t93 = t66 * t67;
t28 = t66 * t72;
t73 = (mrSges(5,1) * t65 + mrSges(5,2) * t67) * t66;
t91 = -mrSges(6,1) * t28 - mrSges(6,2) * t30 - qJD(2) * t73;
t57 = -pkin(3) * t68 - qJ(4) * t66 - pkin(2);
t41 = t57 * qJD(2) + qJD(3);
t82 = qJ(3) * qJD(2);
t56 = qJD(1) * t66 + t68 * t82;
t16 = t65 * t41 + t67 * t56;
t88 = qJ(3) * t68;
t90 = t65 * t57 + t67 * t88;
t86 = qJD(3) * t68;
t85 = qJD(4) * t66;
t83 = t66 * qJD(3);
t81 = qJD(2) * qJD(3);
t80 = pkin(6) * t93;
t78 = t68 * t81;
t32 = t66 * t45;
t24 = qJD(2) * t32;
t33 = t66 * t44;
t25 = qJD(2) * t33;
t11 = -t25 * mrSges(6,1) - t24 * mrSges(6,2);
t15 = t67 * t41 - t56 * t65;
t55 = qJD(1) * t68 - t66 * t82;
t42 = -t65 * t86 - t67 * t85;
t36 = t42 * qJD(2);
t43 = -t65 * t85 + t67 * t86;
t37 = t43 * qJD(2);
t12 = (-pkin(4) * t68 - t80) * qJD(2) + t15;
t14 = -pkin(6) * t79 + t16;
t5 = t12 * t70 - t14 * t69;
t1 = t5 * qJD(5) + t36 * t69 + t37 * t70;
t6 = t12 * t69 + t14 * t70;
t2 = -t6 * qJD(5) + t36 * t70 - t37 * t69;
t77 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t50 = t67 * t57;
t17 = -t80 + t50 + (-qJ(3) * t65 - pkin(4)) * t68;
t20 = -pkin(6) * t94 + t90;
t7 = t17 * t70 - t20 * t69;
t8 = t17 * t69 + t20 * t70;
t76 = t36 * t67 + t37 * t65;
t75 = -t55 * t66 + t56 * t68;
t51 = qJD(4) - t55;
t61 = qJD(5) - t87;
t60 = t63 * qJ(3) * t81;
t54 = (pkin(4) * t65 + qJ(3)) * t66;
t47 = (-t68 * mrSges(5,1) - mrSges(5,3) * t93) * qJD(2);
t46 = (t68 * mrSges(5,2) - mrSges(5,3) * t94) * qJD(2);
t39 = t74 * t66;
t38 = t53 * t66;
t27 = pkin(4) * t79 + t51;
t26 = Ifges(6,4) * t28;
t23 = Ifges(6,5) * t24;
t22 = Ifges(6,6) * t25;
t19 = mrSges(6,1) * t61 - mrSges(6,3) * t30;
t18 = -mrSges(6,2) * t61 - mrSges(6,3) * t28;
t10 = Ifges(6,1) * t30 + Ifges(6,5) * t61 - t26;
t9 = -Ifges(6,2) * t28 + Ifges(6,6) * t61 + t95;
t4 = -t8 * qJD(5) + t42 * t70 - t43 * t69;
t3 = t7 * qJD(5) + t42 * t69 + t43 * t70;
t13 = [-t68 * t11 - t32 * t18 + t33 * t19 + (-t38 * t24 - t39 * t25) * mrSges(6,3) + m(6) * (-t1 * t39 - t2 * t38 - t32 * t6 + t33 * t5) + 0.2e1 * (m(5) * (-t36 * t65 + t37 * t67 - t78) / 0.2e1 + t78 * t99) * t66; t54 * t11 + t61 * (-Ifges(6,5) * t32 + Ifges(6,6) * t33) / 0.2e1 + t43 * t46 + t42 * t47 - t32 * t10 / 0.2e1 + (-Ifges(6,1) * t32 + Ifges(6,4) * t33) * t97 + t27 * (-t33 * mrSges(6,1) - t32 * mrSges(6,2)) - t28 * (-Ifges(6,4) * t32 + Ifges(6,2) * t33) / 0.2e1 + t33 * t9 / 0.2e1 + t3 * t18 + t4 * t19 - t76 * t66 * mrSges(5,3) + (t23 / 0.2e1 - t22 / 0.2e1 + t37 * mrSges(5,2) - t36 * mrSges(5,1) + t77) * t68 + (-t1 * t38 + t2 * t39 + t5 * t32 + t6 * t33) * mrSges(6,3) + (0.2e1 * t100 + ((mrSges(6,1) * t38 - mrSges(6,2) * t39 + t73) * qJD(2) - t91) * t66) * qJD(3) + m(5) * (t37 * t90 + t16 * t43 + t36 * (-t65 * t88 + t50) + t15 * t42 + t60 + t51 * t83) + m(4) * (t60 + (t64 * t82 + t75) * qJD(3)) + m(6) * (t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5 + (qJD(2) * t54 + t27) * t83) + (t8 * mrSges(6,3) - Ifges(6,4) * t39 - Ifges(6,2) * t38 + Ifges(6,6) * t96) * t25 - (-t7 * mrSges(6,3) - Ifges(6,1) * t39 - Ifges(6,4) * t38 + Ifges(6,5) * t96) * t24; t101 * t19 + t102 * t18 + (-t24 * t74 + t25 * t53) * mrSges(6,3) + m(5) * t76 + (t1 * t53 + t101 * t5 + t102 * t6 - t2 * t74) * m(6) + (-t100 - m(4) * t75 + (-t46 * t67 + t47 * t65 - m(5) * (-t15 * t65 + t16 * t67)) * t68 + (-m(5) * t51 + 0.2e1 * t27 * t99 + t91) * t66) * qJD(2); -m(6) * (-t28 * t6 - t30 * t5) + t28 * t18 + t30 * t19 + (m(6) * qJD(3) + t65 * t46 + t67 * t47 + (t15 * t67 + t16 * t65 + qJD(3)) * m(5)) * t84 + t11; -t23 + t22 - t27 * (mrSges(6,1) * t30 - mrSges(6,2) * t28) - t30 * (-Ifges(6,1) * t28 - t95) / 0.2e1 + t9 * t97 - t61 * (-Ifges(6,5) * t28 - Ifges(6,6) * t30) / 0.2e1 - t5 * t18 + t6 * t19 + (-t28 * t5 + t30 * t6) * mrSges(6,3) - t77 + (-Ifges(6,2) * t30 + t10 - t26) * t28 / 0.2e1;];
tauc = t13(:);
