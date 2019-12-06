% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:46
% EndTime: 2019-12-05 17:28:48
% DurationCPUTime: 1.02s
% Computational Cost: add. (1096->164), mult. (2814->263), div. (0->0), fcn. (1899->8), ass. (0->89)
t67 = sin(pkin(9));
t73 = sin(qJ(5));
t70 = cos(pkin(9));
t74 = cos(qJ(5));
t94 = t70 * t74;
t78 = t67 * t73 - t94;
t51 = t78 * qJD(5);
t71 = cos(pkin(8));
t90 = qJD(1) * t71;
t106 = t78 * t90 - t51;
t58 = t67 * t74 + t70 * t73;
t52 = t58 * qJD(5);
t76 = qJD(1) * t58;
t105 = t71 * t76 - t52;
t68 = sin(pkin(8));
t65 = t68 ^ 2;
t66 = t71 ^ 2;
t104 = (t65 + t66) * mrSges(4,3) * qJD(1);
t103 = -m(6) / 0.2e1;
t86 = t68 * qJD(1);
t83 = t67 * t86;
t34 = -t73 * t83 + t86 * t94;
t101 = t34 / 0.2e1;
t100 = -t71 / 0.2e1;
t99 = Ifges(6,4) * t34;
t98 = t67 * t68;
t97 = t67 * t71;
t96 = t68 * t70;
t95 = t70 * t71;
t32 = t68 * t76;
t77 = (mrSges(5,1) * t67 + mrSges(5,2) * t70) * t68;
t93 = -mrSges(6,1) * t32 - mrSges(6,2) * t34 - qJD(1) * t77;
t55 = -cos(pkin(7)) * pkin(1) - pkin(3) * t71 - qJ(4) * t68 - pkin(2);
t30 = t55 * qJD(1) + qJD(3);
t64 = sin(pkin(7)) * pkin(1) + qJ(3);
t62 = qJD(1) * t64;
t42 = qJD(2) * t68 + t62 * t71;
t15 = t67 * t30 + t70 * t42;
t92 = t67 * t55 + t64 * t95;
t89 = qJD(3) * t68;
t88 = qJD(3) * t71;
t87 = qJD(4) * t68;
t85 = qJD(1) * qJD(3);
t84 = pkin(6) * t96;
t82 = t71 * t85;
t36 = t68 * t52;
t27 = qJD(1) * t36;
t37 = t68 * t51;
t28 = qJD(1) * t37;
t13 = -t28 * mrSges(6,1) - t27 * mrSges(6,2);
t14 = t70 * t30 - t42 * t67;
t41 = qJD(2) * t71 - t68 * t62;
t49 = -t67 * t88 - t70 * t87;
t39 = t49 * qJD(1);
t50 = -t67 * t87 + t70 * t88;
t40 = t50 * qJD(1);
t10 = -pkin(6) * t83 + t15;
t9 = (-pkin(4) * t71 - t84) * qJD(1) + t14;
t5 = -t10 * t73 + t74 * t9;
t1 = t5 * qJD(5) + t39 * t73 + t40 * t74;
t6 = t10 * t74 + t73 * t9;
t2 = -t6 * qJD(5) + t39 * t74 - t40 * t73;
t81 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t46 = t70 * t55;
t17 = -t84 + t46 + (-t64 * t67 - pkin(4)) * t71;
t18 = -pkin(6) * t98 + t92;
t7 = t17 * t74 - t18 * t73;
t8 = t17 * t73 + t18 * t74;
t80 = t39 * t70 + t40 * t67;
t79 = -t41 * t68 + t42 * t71;
t38 = qJD(4) - t41;
t63 = qJD(5) - t90;
t60 = t65 * t64 * t85;
t54 = (-t71 * mrSges(5,1) - mrSges(5,3) * t96) * qJD(1);
t53 = (t71 * mrSges(5,2) - mrSges(5,3) * t98) * qJD(1);
t48 = (pkin(4) * t67 + t64) * t68;
t44 = t78 * t68;
t43 = t58 * t68;
t29 = Ifges(6,4) * t32;
t26 = Ifges(6,5) * t27;
t25 = Ifges(6,6) * t28;
t21 = pkin(4) * t83 + t38;
t20 = mrSges(6,1) * t63 - mrSges(6,3) * t34;
t19 = -mrSges(6,2) * t63 - mrSges(6,3) * t32;
t12 = Ifges(6,1) * t34 + Ifges(6,5) * t63 - t29;
t11 = -Ifges(6,2) * t32 + Ifges(6,6) * t63 + t99;
t4 = -t8 * qJD(5) + t49 * t74 - t50 * t73;
t3 = t7 * qJD(5) + t49 * t73 + t50 * t74;
t16 = [t63 * (-Ifges(6,5) * t36 + Ifges(6,6) * t37) / 0.2e1 + t48 * t13 + t50 * t53 + t49 * t54 - t36 * t12 / 0.2e1 + (-Ifges(6,1) * t36 + Ifges(6,4) * t37) * t101 - t32 * (-Ifges(6,4) * t36 + Ifges(6,2) * t37) / 0.2e1 + t21 * (-t37 * mrSges(6,1) - t36 * mrSges(6,2)) + t37 * t11 / 0.2e1 + t3 * t19 + t4 * t20 - t80 * t68 * mrSges(5,3) + (t26 / 0.2e1 - t25 / 0.2e1 + t40 * mrSges(5,2) - t39 * mrSges(5,1) + t81) * t71 + (-t1 * t43 + t2 * t44 + t5 * t36 + t6 * t37) * mrSges(6,3) + (0.2e1 * t104 + ((mrSges(6,1) * t43 - mrSges(6,2) * t44 + t77) * qJD(1) - t93) * t68) * qJD(3) + m(4) * (t60 + (t66 * t62 + t79) * qJD(3)) + m(5) * (t40 * t92 + t15 * t50 + t39 * (-t64 * t97 + t46) + t14 * t49 + t60 + t38 * t89) + m(6) * (t1 * t8 + t2 * t7 + t6 * t3 + t5 * t4 + (qJD(1) * t48 + t21) * t89) + (t8 * mrSges(6,3) - Ifges(6,4) * t44 - Ifges(6,2) * t43 + Ifges(6,6) * t100) * t28 - (-t7 * mrSges(6,3) - Ifges(6,1) * t44 - Ifges(6,4) * t43 + Ifges(6,5) * t100) * t27; -t71 * t13 - t36 * t19 + t37 * t20 + (-t43 * t27 - t44 * t28) * mrSges(6,3) + m(6) * (-t1 * t44 - t2 * t43 - t36 * t6 + t37 * t5) + 0.2e1 * (m(5) * (-t39 * t67 + t40 * t70 - t82) / 0.2e1 + t82 * t103) * t68; t105 * t20 + t106 * t19 + (-t27 * t78 + t28 * t58) * mrSges(6,3) + m(5) * t80 + (t1 * t58 + t105 * t5 + t106 * t6 - t2 * t78) * m(6) + ((-t53 * t70 + t54 * t67) * t71 - m(4) * t79 - m(5) * (-t14 * t97 + t15 * t95) - t104 + (-m(5) * t38 + 0.2e1 * t21 * t103 + t93) * t68) * qJD(1); -m(6) * (-t32 * t6 - t34 * t5) + t32 * t19 + t34 * t20 + (m(6) * qJD(3) + t67 * t53 + t70 * t54 + (t14 * t70 + t15 * t67 + qJD(3)) * m(5)) * t86 + t13; -t26 + t25 - t21 * (mrSges(6,1) * t34 - mrSges(6,2) * t32) - t34 * (-Ifges(6,1) * t32 - t99) / 0.2e1 + t11 * t101 - t63 * (-Ifges(6,5) * t32 - Ifges(6,6) * t34) / 0.2e1 - t5 * t19 + t6 * t20 + (-t32 * t5 + t34 * t6) * mrSges(6,3) - t81 + (-Ifges(6,2) * t34 + t12 - t29) * t32 / 0.2e1;];
tauc = t16(:);
