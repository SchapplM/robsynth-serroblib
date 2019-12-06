% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:37
% EndTime: 2019-12-05 15:12:39
% DurationCPUTime: 0.61s
% Computational Cost: add. (1072->128), mult. (2635->196), div. (0->0), fcn. (1962->8), ass. (0->70)
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t100 = t58 * mrSges(6,1) - t55 * mrSges(6,2) + mrSges(5,1);
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t57 = sin(qJ(3));
t60 = cos(qJ(3));
t44 = -t57 * t53 + t60 * t54;
t52 = qJD(3) + qJD(4);
t67 = (mrSges(6,1) * t55 + mrSges(6,2) * t58) * qJD(5);
t32 = t52 * t67;
t39 = t44 * qJD(1);
t33 = qJD(3) * pkin(3) + t39;
t56 = sin(qJ(4));
t45 = t60 * t53 + t57 * t54;
t40 = t45 * qJD(1);
t59 = cos(qJ(4));
t84 = t59 * t40;
t16 = t56 * t33 + t84;
t41 = t44 * qJD(3);
t36 = qJD(1) * t41;
t42 = t45 * qJD(3);
t37 = qJD(1) * t42;
t7 = t16 * qJD(4) + t56 * t36 + t59 * t37;
t99 = m(6) * t7 + t32;
t14 = t52 * pkin(7) + t16;
t11 = t58 * qJD(2) - t55 * t14;
t12 = t55 * qJD(2) + t58 * t14;
t72 = -t11 * t55 + t12 * t58;
t98 = m(6) * t72;
t86 = t56 * t40;
t15 = t59 * t33 - t86;
t6 = t15 * qJD(4) + t59 * t36 - t56 * t37;
t2 = t11 * qJD(5) + t58 * t6;
t96 = t2 * t58;
t3 = -t12 * qJD(5) - t55 * t6;
t95 = t3 * t55;
t71 = t59 * t44 - t56 * t45;
t94 = t7 * t71;
t93 = Ifges(6,1) * t55;
t92 = Ifges(6,4) * t55;
t88 = t52 * t58;
t87 = t55 * mrSges(6,3);
t82 = Ifges(6,5) * qJD(5);
t81 = Ifges(6,6) * qJD(5);
t80 = qJD(5) * t55;
t79 = qJD(5) * t58;
t77 = t79 / 0.2e1;
t76 = t100 * t52;
t73 = -t11 * t58 - t12 * t55;
t20 = t56 * t44 + t59 * t45;
t46 = qJD(5) * mrSges(6,1) - t52 * t87;
t47 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t88;
t70 = -t46 * t58 - t47 * t55;
t69 = -t55 * t46 + t58 * t47;
t68 = (Ifges(6,2) * t58 + t92) * t52;
t66 = t52 * mrSges(5,2) - t69;
t65 = t73 * qJD(5) - t95;
t62 = t65 + t96;
t13 = -t52 * pkin(4) - t15;
t28 = t68 + t81;
t48 = Ifges(6,4) * t88;
t29 = t52 * t93 + t48 + t82;
t61 = -t6 * mrSges(5,2) + mrSges(6,3) * t96 + t13 * t67 + t29 * t77 + qJD(5) ^ 2 * (Ifges(6,5) * t58 - Ifges(6,6) * t55) / 0.2e1 - (t28 + t68) * t80 / 0.2e1 - t100 * t7 + ((Ifges(6,1) * t58 - t92) * t80 + (0.3e1 * Ifges(6,4) * t58 - 0.2e1 * Ifges(6,2) * t55 + t93) * t77) * t52;
t50 = t56 * pkin(3) + pkin(7);
t18 = t59 * t39 - t86;
t17 = t56 * t39 + t84;
t9 = t20 * qJD(4) + t56 * t41 + t59 * t42;
t8 = t71 * qJD(4) + t59 * t41 - t56 * t42;
t1 = [-t71 * t32 - t76 * t9 + t70 * t20 * qJD(5) + (-t42 * mrSges(4,1) - t41 * mrSges(4,2)) * qJD(3) - t66 * t8 + m(4) * (t36 * t45 - t37 * t44 - t39 * t42 + t40 * t41) + m(5) * (-t15 * t9 + t16 * t8 + t6 * t20 - t94) + m(6) * (t13 * t9 + t62 * t20 + t72 * t8 - t94); m(6) * (t2 * t55 + t3 * t58) + (t98 + (-t55 ^ 2 - t58 ^ 2) * t52 * mrSges(6,3) + t69) * qJD(5); (m(5) * (t56 * t6 - t59 * t7) + ((-m(5) * t15 + m(6) * t13 - t76) * t56 + (m(5) * t16 - t66 + t98) * t59) * qJD(4)) * pkin(3) + t61 - m(6) * (t13 * t17 + t72 * t18) + (t73 * mrSges(6,3) + t70 * t50) * qJD(5) + t66 * t18 + t76 * t17 + m(6) * t62 * t50 - m(5) * (-t15 * t17 + t16 * t18) + (t39 * qJD(3) - t36) * mrSges(4,2) + (t40 * qJD(3) - t37) * mrSges(4,1) - t3 * t87 + t99 * (-t59 * pkin(3) - pkin(4)); (m(6) * (-t11 * t79 - t12 * t80 - t95 + t96) - t47 * t80 - t46 * t79) * pkin(7) + t61 - m(6) * (t13 * t16 + t72 * t15) + t65 * mrSges(6,3) + t76 * t16 + t66 * t15 - t99 * pkin(4); t3 * mrSges(6,1) - t2 * mrSges(6,2) - t11 * t47 + t12 * t46 + ((t82 / 0.2e1 - t13 * mrSges(6,2) - t29 / 0.2e1 - t48 / 0.2e1 + t11 * mrSges(6,3)) * t58 + (-t81 / 0.2e1 - t13 * mrSges(6,1) + t28 / 0.2e1 + t12 * mrSges(6,3) + (t92 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t58) * t52) * t55) * t52;];
tauc = t1(:);
