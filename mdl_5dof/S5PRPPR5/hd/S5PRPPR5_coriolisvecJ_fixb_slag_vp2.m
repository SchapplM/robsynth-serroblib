% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:05
% DurationCPUTime: 0.80s
% Computational Cost: add. (642->131), mult. (1350->189), div. (0->0), fcn. (713->6), ass. (0->72)
t38 = sin(qJ(5));
t72 = mrSges(6,3) * qJD(2);
t28 = qJD(5) * mrSges(6,1) + t38 * t72;
t40 = cos(qJ(5));
t29 = -qJD(5) * mrSges(6,2) - t40 * t72;
t96 = -t38 * t28 + t40 * t29;
t100 = qJD(2) * mrSges(5,2) + t96;
t42 = -pkin(2) - pkin(3);
t41 = cos(qJ(2));
t70 = qJD(1) * t41;
t58 = qJD(3) - t70;
t21 = t42 * qJD(2) + t58;
t39 = sin(qJ(2));
t71 = qJD(1) * t39;
t31 = qJD(2) * qJ(3) + t71;
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t8 = t36 * t21 + t37 * t31;
t6 = -qJD(2) * pkin(6) + t8;
t3 = qJD(4) * t40 - t38 * t6;
t4 = qJD(4) * t38 + t40 * t6;
t56 = -t3 * t38 + t4 * t40;
t98 = m(6) * t56;
t99 = -t96 - t98;
t54 = t40 * mrSges(6,1) - t38 * mrSges(6,2);
t25 = t54 * qJD(2);
t97 = qJD(2) * mrSges(5,1) + t25;
t20 = -t36 * t41 + t37 * t39;
t95 = mrSges(5,1) * t36 + mrSges(5,2) * t37;
t26 = (qJD(3) + t70) * qJD(2);
t64 = qJD(1) * qJD(2);
t61 = t39 * t64;
t10 = t26 * t37 + t36 * t61;
t1 = t3 * qJD(5) + t10 * t40;
t2 = -t4 * qJD(5) - t10 * t38;
t94 = t1 * t40 - t2 * t38;
t7 = t21 * t37 - t31 * t36;
t93 = m(5) * (-t36 * t7 + t37 * t8) + t36 * t25;
t19 = t36 * t39 + t37 * t41;
t9 = t26 * t36 - t37 * t61;
t91 = t19 * t9;
t5 = qJD(2) * pkin(4) - t7;
t88 = t36 * t5;
t87 = t37 * t9;
t83 = Ifges(6,4) * t38;
t82 = Ifges(6,4) * t40;
t81 = Ifges(6,5) * t40;
t80 = Ifges(6,6) * t38;
t79 = t31 * t41;
t73 = t37 * qJ(3) + t36 * t42;
t66 = qJD(5) * t38;
t65 = qJD(5) * t40;
t63 = qJD(2) * qJD(3);
t57 = -t3 * t40 - t38 * t4;
t53 = -mrSges(6,1) * t38 - mrSges(6,2) * t40;
t51 = -qJ(3) * t36 + t37 * t42;
t50 = t38 * (-Ifges(6,1) * t40 + t83);
t49 = t40 * (Ifges(6,2) * t38 - t82);
t48 = qJD(5) * (-t28 * t40 - t29 * t38);
t47 = t53 * qJD(5);
t46 = (-Ifges(6,1) * t38 - t82) * qJD(2);
t45 = (-Ifges(6,2) * t40 - t83) * qJD(2);
t44 = t57 * qJD(5) + t94;
t43 = qJD(2) ^ 2;
t27 = -qJD(2) * pkin(2) + t58;
t22 = pkin(4) - t51;
t18 = qJD(2) * t47;
t17 = Ifges(6,5) * qJD(5) + t46;
t16 = Ifges(6,6) * qJD(5) + t45;
t15 = t19 * qJD(2);
t14 = t20 * qJD(2);
t11 = [t19 * t18 - t97 * t14 + t20 * t48 + t100 * t15 + m(4) * (t26 * t39 + (t79 + (t27 - t70) * t39) * qJD(2)) + m(5) * (t10 * t20 + t14 * t7 + t15 * t8 + t91) + m(6) * (-t14 * t5 + t56 * t15 + t44 * t20 + t91) + ((-mrSges(3,2) + mrSges(4,3)) * t41 + (-mrSges(3,1) - mrSges(4,1)) * t39) * t43; t5 * t47 + qJD(5) ^ 2 * (t80 - t81) / 0.2e1 + t22 * t18 + (m(5) * t73 + mrSges(5,2)) * t10 + (-m(5) * t51 + m(6) * t22 + mrSges(5,1) + t54) * t9 + t95 * t63 + (-t50 - t49) * qJD(2) * qJD(5) + (t45 + t16) * t66 / 0.2e1 - (t46 + t17) * t65 / 0.2e1 + (m(6) * t44 - t28 * t65 - t29 * t66) * (-pkin(6) + t73) + (-pkin(2) * t61 + t26 * qJ(3) - (t27 * t39 + t79) * qJD(1)) * m(4) + (m(5) * t7 - m(6) * t5 - t97) * (t36 * t70 - t37 * t71) + (-m(5) * t8 - t100 - t98) * t19 * qJD(1) + (-t41 * t64 + t26 + t63) * mrSges(4,3) + (t3 * t65 + t4 * t66 - t94) * mrSges(6,3) + (m(6) * (t56 * t37 + t88) + t93 + t96 * t37 + m(4) * t31) * qJD(3); -t37 * t18 + t36 * t48 + m(5) * (t10 * t36 - t87) + m(6) * (t44 * t36 - t87) + (-mrSges(4,3) - t95) * t43 + ((-t31 + t71) * m(4) - m(6) * t88 + t99 * t37 - t93) * qJD(2); m(6) * (t1 * t38 + t2 * t40) + ((t38 ^ 2 + t40 ^ 2) * t72 - t99) * qJD(5); t2 * mrSges(6,1) - t1 * mrSges(6,2) + t4 * t28 - t3 * t29 + (-t5 * t53 + t40 * t17 / 0.2e1 - t38 * t16 / 0.2e1 + (t50 / 0.2e1 + t49 / 0.2e1) * qJD(2) + t57 * mrSges(6,3) + (-t81 / 0.2e1 + t80 / 0.2e1) * qJD(5)) * qJD(2);];
tauc = t11(:);
