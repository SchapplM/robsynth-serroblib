% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:29
% DurationCPUTime: 0.65s
% Computational Cost: add. (850->128), mult. (1641->179), div. (0->0), fcn. (778->6), ass. (0->78)
t53 = cos(qJ(5));
t105 = -t53 / 0.2e1;
t49 = sin(pkin(8));
t52 = sin(qJ(2));
t91 = pkin(1) * qJD(1);
t83 = t52 * t91;
t42 = t49 * t83;
t50 = cos(pkin(8));
t54 = cos(qJ(2));
t82 = t54 * t91;
t34 = t50 * t82 - t42;
t84 = -t34 + qJD(4);
t51 = sin(qJ(5));
t104 = (-Ifges(6,5) * t51 / 0.2e1 + Ifges(6,6) * t105) * qJD(5);
t92 = mrSges(4,1) - mrSges(5,2);
t48 = qJD(1) + qJD(2);
t40 = t48 * pkin(2) + t82;
t15 = t49 * t40 + t50 * t83;
t10 = t48 * qJ(4) + t15;
t21 = t34 * qJD(2);
t11 = t48 * qJD(4) + t21;
t46 = t49 * pkin(2) + qJ(4);
t103 = t84 * t10 + t11 * t46;
t14 = t50 * t40 - t42;
t73 = qJD(4) - t14;
t8 = (-pkin(3) - pkin(7)) * t48 + t73;
t4 = -t51 * qJD(3) + t53 * t8;
t5 = t53 * qJD(3) + t51 * t8;
t75 = t4 * t53 + t5 * t51;
t93 = t50 * t52;
t65 = pkin(1) * (t49 * t54 + t93);
t33 = qJD(2) * t65;
t20 = qJD(1) * t33;
t1 = t4 * qJD(5) + t51 * t20;
t102 = t1 * mrSges(6,3);
t99 = Ifges(6,4) * t51;
t98 = Ifges(6,4) * t53;
t94 = t48 * mrSges(6,3);
t90 = pkin(1) * qJD(2);
t89 = qJD(5) * t5;
t44 = t49 * t52 * pkin(1);
t47 = t54 * pkin(1) + pkin(2);
t78 = t50 * t47 - t44;
t76 = -pkin(3) - t78;
t28 = -pkin(7) + t76;
t88 = qJD(5) * t28;
t81 = -t50 * pkin(2) - pkin(3);
t45 = -pkin(7) + t81;
t87 = qJD(5) * t45;
t85 = qJD(5) * t51;
t74 = -t4 * t51 + t5 * t53;
t2 = t53 * t20 - t89;
t72 = (-t2 - t89) * mrSges(6,3);
t71 = mrSges(3,1) * t52 + mrSges(3,2) * t54;
t70 = mrSges(6,1) * t53 - mrSges(6,2) * t51;
t69 = t51 * mrSges(6,1) + t53 * mrSges(6,2);
t35 = t50 * t54 * t90 - qJD(2) * t44;
t26 = qJD(4) + t35;
t66 = pkin(1) * t93 + t49 * t47;
t31 = qJ(4) + t66;
t68 = t10 * t26 + t11 * t31;
t64 = t51 * (-Ifges(6,2) * t53 - t99);
t63 = t53 * (-Ifges(6,1) * t51 - t98);
t62 = (Ifges(6,1) * t53 - t99) * t48;
t61 = (-Ifges(6,2) * t51 + t98) * t48;
t60 = t70 * qJD(5);
t59 = t71 * t90;
t56 = m(6) * (t74 * qJD(5) + t1 * t51 + t2 * t53);
t24 = Ifges(6,6) * qJD(5) + t61;
t25 = Ifges(6,5) * qJD(5) + t62;
t55 = t4 * mrSges(6,3) * t85 - t21 * mrSges(4,2) - qJD(1) * t59 + t10 * t60 - (t25 + t62) * t85 / 0.2e1 - t92 * t20 + (mrSges(5,3) + t69) * t11 + ((t63 - t64) * t48 + (t24 + t61) * t105 + t104) * qJD(5);
t39 = qJD(5) * mrSges(6,1) - t53 * t94;
t38 = -qJD(5) * mrSges(6,2) - t51 * t94;
t36 = t69 * t48;
t32 = qJD(1) * t65;
t27 = t48 * t60;
t9 = -t48 * pkin(3) + t73;
t3 = [t55 + (-t35 * mrSges(4,2) + t26 * mrSges(5,3) - t92 * t33 - t59) * t48 + m(6) * (t75 * t33 + t68) + t28 * t56 + (t33 * t39 + t38 * t88 + t72) * t53 + (t33 * t38 - t39 * t88 - t102) * t51 + m(4) * (-t14 * t33 + t15 * t35 - t20 * t78 + t21 * t66) + m(5) * (t20 * t76 + t9 * t33 + t68) + t26 * t36 + t31 * t27; t55 + t45 * t56 + (t34 * mrSges(4,2) + t84 * mrSges(5,3) + t92 * t32 + t71 * t91) * t48 + (-t32 * t39 + t38 * t87 + t72) * t53 + (-t32 * t38 - t39 * t87 - t102) * t51 + t84 * t36 + t46 * t27 + (-t75 * t32 + t103) * m(6) + (t20 * t81 - t9 * t32 + t103) * m(5) + ((-t20 * t50 + t21 * t49) * pkin(2) + t14 * t32 - t15 * t34) * m(4); m(6) * (t1 * t53 - t2 * t51) + (-m(6) * t75 - t51 * t38 - t53 * t39 + (-t51 ^ 2 - t53 ^ 2) * t94) * qJD(5); (t53 * t38 - t51 * t39) * qJD(5) + m(5) * t20 + t56 + (-mrSges(5,3) * t48 - t36 + (-m(5) - m(6)) * t10) * t48; t2 * mrSges(6,1) - t1 * mrSges(6,2) - t4 * t38 + t5 * t39 + (-t10 * t70 + t51 * t25 / 0.2e1 + t53 * t24 / 0.2e1 + (-t63 / 0.2e1 + t64 / 0.2e1) * t48 + t74 * mrSges(6,3) + t104) * t48;];
tauc = t3(:);
