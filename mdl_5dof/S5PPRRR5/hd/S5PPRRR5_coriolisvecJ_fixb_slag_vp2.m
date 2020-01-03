% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:39
% EndTime: 2019-12-31 17:35:40
% DurationCPUTime: 0.55s
% Computational Cost: add. (722->114), mult. (1559->180), div. (0->0), fcn. (926->6), ass. (0->63)
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t89 = t45 * mrSges(6,1) - t42 * mrSges(6,2) + mrSges(5,1);
t41 = qJD(3) + qJD(4);
t57 = (mrSges(6,1) * t42 + mrSges(6,2) * t45) * qJD(5);
t25 = t41 * t57;
t47 = cos(qJ(3));
t36 = qJD(3) * pkin(3) + t47 * qJD(2);
t44 = sin(qJ(3));
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t33 = t43 * t47 + t46 * t44;
t56 = t33 * qJD(3);
t71 = qJD(4) * t46;
t72 = qJD(4) * t43;
t7 = t36 * t72 + (t44 * t71 + t56) * qJD(2);
t88 = m(6) * t7 + t25;
t73 = qJD(2) * t44;
t22 = t43 * t36 + t46 * t73;
t16 = t41 * pkin(7) + t22;
t59 = t42 * qJD(1) - t45 * t16;
t9 = -t45 * qJD(1) - t42 * t16;
t64 = -t42 * t9 - t45 * t59;
t87 = m(6) * t64;
t32 = t43 * t44 - t46 * t47;
t55 = t32 * qJD(3);
t6 = t36 * t71 + (-t44 * t72 - t55) * qJD(2);
t2 = qJD(5) * t9 + t45 * t6;
t85 = t2 * t45;
t3 = t59 * qJD(5) - t42 * t6;
t84 = t3 * t42;
t82 = t7 * t32;
t81 = Ifges(6,1) * t42;
t80 = Ifges(6,4) * t42;
t77 = t41 * t45;
t76 = t42 * mrSges(6,3);
t75 = Ifges(6,5) * qJD(5);
t74 = Ifges(6,6) * qJD(5);
t70 = qJD(5) * t42;
t69 = qJD(5) * t45;
t67 = t69 / 0.2e1;
t66 = t89 * t41;
t63 = t42 * t59 - t9 * t45;
t34 = qJD(5) * mrSges(6,1) - t41 * t76;
t35 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t77;
t61 = -t34 * t45 - t35 * t42;
t60 = t42 * t34 - t45 * t35;
t21 = t46 * t36 - t43 * t73;
t58 = (Ifges(6,2) * t45 + t80) * t41;
t54 = t41 * mrSges(5,2) + t60;
t53 = t63 * qJD(5) - t84;
t50 = t53 + t85;
t15 = -t41 * pkin(4) - t21;
t23 = t58 + t74;
t37 = Ifges(6,4) * t77;
t24 = t41 * t81 + t37 + t75;
t49 = -t6 * mrSges(5,2) + mrSges(6,3) * t85 + t15 * t57 + t24 * t67 + qJD(5) ^ 2 * (Ifges(6,5) * t45 - Ifges(6,6) * t42) / 0.2e1 - t89 * t7 - (t23 + t58) * t70 / 0.2e1 + ((Ifges(6,1) * t45 - t80) * t70 + (0.3e1 * Ifges(6,4) * t45 - 0.2e1 * Ifges(6,2) * t42 + t81) * t67) * t41;
t38 = t43 * pkin(3) + pkin(7);
t31 = t32 * qJD(2);
t30 = t33 * qJD(2);
t12 = t33 * qJD(4) + t56;
t11 = -t32 * qJD(4) - t55;
t1 = [m(6) * (-t2 * t42 - t3 * t45) + (-t87 + (t42 ^ 2 + t45 ^ 2) * t41 * mrSges(6,3) + t60) * qJD(5); t32 * t25 + (-t44 * mrSges(4,1) - t47 * mrSges(4,2)) * qJD(3) ^ 2 - t66 * t12 + t61 * t33 * qJD(5) - t54 * t11 + m(5) * (t22 * t11 - t21 * t12 + t6 * t33 + t82) + m(6) * (t64 * t11 + t15 * t12 + t50 * t33 + t82); t49 + (m(5) * (t43 * t6 - t46 * t7) + ((-m(5) * t21 + m(6) * t15 - t66) * t43 + (m(5) * t22 - t54 + t87) * t46) * qJD(4)) * pkin(3) - m(6) * (t15 * t30 - t64 * t31) + (t63 * mrSges(6,3) + t61 * t38) * qJD(5) + m(6) * t50 * t38 - t54 * t31 + t66 * t30 - m(5) * (-t21 * t30 - t22 * t31) - t3 * t76 + t88 * (-t46 * pkin(3) - pkin(4)); t49 + (-t34 * t69 - t35 * t70 + m(6) * (t59 * t70 - t9 * t69 - t84 + t85)) * pkin(7) - m(6) * (t15 * t22 + t64 * t21) + t53 * mrSges(6,3) + t66 * t22 + t54 * t21 - t88 * pkin(4); t3 * mrSges(6,1) - t2 * mrSges(6,2) - t59 * t34 - t9 * t35 + ((t75 / 0.2e1 - t15 * mrSges(6,2) - t24 / 0.2e1 - t37 / 0.2e1 + t9 * mrSges(6,3)) * t45 + (-t74 / 0.2e1 - t15 * mrSges(6,1) + t23 / 0.2e1 - t59 * mrSges(6,3) + (t80 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t45) * t41) * t42) * t41;];
tauc = t1(:);
