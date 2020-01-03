% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:44
% EndTime: 2019-12-31 17:50:46
% DurationCPUTime: 0.63s
% Computational Cost: add. (534->144), mult. (1145->191), div. (0->0), fcn. (452->4), ass. (0->58)
t34 = sin(pkin(7)) * pkin(1) + qJ(3);
t26 = qJD(1) * t34;
t71 = (m(4) + m(5)) * t26;
t70 = Ifges(5,4) + Ifges(6,4);
t68 = qJD(1) / 0.2e1;
t67 = qJD(4) / 0.2e1;
t40 = sin(qJ(4));
t61 = qJD(1) * t40;
t27 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t61;
t28 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t61;
t66 = t27 + t28;
t41 = cos(qJ(4));
t60 = qJD(1) * t41;
t29 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t60;
t30 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t60;
t65 = -t29 - t30;
t33 = -cos(pkin(7)) * pkin(1) - pkin(2) - pkin(6);
t63 = qJ(5) - t33;
t62 = mrSges(4,3) * qJD(1);
t59 = qJD(4) * t40;
t58 = qJD(4) * t41;
t57 = t40 * qJD(2);
t56 = t40 * qJD(5);
t55 = t41 * qJD(2);
t54 = t41 * qJD(5);
t53 = qJ(5) * qJD(1);
t52 = qJD(1) * qJD(3);
t51 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t50 = 0.3e1 / 0.2e1 * Ifges(6,4) + 0.3e1 / 0.2e1 * Ifges(5,4);
t49 = -Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1;
t48 = -Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t47 = m(5) * t33 - mrSges(5,3);
t25 = t40 * pkin(4) + t34;
t14 = t25 * qJD(1) + qJD(5);
t22 = (t40 * mrSges(6,1) + t41 * mrSges(6,2)) * qJD(1);
t46 = -m(6) * t14 - t22;
t16 = t63 * t41;
t21 = t33 * qJD(1) + qJD(3);
t10 = t41 * t21 - t57;
t31 = pkin(4) * t58 + qJD(3);
t6 = -qJD(4) * t57 + t21 * t58;
t11 = t40 * t21 + t55;
t4 = -t41 * t53 + t10;
t3 = qJD(4) * pkin(4) + t4;
t5 = -t40 * t53 + t11;
t44 = m(6) * (-t3 * t41 - t40 * t5);
t43 = t26 * mrSges(5,2) + t14 * mrSges(6,2) - t3 * mrSges(6,3) + ((Ifges(5,1) + Ifges(6,1)) * t41 - t70 * t40) * t68 + (Ifges(6,5) + Ifges(5,5)) * t67;
t42 = -t26 * mrSges(5,1) - t14 * mrSges(6,1) + t5 * mrSges(6,3) + (t70 * t41 + (-Ifges(5,2) - Ifges(6,2)) * t40) * t68 + (Ifges(6,6) + Ifges(5,6)) * t67;
t32 = qJD(1) * mrSges(6,1) * t58;
t24 = t31 * qJD(1);
t23 = qJD(1) * (t40 * mrSges(5,1) + t41 * mrSges(5,2));
t15 = t63 * t40;
t9 = -qJD(4) * t16 - t56;
t8 = t63 * t59 - t54;
t7 = t11 * qJD(4);
t2 = -qJD(1) * t54 + (-t55 + (-t21 + t53) * t40) * qJD(4);
t1 = (-qJ(5) * t58 - t56) * qJD(1) + t6;
t12 = [t25 * t32 + t9 * t27 + t8 * t29 + t31 * t22 + m(6) * (-t1 * t15 + t14 * t31 - t2 * t16 + t24 * t25 + t3 * t8 + t5 * t9) + (t23 + 0.2e1 * t62 + 0.2e1 * t71) * qJD(3) + (mrSges(5,2) * t52 + t24 * mrSges(6,2) - t2 * mrSges(6,3) - t47 * t7 + (t33 * t28 + t48 * qJD(4) + t47 * t11 + (t34 * mrSges(5,1) + t15 * mrSges(6,3) - t50 * t41) * qJD(1) - t42) * qJD(4)) * t41 + (mrSges(5,1) * t52 + t24 * mrSges(6,1) - t1 * mrSges(6,3) + t47 * t6 + (-t33 * t30 + t49 * qJD(4) - t47 * t10 + (-t16 * mrSges(6,3) - t34 * mrSges(5,2) - t25 * mrSges(6,2) + t50 * t40 + (-0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t41) * qJD(1) - t43) * qJD(4)) * t40; m(5) * (t7 * t40 + t6 * t41) + m(6) * (t1 * t41 - t2 * t40) + (t65 * t41 - t66 * t40 + m(5) * (-t10 * t41 - t11 * t40) + t44 + (mrSges(5,3) + mrSges(6,3)) * qJD(1) * (-t40 ^ 2 - t41 ^ 2)) * qJD(4); (t65 * t40 + t66 * t41) * qJD(4) + m(5) * (t6 * t40 - t7 * t41 + (-t10 * t40 + t11 * t41) * qJD(4)) + m(6) * (t1 * t40 + t2 * t41 + (-t3 * t40 + t41 * t5) * qJD(4)) + (-t23 + t46 - t62 - t71) * qJD(1); -t7 * mrSges(5,1) - t6 * mrSges(5,2) - t1 * mrSges(6,2) - t10 * t28 + t11 * t30 - t4 * t27 + (m(6) * pkin(4) + mrSges(6,1)) * t2 + (t29 - m(6) * (-t3 + t4)) * t5 + ((-t10 * mrSges(5,3) - t51 * t61 + t43) * t40 + (t11 * mrSges(5,3) + t51 * t60 + t46 * pkin(4) + (Ifges(6,1) / 0.2e1 - Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t61 + t42) * t41 + (t48 * t41 + (pkin(4) * mrSges(6,3) + t49) * t40) * qJD(4)) * qJD(1); t32 + m(6) * t24 + (-mrSges(6,2) * t59 + t40 * t27 + t41 * t29 - t44) * qJD(1);];
tauc = t12(:);
