% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:39
% DurationCPUTime: 0.48s
% Computational Cost: add. (330->114), mult. (883->163), div. (0->0), fcn. (362->4), ass. (0->49)
t59 = mrSges(4,1) + mrSges(5,1);
t54 = mrSges(5,2) + mrSges(4,3);
t24 = sin(pkin(6)) * pkin(1) + pkin(5);
t18 = t24 * qJD(1);
t32 = sin(qJ(3));
t33 = cos(qJ(3));
t9 = t32 * qJD(2) + t33 * t18;
t3 = qJD(3) * qJ(4) + t9;
t58 = m(5) * t3;
t57 = -Ifges(4,6) / 0.2e1;
t44 = qJD(1) * t32;
t37 = t44 / 0.2e1;
t43 = qJD(1) * t33;
t56 = mrSges(5,2) * t43;
t55 = qJD(3) / 0.2e1;
t42 = t33 * qJD(2);
t49 = t32 * t18;
t8 = t42 - t49;
t45 = -t8 + qJD(4);
t53 = 0.2e1 * t24;
t52 = m(4) / 0.2e1;
t51 = m(5) / 0.2e1;
t7 = qJD(3) * t9;
t50 = t7 * t33;
t48 = t59 * qJD(3) - t54 * t44;
t22 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t43;
t23 = qJD(3) * mrSges(5,3) + t56;
t47 = t22 + t23;
t41 = 0.3e1 / 0.2e1 * Ifges(4,4) - 0.3e1 / 0.2e1 * Ifges(5,5);
t40 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t39 = Ifges(5,6) / 0.2e1 + t57;
t38 = -cos(pkin(6)) * pkin(1) - pkin(2);
t36 = pkin(3) * t32 - qJ(4) * t33;
t15 = -t33 * pkin(3) - t32 * qJ(4) + t38;
t10 = t36 * qJD(3) - t32 * qJD(4);
t2 = -qJD(3) * pkin(3) + t45;
t21 = t38 * qJD(1);
t27 = Ifges(4,4) * t43;
t4 = t15 * qJD(1);
t35 = t2 * mrSges(5,2) + t21 * mrSges(4,2) + (Ifges(5,1) * t32 - Ifges(5,5) * t33) * qJD(1) / 0.2e1 + Ifges(4,1) * t37 + t27 / 0.2e1 - t4 * mrSges(5,3) - t8 * mrSges(4,3) + (Ifges(5,4) + Ifges(4,5)) * t55;
t26 = Ifges(5,5) * t44;
t34 = t21 * mrSges(4,1) + t4 * mrSges(5,1) + Ifges(5,6) * t55 - Ifges(5,3) * t43 / 0.2e1 + t26 / 0.2e1 + qJD(3) * t57 - (Ifges(4,4) * t32 + Ifges(4,2) * t33) * qJD(1) / 0.2e1 - t3 * mrSges(5,2) - t9 * mrSges(4,3);
t25 = qJD(3) * t42;
t17 = t36 * qJD(1);
t16 = (-t33 * mrSges(5,1) - t32 * mrSges(5,3)) * qJD(1);
t6 = -qJD(3) * t49 + t25;
t5 = t10 * qJD(1);
t1 = t25 + (qJD(4) - t49) * qJD(3);
t11 = [m(5) * (t4 * t10 + t5 * t15) + t10 * t16 + (-t5 * mrSges(5,1) + t1 * mrSges(5,2) + t6 * mrSges(4,3) + (t1 * t51 + t6 * t52) * t53 + (t40 * qJD(3) + (-m(4) * t8 + m(5) * t2 - t48) * t24 + (t38 * mrSges(4,2) - t15 * mrSges(5,3) + t41 * t33) * qJD(1) + t35) * qJD(3)) * t33 + (-t5 * mrSges(5,3) + ((t52 + t51) * t53 + t54) * t7 + (t39 * qJD(3) + (-m(4) * t9 - t47 - t58) * t24 + (t15 * mrSges(5,1) + t38 * mrSges(4,1) - t41 * t32 + (-0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t33) * qJD(1) + t34) * qJD(3)) * t32; m(4) * (t6 * t32 - t50) + m(5) * (t1 * t32 - t50) + (t47 * t33 - t48 * t32 + m(4) * (-t32 * t8 + t33 * t9) + m(5) * (t2 * t32 + t3 * t33) + t54 * qJD(1) * (-t32 ^ 2 - t33 ^ 2)) * qJD(3); -t6 * mrSges(4,2) + t1 * mrSges(5,3) - t17 * t16 - t8 * t22 + t48 * t9 - t59 * t7 + t45 * t23 + ((-t27 / 0.2e1 + Ifges(5,5) * t43 / 0.2e1 + (-pkin(3) * mrSges(5,2) + t40) * qJD(3) - t35) * t33 + (-t26 / 0.2e1 + Ifges(4,4) * t37 + (-qJ(4) * mrSges(5,2) + t39) * qJD(3) + (-Ifges(5,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1) * t43 - t34) * t32) * qJD(1) + (-t7 * pkin(3) + t1 * qJ(4) - t4 * t17 - t2 * t9 + t45 * t3) * m(5); t16 * t44 + 0.2e1 * (t7 / 0.2e1 + t4 * t37) * m(5) + (-t23 + t56 - t58) * qJD(3);];
tauc = t11(:);
