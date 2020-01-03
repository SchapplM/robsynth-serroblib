% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:13
% EndTime: 2019-12-31 16:30:15
% DurationCPUTime: 0.72s
% Computational Cost: add. (333->130), mult. (940->185), div. (0->0), fcn. (407->4), ass. (0->60)
t74 = -mrSges(4,1) - mrSges(5,1);
t73 = mrSges(5,2) + mrSges(4,3);
t31 = sin(qJ(2));
t48 = t31 * qJD(1);
t24 = qJD(2) * pkin(5) + t48;
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t72 = (t30 ^ 2 + t32 ^ 2) * t24;
t52 = qJD(2) * t30;
t54 = t74 * qJD(3) + t73 * t52;
t56 = t30 * t24;
t7 = -qJD(3) * pkin(3) + qJD(4) + t56;
t40 = m(5) * t7 + t54;
t8 = qJD(3) * qJ(4) + t32 * t24;
t70 = m(5) * t8;
t69 = -Ifges(4,1) / 0.2e1;
t68 = -Ifges(5,5) * t52 / 0.2e1;
t50 = qJD(2) * t32;
t67 = -Ifges(4,4) * t50 / 0.2e1;
t66 = t50 / 0.2e1;
t65 = -qJD(3) / 0.2e1;
t22 = mrSges(5,2) * t50 + qJD(3) * mrSges(5,3);
t63 = -t22 - t70;
t62 = 0.2e1 * pkin(5);
t61 = m(4) / 0.2e1;
t60 = m(5) / 0.2e1;
t18 = -t32 * pkin(3) - t30 * qJ(4) - pkin(2);
t33 = cos(qJ(2));
t47 = t33 * qJD(1);
t5 = t18 * qJD(2) - t47;
t59 = m(5) * t5;
t42 = qJD(2) * t47;
t49 = qJD(3) * t32;
t4 = t24 * t49 + t30 * t42;
t58 = t30 * t4;
t15 = (-t32 * mrSges(5,1) - t30 * mrSges(5,3)) * qJD(2);
t55 = t15 + (-t32 * mrSges(4,1) + t30 * mrSges(4,2)) * qJD(2);
t21 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t50;
t53 = t21 + t22;
t51 = qJD(2) * t31;
t46 = qJD(2) * qJD(3);
t45 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t44 = 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4);
t43 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t41 = t52 / 0.2e1;
t39 = -t53 - t70;
t38 = pkin(3) * t30 - qJ(4) * t32;
t25 = -qJD(2) * pkin(2) - t47;
t37 = t8 * mrSges(5,2) + Ifges(4,6) * qJD(3) / 0.2e1 + (Ifges(4,4) * t30 + Ifges(4,2) * t32) * qJD(2) / 0.2e1 + Ifges(5,6) * t65 + Ifges(5,3) * t66 + t68 - t25 * mrSges(4,1) - t5 * mrSges(5,1);
t6 = t38 * qJD(3) - t30 * qJD(4);
t36 = t5 * mrSges(5,3) - (Ifges(5,1) * t30 - Ifges(5,5) * t32) * qJD(2) / 0.2e1 + t52 * t69 + t67 - t25 * mrSges(4,2) - t7 * mrSges(5,2) + (Ifges(5,4) + Ifges(4,5)) * t65;
t35 = t40 * t32 + (-t21 + t63) * t30;
t34 = qJD(2) ^ 2;
t23 = t32 * t42;
t14 = (mrSges(4,1) * t30 + mrSges(4,2) * t32) * t46;
t13 = (mrSges(5,1) * t30 - mrSges(5,3) * t32) * t46;
t3 = -qJD(3) * t56 + t23;
t2 = t23 + (qJD(4) - t56) * qJD(3);
t1 = (t6 + t48) * qJD(2);
t9 = [(-t34 * mrSges(3,2) - t14 - t13 - m(5) * t1 + (t53 * t32 + t54 * t30 + m(5) * (t30 * t7 + t32 * t8) + m(4) * t72) * qJD(2)) * t33 + (-t34 * mrSges(3,1) + t55 * qJD(2) + m(4) * (qJD(2) * t25 + t3 * t32 - t42 + t58) + m(5) * (qJD(2) * t5 + t2 * t32 + t58) + t35 * qJD(3)) * t31; m(5) * (t1 * t18 + t5 * t6) - pkin(2) * t14 + t6 * t15 + t18 * t13 + ((-t55 - t59) * t31 + (-pkin(2) * t51 - t25 * t31 - t33 * t72) * m(4)) * qJD(1) + (-t1 * mrSges(5,1) + t2 * mrSges(5,2) + t3 * mrSges(4,3) + (t2 * t60 + t3 * t61) * t62 + (-mrSges(4,1) * t51 + t39 * t33) * qJD(1) + (t40 * pkin(5) + t45 * qJD(3) - t44 * t50 - t36) * qJD(3)) * t32 + (-t1 * mrSges(5,3) + ((t60 + t61) * t62 + t73) * t4 + (mrSges(4,2) * t51 - t40 * t33) * qJD(1) + (t43 * qJD(3) + t39 * pkin(5) + (t44 * t30 + (-0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t32) * qJD(2) - t37) * qJD(3)) * t30; -t3 * mrSges(4,2) + t2 * mrSges(5,3) + qJD(4) * t22 + t74 * t4 + m(5) * (-t4 * pkin(3) + t2 * qJ(4) + t8 * qJD(4)) - t35 * t24 + ((Ifges(5,5) * t66 + t36 + t67) * t32 + (t68 + Ifges(4,4) * t41 + (-Ifges(5,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t69) * t50 + t37) * t30 + ((-pkin(3) * mrSges(5,2) + t45) * t32 + (-qJ(4) * mrSges(5,2) + t43) * t30) * qJD(3) + (-t15 - t59) * t38) * qJD(2); (mrSges(5,2) * t49 + t30 * t15) * qJD(2) + 0.2e1 * (t4 / 0.2e1 + t5 * t41) * m(5) + t63 * qJD(3);];
tauc = t9(:);
