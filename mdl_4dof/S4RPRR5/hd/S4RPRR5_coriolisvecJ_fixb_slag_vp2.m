% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:33
% EndTime: 2019-12-31 16:51:34
% DurationCPUTime: 0.43s
% Computational Cost: add. (604->111), mult. (1077->165), div. (0->0), fcn. (409->4), ass. (0->66)
t32 = -qJD(1) + qJD(3);
t83 = -t32 / 0.2e1;
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t39 = -pkin(1) - pkin(2);
t63 = t38 * qJ(2) + t36 * t39;
t57 = qJ(2) * qJD(3);
t59 = t38 * qJD(2);
t28 = t39 * qJD(1) + qJD(2);
t61 = qJD(3) * t28;
t3 = t38 * t61 + (-t36 * t57 + t59) * qJD(1);
t35 = sin(qJ(4));
t37 = cos(qJ(4));
t58 = qJ(2) * qJD(1);
t15 = t36 * t28 + t38 * t58;
t6 = t32 * pkin(6) + t15;
t62 = qJD(4) * t6;
t1 = t37 * t3 - t35 * t62;
t2 = -t35 * t3 - t37 * t62;
t48 = t1 * t37 - t2 * t35;
t82 = (m(3) * qJ(2) + mrSges(3,3)) * qJD(1);
t56 = (t35 ^ 2 + t37 ^ 2) * t6;
t47 = mrSges(5,1) * t35 + mrSges(5,2) * t37;
t17 = t47 * qJD(4);
t25 = -t37 * mrSges(5,1) + t35 * mrSges(5,2);
t60 = t36 * qJD(2);
t4 = t36 * t61 + (t38 * t57 + t60) * qJD(1);
t14 = t38 * t28 - t36 * t58;
t5 = -t32 * pkin(3) - t14;
t81 = -t3 * mrSges(4,2) + t5 * t17 + (t25 - mrSges(4,1)) * t4;
t66 = t32 * t35;
t23 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t66;
t65 = t32 * t37;
t24 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t65;
t80 = -(t32 * mrSges(4,2) + t35 * t23 - t37 * t24) * t38 + m(4) * (-t14 * t36 + t15 * t38) + m(5) * (t36 * t5 + t38 * t56);
t79 = -t35 / 0.2e1;
t78 = t35 / 0.2e1;
t77 = -t37 / 0.2e1;
t72 = Ifges(5,4) * t35;
t71 = Ifges(5,4) * t37;
t70 = Ifges(5,5) * t37;
t69 = Ifges(5,2) * t37;
t68 = Ifges(5,6) * t35;
t67 = t32 * (t69 + t72);
t16 = t25 * t32;
t54 = t32 * mrSges(4,1) - t16;
t19 = (-Ifges(5,2) * t35 + t71) * qJD(4);
t53 = -t1 * mrSges(5,3) + t19 * t83;
t46 = Ifges(5,1) * t37 - t72;
t20 = t46 * qJD(4);
t52 = t2 * mrSges(5,3) + t20 * t83;
t7 = Ifges(5,6) * qJD(4) + t67;
t51 = -t7 / 0.2e1 - t67 / 0.2e1;
t29 = Ifges(5,4) * t65;
t8 = Ifges(5,1) * t66 + Ifges(5,5) * qJD(4) + t29;
t50 = -t8 / 0.2e1 + (Ifges(5,1) * t35 + t71) * t83;
t45 = t37 * t23 + t35 * t24;
t44 = -t36 * qJ(2) + t38 * t39;
t43 = t19 * t77 + t20 * t79;
t22 = -pkin(6) + t63;
t21 = pkin(3) - t44;
t18 = qJD(4) * (-t68 + t70);
t13 = t63 * qJD(3) + t60;
t12 = t44 * qJD(3) + t59;
t9 = t32 * t17;
t10 = [t13 * t16 + t21 * t9 + 0.2e1 * qJD(2) * t82 + (t12 * t24 + t53) * t37 + (-t12 * t23 + t52) * t35 + m(5) * (t12 * t56 + t5 * t13 + t4 * t21 + t48 * t22) + m(4) * (t15 * t12 - t14 * t13 + t3 * t63 - t4 * t44) + (-t13 * mrSges(4,1) - t12 * mrSges(4,2) + t43) * t32 + (-t18 / 0.2e1 + (-t22 * t23 + t50) * t37 + (-t22 * t24 - t51) * t35) * qJD(4) - t81; t80 * qJD(3) + (m(4) * t3 + m(5) * t48 - t54 * qJD(3) - t45 * qJD(4)) * t36 + (t36 * t54 - t80 - t82) * qJD(1) + (-t9 + (-m(4) - m(5)) * t4) * t38; -pkin(3) * t9 - t15 * t16 + (-t14 * t24 - t53) * t37 + (t14 * t23 - t52) * t35 + (t15 * mrSges(4,1) + t14 * mrSges(4,2) - t43) * t32 + (t18 / 0.2e1 + (-pkin(6) * t23 - t50) * t37 + (-pkin(6) * t24 + t51) * t35) * qJD(4) + (-t4 * pkin(3) + pkin(6) * t48 - t14 * t56 - t5 * t15) * m(5) + t81; t2 * mrSges(5,1) - t1 * mrSges(5,2) + t45 * t6 + (-t5 * t47 + t7 * t78 + (t46 * t79 + t69 * t78) * t32 + (t70 / 0.2e1 - t68 / 0.2e1) * qJD(4) + (t8 + t29) * t77) * t32;];
tauc = t10(:);
