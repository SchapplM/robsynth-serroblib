% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:47
% EndTime: 2019-12-05 15:04:50
% DurationCPUTime: 0.60s
% Computational Cost: add. (632->137), mult. (1661->215), div. (0->0), fcn. (1172->8), ass. (0->68)
t78 = -Ifges(6,1) / 0.2e1;
t52 = cos(qJ(5));
t62 = qJD(3) * t52;
t44 = Ifges(6,4) * t62;
t77 = -t44 / 0.2e1;
t51 = sin(qJ(3));
t53 = cos(qJ(3));
t47 = sin(pkin(8));
t65 = qJD(1) * t47;
t35 = t53 * qJD(2) - t51 * t65;
t46 = sin(pkin(9));
t48 = cos(pkin(9));
t33 = t46 * t51 - t48 * t53;
t49 = cos(pkin(8));
t41 = -t49 * qJD(1) + qJD(4);
t50 = sin(qJ(5));
t36 = t51 * qJD(2) + t53 * t65;
t23 = t48 * t36;
t25 = qJD(3) * pkin(3) + t35;
t12 = t46 * t25 + t23;
t8 = qJD(3) * pkin(6) + t12;
t5 = t52 * t41 - t50 * t8;
t76 = t5 * qJD(5);
t6 = t50 * t41 + t52 * t8;
t75 = t52 * t6;
t34 = t46 * t53 + t48 * t51;
t21 = t34 * t47;
t28 = t35 * qJD(3);
t29 = t36 * qJD(3);
t9 = t46 * t28 + t48 * t29;
t74 = t9 * t21;
t73 = t9 * t33;
t72 = Ifges(6,4) * t50;
t70 = t46 * t36;
t67 = Ifges(6,5) * qJD(5);
t66 = Ifges(6,6) * qJD(5);
t64 = qJD(3) * t47;
t63 = qJD(3) * t50;
t60 = t67 / 0.2e1;
t59 = -t66 / 0.2e1;
t58 = -t5 * t50 + t75;
t22 = t33 * t47;
t16 = -t52 * t22 - t49 * t50;
t15 = t50 * t22 - t49 * t52;
t11 = t48 * t25 - t70;
t39 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t63;
t40 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t62;
t57 = -t50 * t39 + t52 * t40;
t7 = -qJD(3) * pkin(4) - t11;
t56 = t5 * mrSges(6,3) + t63 * t78 + t77 - t67 / 0.2e1 - t7 * mrSges(6,2);
t55 = t6 * mrSges(6,3) + t66 / 0.2e1 + (Ifges(6,2) * t52 + t72) * qJD(3) / 0.2e1 - t7 * mrSges(6,1);
t54 = qJD(3) ^ 2;
t43 = -t48 * pkin(3) - pkin(4);
t42 = t46 * pkin(3) + pkin(6);
t37 = (-t52 * mrSges(6,1) + t50 * mrSges(6,2)) * qJD(3);
t32 = (mrSges(6,1) * t50 + mrSges(6,2) * t52) * qJD(5) * qJD(3);
t27 = t33 * qJD(3);
t26 = t34 * qJD(3);
t18 = t34 * t64;
t17 = t33 * t64;
t14 = t48 * t35 - t70;
t13 = t46 * t35 + t23;
t10 = t48 * t28 - t46 * t29;
t4 = -t16 * qJD(5) + t50 * t18;
t3 = t15 * qJD(5) - t52 * t18;
t2 = -t6 * qJD(5) - t50 * t10;
t1 = t52 * t10 + t76;
t19 = [-t17 * t37 + t21 * t32 + t3 * t40 + t4 * t39 + m(5) * (-t10 * t22 + t11 * t17 - t12 * t18 + t74) + m(6) * (t1 * t16 + t2 * t15 - t7 * t17 + t6 * t3 + t5 * t4 + t74) + (t17 * mrSges(5,1) + t18 * mrSges(5,2) + (-t15 * t52 - t16 * t50) * qJD(5) * mrSges(6,3)) * qJD(3) + (m(4) * (t28 * t53 + t29 * t51 + (-t35 * t53 - t36 * t51) * qJD(3)) + (-t53 * mrSges(4,1) + t51 * mrSges(4,2)) * t54) * t47; t33 * t32 + (-t51 * mrSges(4,1) - t53 * mrSges(4,2)) * t54 + (-qJD(3) * mrSges(5,1) + t37) * t26 + (-t39 * t52 - t40 * t50) * t34 * qJD(5) - (-qJD(3) * mrSges(5,2) + t57) * t27 + m(4) * (t28 * t51 - t29 * t53 + (-t35 * t51 + t36 * t53) * qJD(3)) + m(5) * (t10 * t34 - t11 * t26 - t12 * t27 + t73) + m(6) * (t7 * t26 + t73 - t58 * t27 + (t1 * t52 - t2 * t50 + (-t5 * t52 - t6 * t50) * qJD(5)) * t34); -t29 * mrSges(4,1) - t9 * mrSges(5,1) - t28 * mrSges(4,2) - t10 * mrSges(5,2) - t13 * t37 + t43 * t32 + (t36 * mrSges(4,1) + t13 * mrSges(5,1) + t35 * mrSges(4,2) + t14 * mrSges(5,2)) * qJD(3) + (-t9 * mrSges(6,1) + t1 * mrSges(6,3) - t14 * t40 + (-t42 * t39 + 0.3e1 / 0.2e1 * t44 + t60 - t56) * qJD(5)) * t52 + m(6) * (t9 * t43 + (t1 - t76) * t42 * t52) - m(6) * (t7 * t13 + t14 * t75) + (t9 * mrSges(6,2) + (t59 + (-m(6) * t6 - t40) * t42 + (-0.3e1 / 0.2e1 * t72 + (0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(6,2)) * t52) * qJD(3) - t55) * qJD(5) + (-t42 * m(6) - mrSges(6,3)) * t2 + (t5 * m(6) + t39) * t14) * t50 + (t11 * t13 - t12 * t14 + (t10 * t46 - t48 * t9) * pkin(3)) * m(5); m(6) * (t1 * t50 + t2 * t52) + (m(6) * t58 + (-t50 ^ 2 - t52 ^ 2) * qJD(3) * mrSges(6,3) + t57) * qJD(5); t2 * mrSges(6,1) - t1 * mrSges(6,2) + t6 * t39 - t5 * t40 + ((t60 + t77 + t56) * t52 + (t59 + (t72 / 0.2e1 + (t78 + Ifges(6,2) / 0.2e1) * t52) * qJD(3) + t55) * t50) * qJD(3);];
tauc = t19(:);
