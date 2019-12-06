% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:10:36
% DurationCPUTime: 0.94s
% Computational Cost: add. (651->171), mult. (1792->253), div. (0->0), fcn. (1050->6), ass. (0->76)
t90 = mrSges(5,1) + mrSges(6,1);
t87 = mrSges(6,2) + mrSges(5,3);
t45 = sin(qJ(3));
t47 = cos(qJ(3));
t42 = sin(pkin(8));
t72 = qJD(1) * t42;
t26 = qJD(2) * t47 - t45 * t72;
t16 = t26 * qJD(3);
t44 = sin(qJ(4));
t69 = qJD(3) * t44;
t56 = t69 / 0.2e1;
t89 = -qJD(4) / 0.2e1;
t88 = qJD(4) / 0.2e1;
t27 = qJD(2) * t45 + t47 * t72;
t15 = qJD(3) * pkin(6) + t27;
t46 = cos(qJ(4));
t43 = cos(pkin(8));
t71 = qJD(1) * t43;
t9 = -t15 * t44 - t46 * t71;
t73 = -t9 + qJD(5);
t85 = 2 * m(5);
t84 = 2 * m(6);
t83 = pkin(6) / 0.2e1;
t38 = t44 * t71;
t65 = qJD(4) * t46;
t3 = -qJD(4) * t38 + t15 * t65 + t16 * t44;
t82 = t3 / 0.2e1;
t81 = -t26 / 0.2e1;
t80 = -t27 / 0.2e1;
t77 = t42 * t47;
t18 = t43 * t46 + t44 * t77;
t79 = t18 * t3;
t78 = t3 * t44;
t28 = (-t46 * mrSges(6,1) - t44 * mrSges(6,3)) * qJD(3);
t76 = t28 + (-t46 * mrSges(5,1) + t44 * mrSges(5,2)) * qJD(3);
t75 = t90 * qJD(4) - t87 * t69;
t68 = qJD(3) * t46;
t35 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t68;
t36 = mrSges(6,2) * t68 + qJD(4) * mrSges(6,3);
t74 = t35 + t36;
t67 = qJD(3) * t47;
t66 = qJD(4) * t44;
t17 = t27 * qJD(3);
t63 = qJD(3) * qJD(4);
t62 = t46 * t77;
t61 = pkin(6) * t82;
t60 = qJD(3) * t42 * t45;
t59 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t58 = -0.3e1 / 0.2e1 * Ifges(6,5) + 0.3e1 / 0.2e1 * Ifges(5,4);
t57 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t10 = t15 * t46 - t38;
t24 = (mrSges(6,1) * t44 - mrSges(6,3) * t46) * t63;
t25 = (mrSges(5,1) * t44 + mrSges(5,2) * t46) * t63;
t48 = qJD(3) ^ 2;
t54 = t48 * mrSges(4,2) + t24 + t25;
t53 = pkin(4) * t44 - qJ(5) * t46;
t31 = -pkin(4) * t46 - qJ(5) * t44 - pkin(3);
t51 = -t48 * mrSges(4,1) + qJD(3) * t76;
t13 = qJD(4) * t53 - qJD(5) * t44;
t11 = qJD(3) * t31 - t26;
t14 = -qJD(3) * pkin(3) - t26;
t41 = Ifges(5,4) * t68;
t5 = -qJD(4) * pkin(4) + t73;
t50 = t14 * mrSges(5,2) + t5 * mrSges(6,2) + (Ifges(6,1) * t44 - Ifges(6,5) * t46) * qJD(3) / 0.2e1 + Ifges(5,1) * t56 + t41 / 0.2e1 - t11 * mrSges(6,3) - t9 * mrSges(5,3) + (Ifges(6,4) + Ifges(5,5)) * t88;
t40 = Ifges(6,5) * t69;
t6 = qJD(4) * qJ(5) + t10;
t49 = t11 * mrSges(6,1) + t14 * mrSges(5,1) + Ifges(6,6) * t88 - Ifges(6,3) * t68 / 0.2e1 + t40 / 0.2e1 + Ifges(5,6) * t89 - (Ifges(5,4) * t44 + Ifges(5,2) * t46) * qJD(3) / 0.2e1 - t10 * mrSges(5,3) - t6 * mrSges(6,2);
t30 = t53 * qJD(3);
t19 = -t43 * t44 + t62;
t12 = t46 * t16;
t8 = -qJD(4) * t18 - t46 * t60;
t7 = -qJD(4) * t62 + (qJD(4) * t43 + t60) * t44;
t4 = (t13 + t27) * qJD(3);
t2 = qJD(4) * t9 + t12;
t1 = t12 + (qJD(5) + t9) * qJD(4);
t20 = [t74 * t8 + t75 * t7 + m(5) * (t10 * t8 + t19 * t2 + t7 * t9 + t79) + m(6) * (t1 * t19 - t5 * t7 + t6 * t8 + t79) + t87 * t63 * (t18 * t46 - t19 * t44) + (t54 * t45 + m(5) * (t14 * t67 + t17 * t45) + m(6) * (t11 * t67 + t4 * t45) + t51 * t47) * t42; ((-t44 * t75 + t74 * t46) * qJD(3) + m(5) * (t10 * t68 - t69 * t9 - t17) + m(6) * (t5 * t69 + t6 * t68 - t4) - t54) * t47 + ((-t44 * t74 - t46 * t75) * qJD(4) + m(5) * (qJD(3) * t14 - t10 * t66 + t2 * t46 - t65 * t9 + t78) + m(6) * (qJD(3) * t11 + t1 * t46 + t5 * t65 - t6 * t66 + t78) + t51) * t45; -t17 * mrSges(4,1) - pkin(3) * t25 + t13 * t28 + t31 * t24 + (qJD(3) * mrSges(4,1) - t76) * t27 + (t31 * t4 / 0.2e1 + (t80 + t13 / 0.2e1) * t11) * t84 + (t14 * t80 - pkin(3) * t17 / 0.2e1) * t85 + (-t17 * mrSges(5,1) - t4 * mrSges(6,1) + t1 * mrSges(6,2) + t2 * mrSges(5,3) - t74 * t26 + (t1 * t83 + t6 * t81) * t84 + (t10 * t81 + t2 * t83) * t85 + (t59 * qJD(4) + t58 * t68 + (-m(5) * t9 + m(6) * t5 - t75) * pkin(6) + t50) * qJD(4)) * t46 + (t17 * mrSges(5,2) - t4 * mrSges(6,3) + t87 * t3 + t75 * t26 + (t5 * t81 + t61) * t84 + (t26 * t9 / 0.2e1 + t61) * t85 + (t57 * qJD(4) + (-m(5) * t10 - m(6) * t6 - t74) * pkin(6) + (-t58 * t44 + (-0.3e1 / 0.2e1 * Ifges(6,3) - 0.3e1 / 0.2e1 * Ifges(5,2) + 0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(5,1)) * t46) * qJD(3) + t49) * qJD(4)) * t44; -t2 * mrSges(5,2) + t1 * mrSges(6,3) - t30 * t28 - t9 * t35 + t73 * t36 - t90 * t3 + t75 * t10 + ((-t41 / 0.2e1 + Ifges(6,5) * t68 / 0.2e1 + (-pkin(4) * mrSges(6,2) + t59) * qJD(4) - t50) * t46 + (-t40 / 0.2e1 + Ifges(5,4) * t56 + (-qJ(5) * mrSges(6,2) + t57) * qJD(4) + (-Ifges(6,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t68 - t49) * t44) * qJD(3) + (-t3 * pkin(4) + t1 * qJ(5) - t5 * t10 - t11 * t30 + t73 * t6) * m(6); -qJD(4) * t36 + (mrSges(6,2) * t65 + t44 * t28) * qJD(3) + (t11 * t56 + t6 * t89 + t82) * t84;];
tauc = t20(:);
