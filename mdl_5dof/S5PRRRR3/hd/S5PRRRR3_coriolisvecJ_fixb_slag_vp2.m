% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:04
% EndTime: 2019-12-05 17:06:07
% DurationCPUTime: 0.74s
% Computational Cost: add. (1104->135), mult. (2130->206), div. (0->0), fcn. (1036->6), ass. (0->74)
t45 = sin(qJ(5));
t48 = cos(qJ(5));
t103 = mrSges(6,1) * t48 - mrSges(6,2) * t45 + mrSges(5,1);
t44 = qJD(2) + qJD(3);
t43 = qJD(4) + t44;
t59 = (mrSges(6,1) * t45 + mrSges(6,2) * t48) * qJD(5);
t23 = t43 * t59;
t50 = cos(qJ(3));
t84 = pkin(2) * qJD(2);
t36 = t44 * pkin(3) + t50 * t84;
t47 = sin(qJ(3));
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t86 = t47 * t49;
t63 = t46 * t50 + t86;
t79 = qJD(4) * t49;
t52 = (t63 * qJD(3) + t47 * t79) * pkin(2);
t80 = qJD(4) * t46;
t7 = qJD(2) * t52 + t36 * t80;
t102 = m(6) * t7 + t23;
t72 = t47 * t84;
t20 = t36 * t46 + t49 * t72;
t14 = pkin(8) * t43 + t20;
t10 = qJD(1) * t45 + t14 * t48;
t9 = qJD(1) * t48 - t14 * t45;
t67 = t10 * t48 - t45 * t9;
t101 = m(6) * t67;
t100 = mrSges(4,1) * t47 + mrSges(4,2) * t50;
t87 = t46 * t47;
t62 = t49 * t50 - t87;
t53 = (t62 * qJD(3) - t47 * t80) * pkin(2);
t6 = qJD(2) * t53 + t36 * t79;
t83 = qJD(5) * t9;
t2 = t48 * t6 + t83;
t98 = t2 * t48;
t78 = qJD(5) * t10;
t3 = -t45 * t6 - t78;
t97 = t3 * t45;
t93 = Ifges(6,1) * t45;
t92 = Ifges(6,4) * t45;
t89 = t43 * t48;
t88 = t45 * mrSges(6,3);
t42 = pkin(2) * t50 + pkin(3);
t85 = pkin(2) * t86 + t46 * t42;
t82 = Ifges(6,5) * qJD(5);
t81 = Ifges(6,6) * qJD(5);
t30 = pkin(8) + t85;
t77 = qJD(5) * t30;
t76 = qJD(5) * t45;
t75 = qJD(5) * t48;
t70 = t75 / 0.2e1;
t69 = t103 * t43;
t66 = -t10 * t45 - t48 * t9;
t34 = qJD(5) * mrSges(6,1) - t43 * t88;
t35 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t89;
t64 = -t45 * t34 + t48 * t35;
t61 = -pkin(2) * t87 + t42 * t49;
t60 = (Ifges(6,2) * t48 + t92) * t43;
t19 = t36 * t49 - t46 * t72;
t58 = t43 * mrSges(5,2) - t64;
t57 = t66 * qJD(5) - t97;
t13 = -pkin(4) * t43 - t19;
t21 = t60 + t81;
t37 = Ifges(6,4) * t89;
t22 = t43 * t93 + t37 + t82;
t54 = mrSges(6,3) * t98 + t13 * t59 + t22 * t70 + qJD(5) ^ 2 * (Ifges(6,5) * t48 - Ifges(6,6) * t45) / 0.2e1 - (t21 + t60) * t76 / 0.2e1 - t103 * t7 + ((Ifges(6,1) * t48 - t92) * t76 + (0.3e1 * Ifges(6,4) * t48 - 0.2e1 * Ifges(6,2) * t45 + t93) * t70) * t43;
t51 = -t6 * mrSges(5,2) + t54;
t40 = pkin(3) * t46 + pkin(8);
t32 = t62 * t84;
t31 = t63 * t84;
t29 = -pkin(4) - t61;
t12 = t42 * t80 + t52;
t11 = t42 * t79 + t53;
t1 = [m(6) * (t2 * t45 + t3 * t48) + (t101 + (-t45 ^ 2 - t48 ^ 2) * t43 * mrSges(6,3) + t64) * qJD(5); -t69 * t12 + m(6) * (t12 * t13 + t29 * t7) + m(5) * (t20 * t11 - t19 * t12 + t6 * t85 - t7 * t61) + (-t11 * t43 - t6) * mrSges(5,2) + t54 + (m(6) * (t10 * t11 + t2 * t30 - t9 * t77) + t11 * t35 - mrSges(6,3) * t83 - t34 * t77) * t48 + (m(6) * (-t10 * t77 - t11 * t9 - t3 * t30) - t11 * t34 - t35 * t77 + (-t3 - t78) * mrSges(6,3)) * t45 + t29 * t23 + t100 * pkin(2) * qJD(3) * (-qJD(2) - t44); ((-t48 * t34 - t45 * t35) * t40 + t66 * mrSges(6,3)) * qJD(5) + t58 * t32 + t69 * t31 + m(6) * (t57 + t98) * t40 + t100 * t84 * (-qJD(3) + t44) - m(5) * (-t19 * t31 + t20 * t32) + t51 + (m(5) * (t46 * t6 - t49 * t7) + ((-m(5) * t19 + m(6) * t13 - t69) * t46 + (m(5) * t20 + t101 - t58) * t49) * qJD(4)) * pkin(3) - t3 * t88 - m(6) * (t13 * t31 + t67 * t32) + t102 * (-pkin(3) * t49 - pkin(4)); t57 * mrSges(6,3) + t69 * t20 + t58 * t19 + t51 + (m(6) * (-t10 * t76 - t9 * t75 - t97 + t98) - t34 * t75 - t35 * t76) * pkin(8) - m(6) * (t13 * t20 + t67 * t19) - t102 * pkin(4); t3 * mrSges(6,1) - t2 * mrSges(6,2) + t10 * t34 - t9 * t35 + ((t82 / 0.2e1 - t13 * mrSges(6,2) - t22 / 0.2e1 - t37 / 0.2e1 + t9 * mrSges(6,3)) * t48 + (-t81 / 0.2e1 - t13 * mrSges(6,1) + t21 / 0.2e1 + t10 * mrSges(6,3) + (t92 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t48) * t43) * t45) * t43;];
tauc = t1(:);
