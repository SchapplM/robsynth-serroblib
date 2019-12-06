% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:43
% EndTime: 2019-12-05 17:04:46
% DurationCPUTime: 0.70s
% Computational Cost: add. (1064->126), mult. (2096->192), div. (0->0), fcn. (1026->6), ass. (0->73)
t43 = sin(qJ(4));
t78 = qJD(4) * t43;
t42 = sin(qJ(5));
t45 = cos(qJ(5));
t100 = mrSges(6,1) * t45 - mrSges(6,2) * t42 + mrSges(5,1);
t46 = cos(qJ(4));
t47 = cos(qJ(3));
t73 = qJD(2) + qJD(3);
t68 = t73 * pkin(3);
t83 = pkin(2) * qJD(2);
t53 = t47 * t83 + t68;
t44 = sin(qJ(3));
t72 = t44 * t83;
t19 = t43 * t53 + t46 * t72;
t41 = qJD(4) + t73;
t13 = t41 * pkin(6) + t19;
t8 = qJD(1) * t45 - t13 * t42;
t9 = qJD(1) * t42 + t13 * t45;
t66 = -t42 * t8 + t45 * t9;
t99 = m(6) * t66;
t98 = (mrSges(4,1) * t44 + mrSges(4,2) * t47) * pkin(2);
t86 = t43 * t44;
t60 = t46 * t47 - t86;
t55 = t60 * qJD(3);
t59 = qJD(4) * t68;
t6 = t46 * t59 + (t60 * qJD(4) + t55) * t83;
t82 = qJD(5) * t8;
t2 = t45 * t6 + t82;
t97 = t2 * t45;
t81 = qJD(5) * t9;
t3 = -t42 * t6 - t81;
t96 = t3 * t42;
t85 = t44 * t46;
t61 = t43 * t47 + t85;
t56 = t61 * qJD(3);
t7 = t43 * t59 + (t61 * qJD(4) + t56) * t83;
t93 = t46 * t7;
t92 = Ifges(6,4) * t42;
t18 = t43 * t72 - t46 * t53;
t29 = t61 * t83;
t90 = t18 * t29;
t88 = t41 * t42;
t87 = t41 * t45;
t40 = pkin(2) * t47 + pkin(3);
t84 = pkin(2) * t85 + t43 * t40;
t80 = Ifges(6,5) * qJD(5);
t79 = Ifges(6,6) * qJD(5);
t77 = qJD(4) * t46;
t28 = pkin(6) + t84;
t76 = qJD(5) * t28;
t75 = qJD(5) * t42;
t74 = qJD(5) * t45;
t70 = t74 / 0.2e1;
t69 = t100 * t41;
t12 = t40 * t78 + (t44 * t77 + t56) * pkin(2);
t31 = pkin(2) * t86 - t40 * t46;
t65 = t18 * t12 + t7 * t31;
t34 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t88;
t35 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t87;
t62 = -t42 * t34 + t45 * t35;
t58 = (Ifges(6,2) * t45 + t92) * t41;
t57 = (mrSges(6,1) * t42 + mrSges(6,2) * t45) * qJD(5);
t54 = t41 * mrSges(5,2) - t62;
t20 = t58 + t79;
t36 = Ifges(6,4) * t87;
t21 = Ifges(6,1) * t88 + t36 + t80;
t50 = mrSges(6,3) * t97 + t18 * t57 + t21 * t70 + qJD(5) ^ 2 * (Ifges(6,5) * t45 - Ifges(6,6) * t42) / 0.2e1 - (t20 + t58) * t75 / 0.2e1 - t100 * t7 + ((Ifges(6,1) * t45 - t92) * t75 + (0.3e1 * Ifges(6,4) * t45 + (Ifges(6,1) - 0.2e1 * Ifges(6,2)) * t42) * t70) * t41;
t49 = -t34 * t74 - t35 * t75 + m(6) * (-t8 * t74 - t9 * t75 - t96 + t97);
t48 = (-t96 + (-t42 * t9 - t45 * t8) * qJD(5)) * mrSges(6,3) - t6 * mrSges(5,2) + t50;
t30 = t60 * t83;
t22 = t41 * t57;
t11 = t40 * t77 + (-t44 * t78 + t55) * pkin(2);
t1 = [m(6) * (t2 * t42 + t3 * t45) + (t99 + (-t42 ^ 2 - t45 ^ 2) * t41 * mrSges(6,3) + t62) * qJD(5); -t69 * t12 + m(6) * t65 + m(5) * (t19 * t11 + t6 * t84 + t65) + (-t11 * t41 - t6) * mrSges(5,2) + t31 * t22 + t50 + (-mrSges(6,3) * t82 - t34 * t76 + m(6) * (t11 * t9 + t2 * t28 - t8 * t76) + t11 * t35) * t45 + (-t35 * t76 + m(6) * (-t11 * t8 - t28 * t3 - t9 * t76) - t11 * t34 + (-t3 - t81) * mrSges(6,3)) * t42 + (-(2 * qJD(2)) - qJD(3)) * qJD(3) * t98; t54 * t30 + t69 * t29 + (qJD(2) ^ 2) * t98 + t48 - m(5) * (t19 * t30 + t90) - m(6) * (t66 * t30 + t90) + t49 * (pkin(3) * t43 + pkin(6)) + (-m(6) * t93 + m(5) * (t43 * t6 - t93) + (-t69 + (m(6) + m(5)) * t18) * t78 + (-t22 + (m(5) * t19 - t54 + t99) * qJD(4)) * t46) * pkin(3); t48 + t69 * t19 + (-m(6) * (t19 - t66) - t54) * t18 + t49 * pkin(6); t3 * mrSges(6,1) - t2 * mrSges(6,2) + t9 * t34 - t8 * t35 + ((t80 / 0.2e1 - t18 * mrSges(6,2) - t21 / 0.2e1 - t36 / 0.2e1 + t8 * mrSges(6,3)) * t45 + (-t79 / 0.2e1 - t18 * mrSges(6,1) + t20 / 0.2e1 + t9 * mrSges(6,3) + (t92 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t45) * t41) * t42) * t41;];
tauc = t1(:);
