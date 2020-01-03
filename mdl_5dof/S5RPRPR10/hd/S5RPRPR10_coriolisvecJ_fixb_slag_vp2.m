% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:54
% EndTime: 2019-12-31 18:25:57
% DurationCPUTime: 0.90s
% Computational Cost: add. (1562->159), mult. (2662->239), div. (0->0), fcn. (1238->6), ass. (0->84)
t82 = qJD(1) - qJD(3);
t114 = -qJD(5) / 0.2e1;
t63 = sin(qJ(5));
t100 = Ifges(6,4) * t63;
t65 = cos(qJ(5));
t113 = (Ifges(6,1) * t65 - t100) * t114;
t61 = sin(pkin(8));
t62 = cos(pkin(8));
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t42 = t61 * t66 + t62 * t64;
t95 = t82 * t42;
t41 = t61 * t64 - t62 * t66;
t112 = t82 * t41;
t67 = -pkin(1) - pkin(2);
t80 = -qJ(2) * t64 + t66 * t67;
t55 = qJD(1) * t67 + qJD(2);
t84 = qJ(2) * qJD(1);
t32 = t66 * t55 - t64 * t84;
t21 = -pkin(3) * t82 + t32;
t33 = t55 * t64 + t66 * t84;
t23 = t62 * t33;
t10 = t61 * t21 + t23;
t8 = -pkin(7) * t82 + t10;
t3 = qJD(4) * t65 - t63 * t8;
t4 = qJD(4) * t63 + t65 * t8;
t76 = -t3 * t63 + t4 * t65;
t83 = qJ(2) * qJD(3);
t86 = qJD(3) * t55;
t87 = qJD(2) * t66;
t18 = t66 * t86 + (-t64 * t83 + t87) * qJD(1);
t88 = qJD(2) * t64;
t19 = -t64 * t86 + (-t66 * t83 - t88) * qJD(1);
t6 = t18 * t62 + t19 * t61;
t1 = qJD(5) * t3 + t6 * t65;
t2 = -qJD(5) * t4 - t6 * t63;
t111 = t1 * t65 - t2 * t63;
t43 = (mrSges(6,1) * t63 + mrSges(6,2) * t65) * qJD(5);
t5 = t18 * t61 - t62 * t19;
t52 = -mrSges(6,1) * t65 + mrSges(6,2) * t63;
t98 = t33 * t61;
t9 = t21 * t62 - t98;
t7 = pkin(4) * t82 - t9;
t110 = t19 * mrSges(4,1) - t18 * mrSges(4,2) - t6 * mrSges(5,2) + t7 * t43 + (t52 - mrSges(5,1)) * t5;
t106 = -t65 / 0.2e1;
t107 = t63 / 0.2e1;
t53 = Ifges(6,2) * t65 + t100;
t89 = Ifges(6,6) * qJD(5);
t25 = -t53 * t82 + t89;
t96 = t82 * t65;
t56 = Ifges(6,4) * t96;
t90 = Ifges(6,5) * qJD(5);
t97 = t82 * t63;
t26 = -Ifges(6,1) * t97 - t56 + t90;
t77 = t3 * t65 + t4 * t63;
t109 = mrSges(6,3) * t77 + t26 * t106 + t25 * t107 + (Ifges(6,5) * t65 - Ifges(6,6) * t63) * t114;
t99 = Ifges(6,4) * t65;
t45 = (-Ifges(6,2) * t63 + t99) * qJD(5);
t108 = qJD(5) * ((Ifges(6,1) * t63 + t99) * t106 + t107 * t53) + t45 * t106 + t63 * t113;
t101 = t5 * t41;
t47 = -pkin(3) + t80;
t51 = qJ(2) * t66 + t64 * t67;
t93 = t61 * t47 + t62 * t51;
t85 = qJD(5) * t42;
t81 = m(3) * qJ(2) + mrSges(3,3);
t79 = t82 * t45 / 0.2e1 - t1 * mrSges(6,3);
t78 = -t2 * mrSges(6,3) + t113 * t82;
t74 = t47 * t62 - t51 * t61;
t48 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t97;
t49 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t96;
t69 = -m(6) * t77 - t65 * t48 - t63 * t49;
t58 = -pkin(3) * t62 - pkin(4);
t57 = pkin(3) * t61 + pkin(7);
t34 = t52 * t82;
t31 = -qJD(3) * t51 - t88;
t30 = t80 * qJD(3) + t87;
t27 = t82 * t43;
t16 = -pkin(7) + t93;
t15 = pkin(4) - t74;
t14 = t32 * t62 - t98;
t13 = t32 * t61 + t23;
t12 = t30 * t62 + t31 * t61;
t11 = t30 * t61 - t62 * t31;
t17 = [-t11 * t34 - t15 * t27 + 0.2e1 * t81 * qJD(2) * qJD(1) + (t12 * t49 + t79) * t65 + (-t12 * t48 - t78) * t63 + m(6) * (t11 * t7 + t111 * t16 + t76 * t12 + t15 * t5) + m(4) * (t18 * t51 + t19 * t80 + t33 * t30 + t32 * t31) + m(5) * (t10 * t12 - t9 * t11 - t5 * t74 + t6 * t93) + (t16 * t69 + t109) * qJD(5) - (t31 * mrSges(4,1) - t11 * mrSges(5,1) - t30 * mrSges(4,2) - t12 * mrSges(5,2) + t108) * t82 - t110; -t41 * t27 - t81 * qJD(1) ^ 2 + t95 * t34 + (t112 * t65 - t63 * t85) * t49 + (-t112 * t63 - t65 * t85) * t48 - (mrSges(5,1) * t95 - mrSges(5,2) * t112 + (mrSges(4,1) * t64 + mrSges(4,2) * t66) * t82) * t82 + (t101 + (-qJD(5) * t77 + t111) * t42 + t112 * t76 - t95 * t7) * m(6) + (t112 * t10 + t6 * t42 + t95 * t9 + t101) * m(5) + (t18 * t64 + t19 * t66 - t82 * (-t32 * t64 + t33 * t66)) * m(4); t13 * t34 - t58 * t27 + (-t14 * t49 - t79) * t65 + (t14 * t48 + t78) * t63 + (t57 * t69 - t109) * qJD(5) - (t33 * mrSges(4,1) + t13 * mrSges(5,1) + t32 * mrSges(4,2) + t14 * mrSges(5,2) - t108) * t82 + (t111 * t57 - t7 * t13 - t76 * t14 + t5 * t58) * m(6) + (-t10 * t14 + t9 * t13 + (-t5 * t62 + t6 * t61) * pkin(3)) * m(5) + t110; m(6) * (t1 * t63 + t2 * t65) + (m(6) * t76 + t65 * t49 - t63 * t48 - (-t63 ^ 2 - t65 ^ 2) * t82 * mrSges(6,3)) * qJD(5); t2 * mrSges(6,1) - t1 * mrSges(6,2) - t3 * t49 + t4 * t48 - ((t90 / 0.2e1 - t7 * mrSges(6,2) - t26 / 0.2e1 + t56 / 0.2e1 + t3 * mrSges(6,3)) * t65 + (-t89 / 0.2e1 - t7 * mrSges(6,1) + t25 / 0.2e1 + t4 * mrSges(6,3) - (t100 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t65) * t82) * t63) * t82;];
tauc = t17(:);
