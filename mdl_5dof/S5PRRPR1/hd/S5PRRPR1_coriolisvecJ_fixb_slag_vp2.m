% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:34
% EndTime: 2019-12-05 16:15:37
% DurationCPUTime: 0.81s
% Computational Cost: add. (1099->143), mult. (2015->211), div. (0->0), fcn. (1219->6), ass. (0->76)
t77 = cos(qJ(3));
t92 = qJD(2) * pkin(2);
t81 = -t77 * t92 + qJD(4);
t72 = sin(pkin(9));
t73 = cos(pkin(9));
t74 = sin(qJ(5));
t76 = cos(qJ(5));
t55 = -t74 * t72 + t76 * t73;
t49 = t55 * qJD(5);
t111 = t49 / 0.2e1;
t56 = t76 * t72 + t74 * t73;
t50 = t56 * qJD(5);
t110 = -t50 / 0.2e1;
t71 = qJD(2) + qJD(3);
t37 = t71 * t49;
t109 = t55 * t37;
t38 = t71 * t50;
t108 = t56 * t38;
t85 = qJD(3) * t92;
t82 = t77 * t85;
t54 = t71 * qJD(4) + t82;
t93 = t72 ^ 2 + t73 ^ 2;
t84 = t93 * t54;
t60 = (-pkin(7) - qJ(4)) * t72;
t67 = t73 * pkin(7);
t61 = t73 * qJ(4) + t67;
t30 = t76 * t60 - t74 * t61;
t107 = t30 * qJD(5) + t81 * t55;
t31 = t74 * t60 + t76 * t61;
t106 = -t31 * qJD(5) - t81 * t56;
t104 = t93 * mrSges(5,3);
t105 = t71 * t104;
t43 = t55 * t71;
t44 = t56 * t71;
t80 = -t73 * mrSges(5,1) + t72 * mrSges(5,2);
t103 = -t43 * mrSges(6,1) + t44 * mrSges(6,2) + t80 * t71;
t75 = sin(qJ(3));
t88 = t75 * t92;
t58 = t71 * qJ(4) + t88;
t66 = t73 * qJD(1);
t28 = t66 + (-pkin(7) * t71 - t58) * t72;
t41 = t72 * qJD(1) + t73 * t58;
t29 = t71 * t67 + t41;
t10 = t74 * t28 + t76 * t29;
t3 = -t10 * qJD(5) - t56 * t54;
t9 = t76 * t28 - t74 * t29;
t102 = -t3 * t56 - t9 * t49;
t100 = t44 / 0.2e1;
t98 = t75 * pkin(2);
t97 = t77 * pkin(2);
t94 = Ifges(6,4) * t44;
t91 = qJD(3) * t77;
t86 = mrSges(4,2) * t91;
t64 = -t73 * pkin(4) - pkin(3);
t18 = t38 * mrSges(6,1) + t37 * mrSges(6,2);
t83 = t75 * t85;
t79 = -(-t72 * t58 + t66) * t72 + t41 * t73;
t63 = qJ(4) + t98;
t52 = (-pkin(7) - t63) * t72;
t53 = t73 * t63 + t67;
t26 = t76 * t52 - t74 * t53;
t27 = t74 * t52 + t76 * t53;
t19 = Ifges(6,2) * t43 + Ifges(6,6) * qJD(5) + t94;
t2 = t9 * qJD(5) + t55 * t54;
t39 = Ifges(6,4) * t43;
t20 = Ifges(6,1) * t44 + Ifges(6,5) * qJD(5) + t39;
t42 = t64 * t71 + t81;
t78 = t20 * t111 + t42 * (t50 * mrSges(6,1) + t49 * mrSges(6,2)) + qJD(5) * (Ifges(6,5) * t49 - Ifges(6,6) * t50) / 0.2e1 + t19 * t110 + mrSges(5,3) * t84 + (-t55 * mrSges(6,1) + t56 * mrSges(6,2) + t80) * t83 + (t43 * t110 - t55 * t38) * Ifges(6,2) + (t49 * t100 + t56 * t37) * Ifges(6,1) + (-t10 * t50 + t2 * t55) * mrSges(6,3) + (-t50 * t100 + t43 * t111 - t108 + t109) * Ifges(6,4);
t62 = pkin(2) * t91 + qJD(4);
t59 = t64 - t97;
t57 = -t71 * pkin(3) + t81;
t33 = qJD(5) * mrSges(6,1) - t44 * mrSges(6,3);
t32 = -qJD(5) * mrSges(6,2) + t43 * mrSges(6,3);
t6 = -t27 * qJD(5) - t56 * t62;
t5 = t26 * qJD(5) + t55 * t62;
t1 = [t49 * t32 - t50 * t33 + m(6) * (t10 * t49 + t2 * t56 + t3 * t55 - t9 * t50) + (-t108 - t109) * mrSges(6,3); t78 + m(6) * (t10 * t5 + t2 * t27 + t3 * t26 + t9 * t6) + m(5) * (t79 * t62 + t63 * t84) + t59 * t18 + t5 * t32 + t6 * t33 - mrSges(4,2) * t82 - pkin(2) * t71 * t86 - mrSges(4,1) * t83 + t62 * t105 + (m(6) * (qJD(2) * t59 + t42) + m(5) * (t57 + (-pkin(3) - t97) * qJD(2)) - t71 * mrSges(4,1) + t103) * qJD(3) * t98 + (-t26 * t37 - t27 * t38 + t102) * mrSges(6,3); t78 + ((mrSges(4,1) * t75 + mrSges(4,2) * t77) * t92 + t81 * t104) * t71 + (-t86 + (-mrSges(4,1) * qJD(3) - t103) * t75) * t92 + t106 * t33 + t107 * t32 + (-t30 * t37 - t31 * t38 + t102) * mrSges(6,3) + t64 * t18 + (t107 * t10 + t106 * t9 + t2 * t31 + t3 * t30 - t42 * t88 + t64 * t83) * m(6) + (-pkin(3) * t83 + qJ(4) * t84 + t79 * qJD(4) - (t57 * t75 + t79 * t77) * t92) * m(5); -t43 * t32 + t44 * t33 + (m(5) + m(6)) * t83 - m(6) * (t10 * t43 - t9 * t44) + t18 + (-m(5) * t79 - t105) * t71; Ifges(6,5) * t37 - Ifges(6,6) * t38 - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t42 * (t44 * mrSges(6,1) + t43 * mrSges(6,2)) - t44 * (Ifges(6,1) * t43 - t94) / 0.2e1 + t19 * t100 - qJD(5) * (Ifges(6,5) * t43 - Ifges(6,6) * t44) / 0.2e1 - t9 * t32 + t10 * t33 + (t10 * t44 + t9 * t43) * mrSges(6,3) - (-Ifges(6,2) * t44 + t20 + t39) * t43 / 0.2e1;];
tauc = t1(:);
