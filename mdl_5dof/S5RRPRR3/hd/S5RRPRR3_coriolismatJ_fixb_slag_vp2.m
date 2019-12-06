% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:21
% EndTime: 2019-12-05 18:30:23
% DurationCPUTime: 0.80s
% Computational Cost: add. (2582->128), mult. (5384->181), div. (0->0), fcn. (4320->8), ass. (0->93)
t73 = cos(qJ(5));
t137 = (Ifges(6,1) - Ifges(6,2)) * t73;
t70 = sin(qJ(5));
t141 = -Ifges(6,4) * t70 + t137;
t66 = t70 ^ 2;
t67 = t73 ^ 2;
t116 = t66 + t67;
t133 = t116 * mrSges(6,3);
t138 = -mrSges(5,2) + t133;
t61 = -mrSges(6,1) * t73 + mrSges(6,2) * t70;
t135 = -mrSges(5,1) + t61;
t134 = (-t67 + t66) * Ifges(6,4) - t70 * t137;
t69 = cos(pkin(9));
t72 = sin(qJ(2));
t120 = t69 * t72;
t68 = sin(pkin(9));
t75 = cos(qJ(2));
t54 = (-t68 * t75 - t120) * pkin(1);
t121 = t68 * t72;
t55 = (t69 * t75 - t121) * pkin(1);
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t38 = t54 * t71 + t55 * t74;
t132 = t138 * t38 + t54 * mrSges(4,1) - t55 * mrSges(4,2) + (-t72 * mrSges(3,1) - t75 * mrSges(3,2)) * pkin(1);
t131 = m(6) / 0.2e1;
t130 = -t38 / 0.2e1;
t128 = pkin(2) * t68;
t127 = mrSges(6,1) * t70;
t126 = mrSges(6,2) * t73;
t124 = Ifges(6,4) * t73;
t37 = -t74 * t54 + t55 * t71;
t123 = t37 * mrSges(5,1);
t122 = t37 * t61;
t64 = pkin(1) * t75 + pkin(2);
t100 = -pkin(1) * t121 + t69 * t64;
t105 = t116 * t38;
t48 = pkin(3) + t100;
t52 = pkin(1) * t120 + t64 * t68;
t33 = t48 * t74 - t52 * t71;
t27 = -pkin(4) - t33;
t34 = t48 * t71 + t52 * t74;
t28 = pkin(8) + t34;
t3 = t135 * t37 + m(6) * (t28 * t105 + t27 * t37) + m(5) * (-t33 * t37 + t34 * t38) + m(4) * (t100 * t54 + t52 * t55) + t132;
t115 = t3 * qJD(1);
t106 = t116 * t33;
t79 = -t33 * mrSges(5,2) + mrSges(6,3) * t106 + t135 * t34;
t4 = m(6) * (t28 * t106 + t27 * t34) + t79;
t114 = t4 * qJD(1);
t101 = t124 * t73 + t141 * t70;
t94 = t126 + t127;
t89 = t27 * t94;
t15 = t89 + t101;
t113 = t15 * qJD(1);
t63 = pkin(2) * t69 + pkin(3);
t51 = -t71 * t128 + t63 * t74;
t109 = -t33 / 0.2e1 - t51 / 0.2e1;
t53 = t74 * t128 + t63 * t71;
t108 = t34 / 0.2e1 + t53 / 0.2e1;
t50 = pkin(8) + t53;
t104 = t116 * t50;
t103 = t116 * t51;
t102 = Ifges(6,5) * t73 - Ifges(6,6) * t70;
t90 = pkin(4) * t94;
t96 = -t90 / 0.2e1 + t101;
t49 = -pkin(4) - t51;
t77 = (t28 * t103 + t33 * t104 + t27 * t53 + t34 * t49) * t131;
t80 = m(6) * (-pkin(4) * t37 + pkin(8) * t105);
t2 = t77 - t80 / 0.2e1 + t135 * (-t37 / 0.2e1 + t108) + t138 * (t130 - t109);
t78 = -t51 * mrSges(5,2) + mrSges(6,3) * t103 + t135 * t53;
t6 = m(6) * (t50 * t103 + t49 * t53) + t78;
t93 = t2 * qJD(1) + t6 * qJD(2);
t88 = t49 * t94;
t16 = t88 + t101;
t76 = -t88 / 0.2e1 + t134;
t82 = -t89 / 0.2e1;
t91 = -t126 / 0.2e1 - t127 / 0.2e1;
t84 = t91 * t38;
t8 = t82 + t84 + t76;
t92 = -t8 * qJD(1) + t16 * qJD(2);
t86 = t90 / 0.2e1;
t85 = t91 * t33;
t83 = t91 * t51;
t10 = t82 + t85 + t86 + t134;
t12 = t86 + t83 + t76;
t17 = (mrSges(6,2) * pkin(4) - t124) * t73 + (mrSges(6,1) * pkin(4) - t141) * t70;
t81 = t10 * qJD(1) + t12 * qJD(2) + t17 * qJD(4);
t39 = t88 / 0.2e1;
t18 = t89 / 0.2e1;
t13 = t39 + t83 + t96;
t11 = t18 + t85 + t96;
t9 = t18 + t39 + t84 + t101;
t1 = t77 + t80 / 0.2e1 + t122 / 0.2e1 - t123 / 0.2e1 + t38 * t133 / 0.2e1 + t135 * t108 + (t33 + t51) * mrSges(6,3) * (t67 / 0.2e1 + t66 / 0.2e1) + (t130 + t109) * mrSges(5,2);
t5 = [qJD(2) * t3 + qJD(4) * t4 + qJD(5) * t15, t1 * qJD(4) + t9 * qJD(5) + t115 + (t122 - t123 + 0.2e1 * (t38 * t104 + t37 * t49) * t131 + m(5) * (-t37 * t51 + t38 * t53) + m(4) * (t54 * t69 + t55 * t68) * pkin(2) + t132) * qJD(2), 0, t114 + t1 * qJD(2) + (m(6) * (-pkin(4) * t34 + pkin(8) * t106) + t79) * qJD(4) + t11 * qJD(5), t113 + t9 * qJD(2) + t11 * qJD(4) + (t61 * t28 + t102) * qJD(5); qJD(4) * t2 - qJD(5) * t8 - t115, qJD(4) * t6 + qJD(5) * t16, 0, (m(6) * (-pkin(4) * t53 + pkin(8) * t103) + t78) * qJD(4) + t13 * qJD(5) + t93, t13 * qJD(4) + (t61 * t50 + t102) * qJD(5) + t92; 0, 0, 0, 0, -t94 * qJD(5); -qJD(2) * t2 - qJD(5) * t10 - t114, -qJD(5) * t12 - t93, 0, -t17 * qJD(5), (t61 * pkin(8) + t102) * qJD(5) - t81; qJD(2) * t8 + qJD(4) * t10 - t113, qJD(4) * t12 - t92, 0, t81, 0;];
Cq = t5;
