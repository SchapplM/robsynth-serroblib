% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR8
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:21
% DurationCPUTime: 0.80s
% Computational Cost: add. (1416->160), mult. (3110->246), div. (0->0), fcn. (2098->8), ass. (0->86)
t102 = pkin(2) * qJD(3);
t69 = sin(pkin(9));
t72 = sin(qJ(3));
t108 = t69 * t72;
t73 = sin(qJ(2));
t75 = cos(qJ(3));
t76 = cos(qJ(2));
t57 = t72 * t76 + t73 * t75;
t54 = t57 * qJD(1);
t56 = -t72 * t73 + t75 * t76;
t55 = t56 * qJD(1);
t70 = cos(pkin(9));
t105 = -t54 * t69 + t55 * t70 - (t70 * t75 - t108) * t102;
t107 = t70 * t72;
t104 = -t70 * t54 - t55 * t69 + (t69 * t75 + t107) * t102;
t62 = qJD(2) * pkin(2) + qJD(1) * t76;
t99 = qJD(1) * t73;
t40 = t75 * t62 - t72 * t99;
t68 = qJD(2) + qJD(3);
t33 = pkin(3) * t68 + t40;
t41 = t62 * t72 + t75 * t99;
t36 = t70 * t41;
t16 = t69 * t33 + t36;
t14 = pkin(7) * t68 + t16;
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t11 = qJD(4) * t74 - t14 * t71;
t12 = qJD(4) * t71 + t14 * t74;
t89 = -t11 * t71 + t12 * t74;
t67 = pkin(2) * t75 + pkin(3);
t120 = pkin(2) * t108 - t67 * t70;
t83 = t56 * qJD(2);
t97 = qJD(3) * t75;
t98 = qJD(3) * t72;
t21 = t62 * t97 + (-t73 * t98 + t83) * qJD(1);
t84 = t57 * qJD(2);
t22 = -t62 * t98 + (-t73 * t97 - t84) * qJD(1);
t7 = t21 * t70 + t22 * t69;
t2 = qJD(5) * t11 + t7 * t74;
t119 = t2 * t74;
t26 = -t56 * t70 + t57 * t69;
t6 = t21 * t69 - t70 * t22;
t118 = t26 * t6;
t117 = Ifges(6,1) * t71;
t116 = Ifges(6,4) * t71;
t114 = t11 * mrSges(6,3);
t111 = t41 * t69;
t109 = t68 * t74;
t106 = t71 * mrSges(6,3);
t103 = pkin(2) * t107 + t69 * t67;
t101 = Ifges(6,5) * qJD(5);
t100 = Ifges(6,6) * qJD(5);
t96 = qJD(5) * t12;
t95 = qJD(5) * t71;
t93 = qJD(5) * t74 / 0.2e1;
t91 = -mrSges(6,1) * t74 + mrSges(6,2) * t71;
t90 = -t11 * t74 - t12 * t71;
t15 = t33 * t70 - t111;
t58 = qJD(5) * mrSges(6,1) - t68 * t106;
t59 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t109;
t88 = -t58 * t74 - t59 * t71;
t87 = t71 * t58 - t74 * t59;
t86 = (Ifges(6,2) * t74 + t116) * t68;
t85 = (mrSges(6,1) * t71 + mrSges(6,2) * t74) * qJD(5);
t3 = -t7 * t71 - t96;
t80 = t90 * qJD(5) - t3 * t71 + t119;
t79 = m(6) * t80;
t13 = -pkin(4) * t68 - t15;
t42 = t86 + t100;
t63 = Ifges(6,4) * t109;
t43 = t68 * t117 + t101 + t63;
t78 = t22 * mrSges(4,1) - t21 * mrSges(4,2) - t7 * mrSges(5,2) + mrSges(6,3) * t119 + t13 * t85 + t43 * t93 + qJD(5) ^ 2 * (Ifges(6,5) * t74 - Ifges(6,6) * t71) / 0.2e1 - (t42 + t86) * t95 / 0.2e1 + (-mrSges(5,1) + t91) * t6 + ((Ifges(6,1) * t74 - t116) * t95 + (0.3e1 * Ifges(6,4) * t74 - 0.2e1 * Ifges(6,2) * t71 + t117) * t93) * t68;
t66 = -pkin(3) * t70 - pkin(4);
t65 = pkin(3) * t69 + pkin(7);
t52 = t91 * t68;
t49 = pkin(7) + t103;
t48 = -pkin(4) + t120;
t44 = t68 * t85;
t29 = -t57 * qJD(3) - t84;
t28 = t56 * qJD(3) + t83;
t27 = t56 * t69 + t57 * t70;
t18 = t40 * t70 - t111;
t17 = t40 * t69 + t36;
t10 = t28 * t70 + t29 * t69;
t9 = t28 * t69 - t70 * t29;
t1 = [t26 * t44 + t9 * t52 + (-mrSges(3,1) * t73 - mrSges(3,2) * t76) * qJD(2) ^ 2 - t87 * t10 + t88 * t27 * qJD(5) + m(4) * (t21 * t57 + t22 * t56 + t28 * t41 + t29 * t40) + m(5) * (t10 * t16 - t15 * t9 + t27 * t7 + t118) + m(6) * (t89 * t10 + t13 * t9 + t80 * t27 + t118) + (mrSges(4,1) * t29 - mrSges(5,1) * t9 - mrSges(4,2) * t28 - mrSges(5,2) * t10) * t68; (t105 * mrSges(5,2) + (-pkin(2) * t97 + t55) * mrSges(4,2) - t104 * mrSges(5,1) + (-pkin(2) * t98 + t54) * mrSges(4,1)) * t68 + (-qJD(5) * t49 * t59 + t105 * t58 + (-t3 - t96) * mrSges(6,3)) * t71 + (-t105 * t59 + (-t49 * t58 - t114) * qJD(5)) * t74 + t104 * t52 + t48 * t44 + t78 + t49 * t79 + (t104 * t13 - t105 * t89 + t6 * t48) * m(6) + (t40 * t54 - t41 * t55 + (t21 * t72 + t22 * t75 - t40 * t98 + t41 * t97) * pkin(2)) * m(4) + (t7 * t103 - t104 * t15 - t105 * t16 + t120 * t6) * m(5); (t90 * mrSges(6,3) + t88 * t65) * qJD(5) + t87 * t18 + t65 * t79 + (mrSges(4,1) * t41 + mrSges(5,1) * t17 + mrSges(4,2) * t40 + mrSges(5,2) * t18) * t68 - t3 * t106 + t66 * t44 - t17 * t52 + t78 + (-t13 * t17 - t89 * t18 + t6 * t66) * m(6) + ((-t6 * t70 + t69 * t7) * pkin(3) + t15 * t17 - t16 * t18) * m(5); m(6) * (t2 * t71 + t3 * t74) + (m(6) * t89 + (-t71 ^ 2 - t74 ^ 2) * t68 * mrSges(6,3) - t87) * qJD(5); t3 * mrSges(6,1) - t2 * mrSges(6,2) - t11 * t59 + t12 * t58 + ((t101 / 0.2e1 - t13 * mrSges(6,2) - t43 / 0.2e1 - t63 / 0.2e1 + t114) * t74 + (-t100 / 0.2e1 - t13 * mrSges(6,1) + t42 / 0.2e1 + t12 * mrSges(6,3) + (t116 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t74) * t68) * t71) * t68;];
tauc = t1(:);
