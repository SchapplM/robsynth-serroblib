% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:13
% EndTime: 2022-01-23 09:12:16
% DurationCPUTime: 1.28s
% Computational Cost: add. (909->176), mult. (2294->256), div. (0->0), fcn. (1297->6), ass. (0->99)
t125 = -Ifges(5,4) - Ifges(6,4);
t124 = Ifges(5,1) + Ifges(6,1);
t101 = Ifges(6,5) + Ifges(5,5);
t123 = -Ifges(5,2) - Ifges(6,2);
t100 = Ifges(5,6) + Ifges(6,6);
t62 = cos(qJ(4));
t122 = t100 * t62;
t121 = t125 * t62;
t61 = sin(qJ(4));
t120 = t125 * t61;
t57 = sin(pkin(8));
t53 = t57 ^ 2;
t59 = cos(pkin(8));
t54 = t59 ^ 2;
t119 = (t53 + t54) * mrSges(4,3);
t41 = -cos(pkin(7)) * pkin(1) - pkin(3) * t59 - pkin(6) * t57 - pkin(2);
t103 = t59 * t62;
t51 = sin(pkin(7)) * pkin(1) + qJ(3);
t44 = t51 * t103;
t17 = t61 * t41 + t44;
t50 = -qJD(1) * t59 + qJD(4);
t87 = t57 * qJD(1);
t78 = t62 * t87;
t35 = mrSges(6,1) * t50 - mrSges(6,3) * t78;
t36 = mrSges(5,1) * t50 - mrSges(5,3) * t78;
t97 = t35 + t36;
t86 = t61 * qJD(1);
t80 = t57 * t86;
t33 = -mrSges(6,2) * t50 - mrSges(6,3) * t80;
t34 = -mrSges(5,2) * t50 - mrSges(5,3) * t80;
t98 = t33 + t34;
t118 = t97 * t61 - t98 * t62;
t110 = mrSges(6,2) * t61;
t111 = mrSges(6,1) * t62;
t113 = -t62 / 0.2e1;
t46 = t51 * qJD(1);
t52 = t59 * qJD(2);
t18 = qJD(5) - t52 + (pkin(4) * t86 + t46) * t57;
t29 = t46 * t57 - t52;
t25 = t41 * qJD(1) + qJD(3);
t30 = qJD(2) * t57 + t46 * t59;
t10 = t62 * t25 - t30 * t61;
t11 = t25 * t61 + t30 * t62;
t67 = t10 * t61 - t11 * t62;
t69 = mrSges(5,1) * t62 - mrSges(5,2) * t61;
t76 = qJ(5) * t87;
t8 = -t62 * t76 + t10;
t3 = pkin(4) * t50 + t8;
t9 = -t61 * t76 + t11;
t71 = t3 * t61 - t9 * t62;
t117 = t67 * mrSges(5,3) + t71 * mrSges(6,3) + t29 * t69 + t18 * (-t110 + t111) + (-(t124 * t62 + t120) * t61 / 0.2e1 + (t123 * t61 - t121) * t113) * t87 + (-t101 * t61 - t122 / 0.2e1 + t100 * t113) * t50;
t90 = qJD(4) * t62;
t42 = (pkin(4) * t90 + qJD(3)) * t57;
t40 = qJD(1) * t42;
t112 = m(6) * t40;
t105 = t29 * t57;
t104 = t59 * t61;
t102 = t62 * t57;
t92 = qJD(3) * t59;
t99 = t41 * t90 + t62 * t92;
t94 = qJ(5) * t57;
t93 = qJD(3) * t57;
t91 = qJD(4) * t61;
t89 = qJD(5) * t61;
t88 = qJD(5) * t62;
t85 = qJ(5) * qJD(4);
t84 = qJD(1) * qJD(3);
t83 = qJD(1) * qJD(4);
t82 = t51 * t104;
t81 = t62 * t94;
t79 = t61 * t92;
t77 = 0.3e1 / 0.2e1 * Ifges(5,4) + 0.3e1 / 0.2e1 * Ifges(6,4);
t75 = t59 * t84;
t74 = t57 * t83;
t73 = t53 * t51 * t84 + t29 * t93;
t4 = t25 * t90 - t30 * t91 + t62 * t75;
t1 = (-t62 * t85 - t89) * t87 + t4;
t66 = t11 * qJD(4);
t2 = -t66 + (-t79 + (t61 * t85 - t88) * t57) * qJD(1);
t72 = -t1 * t61 - t2 * t62;
t5 = -t61 * t75 - t66;
t70 = -t4 * t61 - t5 * t62;
t68 = mrSges(6,1) * t61 + mrSges(6,2) * t62;
t65 = t5 * mrSges(5,1) - t4 * mrSges(5,2) - t1 * mrSges(6,2);
t45 = t74 * t111;
t39 = (pkin(4) * t61 + t51) * t57;
t38 = (mrSges(5,1) * t61 + mrSges(5,2) * t62) * t87;
t37 = t68 * t87;
t32 = t62 * t41;
t27 = t69 * t74;
t26 = -t74 * t110 + t45;
t16 = t32 - t82;
t15 = -t61 * t94 + t17;
t14 = -t17 * qJD(4) - t79;
t13 = -qJD(4) * t82 + t99;
t12 = -t81 + t32 + (-t51 * t61 - pkin(4)) * t59;
t7 = -t79 - t57 * t88 + (-t44 + (-t41 + t94) * t61) * qJD(4);
t6 = -t57 * t89 + (-t81 - t82) * qJD(4) + t99;
t19 = [t13 * t34 + t14 * t36 + t39 * t26 + t6 * t33 + t7 * t35 + t42 * t37 + (-t2 * mrSges(6,1) - t65) * t59 + (t70 * mrSges(5,3) + t72 * mrSges(6,3) + qJD(3) * t38 + t51 * t27 + t40 * t68) * t57 + m(4) * (t30 * t92 + t73) + m(5) * (t10 * t14 + t11 * t13 + t16 * t5 + t17 * t4 + t73) + m(6) * (t1 * t15 + t12 * t2 + t18 * t42 + t3 * t7 + t39 * t40 + t6 * t9) + t117 * t57 * qJD(4) + ((m(4) * t51 * t54 + 0.2e1 * t119) * qJD(3) + ((mrSges(5,2) * t93 + (-t17 * mrSges(5,3) - t15 * mrSges(6,3) + t100 * t59 - t77 * t102) * qJD(4)) * t62 + (mrSges(5,1) * t93 + (t12 * mrSges(6,3) + t16 * mrSges(5,3) + t77 * t57 * t61 + t101 * t59 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(6,2)) * t102) * qJD(4)) * t61) * t57) * qJD(1); (-t26 - t27 - t112) * t59 + ((-t98 * t61 - t97 * t62) * qJD(4) + m(5) * (-t10 * t90 - t11 * t91 + t4 * t62 - t5 * t61 - t75) + m(6) * (t1 * t62 - t2 * t61 - t3 * t90 - t9 * t91)) * t57 + (mrSges(5,3) + mrSges(6,3)) * t53 * t83 * (-t61 ^ 2 - t62 ^ 2); -t118 * qJD(4) + m(5) * (-t67 * qJD(4) - t70) + m(6) * (-t71 * qJD(4) - t72) + ((-t37 - t38) * t57 + t118 * t59 - m(4) * (t30 * t59 + t105) - m(5) * (-t10 * t104 + t11 * t103 + t105) - m(6) * (t9 * t103 - t3 * t104 + t18 * t57) - qJD(1) * t119) * qJD(1); -t10 * t34 + t11 * t36 - t8 * t33 + (m(6) * pkin(4) + mrSges(6,1)) * t2 + (t35 - m(6) * (-t3 + t8)) * t9 + (((t123 * t62 + t120) * t61 / 0.2e1 + (-t124 * t61 + t121) * t113) * t87 + (-m(6) * t18 - t37) * t62 * pkin(4) + (-t122 + (mrSges(6,3) * pkin(4) - t101) * t61) * qJD(4) - t117) * t87 + t65; t45 + t112 + (-mrSges(6,2) * t91 + t62 * t35 + t61 * t33 - m(6) * (-t3 * t62 - t61 * t9)) * t87;];
tauc = t19(:);
