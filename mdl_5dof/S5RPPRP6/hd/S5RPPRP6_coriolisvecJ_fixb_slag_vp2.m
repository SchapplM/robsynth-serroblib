% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:55:01
% DurationCPUTime: 1.26s
% Computational Cost: add. (1040->189), mult. (2393->237), div. (0->0), fcn. (1426->4), ass. (0->90)
t113 = Ifges(5,1) + Ifges(6,1);
t112 = Ifges(5,5) + Ifges(6,4);
t63 = sin(pkin(7));
t99 = sin(qJ(4));
t79 = t99 * t63;
t64 = cos(pkin(7));
t66 = cos(qJ(4));
t93 = t66 * t64;
t41 = t79 - t93;
t114 = -mrSges(6,1) - mrSges(5,1);
t90 = mrSges(6,2) + mrSges(5,3);
t89 = Ifges(5,4) - Ifges(6,5);
t111 = Ifges(6,6) - Ifges(5,6);
t42 = t66 * t63 + t64 * t99;
t36 = t42 * qJD(1);
t35 = Ifges(5,4) * t36;
t75 = qJD(1) * t79;
t37 = qJD(1) * t93 - t75;
t96 = Ifges(6,5) * t36;
t110 = t112 * qJD(4) + t113 * t37 - t35 + t96;
t109 = t41 * qJD(3);
t85 = t63 ^ 2 + t64 ^ 2;
t76 = qJD(1) * t85;
t65 = -pkin(1) - qJ(3);
t49 = qJD(1) * t65 + qJD(2);
t77 = -pkin(6) * qJD(1) + t49;
t33 = t77 * t64;
t32 = t77 * t63;
t80 = t99 * t32;
t11 = t66 * t33 - t80;
t108 = -t11 + qJD(5);
t107 = -t36 / 0.2e1;
t106 = t36 / 0.2e1;
t104 = t37 / 0.2e1;
t100 = -pkin(6) + t65;
t44 = t100 * t63;
t45 = t100 * t64;
t21 = t44 * t99 - t66 * t45;
t12 = t66 * t32 + t33 * t99;
t3 = -t109 * qJD(1) + qJD(4) * t12;
t102 = t21 * t3;
t101 = t3 * t41;
t59 = t63 * pkin(3);
t98 = (m(3) * qJ(2));
t97 = Ifges(5,4) * t37;
t95 = t36 * mrSges(5,3);
t94 = t37 * mrSges(5,3);
t92 = -qJD(4) / 0.2e1;
t91 = qJD(4) / 0.2e1;
t27 = -t36 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t88 = -qJD(4) * mrSges(5,2) + t27 - t95;
t87 = -t37 * mrSges(6,2) - t114 * qJD(4) - t94;
t73 = mrSges(4,1) * t63 + mrSges(4,2) * t64;
t86 = mrSges(5,1) * t36 + mrSges(5,2) * t37 + qJD(1) * t73;
t55 = qJ(2) + t59;
t84 = qJD(4) * t66;
t62 = qJD(1) * qJ(2);
t58 = qJD(3) + t62;
t83 = qJD(1) * qJD(2);
t82 = t64 * t84;
t46 = qJD(1) * t59 + t58;
t78 = qJD(4) * t99;
t22 = t66 * t44 + t45 * t99;
t71 = t42 * qJD(3);
t68 = -qJD(1) * t71 + t33 * t84;
t1 = (-t80 + qJD(5)) * qJD(4) + t68;
t38 = -t63 * t84 - t64 * t78;
t39 = -t63 * t78 + t82;
t6 = -qJD(4) * pkin(4) + t108;
t9 = qJD(4) * qJ(5) + t12;
t70 = t1 * t42 - t38 * t6 + t39 * t9 + t101;
t2 = -t32 * t78 + t68;
t69 = t11 * t38 + t12 * t39 + t2 * t42 + t101;
t67 = qJD(1) ^ 2;
t34 = Ifges(6,5) * t37;
t31 = qJD(1) * t82 - qJD(4) * t75;
t30 = qJD(4) * t36;
t29 = t31 * mrSges(6,1);
t28 = t30 * mrSges(5,2);
t19 = mrSges(6,1) * t36 - mrSges(6,3) * t37;
t18 = pkin(4) * t37 + qJ(5) * t36;
t17 = pkin(4) * t42 + qJ(5) * t41 + t55;
t14 = -Ifges(5,2) * t36 + Ifges(5,6) * qJD(4) + t97;
t13 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t36 + t34;
t10 = pkin(4) * t36 - qJ(5) * t37 + t46;
t8 = qJD(4) * t22 - t109;
t7 = -qJD(4) * t21 - t71;
t5 = pkin(4) * t39 - qJ(5) * t38 + qJD(5) * t41 + qJD(2);
t4 = pkin(4) * t31 + qJ(5) * t30 - qJD(5) * t37 + t83;
t15 = [t17 * t29 - t55 * t28 + m(6) * (t1 * t22 + t10 * t5 + t17 * t4 + t6 * t8 + t7 * t9 + t102) + m(5) * (-t11 * t8 + t12 * t7 + t2 * t22 + t102) - t87 * t8 + t88 * t7 - (-mrSges(6,3) * t17 - t113 * t41 + t90 * t21 - t89 * t42) * t30 + (mrSges(5,1) * t55 + (Ifges(5,2) + Ifges(6,3)) * t42 + t89 * t41 - t90 * t22) * t31 - t69 * mrSges(5,3) - t70 * mrSges(6,2) + t4 * (mrSges(6,1) * t42 + mrSges(6,3) * t41) + t5 * t19 + (0.2e1 * mrSges(4,3) * t76 + m(4) * (-t49 * t85 - t65 * t76)) * qJD(3) + ((mrSges(5,1) * t42 - mrSges(5,2) * t41 + (2 * mrSges(3,3)) + t73 + (2 * t98)) * qJD(1) + t86 + m(4) * (t58 + t62) + m(5) * (qJD(1) * t55 + t46)) * qJD(2) + (t46 * mrSges(5,1) + t10 * mrSges(6,1) + t13 / 0.2e1 - t14 / 0.2e1 + Ifges(6,3) * t106 - Ifges(5,2) * t107 + t111 * t91 - t89 * t104) * t39 + (t46 * mrSges(5,2) - t10 * mrSges(6,3) + Ifges(5,4) * t107 + Ifges(6,5) * t106 + t112 * t91 + t113 * t104 + t110 / 0.2e1) * t38; (-mrSges(3,3) - t98) * t67 + t88 * t39 + t87 * t38 + m(5) * t69 + m(6) * t70 + (-m(5) * t46 - m(6) * t10 - t19 + (-qJD(3) * t85 - t58) * m(4) - t86) * qJD(1) + t90 * (-t41 * t30 - t42 * t31); t31 * mrSges(5,1) + t30 * mrSges(6,3) - t28 + t29 + t87 * t37 + t88 * t36 + (m(4) + m(5)) * t83 - t85 * t67 * mrSges(4,3) - m(5) * (-t11 * t37 - t12 * t36) + m(4) * t49 * t76 + (t36 * t9 - t37 * t6 + t4) * m(6); -t46 * (mrSges(5,1) * t37 - mrSges(5,2) * t36) + (Ifges(6,3) * t37 - t96) * t107 - t10 * (t37 * mrSges(6,1) + t36 * mrSges(6,3)) + t14 * t104 + qJD(5) * t27 - t18 * t19 - t2 * mrSges(5,2) + t1 * mrSges(6,3) + t111 * t31 - t112 * t30 + t114 * t3 + (t87 + t94) * t12 + (-t88 - t95) * t11 + (pkin(4) * t30 - qJ(5) * t31 + t36 * t6 + t37 * t9) * mrSges(6,2) + (t111 * t37 - t112 * t36) * t92 + (-t3 * pkin(4) + t1 * qJ(5) - t10 * t18 + t108 * t9 - t6 * t12) * m(6) + (-Ifges(5,2) * t37 + t110 - t35) * t106 - (-t113 * t36 + t13 + t34 - t97) * t37 / 0.2e1; -t30 * mrSges(6,2) - qJD(4) * t27 + t37 * t19 + 0.2e1 * (t3 / 0.2e1 + t9 * t92 + t10 * t104) * m(6);];
tauc = t15(:);
