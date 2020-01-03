% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:03
% EndTime: 2020-01-03 11:52:05
% DurationCPUTime: 0.80s
% Computational Cost: add. (1854->136), mult. (3992->202), div. (0->0), fcn. (2164->8), ass. (0->78)
t62 = sin(qJ(5));
t65 = cos(qJ(5));
t114 = mrSges(6,1) * t65 - mrSges(6,2) * t62 + mrSges(5,1);
t59 = qJD(1) + qJD(3);
t58 = qJD(4) + t59;
t75 = (mrSges(6,1) * t62 + mrSges(6,2) * t65) * qJD(5);
t36 = t58 * t75;
t55 = cos(pkin(9)) * pkin(1) + pkin(2);
t52 = t55 * qJD(1);
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t109 = pkin(1) * sin(pkin(9));
t90 = qJD(1) * t109;
t34 = t67 * t52 - t64 * t90;
t24 = pkin(3) * t59 + t34;
t63 = sin(qJ(4));
t35 = t64 * t52 + t67 * t90;
t66 = cos(qJ(4));
t97 = t66 * t35;
t16 = t24 * t63 + t97;
t29 = t34 * qJD(3);
t30 = t35 * qJD(3);
t7 = qJD(4) * t16 + t29 * t63 + t66 * t30;
t113 = m(6) * t7 + t36;
t14 = pkin(8) * t58 + t16;
t11 = qJD(2) * t65 - t14 * t62;
t12 = qJD(2) * t62 + t14 * t65;
t81 = -t11 * t62 + t12 * t65;
t112 = m(6) * t81;
t85 = -t64 * t109 + t67 * t55;
t46 = pkin(3) + t85;
t47 = t67 * t109 + t64 * t55;
t96 = t63 * t46 + t66 * t47;
t98 = t62 * mrSges(6,3);
t49 = qJD(5) * mrSges(6,1) - t58 * t98;
t100 = t58 * t65;
t50 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t100;
t78 = -t62 * t49 + t65 * t50;
t74 = t58 * mrSges(5,2) - t78;
t111 = m(5) * t16 + t112 - t74;
t101 = t35 * t63;
t15 = t24 * t66 - t101;
t6 = qJD(4) * t15 + t29 * t66 - t30 * t63;
t2 = qJD(5) * t11 + t6 * t65;
t108 = t2 * t65;
t3 = -qJD(5) * t12 - t6 * t62;
t107 = t3 * t62;
t106 = Ifges(6,1) * t62;
t105 = Ifges(6,4) * t62;
t99 = t59 * mrSges(4,1);
t95 = Ifges(6,5) * qJD(5);
t94 = Ifges(6,6) * qJD(5);
t93 = qJD(5) * t62;
t92 = qJD(5) * t65;
t88 = t92 / 0.2e1;
t87 = t114 * t58;
t82 = -t11 * t65 - t12 * t62;
t80 = t46 * t66 - t47 * t63;
t79 = -t65 * t49 - t62 * t50;
t77 = (Ifges(6,2) * t65 + t105) * t58;
t76 = t82 * mrSges(6,3);
t73 = t82 * qJD(5) - t107;
t70 = m(6) * (t73 + t108);
t13 = -pkin(4) * t58 - t15;
t32 = t77 + t94;
t54 = Ifges(6,4) * t100;
t33 = t58 * t106 + t54 + t95;
t69 = -t6 * mrSges(5,2) + mrSges(6,3) * t108 + t13 * t75 + t33 * t88 + qJD(5) ^ 2 * (Ifges(6,5) * t65 - Ifges(6,6) * t62) / 0.2e1 - (t32 + t77) * t93 / 0.2e1 - t114 * t7 + ((Ifges(6,1) * t65 - t105) * t93 + (0.3e1 * Ifges(6,4) * t65 - 0.2e1 * Ifges(6,2) * t62 + t106) * t88) * t58;
t68 = -t30 * mrSges(4,1) - t3 * t98 + t69;
t56 = pkin(3) * t63 + pkin(8);
t44 = t47 * qJD(3);
t43 = t85 * qJD(3);
t20 = pkin(8) + t96;
t19 = -pkin(4) - t80;
t18 = t34 * t66 - t101;
t17 = t34 * t63 + t97;
t9 = t96 * qJD(4) + t43 * t63 + t66 * t44;
t1 = [m(6) * (t13 * t9 + t19 * t7) + (t79 * t20 + t76) * qJD(5) - t87 * t9 + t20 * t70 + m(4) * (t29 * t47 - t30 * t85 - t34 * t44 + t35 * t43) + m(5) * (-t15 * t9 + t6 * t96 - t7 * t80) + (-t43 * t59 - t29) * mrSges(4,2) - t44 * t99 + t19 * t36 + t68 + t111 * (t80 * qJD(4) + t43 * t66 - t44 * t63); m(6) * (t2 * t62 + t3 * t65) + (t112 + (-t62 ^ 2 - t65 ^ 2) * t58 * mrSges(6,3) + t78) * qJD(5); (t79 * t56 + t76) * qJD(5) + t74 * t18 + t87 * t17 + t56 * t70 - m(5) * (-t15 * t17 + t16 * t18) + (t34 * t59 - t29) * mrSges(4,2) + t35 * t99 + (m(5) * (t6 * t63 - t66 * t7) + ((-m(5) * t15 + m(6) * t13 - t87) * t63 + t111 * t66) * qJD(4)) * pkin(3) + t68 - m(6) * (t13 * t17 + t81 * t18) + t113 * (-pkin(3) * t66 - pkin(4)); t73 * mrSges(6,3) + t87 * t16 + t74 * t15 + (m(6) * (-t11 * t92 - t12 * t93 - t107 + t108) - t50 * t93 - t49 * t92) * pkin(8) + t69 - m(6) * (t13 * t16 + t81 * t15) - t113 * pkin(4); t3 * mrSges(6,1) - t2 * mrSges(6,2) - t11 * t50 + t12 * t49 + ((t95 / 0.2e1 - t13 * mrSges(6,2) - t33 / 0.2e1 - t54 / 0.2e1 + t11 * mrSges(6,3)) * t65 + (-t94 / 0.2e1 - t13 * mrSges(6,1) + t32 / 0.2e1 + t12 * mrSges(6,3) + (t105 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t65) * t58) * t62) * t58;];
tauc = t1(:);
