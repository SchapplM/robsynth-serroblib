% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRP2
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:05
% EndTime: 2019-12-31 17:49:08
% DurationCPUTime: 1.11s
% Computational Cost: add. (995->171), mult. (2541->221), div. (0->0), fcn. (1642->6), ass. (0->81)
t104 = Ifges(5,1) + Ifges(6,1);
t103 = Ifges(6,4) + Ifges(5,5);
t57 = sin(pkin(8));
t59 = cos(pkin(8));
t70 = qJD(1) * (t57 ^ 2 + t59 ^ 2);
t106 = mrSges(4,3) * t70;
t105 = mrSges(6,1) + mrSges(5,1);
t82 = mrSges(6,2) + mrSges(5,3);
t81 = -Ifges(5,4) + Ifges(6,5);
t102 = Ifges(6,6) - Ifges(5,6);
t61 = sin(qJ(4));
t91 = cos(qJ(4));
t73 = t91 * t59;
t69 = qJD(1) * t73;
t75 = qJD(1) * t57;
t39 = t61 * t75 - t69;
t37 = Ifges(5,4) * t39;
t47 = t57 * t91 + t61 * t59;
t40 = t47 * qJD(1);
t89 = Ifges(6,5) * t39;
t101 = t103 * qJD(4) + t104 * t40 - t37 + t89;
t51 = sin(pkin(7)) * pkin(1) + qJ(3);
t48 = t51 * qJD(1);
t54 = t59 * qJD(2);
t76 = pkin(6) * qJD(1);
t23 = t54 + (-t48 - t76) * t57;
t74 = t91 * t23;
t32 = t57 * qJD(2) + t59 * t48;
t24 = t59 * t76 + t32;
t85 = t61 * t24;
t7 = t74 - t85;
t100 = qJD(5) - t7;
t99 = -t39 / 0.2e1;
t98 = t39 / 0.2e1;
t96 = t40 / 0.2e1;
t63 = t47 * qJD(3);
t8 = t61 * t23 + t24 * t91;
t3 = qJD(1) * t63 + qJD(4) * t8;
t92 = pkin(6) + t51;
t43 = t92 * t57;
t44 = t92 * t59;
t65 = -t43 * t91 - t61 * t44;
t94 = t65 * t3;
t64 = -t61 * t57 + t73;
t93 = t3 * t64;
t90 = Ifges(5,4) * t40;
t87 = t39 * mrSges(5,3);
t86 = t40 * mrSges(5,3);
t84 = -qJD(4) / 0.2e1;
t83 = qJD(4) / 0.2e1;
t80 = qJD(3) * t69 + qJD(4) * t74;
t28 = -mrSges(6,2) * t39 + qJD(4) * mrSges(6,3);
t79 = -qJD(4) * mrSges(5,2) + t28 - t87;
t78 = -mrSges(6,2) * t40 + t105 * qJD(4) - t86;
t72 = qJD(3) * t75;
t41 = t64 * qJD(4);
t33 = qJD(1) * t41;
t42 = t47 * qJD(4);
t34 = qJD(1) * t42;
t71 = t34 * mrSges(5,1) + t33 * mrSges(5,2);
t68 = -(-t48 * t57 + t54) * t57 + t32 * t59;
t66 = -cos(pkin(7)) * pkin(1) - pkin(3) * t59 - pkin(2);
t19 = -t61 * t43 + t44 * t91;
t38 = qJD(1) * t66 + qJD(3);
t36 = Ifges(6,5) * t40;
t30 = t34 * mrSges(6,1);
t21 = mrSges(6,1) * t39 - mrSges(6,3) * t40;
t20 = pkin(4) * t40 + qJ(5) * t39;
t15 = -Ifges(5,2) * t39 + Ifges(5,6) * qJD(4) + t90;
t14 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t39 + t36;
t13 = -pkin(4) * t64 - qJ(5) * t47 + t66;
t12 = pkin(4) * t42 - qJ(5) * t41 - qJD(5) * t47;
t11 = pkin(4) * t39 - qJ(5) * t40 + t38;
t10 = qJD(4) * t19 + t63;
t9 = qJD(3) * t64 + qJD(4) * t65;
t6 = pkin(4) * t34 - qJ(5) * t33 - qJD(5) * t40;
t5 = qJD(4) * qJ(5) + t8;
t4 = -qJD(4) * pkin(4) + t100;
t2 = (-qJD(4) * t24 - t72) * t61 + t80;
t1 = -t61 * t72 + (qJD(5) - t85) * qJD(4) + t80;
t16 = [t13 * t30 + m(6) * (t1 * t19 + t10 * t4 + t11 * t12 + t13 * t6 + t5 * t9 - t94) + m(5) * (-t10 * t7 + t19 * t2 + t8 * t9 - t94) + (-mrSges(6,3) * t13 + t104 * t47 - t64 * t81 - t65 * t82) * t33 + (-t6 * mrSges(6,3) + t82 * t3 + t81 * t34) * t47 - t78 * t10 + t79 * t9 + (m(4) * (t51 * t70 + t68) + 0.2e1 * t106) * qJD(3) + t66 * t71 + t12 * t21 - (t6 * mrSges(6,1) - t1 * mrSges(6,2) - t2 * mrSges(5,3) + (Ifges(5,2) + Ifges(6,3)) * t34) * t64 - t82 * t19 * t34 + (-t8 * mrSges(5,3) - t5 * mrSges(6,2) + t11 * mrSges(6,1) + t14 / 0.2e1 + t38 * mrSges(5,1) - t15 / 0.2e1 + Ifges(6,3) * t98 - Ifges(5,2) * t99 + t81 * t96 + t102 * t83) * t42 + (t101 / 0.2e1 + t38 * mrSges(5,2) + t4 * mrSges(6,2) - t7 * mrSges(5,3) - t11 * mrSges(6,3) + Ifges(5,4) * t99 + Ifges(6,5) * t98 + t103 * t83 + t104 * t96) * t41; -t78 * t42 + t79 * t41 + m(5) * (t2 * t47 + t41 * t8 - t42 * t7 - t93) + m(6) * (t1 * t47 + t4 * t42 + t41 * t5 - t93) + t82 * (-t33 * t64 - t47 * t34); -t33 * mrSges(6,3) + t30 + t78 * t40 + t79 * t39 - m(5) * (-t39 * t8 - t40 * t7) + t71 + (-m(4) * t68 - t106) * qJD(1) + (t5 * t39 - t4 * t40 + t6) * m(6); (Ifges(6,3) * t40 - t89) * t99 - t11 * (t40 * mrSges(6,1) + t39 * mrSges(6,3)) - t38 * (mrSges(5,1) * t40 - mrSges(5,2) * t39) + t15 * t96 - t20 * t21 + qJD(5) * t28 + t1 * mrSges(6,3) - t2 * mrSges(5,2) + t102 * t34 + t103 * t33 - t105 * t3 + (t78 + t86) * t8 + (-t79 - t87) * t7 + (-pkin(4) * t33 - qJ(5) * t34 + t39 * t4 + t40 * t5) * mrSges(6,2) + (t102 * t40 - t103 * t39) * t84 + (-t3 * pkin(4) + t1 * qJ(5) + t100 * t5 - t11 * t20 - t4 * t8) * m(6) + (-Ifges(5,2) * t40 + t101 - t37) * t98 - (-t104 * t39 + t14 + t36 - t90) * t40 / 0.2e1; t33 * mrSges(6,2) - qJD(4) * t28 + t40 * t21 + 0.2e1 * (t3 / 0.2e1 + t5 * t84 + t11 * t96) * m(6);];
tauc = t16(:);
