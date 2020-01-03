% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:08
% EndTime: 2019-12-31 16:56:12
% DurationCPUTime: 1.21s
% Computational Cost: add. (848->188), mult. (1969->283), div. (0->0), fcn. (970->4), ass. (0->99)
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t79 = qJD(3) * t46;
t47 = cos(qJ(3));
t82 = qJD(1) * t47;
t34 = -t44 * t82 + t79;
t32 = Ifges(5,4) * t34;
t35 = qJD(3) * t44 + t46 * t82;
t45 = sin(qJ(3));
t74 = t45 * qJD(1);
t43 = qJD(4) + t74;
t10 = Ifges(5,1) * t35 + Ifges(5,5) * t43 + t32;
t104 = t46 / 0.2e1;
t105 = -t44 / 0.2e1;
t108 = t35 / 0.2e1;
t48 = -pkin(1) - pkin(5);
t41 = t48 * qJD(1) + qJD(2);
t89 = t41 * t47;
t31 = -qJD(3) * pkin(3) - t89;
t38 = pkin(3) * t45 - pkin(6) * t47 + qJ(2);
t29 = t38 * qJD(1);
t90 = t41 * t45;
t30 = qJD(3) * pkin(6) + t90;
t11 = t29 * t46 - t30 * t44;
t12 = t29 * t44 + t30 * t46;
t57 = t11 * t46 + t12 * t44;
t96 = Ifges(5,4) * t46;
t61 = -Ifges(5,2) * t44 + t96;
t97 = Ifges(5,4) * t44;
t63 = Ifges(5,1) * t46 - t97;
t64 = mrSges(5,1) * t44 + mrSges(5,2) * t46;
t98 = Ifges(5,4) * t35;
t9 = Ifges(5,2) * t34 + Ifges(5,6) * t43 + t98;
t92 = Ifges(5,6) * t44;
t94 = Ifges(5,5) * t46;
t123 = -t57 * mrSges(5,3) + t10 * t104 + t9 * t105 + (-t92 + t94) * t43 / 0.2e1 + t63 * t108 + t61 * t34 / 0.2e1 + t31 * t64;
t100 = Ifges(4,4) * t45;
t118 = qJD(1) / 0.2e1;
t84 = Ifges(4,5) * qJD(3);
t122 = t84 / 0.2e1 + (t47 * Ifges(4,1) - t100) * t118 + t123;
t87 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t34 - mrSges(5,2) * t35 - mrSges(4,3) * t82;
t121 = m(5) * t31 - t87;
t120 = qJD(1) * ((m(3) + m(4)) * qJ(2) + mrSges(3,3));
t67 = pkin(3) * t47 + pkin(6) * t45;
t33 = t67 * qJD(3) + qJD(2);
t26 = t33 * qJD(1);
t78 = qJD(3) * t47;
t71 = t41 * t78;
t1 = t11 * qJD(4) + t26 * t44 + t46 * t71;
t2 = -t12 * qJD(4) + t26 * t46 - t44 * t71;
t66 = t1 * t46 - t2 * t44;
t72 = qJD(3) * qJD(4);
t75 = qJD(4) * t47;
t20 = t46 * t72 + (-t44 * t75 - t45 * t79) * qJD(1);
t80 = qJD(3) * t45;
t21 = -t44 * t72 + (t44 * t80 - t46 * t75) * qJD(1);
t117 = t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t20 + Ifges(5,6) * t21;
t22 = -mrSges(5,2) * t43 + mrSges(5,3) * t34;
t23 = mrSges(5,1) * t43 - mrSges(5,3) * t35;
t116 = -m(5) * t57 - t44 * t22 - t46 * t23;
t113 = t20 / 0.2e1;
t112 = t21 / 0.2e1;
t111 = -t34 / 0.2e1;
t109 = -t35 / 0.2e1;
t107 = -t43 / 0.2e1;
t99 = Ifges(4,4) * t47;
t88 = t45 * t48;
t86 = qJ(2) * mrSges(4,1);
t85 = qJ(2) * mrSges(4,2);
t83 = Ifges(4,6) * qJD(3);
t81 = qJD(3) * mrSges(4,2);
t77 = qJD(4) * t44;
t76 = qJD(4) * t46;
t73 = qJD(1) * qJD(2);
t70 = t48 * t78;
t69 = qJD(1) * t78;
t65 = mrSges(5,1) * t46 - mrSges(5,2) * t44;
t62 = Ifges(5,1) * t44 + t96;
t60 = Ifges(5,2) * t46 + t97;
t58 = Ifges(5,5) * t44 + Ifges(5,6) * t46;
t56 = t11 * t44 - t12 * t46;
t13 = mrSges(5,1) * t69 - mrSges(5,3) * t20;
t14 = -mrSges(5,2) * t69 + mrSges(5,3) * t21;
t55 = -t44 * t13 + t46 * t14;
t25 = t38 * t44 + t46 * t88;
t24 = t38 * t46 - t44 * t88;
t51 = -Ifges(5,3) * t43 - Ifges(5,6) * t34 - Ifges(5,5) * t35 - t11 * mrSges(5,1) + t12 * mrSges(5,2) + t83 / 0.2e1 + (-Ifges(4,2) * t45 + t99) * t118;
t42 = Ifges(5,3) * t69;
t39 = -mrSges(4,3) * t74 - t81;
t37 = t67 * qJD(1);
t36 = qJD(1) * (t45 * mrSges(4,1) + t47 * mrSges(4,2));
t16 = t37 * t44 + t46 * t89;
t15 = t37 * t46 - t44 * t89;
t7 = -t25 * qJD(4) + t33 * t46 - t44 * t70;
t6 = t24 * qJD(4) + t33 * t44 + t46 * t70;
t5 = -mrSges(5,1) * t21 + mrSges(5,2) * t20;
t4 = Ifges(5,1) * t20 + Ifges(5,4) * t21 + Ifges(5,5) * t69;
t3 = Ifges(5,4) * t20 + Ifges(5,2) * t21 + Ifges(5,6) * t69;
t8 = [m(5) * (t1 * t25 + t11 * t7 + t12 * t6 + t2 * t24) + t6 * t22 + t7 * t23 + t24 * t13 + t25 * t14 + (t36 + 0.2e1 * t120) * qJD(2) + (mrSges(4,1) * t73 + t42 / 0.2e1 + (-t84 / 0.2e1 + (-0.2e1 * t85 + 0.3e1 / 0.2e1 * t100) * qJD(1) + t121 * t48 - t122) * qJD(3) + t117) * t45 + (mrSges(4,2) * t73 + t63 * t113 + t61 * t112 + t4 * t104 + t3 * t105 - t48 * t5 + (-t1 * t44 - t2 * t46) * mrSges(5,3) + (-t46 * t9 / 0.2e1 + t10 * t105 + t58 * t107 + t60 * t111 + t62 * t109 + t31 * t65 + t56 * mrSges(5,3)) * qJD(4) + (t48 * t39 - t83 / 0.2e1 + (-m(5) * t48 + t64) * t90 + ((-0.3e1 / 0.2e1 * Ifges(4,4) + t94 / 0.2e1 - t92 / 0.2e1) * t47 + 0.2e1 * t86 + (Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(4,1)) * t45) * qJD(1) - t51) * qJD(3)) * t47; (-t5 + (-m(5) * t56 + t46 * t22 - t44 * t23 + t39) * qJD(3)) * t47 + (m(5) * (-t11 * t76 - t12 * t77 + t66 - t71) - t22 * t77 - t23 * t76 + t55 + t121 * qJD(3)) * t45 + (t116 - t36 - t120) * qJD(1); t62 * t113 + t60 * t112 + t3 * t104 + t44 * t4 / 0.2e1 - t16 * t22 - t15 * t23 - pkin(3) * t5 + t55 * pkin(6) + t66 * mrSges(5,3) + ((-t39 - t81) * t47 + ((-mrSges(4,1) - t65) * qJD(3) + t87) * t45) * t41 + (pkin(6) * t116 + t123) * qJD(4) + (((-t86 + t99 / 0.2e1) * qJD(1) + t51) * t47 + ((t85 - t100 / 0.2e1 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t47) * qJD(1) + t122) * t45 + (-Ifges(4,5) * t45 / 0.2e1 + (-Ifges(4,6) / 0.2e1 + t58 / 0.2e1) * t47) * qJD(3)) * qJD(1) + (-pkin(3) * t41 * t80 + pkin(6) * t66 - t11 * t15 - t12 * t16 - t31 * t90) * m(5); t42 - t31 * (mrSges(5,1) * t35 + mrSges(5,2) * t34) + (Ifges(5,1) * t34 - t98) * t109 + t9 * t108 + (Ifges(5,5) * t34 - Ifges(5,6) * t35) * t107 - t11 * t22 + t12 * t23 + (t11 * t34 + t12 * t35) * mrSges(5,3) + (-Ifges(5,2) * t35 + t10 + t32) * t111 + t117;];
tauc = t8(:);
