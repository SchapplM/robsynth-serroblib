% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:09
% EndTime: 2019-12-31 17:22:10
% DurationCPUTime: 0.73s
% Computational Cost: add. (1193->129), mult. (2892->182), div. (0->0), fcn. (1903->6), ass. (0->88)
t64 = cos(qJ(3));
t60 = sin(qJ(4));
t58 = t60 ^ 2;
t63 = cos(qJ(4));
t59 = t63 ^ 2;
t96 = t58 + t59;
t134 = t96 * t64;
t61 = sin(qJ(3));
t62 = sin(qJ(2));
t102 = t61 * t62;
t65 = cos(qJ(2));
t55 = t65 * pkin(1) + pkin(2);
t36 = -pkin(1) * t102 + t64 * t55;
t31 = -pkin(3) - t36;
t100 = t62 * t64;
t37 = pkin(1) * t100 + t61 * t55;
t32 = pkin(7) + t37;
t137 = t134 * t32 + t31 * t61;
t136 = -Ifges(5,4) * t60 + (Ifges(5,1) - Ifges(5,2)) * t63;
t126 = t96 * mrSges(5,3);
t132 = -mrSges(4,2) + t126;
t117 = t64 * pkin(2);
t54 = -pkin(3) - t117;
t118 = t61 * pkin(2);
t53 = pkin(7) + t118;
t89 = t96 * t53;
t131 = t36 * t89 + t54 * t37;
t46 = -t63 * mrSges(5,1) + t60 * mrSges(5,2);
t128 = -mrSges(4,1) + t46;
t40 = (t64 * t65 - t102) * pkin(1);
t125 = t132 * t40 + (-t62 * mrSges(3,1) - t65 * mrSges(3,2)) * pkin(1);
t124 = m(5) / 0.2e1;
t123 = pkin(3) / 0.2e1;
t122 = -t31 / 0.2e1;
t121 = -t40 / 0.2e1;
t120 = -t54 / 0.2e1;
t114 = t36 * mrSges(4,2);
t112 = t37 * mrSges(4,1);
t111 = t37 * t46;
t39 = (t61 * t65 + t100) * pkin(1);
t110 = t39 * mrSges(4,1);
t109 = t39 * t46;
t105 = t60 * mrSges(5,1);
t104 = t61 * mrSges(4,1);
t103 = t61 * t46;
t99 = t63 * mrSges(5,2);
t98 = t64 * mrSges(4,2);
t91 = t96 * t36;
t67 = mrSges(5,3) * t91 + t111 - t112 - t114;
t3 = m(5) * (t31 * t37 + t32 * t91) + t67;
t95 = t3 * qJD(1);
t90 = t96 * t40;
t4 = t128 * t39 + m(5) * (t31 * t39 + t32 * t90) + m(4) * (-t36 * t39 + t37 * t40) + t125;
t94 = t4 * qJD(1);
t81 = t99 + t105;
t73 = t31 * t81;
t57 = Ifges(5,4) * t63;
t86 = t136 * t60 + t63 * t57;
t12 = t73 + t86;
t93 = t12 * qJD(1);
t87 = Ifges(5,5) * t63 - Ifges(5,6) * t60;
t85 = -t36 / 0.2e1 + t122 + t123;
t84 = mrSges(5,3) * (t59 / 0.2e1 + t58 / 0.2e1);
t83 = t121 + t122 + t120;
t74 = pkin(3) * t81;
t82 = -t74 / 0.2e1 + t86;
t66 = mrSges(5,3) * t134 + t103 - t104 - t98;
t13 = (m(5) * (t134 * t53 + t54 * t61) + t66) * pkin(2);
t68 = m(5) * (-pkin(3) * t39 + pkin(7) * t90);
t2 = (t137 * pkin(2) + t131) * t124 - t68 / 0.2e1 + t128 * (t37 / 0.2e1 + t118 / 0.2e1 - t39 / 0.2e1) + t132 * (t117 / 0.2e1 + t36 / 0.2e1 + t121);
t80 = t2 * qJD(1) + t13 * qJD(2);
t72 = t54 * t81;
t14 = t72 + t86;
t6 = (t83 * mrSges(5,2) - t57) * t63 + (t83 * mrSges(5,1) - t136) * t60;
t79 = -t6 * qJD(1) + t14 * qJD(2);
t78 = -t117 / 0.2e1 + t123 + t120;
t75 = -t99 / 0.2e1 - t105 / 0.2e1;
t10 = (t78 * mrSges(5,2) - t57) * t63 + (t78 * mrSges(5,1) - t136) * t60;
t15 = -t74 + t86;
t8 = (t85 * mrSges(5,2) - t57) * t63 + (t85 * mrSges(5,1) - t136) * t60;
t69 = t8 * qJD(1) + t10 * qJD(2) - t15 * qJD(3);
t35 = t72 / 0.2e1;
t18 = t73 / 0.2e1;
t11 = t75 * t117 + t35 + t82;
t9 = t75 * t36 + t18 + t82;
t7 = t75 * t40 + t18 + t35 + t86;
t1 = t131 * t124 + t111 / 0.2e1 - t114 / 0.2e1 - t112 / 0.2e1 + t68 / 0.2e1 + t109 / 0.2e1 + mrSges(4,2) * t121 - t110 / 0.2e1 + t36 * t84 + (t137 * t124 + t103 / 0.2e1 - t98 / 0.2e1 - t104 / 0.2e1 + t64 * t84) * pkin(2) + t40 * t126 / 0.2e1;
t5 = [t4 * qJD(2) + t3 * qJD(3) + t12 * qJD(4), t1 * qJD(3) + t7 * qJD(4) + t94 + (t109 - t110 + 0.2e1 * (t54 * t39 + t40 * t89) * t124 + m(4) * (-t39 * t64 + t40 * t61) * pkin(2) + t125) * qJD(2), t95 + t1 * qJD(2) + (m(5) * (-pkin(3) * t37 + pkin(7) * t91) + t67) * qJD(3) + t9 * qJD(4), t93 + t7 * qJD(2) + t9 * qJD(3) + (t46 * t32 + t87) * qJD(4); t2 * qJD(3) - t6 * qJD(4) - t94, t13 * qJD(3) + t14 * qJD(4), t11 * qJD(4) + (m(5) * (-pkin(3) * t61 + pkin(7) * t134) + t66) * qJD(3) * pkin(2) + t80, t11 * qJD(3) + (t46 * t53 + t87) * qJD(4) + t79; -t2 * qJD(2) - t8 * qJD(4) - t95, -t10 * qJD(4) - t80, t15 * qJD(4), (t46 * pkin(7) + t87) * qJD(4) - t69; t6 * qJD(2) + t8 * qJD(3) - t93, t10 * qJD(3) - t79, t69, 0;];
Cq = t5;
