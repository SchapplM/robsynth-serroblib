% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:13
% EndTime: 2019-07-18 13:30:14
% DurationCPUTime: 0.72s
% Computational Cost: add. (1096->115), mult. (2793->156), div. (0->0), fcn. (1846->6), ass. (0->84)
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t125 = (-Ifges(6,2) + Ifges(6,1)) * t60;
t68 = Ifges(6,4) * t57 - t125;
t136 = t68 * t57;
t122 = m(6) / 0.2e1;
t58 = sin(qJ(4));
t110 = cos(qJ(3));
t79 = t110 * pkin(2) + pkin(3);
t59 = sin(qJ(3));
t61 = cos(qJ(4));
t98 = t59 * t61;
t36 = pkin(2) * t98 + t58 * t79;
t32 = pkin(6) + t36;
t55 = t57 ^ 2;
t56 = t60 ^ 2;
t95 = t56 + t55;
t130 = t32 * t95;
t72 = -t36 + t130;
t135 = -mrSges(5,2) / 0.2e1 + t72 * t122;
t124 = t95 * mrSges(6,3);
t77 = -mrSges(5,2) + t124;
t87 = m(6) * t95;
t80 = pkin(6) * t87;
t132 = t80 + t77;
t45 = -mrSges(6,1) * t60 + mrSges(6,2) * t57;
t96 = t45 - mrSges(5,1);
t128 = (-t56 / 0.2e1 - t55 / 0.2e1) * mrSges(6,3);
t63 = -Ifges(6,4) * t56 + t136;
t127 = qJD(4) * t63;
t126 = t63 * qJD(5);
t89 = t45 / 0.2e1 - mrSges(5,1) / 0.2e1;
t100 = t58 * t59;
t40 = (t110 * t61 - t100) * pkin(2);
t123 = t77 * t40 + (-t59 * mrSges(4,1) - t110 * mrSges(4,2)) * pkin(2);
t119 = mrSges(5,2) / 0.2e1;
t35 = pkin(2) * t100 - t61 * t79;
t102 = t57 * mrSges(6,1);
t97 = t60 * mrSges(6,2);
t75 = t97 + t102;
t20 = t35 * t75;
t118 = -t20 / 0.2e1;
t117 = t35 / 0.2e1;
t116 = -t40 / 0.2e1;
t112 = pkin(3) * t61;
t111 = t58 * pkin(3);
t108 = Ifges(6,4) * t60;
t39 = (t110 * t58 + t98) * pkin(2);
t106 = t35 * t39;
t105 = t39 * mrSges(5,1);
t104 = t39 * t45;
t103 = t39 * t61;
t85 = t96 * t36;
t3 = t85 + (-m(6) * t72 - t77) * t35;
t94 = t3 * qJD(2);
t4 = t96 * t39 + m(6) * (t40 * t130 + t106) + m(5) * (t36 * t40 + t106) + t123;
t93 = t4 * qJD(2);
t12 = -t20 + t63;
t92 = t12 * qJD(2);
t84 = t96 * t58;
t53 = pkin(6) + t111;
t83 = t95 * t53;
t82 = Ifges(6,5) * t60 - Ifges(6,6) * t57;
t81 = t108 * t60 - t136;
t66 = t75 * t112;
t78 = -t66 / 0.2e1 + t81;
t13 = pkin(3) * t84 + (t77 + (-t111 + t83) * m(6)) * t112;
t67 = t119 - t53 * t87 / 0.2e1;
t2 = (t119 - t80 / 0.2e1) * t40 + t67 * t35 + ((t35 * t122 + t89) * t58 + t135 * t61) * pkin(3) + (t112 / 0.2e1 - t35 / 0.2e1 + t116) * t124 + t96 * (t36 / 0.2e1 - t39 / 0.2e1);
t74 = t2 * qJD(2) + t13 * qJD(3);
t62 = (-t56 + t55) * Ifges(6,4) - t57 * t125;
t14 = t66 + t62;
t70 = t97 / 0.2e1 + t102 / 0.2e1;
t65 = t70 * t40;
t6 = t66 / 0.2e1 + t118 - t65 + t62;
t73 = -t6 * qJD(2) - t14 * qJD(3);
t8 = t118 + (mrSges(6,2) * t117 - t108) * t60 + (mrSges(6,1) * t117 + t68) * t57;
t64 = t8 * qJD(2) + qJD(3) * t63 + t127;
t18 = t20 / 0.2e1;
t11 = -t70 * t112 + t78;
t9 = t70 * t35 + t18 + t81;
t7 = t18 - t65 + t78;
t1 = t104 / 0.2e1 + mrSges(5,2) * t116 - t105 / 0.2e1 + t89 * t36 + (t67 + t128) * t35 + ((m(6) * t117 + t89) * t58 + (-t128 + t135) * t61) * pkin(3) + (t80 + t124) * t40 / 0.2e1;
t5 = [0, 0, 0, 0, -t75 * qJD(5); 0, qJD(3) * t4 + qJD(4) * t3 - qJD(5) * t12, t1 * qJD(4) + t7 * qJD(5) + t93 + (t104 - t105 + 0.2e1 * (-pkin(3) * t103 + t40 * t83) * t122 + m(5) * (t40 * t58 - t103) * pkin(3) + t123) * qJD(3), t1 * qJD(3) + t9 * qJD(5) + t94 + (-t132 * t35 + t85) * qJD(4), -t92 + t7 * qJD(3) + t9 * qJD(4) + (t45 * t32 + t82) * qJD(5); 0, qJD(4) * t2 - qJD(5) * t6 - t93, qJD(4) * t13 - qJD(5) * t14, t11 * qJD(5) + (t132 * t61 + t84) * qJD(4) * pkin(3) + t74, t11 * qJD(4) + (t45 * t53 + t82) * qJD(5) + t73; 0, -qJD(3) * t2 - qJD(5) * t8 - t94, -t74 - t126, -t126, (t45 * pkin(6) + t82) * qJD(5) - t64; 0, qJD(3) * t6 + qJD(4) * t8 + t92, -t73 + t127, t64, 0;];
Cq  = t5;
