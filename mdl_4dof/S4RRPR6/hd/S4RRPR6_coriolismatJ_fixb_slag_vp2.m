% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:24
% EndTime: 2019-12-31 17:04:26
% DurationCPUTime: 0.92s
% Computational Cost: add. (2554->101), mult. (5013->147), div. (0->0), fcn. (5437->6), ass. (0->62)
t79 = sin(pkin(7));
t81 = sin(qJ(2));
t83 = cos(qJ(2));
t96 = cos(pkin(7));
t64 = -t79 * t81 + t83 * t96;
t65 = -t79 * t83 - t81 * t96;
t80 = sin(qJ(4));
t82 = cos(qJ(4));
t88 = t82 * t64 + t65 * t80;
t137 = (Ifges(5,1) - Ifges(5,2)) * t88;
t101 = -qJ(3) - pkin(5);
t71 = t101 * t81;
t72 = t101 * t83;
t121 = t96 * t71 + t79 * t72;
t125 = pkin(6) * t65 + t121;
t54 = t79 * t71 - t96 * t72;
t33 = -pkin(6) * t64 - t54;
t128 = t125 * t82 + t33 * t80;
t28 = t125 * t80 - t33 * t82;
t50 = t64 * t80 - t65 * t82;
t4 = -t28 * mrSges(5,1) - t128 * mrSges(5,2) + Ifges(5,5) * t88 - Ifges(5,6) * t50;
t136 = t4 * qJD(4);
t132 = -t65 * mrSges(4,1) + t64 * mrSges(4,2);
t126 = t50 * mrSges(5,1);
t92 = t126 / 0.2e1;
t102 = t50 * mrSges(5,3);
t77 = -pkin(2) * t83 - pkin(1);
t56 = -t64 * pkin(3) + t77;
t89 = t88 * mrSges(5,2) + t126;
t131 = t56 * t89;
t124 = t88 ^ 2;
t118 = t88 / 0.2e1;
t120 = m(5) / 0.2e1;
t119 = m(4) * pkin(2);
t117 = pkin(2) * t79;
t115 = t81 * pkin(2);
t104 = t88 * mrSges(5,3);
t5 = -t50 * t102 - t88 * t104 - m(5) * (-t128 * t50 + t28 * t88) - m(4) * (t121 * t65 + t54 * t64) + (-t64 ^ 2 - t65 ^ 2) * mrSges(4,3);
t100 = qJD(1) * t5;
t57 = -pkin(3) * t65 + t115;
t1 = t131 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t81) * t81 + (-mrSges(4,2) * t115 - Ifges(4,4) * t65) * t65 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t83 + (Ifges(3,1) - Ifges(3,2)) * t81) * t83 + (-t50 ^ 2 + t124) * Ifges(5,4) + (-mrSges(4,1) * t115 + Ifges(4,4) * t64 + (-Ifges(4,1) + Ifges(4,2)) * t65) * t64 + t50 * t137 + (m(4) * t115 + t132) * t77 + (m(5) * t56 - mrSges(5,1) * t88 + mrSges(5,2) * t50) * t57;
t99 = t1 * qJD(1);
t2 = -Ifges(5,4) * t124 + t28 * t102 - t131 + (-t28 * mrSges(5,3) + Ifges(5,4) * t50 - t137) * t50;
t98 = t2 * qJD(1);
t90 = t96 * pkin(2);
t76 = t90 + pkin(3);
t58 = -t117 * t80 + t76 * t82;
t59 = t117 * t82 + t76 * t80;
t93 = t119 / 0.2e1;
t84 = (-t50 * t58 + t59 * t88) * t120 + (t64 * t79 + t65 * t96) * t93;
t85 = t120 * t57 + t81 * t93;
t8 = t132 - t84 + t85 + t89;
t97 = t8 * qJD(1);
t11 = 0.2e1 * t118 * mrSges(5,2) + 0.2e1 * t92;
t95 = t11 * qJD(1);
t14 = -mrSges(5,1) * t59 - mrSges(5,2) * t58;
t94 = t14 * qJD(4);
t3 = (t118 - t88 / 0.2e1) * Ifges(5,5);
t86 = -t3 * qJD(1) - t14 * qJD(2);
t12 = t92 - t126 / 0.2e1;
t10 = t84 + t85;
t6 = [qJD(2) * t1 - qJD(3) * t5 - qJD(4) * t2, t99 + (-t58 * t104 - t59 * t102 + m(5) * (t128 * t59 - t28 * t58) - t54 * mrSges(4,1) + (t121 * t79 - t54 * t96) * t119 - t121 * mrSges(4,2) + Ifges(4,6) * t65 + Ifges(3,5) * t83 - Ifges(3,6) * t81 + Ifges(4,5) * t64 + (-t83 * mrSges(3,1) + t81 * mrSges(3,2)) * pkin(5) + (t117 * t65 - t64 * t90) * mrSges(4,3) + t4) * qJD(2) + t10 * qJD(3) + t136, qJD(2) * t10 + qJD(4) * t12 - t100, t4 * qJD(2) + t12 * qJD(3) + t136 - t98; -qJD(3) * t8 + qJD(4) * t3 - t99, t94, -t97, -t86 + t94; qJD(2) * t8 + qJD(4) * t11 + t100, t97, 0, t95; -qJD(2) * t3 - qJD(3) * t11 + t98, t86, -t95, 0;];
Cq = t6;
