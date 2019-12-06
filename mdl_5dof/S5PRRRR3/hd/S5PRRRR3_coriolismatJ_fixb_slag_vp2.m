% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:08
% EndTime: 2019-12-05 17:06:10
% DurationCPUTime: 0.75s
% Computational Cost: add. (1194->121), mult. (2895->175), div. (0->0), fcn. (1905->6), ass. (0->90)
t62 = cos(qJ(5));
t135 = (Ifges(6,1) - Ifges(6,2)) * t62;
t59 = sin(qJ(5));
t142 = -Ifges(6,4) * t59 + t135;
t57 = t59 ^ 2;
t58 = t62 ^ 2;
t100 = t57 + t58;
t63 = cos(qJ(4));
t140 = t100 * t63;
t60 = sin(qJ(4));
t61 = sin(qJ(3));
t104 = t60 * t61;
t64 = cos(qJ(3));
t55 = pkin(2) * t64 + pkin(3);
t36 = -pkin(2) * t104 + t55 * t63;
t31 = -pkin(4) - t36;
t103 = t61 * t63;
t37 = pkin(2) * t103 + t55 * t60;
t32 = pkin(8) + t37;
t141 = t140 * t32 + t31 * t60;
t129 = t100 * mrSges(6,3);
t137 = -mrSges(5,2) + t129;
t123 = t63 * pkin(3);
t54 = -pkin(4) - t123;
t124 = t60 * pkin(3);
t53 = pkin(8) + t124;
t91 = t100 * t53;
t136 = t36 * t91 + t37 * t54;
t46 = -mrSges(6,1) * t62 + mrSges(6,2) * t59;
t132 = -mrSges(5,1) + t46;
t130 = (-t58 + t57) * Ifges(6,4) - t59 * t135;
t40 = (t63 * t64 - t104) * pkin(2);
t128 = t137 * t40 + (-mrSges(4,1) * t61 - mrSges(4,2) * t64) * pkin(2);
t127 = m(6) / 0.2e1;
t126 = -t40 / 0.2e1;
t119 = Ifges(6,4) * t62;
t117 = t36 * mrSges(5,2);
t115 = t37 * mrSges(5,1);
t114 = t37 * t46;
t39 = (t60 * t64 + t103) * pkin(2);
t112 = t39 * mrSges(5,1);
t111 = t39 * t46;
t108 = t59 * mrSges(6,1);
t106 = t60 * mrSges(5,1);
t105 = t60 * t46;
t102 = t62 * mrSges(6,2);
t101 = t63 * mrSges(5,2);
t93 = t100 * t36;
t67 = mrSges(6,3) * t93 + t114 - t115 - t117;
t3 = m(6) * (t31 * t37 + t32 * t93) + t67;
t99 = t3 * qJD(2);
t92 = t100 * t40;
t4 = t132 * t39 + m(6) * (t31 * t39 + t32 * t92) + m(5) * (-t36 * t39 + t37 * t40) + t128;
t98 = t4 * qJD(2);
t84 = t102 + t108;
t77 = t31 * t84;
t88 = t119 * t62 + t142 * t59;
t12 = t77 + t88;
t97 = t12 * qJD(2);
t89 = Ifges(6,5) * t62 - Ifges(6,6) * t59;
t87 = mrSges(6,3) * (t58 / 0.2e1 + t57 / 0.2e1);
t78 = pkin(4) * t84;
t85 = -t78 / 0.2e1 + t88;
t66 = mrSges(6,3) * t140 - t101 + t105 - t106;
t13 = (m(6) * (t140 * t53 + t54 * t60) + t66) * pkin(3);
t69 = m(6) * (-pkin(4) * t39 + pkin(8) * t92);
t2 = (t141 * pkin(3) + t136) * t127 - t69 / 0.2e1 + t132 * (t124 / 0.2e1 + t37 / 0.2e1 - t39 / 0.2e1) + t137 * (t123 / 0.2e1 + t36 / 0.2e1 + t126);
t83 = t2 * qJD(2) + t13 * qJD(3);
t76 = t54 * t84;
t14 = t76 + t88;
t65 = -t76 / 0.2e1 + t130;
t71 = -t77 / 0.2e1;
t79 = -t102 / 0.2e1 - t108 / 0.2e1;
t72 = t79 * t40;
t6 = t71 + t72 + t65;
t82 = -t6 * qJD(2) + t14 * qJD(3);
t74 = t78 / 0.2e1;
t73 = t79 * t36;
t68 = t79 * t123;
t10 = t74 + t68 + t65;
t15 = (pkin(4) * mrSges(6,2) - t119) * t62 + (pkin(4) * mrSges(6,1) - t142) * t59;
t8 = t71 + t73 + t74 + t130;
t70 = t8 * qJD(2) + t10 * qJD(3) + t15 * qJD(4);
t35 = t76 / 0.2e1;
t18 = t77 / 0.2e1;
t11 = t35 + t68 + t85;
t9 = t18 + t73 + t85;
t7 = t35 + t18 + t72 + t88;
t1 = t136 * t127 + t114 / 0.2e1 - t117 / 0.2e1 - t115 / 0.2e1 + t69 / 0.2e1 + t111 / 0.2e1 + mrSges(5,2) * t126 - t112 / 0.2e1 + t36 * t87 + (t105 / 0.2e1 - t101 / 0.2e1 - t106 / 0.2e1 + t141 * t127 + t63 * t87) * pkin(3) + t40 * t129 / 0.2e1;
t5 = [0, 0, 0, 0, -t84 * qJD(5); 0, qJD(3) * t4 + qJD(4) * t3 + qJD(5) * t12, t1 * qJD(4) + t7 * qJD(5) + t98 + (t111 - t112 + 0.2e1 * (t54 * t39 + t40 * t91) * t127 + m(5) * (-t39 * t63 + t40 * t60) * pkin(3) + t128) * qJD(3), t99 + t1 * qJD(3) + (m(6) * (-pkin(4) * t37 + pkin(8) * t93) + t67) * qJD(4) + t9 * qJD(5), t97 + t7 * qJD(3) + t9 * qJD(4) + (t46 * t32 + t89) * qJD(5); 0, qJD(4) * t2 - qJD(5) * t6 - t98, qJD(4) * t13 + qJD(5) * t14, t11 * qJD(5) + (m(6) * (-pkin(4) * t60 + pkin(8) * t140) + t66) * qJD(4) * pkin(3) + t83, t11 * qJD(4) + (t46 * t53 + t89) * qJD(5) + t82; 0, -qJD(3) * t2 - qJD(5) * t8 - t99, -qJD(5) * t10 - t83, -t15 * qJD(5), (t46 * pkin(8) + t89) * qJD(5) - t70; 0, qJD(3) * t6 + qJD(4) * t8 - t97, qJD(4) * t10 - t82, t70, 0;];
Cq = t5;
