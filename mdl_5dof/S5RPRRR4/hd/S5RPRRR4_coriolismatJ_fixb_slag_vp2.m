% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:22
% EndTime: 2019-12-05 18:14:25
% DurationCPUTime: 0.84s
% Computational Cost: add. (2340->124), mult. (4635->175), div. (0->0), fcn. (3627->8), ass. (0->94)
t68 = cos(qJ(5));
t141 = (Ifges(6,1) - Ifges(6,2)) * t68;
t65 = sin(qJ(5));
t148 = -Ifges(6,4) * t65 + t141;
t62 = t65 ^ 2;
t63 = t68 ^ 2;
t107 = t62 + t63;
t69 = cos(qJ(4));
t146 = t107 * t69;
t130 = pkin(1) * sin(pkin(9));
t58 = cos(pkin(9)) * pkin(1) + pkin(2);
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t42 = t70 * t130 + t67 * t58;
t66 = sin(qJ(4));
t117 = t42 * t66;
t41 = -t67 * t130 + t70 * t58;
t40 = pkin(3) + t41;
t35 = t40 * t69 - t117;
t33 = -pkin(4) - t35;
t108 = t69 * t42;
t36 = t40 * t66 + t108;
t34 = pkin(8) + t36;
t147 = t146 * t34 + t33 * t66;
t134 = t107 * mrSges(6,3);
t143 = -mrSges(5,2) + t134;
t128 = pkin(3) * t69;
t60 = -pkin(4) - t128;
t129 = pkin(3) * t66;
t59 = pkin(8) + t129;
t98 = t107 * t59;
t142 = t35 * t98 + t36 * t60;
t50 = -mrSges(6,1) * t68 + mrSges(6,2) * t65;
t138 = -mrSges(5,1) + t50;
t38 = t41 * t69 - t117;
t137 = -t42 * mrSges(4,1) - t41 * mrSges(4,2) + t143 * t38;
t135 = (-t63 + t62) * Ifges(6,4) - t65 * t141;
t133 = m(6) / 0.2e1;
t132 = -t38 / 0.2e1;
t126 = Ifges(6,4) * t68;
t124 = t35 * mrSges(5,2);
t122 = t36 * mrSges(5,1);
t121 = t36 * t50;
t37 = t41 * t66 + t108;
t119 = t37 * mrSges(5,1);
t118 = t37 * t50;
t114 = t65 * mrSges(6,1);
t113 = t66 * mrSges(5,1);
t112 = t66 * t50;
t111 = t68 * mrSges(6,2);
t109 = t69 * mrSges(5,2);
t99 = t107 * t38;
t3 = t138 * t37 + m(6) * (t33 * t37 + t34 * t99) + m(5) * (-t35 * t37 + t36 * t38) + t137;
t106 = t3 * qJD(1);
t100 = t107 * t35;
t73 = mrSges(6,3) * t100 + t121 - t122 - t124;
t4 = m(6) * (t34 * t100 + t33 * t36) + t73;
t105 = t4 * qJD(1);
t90 = t111 + t114;
t83 = t33 * t90;
t95 = t126 * t68 + t148 * t65;
t11 = t83 + t95;
t104 = t11 * qJD(1);
t96 = Ifges(6,5) * t68 - Ifges(6,6) * t65;
t94 = mrSges(6,3) * (t62 / 0.2e1 + t63 / 0.2e1);
t84 = pkin(4) * t90;
t92 = -t84 / 0.2e1 + t95;
t72 = mrSges(6,3) * t146 - t109 + t112 - t113;
t15 = (m(6) * (t146 * t59 + t60 * t66) + t72) * pkin(3);
t75 = m(6) * (-pkin(4) * t37 + pkin(8) * t99);
t2 = (t147 * pkin(3) + t142) * t133 - t75 / 0.2e1 + t138 * (t129 / 0.2e1 + t36 / 0.2e1 - t37 / 0.2e1) + t143 * (t128 / 0.2e1 + t35 / 0.2e1 + t132);
t89 = t2 * qJD(1) + t15 * qJD(3);
t82 = t60 * t90;
t16 = t82 + t95;
t71 = -t82 / 0.2e1 + t135;
t77 = -t83 / 0.2e1;
t85 = t114 / 0.2e1 + t111 / 0.2e1;
t78 = t85 * t38;
t6 = t77 - t78 + t71;
t88 = -t6 * qJD(1) + t16 * qJD(3);
t80 = t84 / 0.2e1;
t79 = t85 * t35;
t74 = t85 * t128;
t12 = t80 - t74 + t71;
t17 = (pkin(4) * mrSges(6,2) - t126) * t68 + (pkin(4) * mrSges(6,1) - t148) * t65;
t8 = t77 - t79 + t80 + t135;
t76 = t8 * qJD(1) + t12 * qJD(3) + t17 * qJD(4);
t43 = t82 / 0.2e1;
t18 = t83 / 0.2e1;
t13 = t43 - t74 + t92;
t9 = t18 - t79 + t92;
t7 = t43 + t18 - t78 + t95;
t1 = t142 * t133 + t121 / 0.2e1 - t124 / 0.2e1 - t122 / 0.2e1 + t75 / 0.2e1 + t118 / 0.2e1 + mrSges(5,2) * t132 - t119 / 0.2e1 + t35 * t94 + (t112 / 0.2e1 + t147 * t133 - t109 / 0.2e1 - t113 / 0.2e1 + t69 * t94) * pkin(3) + t38 * t134 / 0.2e1;
t5 = [qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t11, 0, t1 * qJD(4) + t7 * qJD(5) + t106 + (t118 - t119 + 0.2e1 * (t60 * t37 + t38 * t98) * t133 + m(5) * (-t37 * t69 + t38 * t66) * pkin(3) + t137) * qJD(3), t105 + t1 * qJD(3) + (m(6) * (-pkin(4) * t36 + pkin(8) * t100) + t73) * qJD(4) + t9 * qJD(5), t104 + t7 * qJD(3) + t9 * qJD(4) + (t50 * t34 + t96) * qJD(5); 0, 0, 0, 0, -t90 * qJD(5); qJD(4) * t2 - qJD(5) * t6 - t106, 0, qJD(4) * t15 + qJD(5) * t16, t13 * qJD(5) + (m(6) * (-pkin(4) * t66 + pkin(8) * t146) + t72) * qJD(4) * pkin(3) + t89, t13 * qJD(4) + (t50 * t59 + t96) * qJD(5) + t88; -qJD(3) * t2 - qJD(5) * t8 - t105, 0, -qJD(5) * t12 - t89, -t17 * qJD(5), (t50 * pkin(8) + t96) * qJD(5) - t76; qJD(3) * t6 + qJD(4) * t8 - t104, 0, qJD(4) * t12 - t88, t76, 0;];
Cq = t5;
