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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:34:30
% EndTime: 2022-01-23 09:34:31
% DurationCPUTime: 0.81s
% Computational Cost: add. (2340->133), mult. (4635->183), div. (0->0), fcn. (3627->8), ass. (0->92)
t66 = sin(qJ(5));
t63 = t66 ^ 2;
t69 = cos(qJ(5));
t64 = t69 ^ 2;
t103 = t63 + t64;
t70 = cos(qJ(4));
t142 = t103 * t70;
t124 = pkin(1) * sin(pkin(9));
t58 = cos(pkin(9)) * pkin(1) + pkin(2);
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t42 = t71 * t124 + t68 * t58;
t67 = sin(qJ(4));
t112 = t42 * t67;
t41 = -t68 * t124 + t71 * t58;
t40 = pkin(3) + t41;
t35 = t40 * t70 - t112;
t33 = -pkin(4) - t35;
t104 = t70 * t42;
t36 = t40 * t67 + t104;
t34 = pkin(8) + t36;
t143 = t142 * t34 + t33 * t67;
t141 = -Ifges(6,4) * t66 + (Ifges(6,1) - Ifges(6,2)) * t69;
t131 = t103 * mrSges(6,3);
t138 = -mrSges(5,2) + t131;
t122 = t70 * pkin(3);
t60 = -pkin(4) - t122;
t123 = t67 * pkin(3);
t59 = pkin(8) + t123;
t96 = t103 * t59;
t137 = t35 * t96 + t36 * t60;
t50 = -mrSges(6,1) * t69 + mrSges(6,2) * t66;
t134 = mrSges(5,1) - t50;
t38 = t41 * t70 - t112;
t133 = -t42 * mrSges(4,1) - t41 * mrSges(4,2) + t138 * t38;
t130 = m(6) / 0.2e1;
t129 = pkin(4) / 0.2e1;
t128 = -t33 / 0.2e1;
t127 = -t38 / 0.2e1;
t126 = -t60 / 0.2e1;
t119 = t35 * mrSges(5,2);
t117 = t36 * mrSges(5,1);
t116 = t36 * t50;
t37 = t41 * t67 + t104;
t114 = t37 * mrSges(5,1);
t113 = t37 * t50;
t109 = t66 * mrSges(6,1);
t108 = t67 * mrSges(5,1);
t107 = t67 * t50;
t106 = t69 * mrSges(6,2);
t105 = t70 * mrSges(5,2);
t97 = t103 * t38;
t3 = -t134 * t37 + m(6) * (t33 * t37 + t34 * t97) + m(5) * (-t35 * t37 + t36 * t38) + t133;
t102 = t3 * qJD(1);
t98 = t103 * t35;
t73 = mrSges(6,3) * t98 + t116 - t117 - t119;
t4 = m(6) * (t33 * t36 + t34 * t98) + t73;
t101 = t4 * qJD(1);
t87 = t106 + t109;
t79 = t33 * t87;
t62 = Ifges(6,4) * t69;
t93 = t141 * t66 + t69 * t62;
t11 = t79 + t93;
t100 = t11 * qJD(1);
t94 = Ifges(6,5) * t69 - Ifges(6,6) * t66;
t92 = -t35 / 0.2e1 + t128 + t129;
t91 = mrSges(6,3) * (t64 / 0.2e1 + t63 / 0.2e1);
t90 = t127 + t128 + t126;
t80 = pkin(4) * t87;
t89 = -t80 / 0.2e1 + t93;
t72 = mrSges(6,3) * t142 - t105 + t107 - t108;
t15 = (m(6) * (t142 * t59 + t60 * t67) + t72) * pkin(3);
t74 = m(6) * (-pkin(4) * t37 + pkin(8) * t97);
t2 = (t143 * pkin(3) + t137) * t130 - t74 / 0.2e1 + t134 * (-t123 / 0.2e1 - t36 / 0.2e1 + t37 / 0.2e1) + t138 * (t122 / 0.2e1 + t35 / 0.2e1 + t127);
t86 = t2 * qJD(1) + t15 * qJD(3);
t78 = t60 * t87;
t16 = t78 + t93;
t6 = (t90 * mrSges(6,2) - t62) * t69 + (t90 * mrSges(6,1) - t141) * t66;
t85 = -t6 * qJD(1) + t16 * qJD(3);
t84 = -t122 / 0.2e1 + t129 + t126;
t81 = -t106 / 0.2e1 - t109 / 0.2e1;
t12 = (t84 * mrSges(6,2) - t62) * t69 + (t84 * mrSges(6,1) - t141) * t66;
t17 = -t80 + t93;
t8 = (t92 * mrSges(6,2) - t62) * t69 + (t92 * mrSges(6,1) - t141) * t66;
t75 = t8 * qJD(1) + t12 * qJD(3) - t17 * qJD(4);
t43 = t78 / 0.2e1;
t18 = t79 / 0.2e1;
t13 = t81 * t122 + t43 + t89;
t9 = t81 * t35 + t18 + t89;
t7 = t81 * t38 + t18 + t43 + t93;
t1 = t116 / 0.2e1 + t137 * t130 - t119 / 0.2e1 - t117 / 0.2e1 + t113 / 0.2e1 + t74 / 0.2e1 + mrSges(5,2) * t127 - t114 / 0.2e1 + t35 * t91 + (t143 * t130 - t105 / 0.2e1 - t108 / 0.2e1 + t107 / 0.2e1 + t70 * t91) * pkin(3) + t38 * t131 / 0.2e1;
t5 = [qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t11, 0, t1 * qJD(4) + t7 * qJD(5) + t102 + (t113 - t114 + 0.2e1 * (t60 * t37 + t38 * t96) * t130 + m(5) * (-t37 * t70 + t38 * t67) * pkin(3) + t133) * qJD(3), t101 + t1 * qJD(3) + (m(6) * (-pkin(4) * t36 + pkin(8) * t98) + t73) * qJD(4) + t9 * qJD(5), t100 + t7 * qJD(3) + t9 * qJD(4) + (t50 * t34 + t94) * qJD(5); 0, 0, 0, 0, -t87 * qJD(5); qJD(4) * t2 - qJD(5) * t6 - t102, 0, qJD(4) * t15 + qJD(5) * t16, t13 * qJD(5) + (m(6) * (-pkin(4) * t67 + pkin(8) * t142) + t72) * qJD(4) * pkin(3) + t86, t13 * qJD(4) + (t50 * t59 + t94) * qJD(5) + t85; -qJD(3) * t2 - qJD(5) * t8 - t101, 0, -qJD(5) * t12 - t86, t17 * qJD(5), (t50 * pkin(8) + t94) * qJD(5) - t75; qJD(3) * t6 + qJD(4) * t8 - t100, 0, qJD(4) * t12 - t85, t75, 0;];
Cq = t5;
