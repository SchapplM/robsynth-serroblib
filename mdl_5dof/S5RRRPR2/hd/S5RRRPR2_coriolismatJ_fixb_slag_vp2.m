% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:47
% EndTime: 2019-12-05 18:40:49
% DurationCPUTime: 1.12s
% Computational Cost: add. (2637->169), mult. (6035->242), div. (0->0), fcn. (4781->8), ass. (0->112)
t90 = sin(qJ(5));
t93 = cos(qJ(5));
t172 = -Ifges(6,4) * t90 + (Ifges(6,1) - Ifges(6,2)) * t93;
t86 = t90 ^ 2;
t87 = t93 ^ 2;
t134 = t86 + t87;
t166 = t134 * mrSges(6,3);
t111 = -mrSges(5,2) + t166;
t75 = -t93 * mrSges(6,1) + t90 * mrSges(6,2);
t135 = -mrSges(5,1) + t75;
t92 = sin(qJ(2));
t94 = cos(qJ(3));
t139 = t92 * t94;
t95 = cos(qJ(2));
t83 = t95 * pkin(1) + pkin(2);
t91 = sin(qJ(3));
t66 = pkin(1) * t139 + t91 * t83;
t89 = cos(pkin(9));
t50 = t89 * t66;
t141 = t91 * t92;
t65 = -pkin(1) * t141 + t94 * t83;
t88 = sin(pkin(9));
t37 = t88 * t65 + t50;
t168 = -t66 * mrSges(4,1) - t65 * mrSges(4,2) + t135 * t37;
t144 = t89 * t91;
t67 = (t88 * t94 + t144) * pkin(2);
t167 = t135 * t67 + (-t91 * mrSges(4,1) - t94 * mrSges(4,2)) * pkin(2);
t70 = (t94 * t95 - t141) * pkin(1);
t147 = t70 * mrSges(4,2);
t69 = (-t91 * t95 - t139) * pkin(1);
t148 = t69 * mrSges(4,1);
t42 = t88 * t69 + t89 * t70;
t164 = t111 * t42 - t147 + t148 + (-t92 * mrSges(3,1) - t95 * mrSges(3,2)) * pkin(1);
t163 = 2 * qJD(3);
t162 = m(5) / 0.2e1;
t161 = m(6) / 0.2e1;
t160 = m(5) * pkin(3);
t146 = t88 * t66;
t60 = pkin(3) + t65;
t35 = t89 * t60 - t146;
t29 = -pkin(4) - t35;
t159 = -t29 / 0.2e1;
t38 = t89 * t65 - t146;
t158 = -t38 / 0.2e1;
t157 = -t42 / 0.2e1;
t145 = t88 * t91;
t152 = t94 * pkin(2);
t82 = pkin(3) + t152;
t63 = -pkin(2) * t145 + t89 * t82;
t54 = -pkin(4) - t63;
t156 = -t54 / 0.2e1;
t68 = (t89 * t94 - t145) * pkin(2);
t155 = -t68 / 0.2e1;
t81 = -t89 * pkin(3) - pkin(4);
t154 = -t81 / 0.2e1;
t41 = -t89 * t69 + t88 * t70;
t150 = t41 * mrSges(5,1);
t149 = t41 * t75;
t143 = t90 * mrSges(6,1);
t138 = t93 * mrSges(6,2);
t36 = t88 * t60 + t50;
t30 = pkin(8) + t36;
t125 = t30 * t134;
t3 = t111 * t38 + m(5) * (-t35 * t37 + t36 * t38) + m(6) * (t38 * t125 + t29 * t37) + t168;
t133 = t3 * qJD(1);
t4 = t135 * t41 + m(5) * (-t35 * t41 + t36 * t42) + m(6) * (t42 * t125 + t29 * t41) + m(4) * (t65 * t69 + t66 * t70) + t164;
t132 = t4 * qJD(1);
t108 = t138 + t143;
t104 = t29 * t108;
t85 = Ifges(6,4) * t93;
t120 = t172 * t90 + t93 * t85;
t15 = t104 + t120;
t131 = t15 * qJD(1);
t130 = t160 / 0.2e1;
t128 = t67 / 0.2e1 + t37 / 0.2e1;
t127 = t155 + t158;
t64 = pkin(2) * t144 + t88 * t82;
t55 = pkin(8) + t64;
t124 = t134 * t55;
t123 = t134 * t68;
t80 = t88 * pkin(3) + pkin(8);
t122 = t134 * t80;
t121 = Ifges(6,5) * t93 - Ifges(6,6) * t90;
t119 = -pkin(2) * t91 / 0.2e1 - t66 / 0.2e1;
t118 = -t152 / 0.2e1 - t65 / 0.2e1;
t117 = t158 + t159 + t154;
t116 = t157 + t159 + t156;
t114 = t155 + t154 + t156;
t102 = t81 * t108;
t112 = t102 / 0.2e1 + t120;
t96 = (t30 * t123 + t38 * t124 + t67 * t29 + t54 * t37) * t161 + (-t67 * t35 + t68 * t36 - t63 * t37 + t64 * t38) * t162;
t97 = -m(6) * (t42 * t122 + t81 * t41) / 0.2e1 - (-t41 * t89 + t42 * t88) * t160 / 0.2e1;
t2 = (t70 / 0.2e1 + t118) * mrSges(4,2) + (-t69 / 0.2e1 + t119) * mrSges(4,1) + t96 + t97 + t135 * (-t41 / 0.2e1 + t128) + t111 * (t157 - t127);
t6 = t111 * t68 + m(6) * (t55 * t123 + t54 * t67) + m(5) * (-t63 * t67 + t64 * t68) + t167;
t107 = t2 * qJD(1) + t6 * qJD(2);
t103 = t54 * t108;
t17 = t103 + t120;
t8 = (t116 * mrSges(6,2) - t85) * t93 + (t116 * mrSges(6,1) - t172) * t90;
t106 = -t8 * qJD(1) + t17 * qJD(2);
t105 = t143 / 0.2e1 + t138 / 0.2e1;
t99 = qJD(3) * t111;
t10 = (t117 * mrSges(6,2) - t85) * t93 + (t117 * mrSges(6,1) - t172) * t90;
t13 = (t114 * mrSges(6,2) - t85) * t93 + (t114 * mrSges(6,1) - t172) * t90;
t18 = t102 + t120;
t98 = t10 * qJD(1) + t13 * qJD(2) - t18 * qJD(3);
t43 = t103 / 0.2e1;
t19 = t104 / 0.2e1;
t14 = -t105 * t68 + t112 + t43;
t11 = -t105 * t38 + t112 + t19;
t9 = -t105 * t42 + t120 + t19 + t43;
t1 = t149 / 0.2e1 - t150 / 0.2e1 + t148 / 0.2e1 - t147 / 0.2e1 + t118 * mrSges(4,2) + t119 * mrSges(4,1) + t96 - t97 + t42 * t166 / 0.2e1 + t135 * t128 + (t38 + t68) * mrSges(6,3) * (t87 / 0.2e1 + t86 / 0.2e1) + (t157 + t127) * mrSges(5,2);
t5 = [qJD(2) * t4 + qJD(3) * t3 + qJD(5) * t15, t1 * qJD(3) + t9 * qJD(5) + t132 + (t149 - t150 + 0.2e1 * (-t41 * t63 + t42 * t64) * t162 + 0.2e1 * (t42 * t124 + t54 * t41) * t161 + m(4) * (t69 * t94 + t70 * t91) * pkin(2) + t164) * qJD(2), t133 + t1 * qJD(2) + t168 * qJD(3) + t11 * qJD(5) + t38 * t99 + ((t38 * t122 + t81 * t37) * t161 + (-t37 * t89 + t38 * t88) * t130) * t163, 0, t131 + t9 * qJD(2) + t11 * qJD(3) + (t75 * t30 + t121) * qJD(5); qJD(3) * t2 - qJD(5) * t8 - t132, qJD(3) * t6 + qJD(5) * t17, t167 * qJD(3) + t14 * qJD(5) + t68 * t99 + ((t68 * t122 + t81 * t67) * t161 + (-t67 * t89 + t68 * t88) * t130) * t163 + t107, 0, t14 * qJD(3) + (t75 * t55 + t121) * qJD(5) + t106; -qJD(2) * t2 - qJD(5) * t10 - t133, -qJD(5) * t13 - t107, t18 * qJD(5), 0, (t75 * t80 + t121) * qJD(5) - t98; 0, 0, 0, 0, -t108 * qJD(5); qJD(2) * t8 + qJD(3) * t10 - t131, qJD(3) * t13 - t106, t98, 0, 0;];
Cq = t5;
