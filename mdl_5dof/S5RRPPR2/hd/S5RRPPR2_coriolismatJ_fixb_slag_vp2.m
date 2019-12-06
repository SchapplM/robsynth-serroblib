% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:08
% EndTime: 2019-12-05 18:20:10
% DurationCPUTime: 0.91s
% Computational Cost: add. (2988->171), mult. (6286->242), div. (0->0), fcn. (5136->8), ass. (0->119)
t102 = sin(qJ(5));
t178 = t102 / 0.2e1;
t177 = qJD(1) + qJD(2);
t104 = cos(qJ(5));
t145 = t104 * mrSges(6,2);
t149 = t102 * mrSges(6,1);
t128 = t145 + t149;
t98 = sin(pkin(9));
t96 = t98 ^ 2;
t100 = cos(pkin(9));
t97 = t100 ^ 2;
t176 = t96 * t128 + (t96 + t97) * mrSges(5,3);
t101 = cos(pkin(8));
t103 = sin(qJ(2));
t138 = t101 * t103;
t105 = cos(qJ(2));
t158 = t105 * pkin(1);
t93 = pkin(2) + t158;
t99 = sin(pkin(8));
t131 = pkin(1) * t138 + t99 * t93;
t75 = qJ(4) + t131;
t67 = t96 * t75;
t92 = t99 * pkin(2) + qJ(4);
t85 = t96 * t92;
t154 = t67 + t85;
t140 = t100 * t102;
t123 = -t100 * pkin(4) - t98 * pkin(7) - pkin(3);
t91 = t99 * t103 * pkin(1);
t132 = t101 * t93 - t91;
t60 = t123 - t132;
t37 = t104 * t60 - t75 * t140;
t139 = t100 * t104;
t38 = t102 * t60 + t75 * t139;
t159 = t101 * pkin(2);
t79 = t123 - t159;
t56 = t104 * t79 - t92 * t140;
t57 = t102 * t79 + t92 * t139;
t175 = ((t38 + t57) * t104 + (-t37 - t56) * t102) * t100 + t154;
t129 = mrSges(6,1) * t104 - mrSges(6,2) * t102;
t160 = t104 / 0.2e1;
t166 = t104 ^ 2;
t167 = t102 ^ 2;
t147 = t102 * t98;
t137 = mrSges(6,3) * t147;
t82 = t100 * mrSges(6,2) - t137;
t143 = t104 * t98;
t136 = mrSges(6,3) * t143;
t83 = -t100 * mrSges(6,1) - t136;
t106 = (t166 / 0.2e1 + t167 / 0.2e1) * t96 * mrSges(6,3) + (t100 * t129 / 0.2e1 + t82 * t178 + t83 * t160) * t98;
t174 = qJD(3) * t106;
t115 = (-t149 / 0.2e1 - t145 / 0.2e1) * t100;
t171 = -t104 * t82 / 0.2e1 + t83 * t178;
t109 = t115 + t171;
t173 = qJD(4) * t109;
t112 = t115 - t171;
t172 = t112 * qJD(4);
t62 = (t166 - 0.1e1 + t167) * t98 * t100;
t161 = m(6) * t62;
t134 = t161 / 0.2e1;
t155 = qJD(4) * t134 - qJD(5) * t106;
t164 = -m(6) / 0.2e1;
t133 = qJD(3) * t164;
t170 = -qJD(5) * t109 + t62 * t133;
t169 = qJD(3) * t134 + qJD(5) * t112;
t80 = (t105 * t99 + t138) * pkin(1);
t81 = t101 * t158 - t91;
t47 = t104 * t80 - t81 * t140;
t48 = t102 * t80 + t81 * t139;
t168 = (-mrSges(4,2) + t176) * t81 + t47 * t83 + t48 * t82 + (-t103 * mrSges(3,1) - t105 * mrSges(3,2)) * pkin(1);
t165 = m(5) / 0.2e1;
t163 = m(6) / 0.2e1;
t24 = (-t100 * t81 - t102 * t47 + t104 * t48) * t98;
t162 = m(6) * t24;
t157 = t97 * t75;
t156 = t97 * t92;
t153 = Ifges(6,4) * t102;
t152 = Ifges(6,4) * t104;
t151 = Ifges(6,5) * t100;
t150 = Ifges(6,6) * t100;
t118 = t82 * t139 - t83 * t140 + t176;
t10 = m(6) * (t67 + (-t102 * t37 + t104 * t38) * t100) + m(5) * (t67 + t157) + t118;
t141 = t10 * qJD(1);
t135 = t162 / 0.2e1;
t14 = m(6) * (t85 + (-t102 * t56 + t104 * t57) * t100) + m(5) * (t85 + t156) + t118;
t108 = ((t75 + t92) * t97 + t154) * t165 + t118;
t111 = (t102 * t48 + t104 * t47) * t163 + t80 * t165;
t5 = t175 * t164 - t108 + t111;
t127 = t5 * qJD(1) - t14 * qJD(2);
t110 = Ifges(6,4) * t143 + (Ifges(6,1) - Ifges(6,2)) * t147 - t150;
t117 = (-Ifges(6,4) * t147 - t151) * t102;
t31 = t37 * t82;
t32 = t38 * t83;
t36 = t37 * t137;
t121 = t96 * t129;
t51 = t75 * t121;
t4 = -t51 + t32 - t36 - t31 + (t117 + (t38 * mrSges(6,3) + t110) * t104) * t98;
t126 = -t4 * qJD(1) - t174;
t125 = t177 * t106;
t124 = t177 * t109;
t122 = -t47 * mrSges(6,1) / 0.2e1 + t48 * mrSges(6,2) / 0.2e1;
t120 = (-Ifges(6,5) * t102 - Ifges(6,6) * t104) * t98;
t52 = t81 * t67;
t84 = -t100 * mrSges(5,1) + t98 * mrSges(5,2);
t3 = (-mrSges(4,1) + t84) * t80 + m(6) * (t37 * t47 + t38 * t48 + t52) + m(4) * (t131 * t81 - t132 * t80) + m(5) * (t81 * t157 + t52 + (-pkin(3) - t132) * t80) + t168;
t119 = -t3 * qJD(1) + t24 * t133;
t41 = t56 * t82;
t42 = t57 * t83;
t49 = t56 * t137;
t63 = t92 * t121;
t107 = -(-t151 + (Ifges(6,1) * t104 - t153) * t98) * t147 / 0.2e1 - (-t150 + (-Ifges(6,2) * t102 + t152) * t98) * t143 / 0.2e1 - t100 * t120 / 0.2e1 + (-t38 / 0.2e1 - t57 / 0.2e1) * t136 + t31 / 0.2e1 - t32 / 0.2e1 + t36 / 0.2e1 + t41 / 0.2e1 - t42 / 0.2e1 + t49 / 0.2e1 + t51 / 0.2e1 + t63 / 0.2e1 + (-t102 * (-Ifges(6,2) * t104 - t153) / 0.2e1 + (-Ifges(6,1) * t102 - t152) * t160) * t96;
t2 = t107 + t122;
t7 = -t41 + t42 - t49 - t63 + (t117 + (t57 * mrSges(6,3) + t110) * t104) * t98;
t116 = t2 * qJD(1) - t7 * qJD(2) - t174;
t64 = t81 * t85;
t35 = 0.2e1 * (qJD(1) / 0.4e1 + qJD(2) / 0.4e1) * t161;
t8 = qJD(2) * t135 + t155;
t6 = t175 * t163 + t108 + t111;
t1 = t107 - t122;
t9 = [qJD(2) * t3 + qJD(4) * t10 - qJD(5) * t4, t6 * qJD(4) + t1 * qJD(5) - t119 + (-t80 * mrSges(4,1) + t80 * t84 + 0.2e1 * (t56 * t47 + t57 * t48 + t64) * t163 + m(4) * (-t101 * t80 + t81 * t99) * pkin(2) + 0.2e1 * (t81 * t156 + t64 + (-pkin(3) - t159) * t80) * t165 + t168) * qJD(2), t8, t6 * qJD(2) + t141 + t169, t1 * qJD(2) + t172 + (-t38 * mrSges(6,1) - t37 * mrSges(6,2) + t120) * qJD(5) + t126; -t5 * qJD(4) + t2 * qJD(5) + t119, qJD(4) * t14 - qJD(5) * t7, -qJD(1) * t162 / 0.2e1 + t155, -t127 + t169, t172 + (-t57 * mrSges(6,1) - t56 * mrSges(6,2) + t120) * qJD(5) + t116; t8, qJD(1) * t135 + t155, 0, t35, -t129 * qJD(5) * t98 - t125; t5 * qJD(2) - t141 + t170, t127 + t170, -t35, 0, -t128 * qJD(5) - t124; -qJD(2) * t2 - t126 + t173, -t116 + t173, t125, t124, 0;];
Cq = t9;
