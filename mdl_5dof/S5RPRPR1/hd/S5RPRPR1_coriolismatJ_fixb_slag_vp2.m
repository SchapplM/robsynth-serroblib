% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:05
% EndTime: 2019-12-05 17:47:09
% DurationCPUTime: 1.36s
% Computational Cost: add. (4888->140), mult. (8290->186), div. (0->0), fcn. (9185->6), ass. (0->87)
t105 = sin(qJ(5));
t107 = cos(qJ(5));
t104 = sin(pkin(8));
t106 = sin(qJ(3));
t108 = cos(qJ(3));
t138 = cos(pkin(8));
t92 = t104 * t106 - t108 * t138;
t93 = -t104 * t108 - t106 * t138;
t118 = t105 * t93 - t107 * t92;
t120 = t105 * t92 + t107 * t93;
t109 = -pkin(1) - pkin(6);
t136 = t106 * t109;
t94 = -t106 * qJ(4) + t136;
t95 = (-qJ(4) + t109) * t108;
t144 = -t104 * t94 + t138 * t95;
t177 = pkin(7) * t92 + t144;
t73 = -t104 * t95 - t138 * t94;
t48 = pkin(7) * t93 - t73;
t196 = -t105 * t48 + t107 * t177;
t35 = t105 * t177 + t107 * t48;
t3 = -t35 * mrSges(6,1) - t196 * mrSges(6,2) + Ifges(6,5) * t120 - Ifges(6,6) * t118;
t208 = t3 * qJD(5);
t101 = t106 * pkin(3) + qJ(2);
t79 = -pkin(4) * t93 + t101;
t207 = m(6) * t79;
t201 = t120 * mrSges(6,1);
t128 = t201 / 0.2e1;
t180 = t120 * mrSges(6,2);
t125 = t118 * mrSges(6,1) + t180;
t163 = Ifges(6,4) * t118;
t168 = -t118 / 0.2e1;
t185 = t120 / 0.2e1;
t193 = t118 / 0.2e1;
t2 = (Ifges(6,2) * t120 + t163) * t168 + (Ifges(6,1) * t120 - t163) * t193 + (0.2e1 * Ifges(6,4) * t120 + (Ifges(6,1) - Ifges(6,2)) * t118) * t185 + t79 * t125;
t182 = t118 * mrSges(6,2);
t203 = t201 - t182;
t202 = m(5) * t101;
t199 = t168 + t193;
t197 = t118 ^ 2 + t120 ^ 2;
t129 = t180 / 0.2e1;
t100 = pkin(3) * t138 + pkin(4);
t158 = pkin(3) * t104;
t81 = t100 * t107 - t105 * t158;
t82 = t100 * t105 + t107 * t158;
t186 = t118 * t82 + t120 * t81;
t178 = t104 * t92 - t138 * t93;
t175 = -m(5) / 0.2e1;
t173 = m(6) / 0.2e1;
t172 = m(5) * pkin(3);
t6 = (t118 * t185 - t120 * t193) * mrSges(6,3);
t166 = t6 * qJD(3);
t157 = t108 * pkin(3);
t127 = t92 ^ 2 + t93 ^ 2;
t9 = t197 * mrSges(6,3) + t127 * mrSges(5,3) + m(6) * (-t118 * t196 + t120 * t35) + m(5) * (t144 * t92 - t73 * t93);
t143 = qJD(1) * t9;
t140 = t108 * mrSges(4,2);
t139 = t2 * qJD(1);
t124 = -t93 * mrSges(5,1) - t92 * mrSges(5,2);
t114 = -t106 * mrSges(4,1) - t124 - t140 + t203;
t22 = mrSges(3,3) + (m(4) + m(3)) * qJ(2) + t202 + t207 - t114;
t137 = qJD(1) * t22;
t131 = t172 / 0.2e1;
t112 = (-t118 * t81 + t120 * t82) * t173 + (t104 * t93 + t138 * t92) * t131;
t80 = -pkin(4) * t92 + t157;
t116 = t108 * t131 + t173 * t80;
t84 = t93 * mrSges(5,2);
t15 = t92 * mrSges(5,1) + t112 - t116 - t125 - t84;
t135 = t15 * qJD(1);
t113 = t127 * t175 - t197 * t173;
t130 = t175 - m(6) / 0.2e1;
t18 = t113 + t130;
t134 = t18 * qJD(1);
t23 = 0.2e1 * t193 * mrSges(6,1) + 0.2e1 * t129;
t133 = t23 * qJD(1);
t27 = -mrSges(6,1) * t82 - mrSges(6,2) * t81;
t132 = t27 * qJD(5);
t1 = t101 * t84 + (-t101 * mrSges(5,1) - Ifges(5,4) * t92) * t92 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t108 + pkin(3) * t124) * t108 + t157 * t202 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t106 + (-Ifges(4,1) + Ifges(4,2)) * t108) * t106 + (Ifges(5,4) * t93 + (-Ifges(5,1) + Ifges(5,2)) * t92) * t93 + t2 + (-t203 + t207) * t80;
t119 = t1 * qJD(1) - t6 * qJD(2);
t117 = t6 * qJD(1);
t13 = t128 - t201 / 0.2e1 + t199 * mrSges(6,2);
t4 = t199 * Ifges(6,6) + (t185 - t120 / 0.2e1) * Ifges(6,5);
t115 = t4 * qJD(1) + t13 * qJD(2) + t27 * qJD(3);
t24 = t129 - t180 / 0.2e1;
t21 = t112 + t116;
t17 = t113 - t130;
t14 = -t182 + 0.2e1 * t128;
t5 = [qJD(2) * t22 + qJD(3) * t1 + qJD(4) * t9 + qJD(5) * t2, qJD(4) * t17 + t137 - t166, (-t109 * t140 - mrSges(4,1) * t136 + m(6) * (t196 * t82 - t35 * t81) + (t104 * t144 + t138 * t73) * t172 - t144 * mrSges(5,2) + t73 * mrSges(5,1) + Ifges(5,5) * t93 + Ifges(5,6) * t92 - Ifges(4,5) * t106 - Ifges(4,6) * t108 + t178 * mrSges(5,3) * pkin(3) - t186 * mrSges(6,3) + t3) * qJD(3) + t21 * qJD(4) + t208 + t119, qJD(2) * t17 + qJD(3) * t21 + qJD(5) * t24 + t143, t3 * qJD(3) + t24 * qJD(4) + t139 + t208; qJD(4) * t18 - t137 - t166, 0, (m(6) * t186 - t178 * t172 + t114) * qJD(3) + t14 * qJD(5) - t117, t134, t14 * qJD(3) + t203 * qJD(5); qJD(4) * t15 + qJD(5) * t4 - t119, qJD(5) * t13 + t117, t132, t135, t115 + t132; -qJD(2) * t18 - qJD(3) * t15 + qJD(5) * t23 - t143, -t134, -t135, 0, t133; -qJD(3) * t4 - qJD(4) * t23 - t139, -t13 * qJD(3), -t115, -t133, 0;];
Cq = t5;
