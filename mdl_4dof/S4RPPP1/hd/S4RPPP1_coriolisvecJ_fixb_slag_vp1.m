% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:19
% EndTime: 2019-03-08 18:26:21
% DurationCPUTime: 1.08s
% Computational Cost: add. (1105->152), mult. (2933->176), div. (0->0), fcn. (2644->6), ass. (0->86)
t92 = sin(pkin(4));
t94 = cos(qJ(1));
t147 = t92 * t94;
t159 = rSges(5,1) + pkin(3);
t132 = cos(pkin(6));
t133 = cos(pkin(4));
t107 = t133 * t132;
t131 = sin(pkin(6));
t93 = sin(qJ(1));
t57 = -t94 * t107 + t93 * t131;
t111 = -t57 * rSges(5,2) + t159 * t147;
t155 = rSges(5,3) + qJ(4);
t106 = t133 * t131;
t58 = t106 * t94 + t132 * t93;
t160 = -t155 * t58 + t111;
t59 = t107 * t93 + t131 * t94;
t60 = -t93 * t106 + t94 * t132;
t115 = t60 * pkin(2) + qJ(3) * t59;
t128 = qJD(3) * t57;
t135 = qJ(2) * t92;
t91 = t94 * pkin(1);
t136 = t93 * t135 + t91;
t129 = qJD(2) * t92;
t80 = t94 * t129;
t157 = qJD(1) * t136 - t80;
t158 = qJD(1) * t115 + t128 + t157;
t156 = pkin(1) * t93;
t148 = t92 * t93;
t103 = -rSges(3,1) * t58 + rSges(3,2) * t57 + rSges(3,3) * t147;
t63 = -t94 * t135 + t156;
t154 = t103 - t63;
t49 = t57 * qJ(3);
t153 = -pkin(2) * t58 - t49;
t79 = t93 * t129;
t137 = -qJD(1) * t63 + t79;
t121 = qJD(1) * t147;
t138 = qJ(2) * t121 + t79;
t48 = qJD(3) * t59;
t40 = t57 * qJD(1);
t41 = t58 * qJD(1);
t99 = -t41 * pkin(2) - qJ(3) * t40 + t48;
t152 = -qJD(1) * t153 - t137 + t138 - t48 + t99;
t22 = rSges(4,1) * t148 - t60 * rSges(4,2) + t59 * rSges(4,3);
t12 = qJD(1) * t22 + t158;
t45 = qJD(4) * t58;
t151 = t45 + t158;
t42 = t59 * qJD(1);
t142 = t42 * qJ(3) + t128;
t43 = t60 * qJD(1);
t20 = t43 * pkin(2) + t142;
t146 = -t20 - t157;
t145 = t59 * rSges(5,2) + t159 * t148 + t155 * t60;
t144 = -t63 + t153;
t120 = qJD(1) * t129;
t130 = qJD(1) * t93;
t143 = qJD(1) * (-pkin(1) * t130 + t138) + t93 * t120;
t141 = t48 + t79;
t73 = t94 * t120;
t139 = -qJD(3) * t40 + t73;
t46 = qJD(4) * t60;
t127 = -pkin(2) - t155;
t125 = -t41 * rSges(3,1) + t40 * rSges(3,2) + rSges(3,3) * t121;
t124 = t60 * rSges(3,1) - t59 * rSges(3,2) + rSges(3,3) * t148;
t123 = rSges(4,1) * t121 + t41 * rSges(4,2) - t40 * rSges(4,3);
t122 = t92 * t130;
t119 = -t42 * rSges(5,2) - t45;
t118 = rSges(4,1) * t147 - rSges(4,3) * t57;
t114 = qJD(1) * t99 + qJD(3) * t42 + t143;
t112 = -t49 - t63;
t109 = -rSges(3,1) * t43 + rSges(3,2) * t42;
t108 = rSges(4,2) * t43 - rSges(4,3) * t42;
t105 = t115 + t136;
t102 = rSges(4,2) * t58 + t118;
t96 = -t40 * rSges(5,2) + t159 * t121 - t155 * t41 + t46;
t18 = qJD(1) * t124 + t157;
t17 = t154 * qJD(1) + t79;
t14 = t73 + (-rSges(3,3) * t122 + t109 - t157) * qJD(1);
t13 = qJD(1) * t125 + t143;
t11 = (t102 + t144) * qJD(1) + t141;
t8 = (-rSges(4,1) * t122 + t108 + t146) * qJD(1) + t139;
t7 = qJD(1) * t123 + t114;
t6 = t145 * qJD(1) + t151;
t5 = t46 + (t144 + t160) * qJD(1) + t141;
t2 = -qJD(4) * t41 + (-t122 * t159 - t155 * t43 + t119 + t146) * qJD(1) + t139;
t1 = qJD(1) * t96 + qJD(4) * t43 + t114;
t3 = [(t1 * (t105 + t145) + (t127 * t58 + t111 + t112) * t2 + (t96 - t46 + (-t156 - t160) * qJD(1) + t152) * t6 + (t127 * t43 + t119 - t142 + t80 + (-t91 + (-qJ(2) - t159) * t148 + t145) * qJD(1) + t151) * t5) * m(5) + (t8 * ((rSges(4,2) - pkin(2)) * t58 + t112 + t118) + t7 * (t105 + t22) + (t123 + (-t102 - t156) * qJD(1) + t152) * t12 + (t108 - t20 + t80 + (-t91 + (-rSges(4,1) - qJ(2)) * t148) * qJD(1) + t12) * t11) * m(4) + (t14 * t154 + t17 * (t109 + t80) + t13 * (t124 + t136) + t18 * (t125 + t138) + (-t17 * t91 + (-t18 * pkin(1) + t17 * (-rSges(3,3) - qJ(2)) * t92) * t93) * qJD(1) - (qJD(1) * t103 + t137 - t17) * t18) * m(3); 0.2e1 * (m(3) * (-t13 * t94 + t14 * t93) / 0.2e1 + m(4) * (-t7 * t94 + t8 * t93) / 0.2e1 + m(5) * (-t1 * t94 + t2 * t93) / 0.2e1) * t92; m(4) * (-t11 * t40 + t12 * t42 + t57 * t7 + t59 * t8) + m(5) * (t1 * t57 + t2 * t59 - t5 * t40 + t6 * t42) + 0.2e1 * (-m(4) * (-t11 * t57 + t12 * t59) / 0.2e1 - m(5) * (-t5 * t57 + t6 * t59) / 0.2e1) * qJD(1); (t1 * t58 + t2 * t60 - t41 * t5 + t43 * t6 - (-t5 * t58 + t6 * t60) * qJD(1)) * m(5);];
tauc  = t3(:);
