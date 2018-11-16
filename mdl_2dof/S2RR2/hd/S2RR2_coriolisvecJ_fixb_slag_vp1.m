% Calculate vector of centrifugal and coriolis load on the joints for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [2x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S2RR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert( isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:48:45
% EndTime: 2018-11-16 16:48:47
% DurationCPUTime: 1.70s
% Computational Cost: add. (1057->177), mult. (2851->261), div. (0->0), fcn. (2266->4), ass. (0->111)
t85 = cos(qJ(2));
t80 = Icges(3,4) * t85;
t83 = sin(qJ(2));
t102 = -Icges(3,2) * t83 + t80;
t86 = cos(qJ(1));
t131 = Icges(3,6) * t86;
t84 = sin(qJ(1));
t41 = -t102 * t84 + t131;
t132 = Icges(3,6) * t84;
t133 = Icges(3,2) * t86;
t138 = Icges(3,4) * t86;
t42 = -t133 * t83 + t138 * t85 + t132;
t171 = Icges(3,1) * t83 + t80;
t51 = t171 * t84;
t52 = t171 * t86;
t166 = (t42 + t52) * t84 + (t41 - t51) * t86;
t134 = Icges(3,2) * t85;
t136 = Icges(3,5) * t86;
t140 = Icges(3,1) * t85;
t139 = Icges(3,4) * t83;
t76 = t84 * t139;
t43 = -t140 * t84 + t136 + t76;
t137 = Icges(3,5) * t84;
t77 = t83 * t138;
t44 = t140 * t86 + t137 - t77;
t94 = (-t133 * t85 + t44 - t77) * t84 + (t134 * t84 + t43 + t76) * t86;
t177 = t166 * t85 + t94 * t83;
t176 = t86 * t84;
t128 = qJD(2) * t84;
t67 = rSges(3,1) * t83 + rSges(3,2) * t85;
t55 = t67 * t128;
t142 = t171 + t102;
t63 = t134 + t139;
t66 = -t139 + t140;
t143 = t63 - t66;
t175 = (t142 * t83 + t143 * t85) * qJD(1);
t174 = 0.2e1 * qJD(2);
t106 = t42 * t83 - t44 * t85;
t172 = t106 * t86;
t57 = t102 * qJD(2);
t58 = t66 * qJD(2);
t61 = Icges(3,5) * t83 + Icges(3,6) * t85;
t168 = -qJD(1) * t61 + t57 * t83 - t58 * t85 + (t171 * t83 + t63 * t85) * qJD(2);
t163 = pkin(1) * qJD(1) ^ 2;
t162 = -rSges(3,3) - pkin(1);
t160 = rSges(3,1) * t85;
t158 = t41 * t83;
t157 = t43 * t85;
t47 = t61 * t84;
t156 = t61 * t86;
t155 = t83 * t84;
t154 = t83 * t86;
t153 = t84 * t85;
t152 = t85 * t86;
t82 = t86 * rSges(3,3);
t62 = Icges(3,5) * t85 - Icges(3,6) * t83;
t39 = Icges(3,3) * t86 - t62 * t84;
t149 = t41 * t155 + t86 * t39;
t148 = t43 * t152 + t84 * t39;
t130 = Icges(3,3) * t84;
t127 = qJD(2) * t86;
t104 = -t171 * t85 + t83 * t63;
t22 = -t104 * t86 + t47;
t126 = t22 * qJD(1);
t125 = t62 * qJD(1);
t124 = rSges(3,1) * t152;
t79 = rSges(3,2) * t154;
t123 = qJD(1) * t79 + t55;
t122 = t67 * t127;
t45 = -rSges(3,1) * t153 + rSges(3,2) * t155 + t82;
t121 = pkin(1) * t86 + t45;
t118 = -t128 / 0.2e1;
t116 = -t127 / 0.2e1;
t115 = t127 / 0.2e1;
t40 = -t131 * t83 + t136 * t85 + t130;
t114 = -t40 - t157;
t109 = -rSges(3,2) * t83 + t160;
t101 = -rSges(3,3) * t84 - t124;
t46 = -t101 - t79;
t23 = t55 + (-pkin(1) * t84 - t46) * qJD(1);
t24 = qJD(1) * t121 - t122;
t108 = t23 * t84 - t24 * t86;
t19 = t41 * t85 + t43 * t83;
t107 = -t157 + t158;
t20 = t42 * t85 + t44 * t83;
t14 = -t153 * t44 + t42 * t155 + t86 * t40;
t54 = t67 * t86;
t13 = -t153 * t43 + t149;
t100 = (t13 * t86 + t14 * t84) * qJD(2);
t15 = -t154 * t41 + t148;
t16 = t84 * t40 - t172;
t99 = (t15 * t86 + t16 * t84) * qJD(2);
t17 = (-t45 * t84 + t46 * t86) * qJD(2);
t93 = -qJD(2) * t156 + (t106 + t39) * qJD(1);
t92 = qJD(2) * t47 + (-t62 * t86 + t107 - t130) * qJD(1);
t91 = t104 * qJD(1) + t62 * qJD(2);
t59 = t109 * qJD(2);
t53 = t67 * t84;
t32 = qJD(1) * t101 + t123;
t31 = -qJD(2) * t54 + (-t109 * t84 + t82) * qJD(1);
t21 = t104 * t84 + t156;
t18 = t21 * qJD(1);
t12 = -t86 * t163 + t59 * t128 + (-t31 + t122) * qJD(1);
t11 = -t84 * t163 - t59 * t127 + (t32 + t55) * qJD(1);
t10 = t168 * t84 + t91 * t86;
t9 = -t168 * t86 + t91 * t84;
t8 = -t106 * qJD(2) + (qJD(1) * t41 - t127 * t63) * t85 + (-qJD(2) * t52 + (-t66 * t84 + t136) * qJD(1)) * t83;
t7 = -qJD(2) * t107 + (t63 * t128 + (-t102 * t86 - t132) * qJD(1)) * t85 + (qJD(2) * t51 + (-t66 * t86 - t137) * qJD(1)) * t83;
t6 = t99 + t126;
t5 = t18 + t100;
t1 = [(t18 + ((t149 + t16 + t172) * t86 + (-t15 + (t114 - t158) * t86 + t14 + t148) * t84) * qJD(2)) * t118 + m(3) * (t11 * t121 + t24 * t123 + t12 * (t162 * t84 - t124 + t79) + t23 * t122) + (t6 - t126 + ((t14 + (-t40 + t158) * t86 - t148) * t86 + (t114 * t84 - t13 + t149) * t84) * qJD(2)) * t116 + (t7 + t10) * t115 + (t5 + t9 + t8) * t128 / 0.2e1 + (-qJD(2) * t104 + t57 * t85 + t58 * t83 + m(3) * ((-t160 * t24 + t162 * t23) * t86 + (t109 * t23 + t162 * t24) * t84) + (t19 + t21) * t118 + (t22 + t20) * t115) * qJD(1); ((-t128 * t156 + t125) * t84 + (-t175 + (t84 * t47 - t177) * qJD(2)) * t86) * t118 + ((t47 * t127 + t125) * t86 + (t175 + (-t156 * t86 + t177) * qJD(2)) * t84) * t116 + (t9 * qJD(1) + (t92 * t176 + t93 * t84 ^ 2 + (-t15 * t84 + t16 * t86) * qJD(1)) * t174) * t84 / 0.2e1 + (t10 * qJD(1) + (t92 * t86 ^ 2 + t93 * t176 + (-t13 * t84 + t14 * t86) * qJD(1)) * t174) * t86 / 0.2e1 + (0.2e1 * t17 * (t86 * t31 - t84 * t32 + (-t45 * t86 - t46 * t84) * qJD(1)) + t108 * t59 + (-t11 * t86 + t12 * t84 + (t23 * t86 + t24 * t84) * qJD(1)) * t67 - (t23 * t54 + t24 * t53) * qJD(1) - (t17 * (-t53 * t84 - t54 * t86) + t108 * t109) * qJD(2)) * m(3) - ((t142 * t85 - t143 * t83) * qJD(1) + (-t166 * t83 + t85 * t94) * qJD(2) + (t5 + t100) * t84) * qJD(1) / 0.2e1 + ((-t19 * qJD(1) + t8) * t84 + (t20 * qJD(1) + t6 + t7 + t99) * t86) * qJD(1) / 0.2e1;];
tauc  = t1(:);
