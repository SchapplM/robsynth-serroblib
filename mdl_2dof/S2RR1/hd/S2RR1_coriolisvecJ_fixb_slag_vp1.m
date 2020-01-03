% Calculate vector of centrifugal and Coriolis load on the joints for
% S2RR1
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S2RR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:07
% EndTime: 2020-01-03 11:19:09
% DurationCPUTime: 1.97s
% Computational Cost: add. (1057->169), mult. (2851->267), div. (0->0), fcn. (2266->4), ass. (0->113)
t77 = sin(qJ(1));
t79 = cos(qJ(1));
t78 = cos(qJ(2));
t137 = Icges(3,4) * t78;
t76 = sin(qJ(2));
t99 = -Icges(3,2) * t76 + t137;
t40 = Icges(3,6) * t79 + t99 * t77;
t135 = Icges(3,5) * t79;
t139 = Icges(3,1) * t78;
t138 = Icges(3,4) * t76;
t73 = t77 * t138;
t42 = t77 * t139 + t135 - t73;
t103 = t40 * t76 - t42 * t78;
t136 = Icges(3,5) * t78;
t98 = -Icges(3,6) * t76 + t136;
t39 = -Icges(3,3) * t77 + t98 * t79;
t170 = t103 - t39;
t169 = t77 * t79;
t104 = rSges(3,1) * t76 + rSges(3,2) * t78;
t123 = t77 * qJD(2);
t121 = t104 * t123;
t128 = qJD(1) * t77;
t146 = t78 * t79;
t148 = t76 * t79;
t110 = rSges(3,1) * t146 - rSges(3,2) * t148;
t154 = rSges(3,3) * t77;
t93 = -t110 + t154;
t23 = pkin(1) * t128 + qJD(1) * t93 + t121;
t147 = t77 * t78;
t122 = rSges(3,1) * t147;
t149 = t76 * t77;
t74 = rSges(3,2) * t149;
t111 = t74 - t122;
t44 = rSges(3,3) * t79 - t111;
t126 = qJD(2) * t79;
t54 = t104 * t126;
t24 = -(pkin(1) * t79 + t44) * qJD(1) - t54;
t51 = t104 * t77;
t52 = t104 * t79;
t168 = t23 * t52 + t24 * t51;
t167 = 0.2e1 * qJD(2);
t164 = -t23 * t77 + t24 * t79;
t133 = Icges(3,2) * t78;
t62 = -t133 - t138;
t64 = -Icges(3,1) * t76 - t137;
t101 = t76 * t62 - t78 * t64;
t163 = t101 * qJD(1) + t98 * qJD(2);
t132 = Icges(3,6) * t77;
t41 = t99 * t79 - t132;
t150 = t41 * t76;
t100 = -t138 + t139;
t43 = -Icges(3,5) * t77 + t100 * t79;
t102 = -t43 * t78 + t150;
t130 = Icges(3,3) * t79;
t60 = -Icges(3,5) * t76 - Icges(3,6) * t78;
t90 = qJD(2) * t60;
t162 = -t79 * t90 + (t98 * t77 + t102 + t130) * qJD(1);
t161 = t170 * qJD(1) - t77 * t90;
t160 = -t40 * t79 + t41 * t77;
t157 = pkin(1) * qJD(1) ^ 2;
t156 = rSges(3,3) + pkin(1);
t155 = rSges(3,1) * t78;
t45 = t60 * t77;
t46 = t60 * t79;
t36 = t79 * t39;
t145 = -qJD(1) / 0.2e1;
t144 = t43 * t147 + t36;
t141 = t62 + t100;
t140 = t99 - t64;
t127 = qJD(1) * t79;
t21 = -t101 * t77 + t46;
t125 = t21 * qJD(1);
t124 = t98 * qJD(1);
t117 = -t126 / 0.2e1;
t114 = t123 / 0.2e1;
t38 = -t76 * t132 + t77 * t136 + t130;
t113 = -t38 - t150;
t109 = t42 * t146 - t40 * t148;
t105 = -rSges(3,2) * t76 + t155;
t19 = -t40 * t78 - t42 * t76;
t20 = -t41 * t78 - t43 * t76;
t96 = t104 * qJD(2);
t13 = -t103 * t77 + t38 * t79;
t14 = -t41 * t149 + t144;
t95 = (t13 * t79 - t14 * t77) * qJD(2);
t15 = -t38 * t77 + t109;
t35 = t43 * t146;
t16 = -t41 * t148 - t77 * t39 + t35;
t94 = (t15 * t79 - t16 * t77) * qJD(2);
t92 = qJD(2) * t64;
t91 = qJD(2) * t62;
t89 = (t77 * t133 - t42 + t73) * t79 - (-t62 * t79 - t43) * t77;
t88 = (t140 * t76 - t141 * t78) * qJD(1);
t17 = (-t77 * t44 + t79 * t93) * qJD(2);
t56 = t99 * qJD(2);
t57 = t100 * qJD(2);
t83 = -qJD(1) * t60 + t56 * t76 - t57 * t78 + (-t62 * t78 - t64 * t76) * qJD(2);
t81 = t160 * t78 + t89 * t76;
t70 = qJD(1) * t74;
t58 = t105 * qJD(2);
t32 = -t77 * t96 + (t105 * t79 - t154) * qJD(1);
t31 = -qJD(1) * t122 + t70 + (-rSges(3,3) * qJD(1) - t96) * t79;
t22 = -t101 * t79 - t45;
t18 = t22 * qJD(1);
t12 = -t79 * t157 - t58 * t123 + (t31 - t54) * qJD(1);
t11 = t77 * t157 - t58 * t126 + (-t32 + t121) * qJD(1);
t10 = -t163 * t79 + t83 * t77;
t9 = t163 * t77 + t83 * t79;
t8 = t102 * qJD(2) - (-t40 * qJD(1) + t79 * t91) * t78 - (t79 * t92 + (-t100 * t77 - t135) * qJD(1)) * t76;
t7 = t103 * qJD(2) - (t41 * qJD(1) + t77 * t91) * t78 - (t43 * qJD(1) + t77 * t92) * t76;
t6 = t18 + t94;
t5 = t95 + t125;
t1 = [(t18 + ((t109 - t14 + t144) * t79 + (-t170 * t77 - t13 - t35) * t77) * qJD(2)) * t117 + (t101 * qJD(2) + t56 * t78 + t57 * t76) * qJD(1) - (t8 + t9) * t123 / 0.2e1 + (t5 - t125 + ((t113 * t79 - t16 + t35) * t79 + (t113 * t77 + t144 - t15 - t36) * t77) * qJD(2)) * t114 + (t12 * (-t156 * t77 + t110) - t23 * t70 + t11 * (-t156 * t79 + t111) + t168 * qJD(2) + ((t23 * t155 + t24 * t156) * t77 + (-t24 * t105 + t23 * t156) * t79) * qJD(1)) * m(3) + (t6 + t10 + t7) * t126 / 0.2e1 + ((t21 + t19) * t77 + (t20 + t22) * t79) * qJD(2) * t145; qJD(1) * (t7 * t79 - t8 * t77 + (-t19 * t77 - t20 * t79) * qJD(1)) / 0.2e1 + ((t45 * t126 - t124) * t79 + (t88 + (-t79 * t46 + t81) * qJD(2)) * t77) * t117 + ((t140 * t78 + t141 * t76) * qJD(1) + (-t160 * t76 + t89 * t78) * qJD(2)) * t145 + ((t46 * t123 + t124) * t77 + (t88 + (-t77 * t45 + t81) * qJD(2)) * t79) * t114 - (qJD(1) * t9 + (t161 * t169 - t162 * t77 ^ 2 + (-t15 * t77 - t16 * t79) * qJD(1)) * t167) * t77 / 0.2e1 + (t10 * qJD(1) + (-t161 * t79 ^ 2 + t162 * t169 + (-t13 * t77 - t14 * t79) * qJD(1)) * t167) * t79 / 0.2e1 - (t5 + t95) * t128 / 0.2e1 - (t6 + t94) * t127 / 0.2e1 + (0.2e1 * t17 * (-t79 * t31 - t77 * t32 + (-t79 * t44 - t77 * t93) * qJD(1)) - t168 * qJD(1) - (t17 * (t51 * t77 + t52 * t79) - t164 * t105) * qJD(2) - t164 * t58 - (t11 * t79 + t12 * t77 - t127 * t23 - t128 * t24) * t104) * m(3);];
tauc = t1(:);
