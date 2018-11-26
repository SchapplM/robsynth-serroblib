% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:44:24
% EndTime: 2018-11-16 16:44:26
% DurationCPUTime: 1.48s
% Computational Cost: add. (1057->168), mult. (2851->261), div. (0->0), fcn. (2266->4), ass. (0->113)
t83 = sin(qJ(1));
t85 = cos(qJ(1));
t179 = t85 * t83;
t82 = sin(qJ(2));
t84 = cos(qJ(2));
t102 = Icges(3,5) * t82 + Icges(3,6) * t84;
t47 = t102 * t83;
t48 = t102 * t85;
t113 = rSges(3,1) * t82 + rSges(3,2) * t84;
t134 = qJD(2) * t85;
t129 = t113 * t134;
t135 = qJD(2) * t83;
t130 = t113 * t135;
t146 = Icges(3,4) * t84;
t105 = -Icges(3,2) * t82 + t146;
t40 = Icges(3,6) * t85 + t105 * t83;
t178 = 0.2e1 * qJD(2);
t143 = Icges(3,5) * t85;
t148 = Icges(3,1) * t84;
t147 = Icges(3,4) * t82;
t75 = t83 * t147;
t42 = t148 * t83 + t143 - t75;
t163 = t42 * t84;
t111 = t40 * t82 - t163;
t138 = Icges(3,3) * t85;
t140 = Icges(3,6) * t83;
t144 = Icges(3,5) * t84;
t38 = -t140 * t82 + t144 * t83 + t138;
t13 = -t111 * t83 + t85 * t38;
t103 = -Icges(3,6) * t82 + t144;
t104 = Icges(3,2) * t84 + t147;
t106 = Icges(3,1) * t82 + t146;
t158 = t84 * t106;
t108 = -t104 * t82 + t158;
t176 = t108 * qJD(1) + t103 * qJD(2);
t145 = Icges(3,5) * t83;
t76 = t85 * t147;
t43 = -t148 * t85 + t145 + t76;
t162 = t43 * t84;
t41 = -t105 * t85 + t140;
t109 = t41 * t82 - t162;
t96 = qJD(2) * t102;
t175 = -t85 * t96 + (-t103 * t83 + t109 - t138) * qJD(1);
t39 = Icges(3,3) * t83 - t103 * t85;
t174 = t83 * t96 + (t111 + t39) * qJD(1);
t173 = (t43 + t76) * t83 + (t42 - t75) * t85;
t170 = pkin(1) * qJD(1) ^ 2;
t169 = rSges(3,3) + pkin(1);
t167 = rSges(3,2) * t82;
t165 = (t40 * t84 + t42 * t82) * t83;
t159 = t83 * t84;
t131 = rSges(3,1) * t159;
t161 = t82 * t83;
t77 = rSges(3,2) * t161;
t117 = -t77 + t131;
t44 = rSges(3,3) * t85 + t117;
t136 = qJD(1) * t85;
t80 = pkin(1) * t136;
t24 = qJD(1) * t44 + t129 + t80;
t164 = t24 * t85;
t160 = t82 * t85;
t81 = t83 * rSges(3,3);
t157 = t84 * t85;
t155 = qJD(1) / 0.2e1;
t154 = t40 * t160 + t83 * t38;
t153 = t41 * t160 + t83 * t39;
t107 = -t147 + t148;
t150 = -t104 + t107;
t149 = -t105 - t106;
t133 = (t104 * t161 - t158 * t83 - t48) * qJD(1);
t132 = t103 * qJD(1);
t45 = -rSges(3,1) * t157 + rSges(3,2) * t160 + t81;
t128 = pkin(1) * t83 + t45;
t125 = -t135 / 0.2e1;
t122 = t134 / 0.2e1;
t33 = t41 * t161;
t121 = t85 * t39 - t33;
t120 = t38 + t162;
t118 = rSges(3,3) * t136 + qJD(1) * t131 + t129;
t114 = rSges(3,1) * t84 - t167;
t23 = qJD(1) * t128 + t130;
t112 = -t23 * t83 - t164;
t53 = t113 * t83;
t14 = -t159 * t43 - t121;
t100 = (-t13 * t85 + t14 * t83) * qJD(2);
t15 = t157 * t42 - t154;
t16 = -t157 * t43 + t153;
t99 = (-t15 * t85 + t16 * t83) * qJD(2);
t17 = (-t44 * t83 + t45 * t85) * qJD(2);
t98 = qJD(2) * t106;
t97 = qJD(2) * t104;
t95 = t40 * t85 + t41 * t83;
t94 = (t149 * t82 + t150 * t84) * qJD(1);
t56 = t105 * qJD(2);
t57 = t107 * qJD(2);
t88 = -qJD(1) * t102 - t56 * t82 + t57 * t84 + (-t104 * t84 - t106 * t82) * qJD(2);
t87 = t173 * t82 + t95 * t84;
t58 = t114 * qJD(2);
t54 = t113 * t85;
t32 = qJD(2) * t53 + (-t114 * t85 + t81) * qJD(1);
t31 = -qJD(1) * t77 + t118;
t22 = t108 * t85 - t47;
t20 = t41 * t84 + t43 * t82;
t18 = t22 * qJD(1);
t12 = t85 * t170 + t58 * t135 + (t31 + t129) * qJD(1);
t11 = -t83 * t170 + t58 * t134 + (-t32 - t130) * qJD(1);
t10 = t176 * t85 + t88 * t83;
t9 = -t176 * t83 + t88 * t85;
t8 = qJD(2) * t109 - (qJD(1) * t40 + t85 * t97) * t84 - (t85 * t98 + (t107 * t83 + t143) * qJD(1)) * t82;
t7 = -t111 * qJD(2) - (qJD(1) * t41 + t83 * t97) * t84 - (t83 * t98 + (-t107 * t85 + t145) * qJD(1)) * t82;
t6 = t18 + t99;
t5 = t100 - t133;
t1 = [(t18 + ((-t33 + t14 + (t39 - t163) * t85 + t154) * t85 + t153 * t83) * qJD(2)) * t122 + (qJD(2) * t108 + t56 * t84 + t57 * t82) * qJD(1) + m(3) * (t12 * t128 + t23 * (t80 + t118) + t11 * (t169 * t85 + t117) - t24 * t130 + (t114 * t164 + (-t167 * t23 - t169 * t24) * t83) * qJD(1)) + (t8 + t9) * t135 / 0.2e1 + (t5 + 0.2e1 * t133 + ((t120 * t85 - t153 + t16) * t85 + (t120 * t83 + t121 + t15) * t83) * qJD(2)) * t125 - (t6 + t10 + t7) * t134 / 0.2e1 + (t165 + (-t20 + t22) * t85) * qJD(2) * t155; ((t134 * t47 + t132) * t85 + (t94 + (-t85 * t48 + t87) * qJD(2)) * t83) * t122 - qJD(1) * ((-t149 * t84 + t150 * t82) * qJD(1) + (-t173 * t84 + t95 * t82) * qJD(2)) / 0.2e1 + ((t135 * t48 - t132) * t83 + (t94 + (-t83 * t47 + t87) * qJD(2)) * t85) * t125 + (t9 * qJD(1) + (-t174 * t179 - t175 * t83 ^ 2 + (t15 * t83 + t16 * t85) * qJD(1)) * t178) * t83 / 0.2e1 - (t10 * qJD(1) + (t174 * t85 ^ 2 + t175 * t179 + (t13 * t83 + t14 * t85) * qJD(1)) * t178) * t85 / 0.2e1 + (t6 + t99) * t136 / 0.2e1 + (0.2e1 * t17 * (t85 * t31 + t83 * t32 + (-t44 * t85 - t45 * t83) * qJD(1)) - t112 * t58 - (-t11 * t85 - t12 * t83 + (-t23 * t85 + t24 * t83) * qJD(1)) * t113 - (t23 * t54 - t24 * t53) * qJD(1) - (t17 * (t53 * t83 + t54 * t85) - t112 * t114) * qJD(2)) * m(3) + (-t7 * t85 + (-t20 * t85 + t165) * qJD(1) + (t8 + t5 + t100) * t83) * t155;];
tauc  = t1(:);
