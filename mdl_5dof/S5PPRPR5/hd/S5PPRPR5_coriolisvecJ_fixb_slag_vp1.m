% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRPR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:29
% DurationCPUTime: 3.21s
% Computational Cost: add. (2920->276), mult. (6943->388), div. (0->0), fcn. (7308->6), ass. (0->149)
t193 = sin(pkin(7));
t194 = cos(pkin(7));
t215 = sin(qJ(3));
t216 = cos(qJ(3));
t104 = -t193 * t215 - t194 * t216;
t105 = -t193 * t216 + t194 * t215;
t124 = sin(qJ(5));
t125 = cos(qJ(5));
t110 = Icges(6,5) * t124 + Icges(6,6) * t125;
t49 = -Icges(6,3) * t104 + t105 * t110;
t199 = t104 * t49;
t123 = Icges(6,4) * t124;
t143 = Icges(6,2) * t125 + t123;
t52 = -Icges(6,6) * t104 + t105 * t143;
t234 = t124 * t52;
t233 = t125 * t52;
t103 = t105 * pkin(3);
t159 = -qJ(4) * t104 - t103;
t178 = qJD(4) * t105;
t232 = qJD(3) * t159 + t178;
t142 = Icges(6,5) * t125 - Icges(6,6) * t124;
t64 = t142 * t104;
t63 = t142 * t105;
t173 = qJD(5) * t125;
t94 = qJD(3) * t105;
t137 = -t104 * t173 - t124 * t94;
t191 = Icges(6,4) * t125;
t113 = Icges(6,2) * t124 - t191;
t114 = Icges(6,1) * t124 + t191;
t179 = t113 - t114;
t115 = -Icges(6,1) * t125 + t123;
t180 = t115 + t143;
t231 = (t124 * t179 - t125 * t180) * qJD(3);
t230 = 0.2e1 * qJD(5);
t99 = t105 * qJ(4);
t78 = -t104 * pkin(3) + t99;
t75 = t105 * rSges(5,2) - t104 * rSges(5,3);
t229 = qJD(3) * t75 + t232;
t176 = qJD(5) * t105;
t154 = rSges(6,1) * t125 - t124 * rSges(6,2);
t175 = qJD(5) * t154;
t95 = qJD(3) * t104;
t228 = t95 * pkin(3) - qJ(4) * t94;
t140 = t125 * t113 + t124 * t115;
t35 = t105 * t140 + t64;
t153 = t124 * rSges(6,1) + rSges(6,2) * t125;
t58 = -t104 * rSges(6,3) + t153 * t105;
t107 = t143 * qJD(5);
t108 = t114 * qJD(5);
t225 = qJD(5) * (t113 * t124 - t115 * t125) - t107 * t125 - t108 * t124;
t184 = t104 * t125;
t185 = t104 * t124;
t60 = -rSges(6,1) * t185 - rSges(6,2) * t184 - t105 * rSges(6,3);
t224 = -(pkin(6) * t105 - t60) * qJD(3) + t105 * t175 + t232;
t220 = pkin(6) * t95;
t217 = -rSges(5,2) + pkin(3);
t214 = rSges(5,3) * t95;
t98 = t104 * qJD(4);
t138 = t228 - t98;
t211 = qJD(3) * t138 - qJD(4) * t95;
t210 = t105 * t115 + t52;
t54 = -Icges(6,6) * t105 - t104 * t143;
t209 = -t104 * t115 + t54;
t183 = t105 * t124;
t190 = Icges(6,5) * t104;
t182 = t105 * t125;
t86 = Icges(6,4) * t182;
t55 = Icges(6,1) * t183 - t190 + t86;
t208 = -Icges(6,2) * t183 + t55 + t86;
t87 = Icges(6,4) * t184;
t57 = -Icges(6,1) * t185 - Icges(6,5) * t105 - t87;
t207 = Icges(6,2) * t185 + t57 - t87;
t206 = -qJD(3) * t78 - t98;
t205 = -t95 * qJ(4) + t178;
t204 = -t95 * rSges(5,2) - t94 * rSges(5,3);
t203 = m(4) * qJD(3);
t202 = rSges(5,3) * t105;
t16 = -t105 * t49 - t52 * t184 - t55 * t185;
t200 = t104 * t16;
t197 = t124 * t95;
t196 = t125 * t94;
t195 = t125 * t95;
t36 = -t104 * t140 + t63;
t181 = t36 * qJD(3);
t177 = qJD(5) * t104;
t174 = qJD(5) * t124;
t172 = t110 * qJD(3);
t171 = rSges(6,3) + pkin(3) + pkin(6);
t51 = -Icges(6,3) * t105 - t104 * t110;
t15 = -t104 * t51 + t54 * t182 + t57 * t183;
t170 = t105 * t51 + t54 * t184 + t57 * t185;
t169 = t104 * t174;
t164 = t177 / 0.2e1;
t163 = t176 / 0.2e1;
t162 = t137 * rSges(6,1) - rSges(6,2) * t196 + t95 * rSges(6,3);
t160 = qJD(2) * t194;
t155 = rSges(4,1) * t104 + rSges(4,2) * t105;
t152 = -rSges(5,2) * t104 - t202;
t122 = qJD(2) * t193;
t18 = t122 + t224;
t145 = -t98 - t160;
t19 = -t104 * t175 + (t104 * pkin(6) - t58 - t78) * qJD(3) + t145;
t151 = t104 * t19 - t105 * t18;
t150 = -t104 * t60 + t105 * t58;
t147 = t124 * t55 + t233;
t23 = -t125 * t55 + t234;
t146 = t124 * t57 + t125 * t54;
t24 = t124 * t54 - t125 * t57;
t136 = t169 - t196;
t135 = t105 * t173 - t197;
t134 = -t105 * t174 - t195;
t14 = t105 * t147 - t199;
t133 = qJD(5) * (-t104 * t14 - t105 * t15);
t132 = qJD(5) * (t105 * t170 - t200);
t131 = t124 * t210 - t125 * t208;
t130 = -t124 * t209 + t125 * t207;
t109 = t153 * qJD(5);
t32 = rSges(6,1) * t135 + rSges(6,2) * t134 - t94 * rSges(6,3);
t48 = -t94 * pkin(3) + t205;
t83 = qJD(4) * t94;
t10 = -t83 + (t104 * t109 - t154 * t94) * qJD(5) + (t94 * pkin(6) - t32 - t48) * qJD(3);
t31 = rSges(6,2) * t169 + t162;
t11 = (-t105 * t109 - t154 * t95) * qJD(5) + (t31 + t220) * qJD(3) + t211;
t129 = t10 * t104 - t105 * t11 + t18 * t95 + t19 * t94;
t128 = -t104 * t31 + t105 * t32 - t58 * t95 - t60 * t94;
t106 = t110 * qJD(5);
t77 = -rSges(4,1) * t105 + rSges(4,2) * t104;
t72 = t154 * t104;
t71 = t154 * t105;
t70 = rSges(4,1) * t95 + rSges(4,2) * t94;
t69 = -rSges(4,1) * t94 + rSges(4,2) * t95;
t62 = qJD(3) * t155 - t160;
t61 = qJD(3) * t77 + t122;
t56 = -t105 * t114 + t190;
t34 = (t152 - t78) * qJD(3) + t145;
t33 = t122 + t229;
t26 = Icges(6,5) * t135 + Icges(6,6) * t134 - Icges(6,3) * t94;
t25 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t95;
t22 = qJD(3) * t204 + t211;
t21 = -t83 + (-rSges(5,2) * t94 + t214 - t48) * qJD(3);
t20 = qJD(5) * t150 + qJD(1);
t13 = -t104 * t106 - t105 * t225 - t140 * t95 + t142 * t94;
t12 = t104 * t225 - t105 * t106 - t140 * t94 - t142 * t95;
t9 = qJD(5) * t146 + t124 * (Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t95) - t125 * (Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t95);
t8 = qJD(5) * t147 + t124 * (Icges(6,4) * t135 + Icges(6,2) * t134 - Icges(6,6) * t94) - t125 * (Icges(6,1) * t135 + Icges(6,4) * t134 - Icges(6,5) * t94);
t7 = t128 * qJD(5);
t6 = t132 - t181;
t5 = -qJD(3) * t35 + t133;
t1 = [m(6) * t7; (t193 * t70 + t194 * t69) * t203 + m(5) * (t193 * t22 - t194 * t21) + m(6) * (-t10 * t194 + t11 * t193); (qJD(5) * t140 + t124 * t107 - t108 * t125) * qJD(3) - t6 * t177 / 0.2e1 + m(4) * (t61 * t70 - t62 * t69 + (-t155 * t69 + t70 * t77) * qJD(3)) - (t155 * t61 - t62 * t77) * t203 + (t11 * (-t103 + t60) - t10 * t99 + (-t11 * pkin(6) - t10 * t153) * t105 + (rSges(6,1) * t197 + rSges(6,2) * t195 - t154 * t176 + t171 * t94 - t205 + t224) * t19 + (qJD(3) * t58 + t162 - t206 + t220 + t228) * t18 + (-t11 * qJ(4) + t10 * t171 + (rSges(6,2) * t174 - pkin(6) * qJD(3) - qJD(4) + t175) * t18) * t104) * m(6) + (t22 * (t159 + t75) + t21 * (t104 * t217 - t202 - t99) + (t217 * t94 - t205 + t214 + t229) * t34 + (-qJD(3) * t152 + t138 + t204 - t206) * t33) * m(5) + (-t181 + (-t200 + (-t199 - t14 + (-t124 * t56 + t233) * t105 + t170) * t105) * qJD(5) + t8 + t13) * t164 + (t9 + t12 + t5 + (-t199 * t104 + (-t16 + (t146 - t49) * t105 + (-t51 + (t55 + t56) * t124) * t104) * t105) * qJD(5) + (t125 * t56 - t23 + t234 + t35) * qJD(3)) * t163 + ((t23 + t35) * t94 - (t24 + t36) * t95) * qJD(5) / 0.2e1; m(5) * (-t104 * t21 + t105 * t22 - t33 * t95 - t34 * t94) - m(6) * t129 + 0.2e1 * (-m(5) * (-t104 * t33 - t105 * t34) / 0.2e1 - m(6) * (-t104 * t18 - t105 * t19) / 0.2e1) * qJD(3); -qJD(3) * (-t104 * t8 - t105 * t9 - t23 * t94 + t24 * t95) / 0.2e1 + ((t63 * t177 + t172) * t104 + (t231 + (-t130 * t105 + (-t64 + t131) * t104) * qJD(5)) * t105) * t164 + ((-t64 * t176 + t172) * t105 + (-t231 + (-t131 * t104 + (t63 + t130) * t105) * qJD(5)) * t104) * t163 + qJD(3) * (-(t180 * t124 + t179 * t125) * qJD(3) + ((-t104 * t210 - t105 * t209) * t125 + (-t104 * t208 - t105 * t207) * t124) * qJD(5)) / 0.2e1 - (-qJD(3) * t13 + (-t104 * (-t104 * t26 - t147 * t95 - t94 * t49) - t105 * (-t104 * t25 - t146 * t95 - t94 * t51) - t14 * t94 + t15 * t95) * t230) * t104 / 0.2e1 - (-qJD(3) * t12 + (-(-t105 * t26 - t147 * t94 + t95 * t49) * t104 - t105 * (-t105 * t25 - t146 * t94 + t95 * t51) - t16 * t94 - t170 * t95) * t230) * t105 / 0.2e1 - (t5 + t133) * t94 / 0.2e1 + (t6 + t132) * t95 / 0.2e1 + (t151 * t109 - t129 * t154 + t20 * t128 + t7 * t150 - (-t18 * t72 - t19 * t71) * qJD(3) - (t20 * (t104 * t72 + t105 * t71) + t151 * t153) * qJD(5)) * m(6);];
tauc = t1(:);
