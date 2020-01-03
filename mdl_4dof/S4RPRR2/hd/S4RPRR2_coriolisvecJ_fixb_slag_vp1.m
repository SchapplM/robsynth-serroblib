% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:11
% DurationCPUTime: 2.58s
% Computational Cost: add. (5467->232), mult. (4244->330), div. (0->0), fcn. (3240->8), ass. (0->149)
t137 = qJD(1) ^ 2;
t245 = 2 * qJD(4);
t131 = qJD(1) + qJD(3);
t132 = qJ(1) + pkin(7);
t128 = qJ(3) + t132;
t124 = sin(t128);
t133 = sin(qJ(4));
t206 = t124 * t133;
t101 = rSges(5,2) * t206;
t125 = cos(t128);
t201 = t125 * rSges(5,3) + t101;
t135 = cos(qJ(4));
t205 = t124 * t135;
t65 = rSges(5,1) * t205 - t201;
t58 = t131 * t65;
t120 = t125 * pkin(6);
t83 = pkin(3) * t124 - t120;
t244 = -t131 * t83 - t58;
t202 = t125 * t135;
t243 = rSges(5,1) * t202 + t124 * rSges(5,3);
t242 = -t125 * pkin(3) - t124 * pkin(6);
t241 = pkin(2) * cos(t132) + cos(qJ(1)) * pkin(1);
t129 = Icges(5,4) * t135;
t166 = -Icges(5,2) * t133 + t129;
t107 = Icges(5,1) * t133 + t129;
t203 = t125 * t133;
t193 = rSges(5,2) * t203;
t66 = -t193 + t243;
t152 = t66 - t242;
t220 = rSges(5,2) * t135;
t109 = rSges(5,1) * t133 + t220;
t197 = qJD(4) * t124;
t80 = t109 * t197;
t239 = -t131 * t152 + t80;
t104 = Icges(5,5) * t135 - Icges(5,6) * t133;
t212 = Icges(5,4) * t133;
t105 = Icges(5,2) * t135 + t212;
t165 = t133 * t105 - t135 * t107;
t238 = t104 * qJD(4) + t131 * t165;
t103 = Icges(5,5) * t133 + Icges(5,6) * t135;
t146 = Icges(5,3) * t131 - qJD(4) * t103;
t155 = t166 * t125;
t62 = Icges(5,6) * t124 + t155;
t215 = t133 * t62;
t108 = Icges(5,1) * t135 - t212;
t156 = t108 * t125;
t64 = Icges(5,5) * t124 + t156;
t168 = -t135 * t64 + t215;
t207 = t124 * t131;
t237 = -t104 * t207 + t125 * t146 + t131 * t168;
t154 = t104 * t125;
t61 = Icges(5,4) * t205 - Icges(5,2) * t206 - Icges(5,6) * t125;
t216 = t133 * t61;
t100 = Icges(5,4) * t206;
t63 = Icges(5,1) * t205 - Icges(5,5) * t125 - t100;
t169 = -t135 * t63 + t216;
t236 = t124 * t146 + (t154 + t169) * t131;
t59 = Icges(5,5) * t205 - Icges(5,6) * t206 - Icges(5,3) * t125;
t24 = -t124 * t169 - t125 * t59;
t224 = -Icges(5,2) * t205 - t100 + t63;
t226 = t107 * t124 + t61;
t234 = -t133 * t224 - t135 * t226;
t231 = t131 / 0.2e1;
t161 = t241 * qJD(1);
t118 = t125 * rSges(4,1);
t82 = -rSges(4,2) * t124 + t118;
t57 = t131 * t82 + t161;
t81 = rSges(4,1) * t124 + rSges(4,2) * t125;
t229 = t57 * t81;
t228 = -t124 * t59 - t63 * t202;
t60 = Icges(5,3) * t124 + t154;
t227 = t124 * t60 + t64 * t202;
t225 = -t107 * t125 - t62;
t223 = -t105 * t125 + t64;
t222 = rSges(5,1) * t135;
t174 = -pkin(2) * sin(t132) - sin(qJ(1)) * pkin(1);
t162 = t174 * qJD(1);
t196 = qJD(4) * t125;
t191 = t109 * t196;
t145 = t162 - t191;
t28 = (-t65 - t83) * t131 + t145;
t219 = t125 * t28;
t217 = t131 * t81;
t209 = t103 * t125;
t47 = -t124 * t165 - t209;
t214 = t47 * t131;
t210 = t103 * t124;
t208 = t104 * t131;
t204 = t125 * t131;
t200 = -t105 + t108;
t199 = t107 + t166;
t195 = qJD(4) * t133;
t192 = qJD(4) * t220;
t194 = t131 * t193 + (rSges(5,1) * t195 + t192) * t124;
t190 = t125 * t195;
t187 = -pkin(3) - t222;
t186 = -t197 / 0.2e1;
t183 = t196 / 0.2e1;
t51 = t64 * t205;
t181 = t125 * t60 - t51;
t160 = rSges(5,3) * t204 + t131 * t101 - t125 * t192;
t44 = (-t131 * t205 - t190) * rSges(5,1) + t160;
t180 = t44 + t58;
t45 = t131 * t243 - t194;
t179 = -t131 * t66 + t45;
t178 = -t59 + t215;
t170 = -rSges(5,2) * t133 + t222;
t35 = t133 * t63 + t135 * t61;
t36 = t133 * t64 + t135 * t62;
t164 = t174 * t137;
t163 = t241 * t137;
t25 = -t206 * t62 - t181;
t158 = (t124 * t25 - t125 * t24) * qJD(4);
t26 = -t203 * t61 - t228;
t27 = -t203 * t62 + t227;
t157 = (t124 * t27 - t125 * t26) * qJD(4);
t151 = -t133 * t223 + t135 * t225;
t150 = t124 * t187 + t120 + t201;
t56 = t162 - t217;
t149 = (-t133 * t199 + t135 * t200) * t131;
t148 = Icges(5,5) * t131 - qJD(4) * t107;
t147 = Icges(5,6) * t131 - qJD(4) * t105;
t48 = -t125 * t165 + t210;
t46 = t48 * t131;
t10 = t46 + t157;
t14 = -qJD(4) * t169 + t133 * (t124 * t148 + t131 * t156) + t135 * (t124 * t147 + t131 * t155);
t15 = -qJD(4) * t168 + t133 * (-t108 * t207 + t125 * t148) + t135 * (t125 * t147 - t166 * t207);
t89 = t166 * qJD(4);
t90 = t108 * qJD(4);
t139 = t103 * t131 - t133 * t89 + t135 * t90 + (-t105 * t135 - t107 * t133) * qJD(4);
t18 = t124 * t238 + t139 * t125;
t19 = t139 * t124 - t125 * t238;
t9 = t158 + t214;
t144 = (t46 + ((t25 - t51 + (t60 + t216) * t125 + t228) * t125 + t227 * t124) * qJD(4)) * t183 + (-qJD(4) * t165 + t133 * t90 + t135 * t89) * t131 + (t9 - t214 + ((t125 * t178 - t227 + t27) * t125 + (t124 * t178 + t181 + t26) * t124) * qJD(4)) * t186 + (t15 + t18) * t197 / 0.2e1 - (t10 + t14 + t19) * t196 / 0.2e1 + ((t35 + t47) * t124 + (t36 + t48) * t125) * qJD(4) * t231;
t29 = t161 - t239;
t97 = pkin(6) * t204;
t138 = t28 * t194 + t29 * (-rSges(5,1) * t190 + t160 + t97) + (t187 * t219 + (t28 * (-rSges(5,3) - pkin(6)) + t29 * t187) * t124) * t131;
t95 = rSges(4,2) * t207;
t91 = t170 * qJD(4);
t78 = t109 * t125;
t77 = t109 * t124;
t68 = rSges(4,1) * t204 - t95;
t50 = -t131 * t68 - t163;
t49 = -t131 * t217 + t164;
t32 = qJD(2) + (t124 * t65 + t125 * t66) * qJD(4);
t23 = -t91 * t196 - t163 + (t131 * t242 - t45 + t80) * t131;
t22 = -t91 * t197 + t164 + (-pkin(3) * t207 - t191 + t44 + t97) * t131;
t11 = (t179 * t124 + t180 * t125) * qJD(4);
t1 = [m(4) * (t50 * (t174 - t81) + t56 * t95 + t49 * (t82 + t241) + (-t118 * t56 - t229) * t131 + (t174 * t57 - t241 * t56) * qJD(1)) + t144 + (-(-t28 + t145 + t244) * t29 + t23 * (t150 + t174) + t22 * (t152 + t241) + (t174 * t29 - t241 * t28) * qJD(1) + t138) * m(5); m(5) * t11; t144 + (-t28 * t239 - t29 * (-t191 + t244) + t23 * t150 + t22 * t152 + t138) * m(5) + (t49 * t82 - t50 * t81 - t56 * t68 - t57 * t217 - (-t56 * t82 - t229) * t131) * m(4); ((t131 * t36 - t14) * t125 + (t131 * t35 + t15) * t124) * t231 + ((-t197 * t209 + t208) * t124 + (t149 + (-t234 * t125 + (t210 + t151) * t124) * qJD(4)) * t125) * t186 + ((-t196 * t210 - t208) * t125 + (t149 + (t151 * t124 + (-t234 + t209) * t125) * qJD(4)) * t124) * t183 - t131 * ((t200 * t133 + t199 * t135) * t131 + ((t124 * t223 - t125 * t224) * t135 + (t124 * t225 + t125 * t226) * t133) * qJD(4)) / 0.2e1 + (t131 * t18 + ((-t124 * t236 + t131 * t27) * t125 + (t124 * t237 + t131 * t26) * t124) * t245) * t124 / 0.2e1 - (t131 * t19 + ((t125 * t236 + t131 * t25) * t125 + (-t125 * t237 + t131 * t24) * t124) * t245) * t125 / 0.2e1 + (t9 + t158) * t207 / 0.2e1 + (t10 + t157) * t204 / 0.2e1 + ((t11 * t66 + t32 * t180 - t28 * t91) * t125 + (t11 * t65 + t32 * t179 - t29 * t91) * t124 + ((-t131 * t29 - t23) * t125 + (t131 * t28 - t22) * t124) * t109 - (t28 * t77 - t29 * t78) * t131 - (t32 * (-t124 * t77 - t125 * t78) + (-t124 * t29 - t219) * t170) * qJD(4)) * m(5);];
tauc = t1(:);
