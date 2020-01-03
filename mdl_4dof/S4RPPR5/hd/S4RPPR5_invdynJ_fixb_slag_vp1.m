% Calculate vector of inverse dynamics joint torques for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:48
% DurationCPUTime: 4.14s
% Computational Cost: add. (3115->332), mult. (7147->425), div. (0->0), fcn. (7346->6), ass. (0->180)
t233 = sin(pkin(6));
t234 = cos(pkin(6));
t249 = sin(qJ(1));
t250 = cos(qJ(1));
t103 = -t233 * t249 - t234 * t250;
t104 = t250 * t233 - t249 * t234;
t147 = sin(qJ(4));
t148 = cos(qJ(4));
t231 = Icges(5,4) * t148;
t176 = -Icges(5,2) * t147 + t231;
t54 = Icges(5,6) * t103 + t104 * t176;
t232 = Icges(5,4) * t147;
t178 = Icges(5,1) * t148 - t232;
t57 = Icges(5,5) * t103 + t104 * t178;
t269 = t147 * t54 - t148 * t57;
t174 = Icges(5,5) * t148 - Icges(5,6) * t147;
t51 = Icges(5,3) * t103 + t104 * t174;
t14 = t103 * t51 - t104 * t269;
t53 = Icges(5,3) * t104 - t103 * t174;
t273 = t104 * t53;
t188 = rSges(5,1) * t147 + rSges(5,2) * t148;
t264 = qJD(4) * t188;
t272 = t104 * t264;
t271 = t264 * t103;
t182 = -t147 * t57 - t148 * t54;
t145 = t250 * pkin(2);
t216 = t250 * pkin(1) + t249 * qJ(2);
t257 = t145 + t216;
t77 = -t103 * rSges(4,1) - t104 * rSges(4,2);
t268 = t77 + t257;
t173 = Icges(5,5) * t147 + Icges(5,6) * t148;
t65 = t173 * t103;
t64 = t173 * t104;
t205 = t249 * pkin(2);
t206 = t249 * pkin(1);
t162 = -t206 - t205;
t267 = t162 + t205;
t141 = t250 * qJ(2);
t120 = t206 - t141;
t138 = qJD(2) * t249;
t201 = qJD(1) * t250;
t219 = qJ(2) * t201 + t138;
t266 = qJD(1) * t120 - t138 + t219;
t225 = t103 * t148;
t226 = t103 * t147;
t62 = -rSges(5,1) * t225 + rSges(5,2) * t226 + t104 * rSges(5,3);
t240 = -t103 * pkin(3) + pkin(5) * t104 + t62;
t265 = t257 + t240;
t93 = t103 * qJD(1);
t94 = t104 * qJD(1);
t263 = t94 * pkin(3) + t93 * pkin(5);
t139 = qJD(2) * t250;
t262 = -qJD(1) * t257 + t139;
t261 = t104 * rSges(4,1) - t103 * rSges(4,2);
t175 = Icges(5,2) * t148 + t232;
t177 = Icges(5,1) * t147 + t231;
t171 = -t147 * t175 + t148 * t177;
t37 = t104 * t171 + t65;
t260 = qJD(1) * t37;
t259 = t104 * pkin(3) + t103 * pkin(5);
t217 = t250 * rSges(3,1) + t249 * rSges(3,3);
t85 = t216 + t217;
t223 = t104 * t148;
t224 = t104 * t147;
t95 = t103 * rSges(5,3);
t60 = rSges(5,1) * t223 - rSges(5,2) * t224 + t95;
t255 = t60 + t259;
t242 = t175 * t104 - t57;
t244 = -t177 * t104 - t54;
t254 = t147 * t242 + t148 * t244;
t149 = qJD(1) ^ 2;
t72 = qJD(4) * t93 + qJDD(4) * t104;
t253 = t72 / 0.2e1;
t73 = qJD(4) * t94 - qJDD(4) * t103;
t252 = t73 / 0.2e1;
t251 = t94 * pkin(5);
t247 = t94 * rSges(5,3);
t59 = Icges(5,5) * t104 - t103 * t178;
t246 = -t103 * t53 - t59 * t223;
t245 = -t59 * t225 + t273;
t56 = Icges(5,6) * t104 - t103 * t176;
t243 = -t177 * t103 + t56;
t241 = t175 * t103 + t59;
t235 = t148 * t94;
t239 = rSges(5,1) * t235 + t93 * rSges(5,3);
t238 = t94 * rSges(4,1) - t93 * rSges(4,2);
t237 = t147 * t56;
t236 = t147 * t94;
t222 = -t175 + t178;
t221 = -t176 - t177;
t220 = qJD(1) * t139 + qJDD(2) * t249;
t215 = qJD(4) * t103;
t214 = qJD(4) * t104;
t123 = -rSges(5,1) * t148 + rSges(5,2) * t147;
t108 = t123 * qJD(4);
t213 = qJD(4) * t108;
t211 = qJD(4) * t147;
t210 = qJD(4) * t148;
t209 = t174 * qJD(1);
t208 = -t250 / 0.2e1;
t207 = t249 / 0.2e1;
t204 = t249 * rSges(3,1);
t202 = t103 * t211;
t200 = qJD(1) * t249;
t199 = -t215 / 0.2e1;
t198 = t215 / 0.2e1;
t197 = -t214 / 0.2e1;
t196 = t214 / 0.2e1;
t195 = -t120 - t205;
t191 = t93 * rSges(4,1) + t94 * rSges(4,2);
t180 = t147 * t59 + t148 * t56;
t166 = t202 + t235;
t167 = t103 * t210 - t236;
t30 = Icges(5,4) * t166 + Icges(5,2) * t167 + Icges(5,6) * t93;
t32 = Icges(5,1) * t166 + Icges(5,4) * t167 + Icges(5,5) * t93;
t152 = qJD(4) * t180 + t147 * t30 - t148 * t32;
t164 = t104 * t211 - t148 * t93;
t165 = t104 * t210 + t147 * t93;
t29 = Icges(5,4) * t164 + Icges(5,2) * t165 + Icges(5,6) * t94;
t31 = Icges(5,1) * t164 + Icges(5,4) * t165 + Icges(5,5) * t94;
t153 = qJD(4) * t182 + t147 * t29 - t148 * t31;
t179 = -t148 * t59 + t237;
t27 = Icges(5,5) * t164 + Icges(5,6) * t165 + Icges(5,3) * t94;
t28 = Icges(5,5) * t166 + Icges(5,6) * t167 + Icges(5,3) * t93;
t190 = -(-t103 * t27 + t104 * t153 - t269 * t93 - t51 * t94) * t103 + t104 * (-t103 * t28 + t104 * t152 + t179 * t93 + t53 * t94);
t189 = -t103 * (t103 * t153 + t104 * t27 + t269 * t94 - t51 * t93) + t104 * (t103 * t152 + t104 * t28 - t179 * t94 + t53 * t93);
t15 = t224 * t56 + t246;
t187 = -t103 * t14 + t104 * t15;
t16 = -t104 * t51 + t225 * t57 - t226 * t54;
t17 = t226 * t56 + t245;
t186 = -t103 * t16 + t104 * t17;
t61 = t104 * t123 - t95;
t168 = -t61 + t195 + t259;
t21 = qJD(1) * t168 + t138 + t271;
t22 = qJD(1) * t265 - t139 + t272;
t185 = -t103 * t21 - t104 * t22;
t33 = rSges(5,1) * t164 + rSges(5,2) * t165 + t247;
t34 = rSges(5,1) * t202 + rSges(5,2) * t167 + t239;
t184 = -t103 * t34 - t104 * t33;
t183 = t103 * t62 + t104 * t61;
t172 = -t147 * t177 - t148 * t175;
t170 = t261 + t195;
t169 = pkin(3) - t123;
t163 = -t145 * t149 + t220;
t161 = qJDD(1) * t216 - qJDD(2) * t250 + (-pkin(1) * t200 + t138 + t219) * qJD(1);
t160 = -t204 - t206;
t126 = rSges(2,1) * t250 - rSges(2,2) * t249;
t122 = rSges(2,1) * t249 + rSges(2,2) * t250;
t157 = t147 * t241 + t148 * t243;
t156 = t141 + t162;
t155 = (t147 * t221 + t148 * t222) * qJD(1);
t151 = qJDD(1) * t145 - t149 * t205 + t161;
t106 = t176 * qJD(4);
t107 = t178 * qJD(4);
t150 = qJD(4) * t172 - t106 * t147 + t107 * t148;
t143 = t250 * rSges(3,3);
t136 = rSges(3,3) * t201;
t121 = t204 - t143;
t105 = t174 * qJD(4);
t92 = qJD(1) * t216 - t139;
t71 = t188 * t103;
t70 = t188 * t104;
t49 = qJD(1) * t170 + t138;
t40 = qJDD(1) * t217 + qJD(1) * (-rSges(3,1) * t200 + t136) + t161;
t39 = -qJD(1) * t92 + (-t120 - t121) * qJDD(1) - t149 * t217 + t220;
t38 = t103 * t171 - t64;
t36 = t38 * qJD(1);
t24 = qJD(1) * t238 + qJDD(1) * t77 + t151;
t23 = (t191 - t92) * qJD(1) + t170 * qJDD(1) + t163;
t20 = qJD(4) * t183 - qJD(3);
t13 = t103 * t150 - t104 * t105 - t171 * t94 - t173 * t93;
t12 = t103 * t105 + t104 * t150 + t171 * t93 - t173 * t94;
t11 = qJD(4) * t179 - t147 * t32 - t148 * t30;
t10 = -t269 * qJD(4) - t147 * t31 - t148 * t29;
t9 = -t104 * t213 + t72 * t188 + t240 * qJDD(1) + (t263 + t34) * qJD(1) + t151;
t8 = -t103 * t213 - t73 * t188 + (t93 * pkin(3) - t251 - t33 - t92) * qJD(1) + t168 * qJDD(1) + t163;
t7 = qJD(4) * t184 - t61 * t72 + t62 * t73 + qJDD(3);
t6 = qJD(4) * t186 + t36;
t5 = qJD(4) * t187 + t260;
t1 = [-m(2) * (-g(1) * t122 + g(2) * t126) + (t171 * qJD(4) + t106 * t148 + t107 * t147) * qJD(1) + (t36 + ((t15 - t16 - t246) * t103 + t245 * t104) * qJD(4)) * t198 + (t38 - t180) * t253 + (t37 - t182) * t252 + (t11 + t13) * t196 + (t10 + t12 + t6) * t199 + (-t260 + ((t16 + (-t179 + t51) * t104) * t104 + (t17 + (t51 - t237) * t103 + t273 - t245) * t103) * qJD(4) + t5) * t197 + (-t20 * (t60 + t61) * t215 - g(1) * (t156 + t255) + (-g(2) + t9) * t265 + (t156 + (rSges(5,3) + pkin(5)) * t103 + t169 * t104) * t8 + (t169 * t93 - t247 - t251 + t262 - t272) * t21 + (-rSges(5,2) * t236 + t271 - t188 * t215 + t21 + t239 + (-t255 + t267) * qJD(1) + t263 + t266) * t22) * m(5) + ((-g(2) + t24) * t268 + (t191 + t262) * t49 + (t49 + t238 + (-t261 + t267) * qJD(1) + t266) * (qJD(1) * t268 - t139) + (-g(1) + t23) * (t156 + t261)) * m(4) + ((-g(1) + t39) * (t141 + t143 + t160) + (-g(2) + t40) * t85 + (t136 + (t121 + t160) * qJD(1) + t266) * (t85 * qJD(1) - t139)) * m(3) + (-t172 + Icges(2,3) + Icges(3,2) + Icges(4,3) + m(2) * (t122 ^ 2 + t126 ^ 2)) * qJDD(1); (-m(3) - m(4) - m(5)) * (g(1) * t249 - g(2) * t250) + 0.2e1 * (t207 * t8 + t208 * t9) * m(5) + 0.2e1 * (t207 * t23 + t208 * t24) * m(4) + 0.2e1 * (t207 * t39 + t208 * t40) * m(3); (t7 + g(3)) * m(5) + (qJDD(3) + g(3)) * m(4); t93 * t6 / 0.2e1 + t104 * (qJD(1) * t13 + t189 * qJD(4) + qJDD(1) * t38 + t16 * t73 + t17 * t72) / 0.2e1 + t186 * t253 + (t16 * t94 + t17 * t93 + t189) * t196 + t94 * t5 / 0.2e1 - t103 * (qJD(1) * t12 + t190 * qJD(4) + qJDD(1) * t37 + t14 * t73 + t15 * t72) / 0.2e1 + t187 * t252 + (t14 * t94 + t15 * t93 + t190) * t199 + qJDD(1) * (t103 * t182 - t104 * t180) / 0.2e1 + qJD(1) * (-t10 * t103 + t104 * t11 - t180 * t93 - t182 * t94) / 0.2e1 + ((t65 * t214 - t209) * t104 + (t155 + (-t254 * t103 + (-t64 + t157) * t104) * qJD(4)) * t103) * t197 + ((t64 * t215 + t209) * t103 + (t155 + (t157 * t104 + (-t254 - t65) * t103) * qJD(4)) * t104) * t198 - qJD(1) * ((t222 * t147 - t221 * t148) * qJD(1) + ((t103 * t242 - t104 * t241) * t148 + (-t103 * t244 + t104 * t243) * t147) * qJD(4)) / 0.2e1 + (-t7 * t183 + t20 * (t61 * t93 - t62 * t94 - t184) + t185 * t108 - (-t103 * t8 - t104 * t9 + t21 * t94 - t22 * t93) * t188 - (-t21 * t70 + t22 * t71) * qJD(1) - (t20 * (t103 * t71 + t104 * t70) + t185 * t123) * qJD(4) - g(1) * t71 - g(2) * t70 - g(3) * t123) * m(5);];
tau = t1;
