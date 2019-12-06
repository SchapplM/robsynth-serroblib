% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:20
% EndTime: 2019-12-05 17:49:28
% DurationCPUTime: 4.02s
% Computational Cost: add. (8428->320), mult. (5174->404), div. (0->0), fcn. (3828->10), ass. (0->193)
t163 = qJD(1) ^ 2;
t155 = pkin(9) + qJ(5);
t149 = sin(t155);
t151 = cos(t155);
t120 = rSges(6,1) * t149 + rSges(6,2) * t151;
t303 = qJD(5) * t120;
t157 = qJ(1) + pkin(8);
t153 = qJ(3) + t157;
t145 = sin(t153);
t160 = -pkin(7) - qJ(4);
t135 = t145 * t160;
t146 = cos(t153);
t245 = qJ(4) * t145;
t159 = cos(pkin(9));
t147 = pkin(4) * t159 + pkin(3);
t276 = pkin(3) - t147;
t238 = t146 * t149;
t129 = rSges(6,2) * t238;
t237 = t146 * t151;
t191 = -rSges(6,1) * t237 - rSges(6,3) * t145;
t79 = -t129 - t191;
t305 = t146 * t276 + t135 + t245 - t79;
t156 = qJD(1) + qJD(3);
t143 = Icges(6,4) * t151;
t195 = -Icges(6,2) * t149 + t143;
t299 = Icges(6,1) * t149 + t143;
t231 = t299 + t195;
t250 = Icges(6,4) * t149;
t110 = Icges(6,2) * t151 + t250;
t113 = Icges(6,1) * t151 - t250;
t232 = t110 - t113;
t304 = (t149 * t231 + t232 * t151) * t156;
t152 = cos(t157);
t162 = cos(qJ(1));
t279 = pkin(1) * t162;
t206 = pkin(2) * t152 + t279;
t295 = t206 * qJD(1);
t150 = sin(t157);
t277 = pkin(2) * t150;
t161 = sin(qJ(1));
t280 = pkin(1) * t161;
t207 = t277 + t280;
t189 = t207 * qJD(1);
t302 = 0.2e1 * qJD(5);
t201 = rSges(4,1) * t145 + rSges(4,2) * t146;
t268 = rSges(4,1) * t146;
t103 = -t145 * rSges(4,2) + t268;
t235 = t156 * t103;
t70 = -t295 - t235;
t301 = t201 * t70;
t74 = Icges(6,4) * t237 - Icges(6,2) * t238 + Icges(6,6) * t145;
t128 = Icges(6,4) * t238;
t76 = Icges(6,1) * t237 + Icges(6,5) * t145 - t128;
t197 = t149 * t74 - t151 * t76;
t300 = t146 * t197;
t90 = t201 * t156;
t144 = t150 * rSges(3,2);
t269 = rSges(3,1) * t152;
t205 = -t269 - t279;
t298 = t144 + t205;
t138 = qJD(4) * t145;
t253 = -t138 - t156 * (-pkin(3) * t145 + qJ(4) * t146);
t171 = t189 + t253;
t226 = qJD(5) * t146;
t263 = rSges(6,2) * t149;
t265 = rSges(6,1) * t151;
t199 = -t263 + t265;
t259 = t146 * rSges(6,3);
t290 = t145 * t199 - t259;
t293 = t120 * t226 + (-(-qJ(4) - t160) * t146 - t276 * t145 + t290) * t156;
t25 = t171 + t293;
t102 = pkin(3) * t146 + t245;
t139 = qJD(4) * t146;
t177 = t139 - t295;
t227 = qJD(5) * t145;
t94 = t120 * t227;
t26 = t94 + (-t102 + t305) * t156 + t177;
t296 = t145 * t26 + t146 * t25;
t294 = t305 * t156 + t94;
t108 = Icges(6,5) * t149 + Icges(6,6) * t151;
t289 = -Icges(6,3) * t156 + qJD(5) * t108;
t288 = -Icges(6,6) * t156 + qJD(5) * t110;
t100 = t113 * qJD(5);
t99 = t195 * qJD(5);
t287 = qJD(5) * (t110 * t151 + t149 * t299) - t100 * t151 - t108 * t156 + t149 * t99;
t286 = -Icges(6,5) * t156 + qJD(5) * t299;
t270 = -Icges(6,2) * t237 - t128 + t76;
t272 = t146 * t299 + t74;
t284 = t149 * t270 + t151 * t272;
t241 = t145 * t149;
t127 = Icges(6,4) * t241;
t240 = t145 * t151;
t75 = -Icges(6,1) * t240 + Icges(6,5) * t146 + t127;
t271 = Icges(6,2) * t240 + t127 + t75;
t73 = Icges(6,6) * t146 - t145 * t195;
t273 = -t145 * t299 + t73;
t283 = -t149 * t271 - t151 * t273;
t278 = pkin(1) * t163;
t109 = Icges(6,5) * t151 - Icges(6,6) * t149;
t71 = Icges(6,3) * t146 - t109 * t145;
t275 = t146 * t71 + t73 * t241;
t274 = t145 * t71 + t75 * t237;
t267 = rSges(5,1) * t159;
t264 = rSges(5,2) * sin(pkin(9));
t261 = rSges(5,3) * t146;
t257 = t149 * t73;
t256 = t151 * t75;
t193 = t149 * t110 - t151 * t299;
t82 = t108 * t145;
t46 = -t146 * t193 + t82;
t255 = t46 * t156;
t254 = rSges(5,3) + qJ(4);
t252 = t156 * t102 - t139;
t244 = t108 * t146;
t243 = t109 * t156;
t88 = t120 * t145;
t242 = t120 * t146;
t239 = t145 * t156;
t236 = t146 * t156;
t234 = t156 * t160;
t233 = t146 * t234 + t147 * t239;
t230 = t129 + t135;
t133 = pkin(3) * t239;
t229 = t133 - t138;
t148 = t161 * t278;
t228 = t163 * t277 + t148;
t134 = t146 * t264;
t223 = -t303 * t146 - t239 * t265;
t222 = t156 * t129 + t303 * t145;
t218 = -pkin(3) - t267;
t217 = -t227 / 0.2e1;
t215 = -t226 / 0.2e1;
t214 = t226 / 0.2e1;
t72 = Icges(6,5) * t237 - Icges(6,6) * t238 + Icges(6,3) * t145;
t212 = -t72 - t256;
t211 = -t147 - t265;
t210 = -qJ(4) * t236 - t138 + t229;
t202 = rSges(3,1) * t150 + rSges(3,2) * t152;
t200 = t264 - t267;
t40 = t149 * t75 + t151 * t73;
t198 = -t256 + t257;
t41 = t149 * t76 + t151 * t74;
t28 = t146 * t72 - t240 * t76 + t74 * t241;
t192 = -rSges(5,3) * t145 - t146 * t267;
t190 = t206 * t163;
t27 = -t240 * t75 + t275;
t186 = (t145 * t28 + t146 * t27) * qJD(5);
t29 = -t238 * t73 + t274;
t30 = t145 * t72 - t300;
t185 = (t145 * t30 + t146 * t29) * qJD(5);
t184 = t113 * t156;
t183 = t195 * t156;
t176 = -t243 * t145 - t289 * t146 + t156 * t197;
t175 = t289 * t145 - t243 * t146 + t156 * t198;
t174 = t109 * qJD(5) + t156 * t193;
t173 = -t190 + (t139 - t252) * t156;
t172 = t295 + t252;
t10 = t185 + t255;
t16 = -qJD(5) * t198 + t149 * (t286 * t145 - t146 * t184) + t151 * (t288 * t145 - t146 * t183);
t17 = -qJD(5) * t197 + t149 * (-t145 * t184 - t286 * t146) + t151 * (-t145 * t183 - t288 * t146);
t20 = t174 * t145 - t287 * t146;
t21 = t287 * t145 + t174 * t146;
t45 = t145 * t193 + t244;
t42 = t45 * t156;
t9 = t42 + t186;
t170 = (t42 + ((t275 + t30 + t300) * t146 + (-t29 + (t212 - t257) * t146 + t28 + t274) * t145) * qJD(5)) * t217 + (-t255 + ((t28 + (-t72 + t257) * t146 - t274) * t146 + (t145 * t212 - t27 + t275) * t145) * qJD(5) + t10) * t215 + (t16 + t21) * t214 + (t17 + t20 + t9) * t227 / 0.2e1 + (-qJD(5) * t193 + t100 * t149 + t151 * t99 + (t40 + t45) * t217 + (t41 + t46) * t214) * t156;
t167 = t290 * t145 + t146 * t79;
t53 = (rSges(6,2) * t241 + t259) * t156 + t223;
t54 = t156 * t191 + t222;
t166 = (-rSges(6,3) * t236 + t53) * t146 + (-t54 + (t146 * t199 - t79) * t156) * t145;
t123 = t145 * t234;
t101 = t199 * qJD(5);
t14 = t101 * t227 + (-t133 - t53 + (qJ(4) * t156 + t303) * t146 + t210 + t233) * t156 + t228;
t15 = -t101 * t226 + (t94 + t123 + t54 + (-t146 * t147 + t102) * t156) * t156 + t173;
t165 = t26 * (-t138 - t223 + t233) + (t14 * t211 + t15 * (rSges(6,3) - t160) + (-t26 * rSges(6,3) - t211 * t25) * t156) * t146 + (-t26 * t156 * t263 + t15 * (-t147 - t199) + (t156 * t25 - t14) * rSges(6,3)) * t145 - t25 * (t123 + t139 + t222);
t114 = t239 * t267;
t115 = t156 * t134;
t31 = (t114 + (-t145 * t264 - t261) * t156 + t210) * t156 + t228;
t32 = (t156 * t192 + t115) * t156 + t173;
t77 = t156 * (t145 * t200 + t261);
t43 = -t77 + t171;
t81 = -t134 - t192;
t44 = (-t102 - t81) * t156 + t177;
t164 = t44 * (t114 + t229) + (-t31 * t254 + t32 * (-pkin(3) + t200)) * t145 + (t31 * t218 + t32 * t254) * t146 + ((t254 * t43 - t264 * t44) * t145 + (-t218 * t43 - t254 * t44) * t146) * t156 - t43 * (t115 + t139);
t132 = rSges(4,2) * t239;
t91 = -rSges(4,1) * t236 + t132;
t78 = t156 * t81;
t69 = t189 + t90;
t61 = t156 * t91 - t190;
t60 = t156 * t90 + t228;
t39 = qJD(5) * t167 + qJD(2);
t11 = t166 * qJD(5);
t1 = [t170 + m(3) * ((t163 * t202 + t148) * t298 + (t162 * t278 + (-0.2e1 * t144 - t205 + t269 + t298) * t163) * (t202 + t280)) + (t14 * (-t206 + t230) - t15 * t207 + (t206 * t25 + t207 * t26) * qJD(1) + t165 - (t26 + t172 - t294) * t25) * m(6) + (t31 * (t134 - t206) - t32 * t207 + (t206 * t43 + t207 * t44) * qJD(1) + t164 - (t44 + t78 + t172) * t43) * m(5) + (t60 * (-t103 - t206) + t61 * (-t201 - t207) + t301 * t156 + t70 * t189 + (t268 * t156 - t132 - t235 - t70) * t69) * m(4); m(6) * t11; t170 + (t14 * t230 + t165 - t26 * (t253 + t293) + t25 * (-t252 + t294)) * m(6) + (t31 * t134 + t164 - t44 * (-t77 + t253) + t43 * (-t78 - t252)) * m(5) + (-(t69 * t103 + t301) * t156 - t60 * t103 - t61 * t201 - t69 * t91 + t70 * t90) * m(4); m(5) * (t32 * t145 + t31 * t146) + m(6) * (t14 * t146 + t15 * t145); t156 * ((t156 * t41 + t16) * t146 + (-t156 * t40 + t17) * t145) / 0.2e1 - t156 * ((-t232 * t149 + t231 * t151) * t156 + ((t145 * t270 + t146 * t271) * t151 + (-t272 * t145 - t273 * t146) * t149) * qJD(5)) / 0.2e1 + ((t82 * t226 + t243) * t146 + (t304 + (t284 * t145 + (-t283 - t244) * t146) * qJD(5)) * t145) * t215 + ((-t227 * t244 + t243) * t145 + (-t304 + (t283 * t146 + (-t284 + t82) * t145) * qJD(5)) * t146) * t217 + (t156 * t20 + ((t175 * t145 + t156 * t30) * t146 + (t176 * t145 - t156 * t29) * t145) * t302) * t145 / 0.2e1 + (t156 * t21 + ((t175 * t146 + t156 * t28) * t146 + (t176 * t146 - t156 * t27) * t145) * t302) * t146 / 0.2e1 - (t9 + t186) * t239 / 0.2e1 + (t10 + t185) * t236 / 0.2e1 + (t11 * t167 + t14 * t88 - t15 * t242 + t39 * t166 - (t242 * t26 - t25 * t88) * t156 - (t39 * (-t145 * t88 - t146 * t242) + t296 * t199) * qJD(5) + (t26 * t236 - t25 * t239) * t120 + t296 * t101) * m(6);];
tauc = t1(:);
