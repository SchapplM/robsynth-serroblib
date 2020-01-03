% Calculate time derivative of joint inertia matrix for
% S5RPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR11_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR11_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR11_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:11
% EndTime: 2019-12-31 18:27:24
% DurationCPUTime: 7.63s
% Computational Cost: add. (8522->485), mult. (12081->704), div. (0->0), fcn. (11256->8), ass. (0->244)
t173 = pkin(8) + qJ(3);
t167 = sin(t173);
t168 = cos(t173);
t278 = Icges(5,5) * t168;
t282 = Icges(4,4) * t168;
t321 = -t278 + t282 + (Icges(4,1) + Icges(5,1)) * t167;
t320 = t321 * qJD(3);
t179 = sin(qJ(5));
t181 = cos(qJ(5));
t137 = t167 * t181 - t168 * t179;
t308 = qJD(3) - qJD(5);
t313 = t308 * t137;
t319 = -t313 / 0.2e1;
t209 = t167 * t179 + t168 * t181;
t92 = t308 * t209;
t318 = -t92 / 0.2e1;
t180 = sin(qJ(1));
t317 = t180 / 0.2e1;
t182 = cos(qJ(1));
t316 = -t182 / 0.2e1;
t315 = -qJD(1) / 0.2e1;
t314 = qJD(1) / 0.2e1;
t312 = Icges(6,5) * t318 + Icges(6,6) * t319;
t219 = -Icges(4,2) * t167 + t282;
t108 = Icges(4,6) * t180 + t182 * t219;
t283 = Icges(4,4) * t167;
t223 = Icges(4,1) * t168 - t283;
t112 = Icges(4,5) * t180 + t182 * t223;
t210 = t108 * t167 - t112 * t168;
t200 = t210 * t180;
t107 = -Icges(4,6) * t182 + t180 * t219;
t111 = -Icges(4,5) * t182 + t180 * t223;
t211 = t107 * t167 - t111 * t168;
t201 = t211 * t182;
t215 = Icges(5,3) * t167 + t278;
t102 = Icges(5,6) * t180 + t182 * t215;
t279 = Icges(5,5) * t167;
t221 = Icges(5,1) * t168 + t279;
t110 = Icges(5,4) * t180 + t182 * t221;
t212 = t102 * t167 + t110 * t168;
t202 = t212 * t180;
t101 = -Icges(5,6) * t182 + t180 * t215;
t109 = -Icges(5,4) * t182 + t180 * t221;
t213 = t101 * t167 + t109 * t168;
t203 = t213 * t182;
t172 = t180 * rSges(5,2);
t269 = t168 * t182;
t270 = t167 * t182;
t119 = rSges(5,1) * t269 + rSges(5,3) * t270 + t172;
t162 = pkin(3) * t269;
t128 = qJ(4) * t270 + t162;
t311 = -t119 - t128;
t171 = t180 * rSges(4,3);
t310 = -rSges(4,2) * t270 + t171;
t178 = -pkin(6) - qJ(2);
t177 = cos(pkin(8));
t164 = pkin(2) * t177 + pkin(1);
t293 = rSges(4,1) * t168;
t233 = -rSges(4,2) * t167 + t293;
t206 = -t164 - t233;
t89 = (rSges(4,3) - t178) * t182 + t206 * t180;
t120 = rSges(4,1) * t269 + t310;
t157 = t182 * t164;
t237 = -t180 * t178 + t157;
t90 = t120 + t237;
t309 = t180 * t90 + t182 * t89;
t216 = Icges(4,5) * t168 - Icges(4,6) * t167;
t103 = -Icges(4,3) * t182 + t180 * t216;
t217 = Icges(5,4) * t168 + Icges(5,6) * t167;
t105 = -Icges(5,2) * t182 + t180 * t217;
t207 = rSges(3,1) * t177 - rSges(3,2) * sin(pkin(8)) + pkin(1);
t287 = rSges(3,3) + qJ(2);
t100 = t180 * t287 + t182 * t207;
t307 = 2 * m(4);
t306 = 2 * m(5);
t305 = 2 * m(6);
t174 = t180 ^ 2;
t175 = t182 ^ 2;
t304 = m(5) / 0.2e1;
t303 = m(6) / 0.2e1;
t302 = -pkin(3) - pkin(4);
t301 = -rSges(5,1) - pkin(3);
t300 = rSges(6,3) + pkin(7);
t150 = rSges(4,1) * t167 + rSges(4,2) * t168;
t299 = m(4) * t150;
t298 = pkin(4) * t167;
t297 = pkin(7) * t180;
t296 = pkin(7) * t182;
t257 = qJD(1) * t180;
t54 = -t137 * t257 + t182 * t92;
t55 = -t182 * t313 - t209 * t257;
t294 = t55 * rSges(6,1) + t54 * rSges(6,2);
t292 = rSges(5,1) * t167;
t125 = t137 * t182;
t126 = t209 * t182;
t67 = Icges(6,5) * t126 + Icges(6,6) * t125 - Icges(6,3) * t180;
t291 = t180 * t67;
t289 = t182 * t67;
t286 = -rSges(5,3) - qJ(4);
t123 = t137 * t180;
t124 = t209 * t180;
t231 = -t124 * rSges(6,1) - t123 * rSges(6,2);
t72 = rSges(6,3) * t182 - t231;
t285 = t180 * t168 * pkin(4) + t296 + t72;
t161 = pkin(4) * t269;
t267 = t126 * rSges(6,1) + t125 * rSges(6,2);
t73 = -rSges(6,3) * t180 + t267;
t284 = t161 - t297 + t73;
t273 = qJ(4) * t167;
t272 = qJ(4) * t168;
t96 = rSges(6,1) * t137 - rSges(6,2) * t209;
t271 = qJD(1) * t96;
t268 = t178 * t182;
t230 = pkin(3) * t168 + t273;
t116 = qJD(3) * t230 - qJD(4) * t168;
t232 = rSges(5,1) * t168 + rSges(5,3) * t167;
t266 = -t232 * qJD(3) - t116;
t127 = t230 * t180;
t265 = t180 * t127 + t182 * t128;
t148 = pkin(3) * t167 - t272;
t149 = -rSges(5,3) * t168 + t292;
t264 = -t148 - t149;
t252 = qJD(3) * t182;
t240 = t168 * t252;
t251 = qJD(4) * t167;
t263 = qJ(4) * t240 + t182 * t251;
t256 = qJD(1) * t182;
t262 = rSges(5,2) * t256 + rSges(5,3) * t240;
t170 = qJD(2) * t182;
t261 = t178 * t257 + t170;
t260 = t174 + t175;
t104 = Icges(4,3) * t180 + t182 * t216;
t259 = qJD(1) * t104;
t106 = Icges(5,2) * t180 + t182 * t217;
t258 = qJD(1) * t106;
t255 = qJD(3) * t167;
t254 = qJD(3) * t168;
t253 = qJD(3) * t180;
t26 = Icges(6,4) * t55 + Icges(6,2) * t54 - Icges(6,6) * t256;
t27 = Icges(6,1) * t55 + Icges(6,4) * t54 - Icges(6,5) * t256;
t49 = Icges(6,4) * t92 + Icges(6,2) * t313;
t50 = Icges(6,1) * t92 + Icges(6,4) * t313;
t69 = Icges(6,4) * t126 + Icges(6,2) * t125 - Icges(6,6) * t180;
t71 = Icges(6,1) * t126 + Icges(6,4) * t125 - Icges(6,5) * t180;
t93 = Icges(6,5) * t137 - Icges(6,6) * t209;
t94 = Icges(6,4) * t137 - Icges(6,2) * t209;
t95 = Icges(6,1) * t137 - Icges(6,4) * t209;
t250 = t125 * t49 / 0.2e1 + t126 * t50 / 0.2e1 + t180 * t312 - t256 * t93 / 0.2e1 + t54 * t94 / 0.2e1 + t55 * t95 / 0.2e1 - t209 * t26 / 0.2e1 + t137 * t27 / 0.2e1 + t69 * t313 / 0.2e1 + t71 * t92 / 0.2e1;
t56 = qJD(1) * t125 + t180 * t92;
t57 = qJD(1) * t126 - t180 * t313;
t189 = Icges(6,4) * t57 + Icges(6,2) * t56 - Icges(6,6) * t257;
t190 = Icges(6,1) * t57 + Icges(6,4) * t56 - Icges(6,5) * t257;
t68 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t182;
t70 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t182;
t249 = -t123 * t49 / 0.2e1 - t124 * t50 / 0.2e1 + t182 * t312 + t257 * t93 / 0.2e1 - t56 * t94 / 0.2e1 - t57 * t95 / 0.2e1 + t209 * t189 / 0.2e1 - t137 * t190 / 0.2e1 + t68 * t319 + t70 * t318;
t248 = -t178 - t300;
t242 = t167 * t253;
t154 = pkin(3) * t242;
t241 = t167 * t252;
t243 = t167 * t257;
t247 = t127 * t256 + t180 * (qJD(1) * t162 + t180 * t251 - t154 + (t167 * t256 + t168 * t253) * qJ(4)) + t182 * (-qJ(4) * t243 + (-t168 * t257 - t241) * pkin(3) + t263);
t169 = qJD(2) * t180;
t245 = t169 + t263;
t244 = t154 + t261;
t239 = t96 + t298;
t51 = rSges(6,1) * t92 + rSges(6,2) * t313;
t238 = t260 * t51;
t98 = t264 * t182;
t236 = -t148 - t239;
t235 = -t57 * rSges(6,1) - t56 * rSges(6,2);
t188 = Icges(6,5) * t57 + Icges(6,6) * t56 - Icges(6,3) * t257;
t19 = t125 * t69 + t126 * t71 - t291;
t66 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t182;
t191 = t125 * t68 + t126 * t70 - t180 * t66;
t25 = Icges(6,5) * t55 + Icges(6,6) * t54 - Icges(6,3) * t256;
t1 = (t125 * t26 + t126 * t27 - t180 * t25 + t54 * t69 + t55 * t71) * t180 - (t125 * t189 + t126 * t190 - t180 * t188 - t256 * t66 + t54 * t68 + t55 * t70) * t182 + (t19 * t182 + (t191 - t289) * t180) * qJD(1);
t18 = t123 * t69 + t124 * t71 + t289;
t192 = t123 * t68 + t124 * t70 + t182 * t66;
t2 = (t123 * t26 + t124 * t27 + t182 * t25 + t56 * t69 + t57 * t71) * t180 - (t123 * t189 + t124 * t190 + t182 * t188 - t257 * t66 + t56 * t68 + t57 * t70) * t182 + (t18 * t182 + (t192 - t291) * t180) * qJD(1);
t234 = t180 * t1 - t182 * t2;
t193 = t168 * t302 - t164 - t273;
t184 = t180 * t193 + t182 * t248;
t35 = t184 + t231;
t36 = t180 * t248 + t128 + t157 + t161 + t267;
t225 = t180 * t36 + t182 * t35;
t187 = t167 * t286 + t168 * t301 - t164;
t185 = t187 * t180;
t62 = (rSges(5,2) - t178) * t182 + t185;
t63 = t237 - t311;
t224 = t180 * t63 + t182 * t62;
t218 = Icges(4,2) * t168 + t283;
t208 = -pkin(4) * t254 - t116 - t51;
t61 = t236 * t182;
t204 = qJD(3) * t150;
t197 = qJD(3) * t218;
t196 = qJD(3) * (-Icges(5,4) * t167 + Icges(5,6) * t168);
t195 = qJD(3) * (-Icges(4,5) * t167 - Icges(4,6) * t168);
t186 = rSges(4,2) * t243 + rSges(4,3) * t256 - t182 * t204;
t99 = -t180 * t207 + t182 * t287;
t139 = t233 * qJD(3);
t129 = t148 * t257;
t118 = -rSges(4,3) * t182 + t180 * t233;
t117 = -rSges(5,2) * t182 + t180 * t232;
t97 = t264 * t180;
t88 = -t100 * qJD(1) + t170;
t87 = qJD(1) * t99 + t169;
t79 = t180 * t196 + t258;
t78 = -qJD(1) * t105 + t182 * t196;
t77 = t180 * t195 + t259;
t76 = -qJD(1) * t103 + t182 * t195;
t60 = t236 * t180;
t59 = t150 * t253 + (t182 * t206 - t171) * qJD(1) + t261;
t58 = t169 + (-t268 + (-t164 - t293) * t180) * qJD(1) + t186;
t47 = qJD(1) * t98 + t180 * t266;
t46 = t149 * t257 + t182 * t266 + t129;
t45 = t180 * t104 - t210 * t182;
t44 = t180 * t103 - t201;
t43 = t180 * t106 + t212 * t182;
t42 = t180 * t105 + t203;
t41 = -t104 * t182 - t200;
t40 = -t103 * t182 - t180 * t211;
t39 = -t106 * t182 + t202;
t38 = -t105 * t182 + t180 * t213;
t37 = t180 * t117 + t119 * t182 + t265;
t34 = -t180 * t72 - t182 * t73;
t33 = (-t251 + (t168 * t286 + t292) * qJD(3)) * t180 + (t182 * t187 - t172) * qJD(1) + t244;
t32 = t301 * t241 + (t185 - t268) * qJD(1) + t245 + t262;
t31 = t137 * t71 - t209 * t69;
t30 = t137 * t70 - t209 * t68;
t29 = -rSges(6,3) * t257 - t235;
t28 = -rSges(6,3) * t256 + t294;
t24 = t125 * t94 + t126 * t95 - t180 * t93;
t23 = t123 * t94 + t124 * t95 + t182 * t93;
t22 = t180 * t285 + t182 * t284 + t265;
t21 = qJD(1) * t61 + t180 * t208;
t20 = t182 * t208 + t239 * t257 + t129;
t15 = (-t251 + (-t272 + t298) * qJD(3)) * t180 + (t180 * t300 + t182 * t193) * qJD(1) + t235 + t244;
t14 = qJD(1) * t184 + t241 * t302 + t245 + t294;
t13 = t182 * t262 + (-t149 * t174 - t175 * t292) * qJD(3) + (t182 * t117 + (t172 + t311) * t180) * qJD(1) + t247;
t12 = t180 * t19 - t182 * t191;
t11 = t18 * t180 - t182 * t192;
t10 = -t180 * t29 - t182 * t28 + (t180 * t73 - t182 * t72) * qJD(1);
t5 = (-pkin(4) * t241 + t28 + (t285 - t296) * qJD(1)) * t182 + (-pkin(4) * t242 + t29 + (-t128 - t284 - t297) * qJD(1)) * t180 + t247;
t3 = [t92 * t95 + t137 * t50 + t313 * t94 - t209 * t49 + (t14 * t36 + t15 * t35) * t305 + (t32 * t63 + t33 * t62) * t306 + (t58 * t90 + t59 * t89) * t307 + 0.2e1 * m(3) * (t100 * t87 + t88 * t99) + (-Icges(5,3) * t168 - t218 + t221 + t223 + t279) * t255 + (t219 - t215 + t321) * t254; m(6) * (qJD(1) * t225 - t14 * t182 + t180 * t15) + m(5) * (qJD(1) * t224 + t180 * t33 - t182 * t32) + m(4) * (t309 * qJD(1) + t180 * t59 - t182 * t58) + m(3) * (t180 * t88 - t182 * t87 + (t100 * t180 + t182 * t99) * qJD(1)); 0; ((t102 * t314 + t108 * t315 + t197 * t317) * t168 + (t320 * t317 + (t110 + t112) * t315) * t167 + t249) * t182 + ((t101 * t314 + t107 * t315 + t197 * t316) * t168 + (t320 * t316 + (t109 + t111) * t315) * t167 + t250) * t180 + m(6) * (t14 * t60 + t15 * t61 + t20 * t35 + t21 * t36) + m(5) * (t32 * t97 + t33 * t98 + t46 * t62 + t47 * t63) + m(4) * ((-t180 * t58 - t182 * t59) * t150 - t309 * t139) + ((t24 / 0.2e1 + t31 / 0.2e1 - t90 * t299 + (-t102 / 0.2e1 + t108 / 0.2e1) * t168 + (t110 / 0.2e1 + t112 / 0.2e1) * t167) * t182 + (t23 / 0.2e1 + t30 / 0.2e1 + t89 * t299 + (-t101 / 0.2e1 + t107 / 0.2e1) * t168 + (t109 / 0.2e1 + t111 / 0.2e1) * t167) * t180) * qJD(1) + (-t203 / 0.2e1 + t201 / 0.2e1 + t202 / 0.2e1 - t200 / 0.2e1 + (t216 + t217) * (t174 / 0.2e1 + t175 / 0.2e1)) * qJD(3); m(5) * (t46 * t180 - t182 * t47 + (t180 * t97 + t182 * t98) * qJD(1)) + m(6) * (t20 * t180 - t182 * t21 + (t180 * t60 + t182 * t61) * qJD(1)); (t20 * t61 + t21 * t60 + t22 * t5) * t305 + (t37 * t13 + t46 * t98 + t47 * t97) * t306 - t182 * ((t182 * t77 + (t41 + t201) * qJD(1)) * t182 + (t40 * qJD(1) + (-t108 * t254 - t112 * t255 + t259) * t180 + (-t76 + (t107 * t168 + t111 * t167) * qJD(3) - t210 * qJD(1)) * t182) * t180) - t182 * ((t182 * t79 + (t39 - t203) * qJD(1)) * t182 + (t38 * qJD(1) + (t102 * t254 - t110 * t255 + t258) * t180 + (-t78 + (-t101 * t168 + t109 * t167) * qJD(3) + t212 * qJD(1)) * t182) * t180) + t180 * ((t180 * t76 + (t44 + t200) * qJD(1)) * t180 + (t45 * qJD(1) + (t107 * t254 + t111 * t255) * t182 + (-t77 + (-t108 * t168 - t112 * t167) * qJD(3) + (t104 - t211) * qJD(1)) * t180) * t182) + t180 * ((t180 * t78 + (t42 - t202) * qJD(1)) * t180 + (t43 * qJD(1) + (-t101 * t254 + t109 * t255) * t182 + (-t79 + (t102 * t168 - t110 * t167) * qJD(3) + (t106 + t213) * qJD(1)) * t180) * t182) + ((t180 * t118 + t120 * t182) * ((qJD(1) * t118 + t186) * t182 + (-t180 * t204 + (-t120 + t310) * qJD(1)) * t180) + t260 * t150 * t139) * t307 + t234 + (t11 + (-t38 - t40) * t182 + (t39 + t41) * t180) * t257 + (t12 + (-t42 - t44) * t182 + (t43 + t45) * t180) * t256; 0.2e1 * (t224 * t304 + t225 * t303) * t254 + 0.2e1 * ((t14 * t180 + t15 * t182 + t256 * t36 - t257 * t35) * t303 + (t180 * t32 + t182 * t33 + t256 * t63 - t257 * t62) * t304) * t167; 0; 0.2e1 * ((t252 * t61 + t253 * t60 - t5) * t303 + (t252 * t98 + t253 * t97 - t13) * t304) * t168 + 0.2e1 * ((qJD(3) * t22 + t180 * t21 + t182 * t20 + t256 * t60 - t257 * t61) * t303 + (qJD(3) * t37 + t180 * t47 + t182 * t46 + t256 * t97 - t257 * t98) * t304) * t167; 0.4e1 * (t304 + t303) * (-0.1e1 + t260) * t167 * t254; (m(6) * (t15 * t96 + t271 * t36 + t35 * t51) - t249) * t182 + (m(6) * (t14 * t96 - t271 * t35 + t36 * t51) - t250) * t180 + ((t24 + t31) * t182 + (t23 + t30) * t180) * t315; 0; m(6) * (t10 * t22 + t34 * t5) + (m(6) * (t20 * t96 + t271 * t60 + t51 * t61) - qJD(1) * t12 + t2) * t182 + (m(6) * (t21 * t96 - t271 * t61 + t51 * t60) - t1 - qJD(1) * t11) * t180; m(6) * (-t10 * t168 + t167 * t238 + (t168 * t260 * t96 + t167 * t34) * qJD(3)); (t34 * t10 + t238 * t96) * t305 + (t180 * t11 + t182 * t12) * qJD(1) + t234;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
