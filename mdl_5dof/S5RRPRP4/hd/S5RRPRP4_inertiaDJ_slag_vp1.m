% Calculate time derivative of joint inertia matrix for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:44
% DurationCPUTime: 4.87s
% Computational Cost: add. (4493->306), mult. (4850->427), div. (0->0), fcn. (3536->6), ass. (0->180)
t163 = sin(qJ(4));
t165 = cos(qJ(4));
t259 = Icges(5,4) * t165;
t126 = -Icges(5,2) * t163 + t259;
t256 = Icges(6,5) * t163;
t260 = Icges(5,4) * t163;
t307 = t256 - t260 + (Icges(5,1) + Icges(6,1)) * t165;
t309 = -t165 * t126 - t307 * t163;
t191 = -Icges(6,3) * t165 + t256;
t194 = Icges(5,2) * t165 + t260;
t308 = (t191 - t194) * t163;
t161 = qJD(1) + qJD(2);
t305 = t161 * t165;
t255 = Icges(6,5) * t165;
t123 = Icges(6,3) * t163 + t255;
t298 = -t123 * t165 - t309;
t192 = Icges(5,5) * t163 + Icges(5,6) * t165;
t193 = Icges(6,4) * t163 - Icges(6,6) * t165;
t195 = Icges(6,1) * t163 - t255;
t304 = t195 * t305 + (-t192 - t193) * qJD(4) + (t298 + t308) * t161;
t303 = rSges(6,1) + pkin(4);
t301 = rSges(6,3) + qJ(5);
t292 = t163 * t303;
t302 = rSges(4,2) - pkin(2);
t124 = Icges(5,5) * t165 - Icges(5,6) * t163;
t125 = Icges(6,4) * t165 + Icges(6,6) * t163;
t300 = t124 + t125;
t162 = qJ(1) + qJ(2);
t158 = sin(t162);
t236 = qJD(4) * t165;
t159 = cos(t162);
t246 = t159 * t163;
t296 = -t158 * t236 - t161 * t246;
t212 = t301 * t165;
t75 = Icges(5,6) * t159 + t158 * t194;
t196 = Icges(5,1) * t163 + t259;
t183 = t196 * t158;
t79 = Icges(5,5) * t159 + t183;
t200 = t163 * t79 + t165 * t75;
t294 = t159 * t200;
t69 = Icges(6,6) * t159 + t158 * t191;
t77 = Icges(6,4) * t159 + t158 * t195;
t204 = t163 * t77 - t165 * t69;
t293 = t159 * t204;
t290 = t303 * t246;
t152 = t159 * rSges(6,2);
t249 = t158 * t163;
t289 = -t303 * t249 - t152;
t287 = -Icges(5,3) * t161 + qJD(4) * t124;
t286 = -Icges(6,2) * t161 + qJD(4) * t125;
t242 = t301 * t163 + t303 * t165;
t210 = t242 * t161;
t261 = qJD(5) * t163 + (t212 - t292) * qJD(4);
t27 = t158 * t210 - t159 * t261;
t28 = t158 * t261 + t159 * t210;
t63 = t242 * t158;
t64 = t242 * t159;
t280 = t158 * (-t161 * t64 + t28) + t159 * (t161 * t63 - t27);
t237 = qJD(4) * t163;
t224 = t158 * t237;
t235 = qJD(5) * t165;
t171 = t158 * t235 - t301 * t224 + t303 * t296;
t273 = -pkin(2) - pkin(7);
t234 = -rSges(6,2) + t273;
t173 = t158 * t234 - t159 * t212;
t247 = t159 * t161;
t243 = qJ(3) * t247 + qJD(3) * t158;
t12 = t161 * t173 - t171 + t243;
t145 = qJD(3) * t159;
t221 = t159 * t236;
t222 = t159 * t237;
t248 = t158 * t165;
t172 = t159 * t235 - t303 * t221 + t301 * (-t161 * t248 - t222);
t13 = t145 + (t234 * t159 + (-qJ(3) - t292) * t158) * t161 - t172;
t147 = t159 * qJ(3);
t35 = t147 + t173 + t290;
t155 = t159 * pkin(2);
t239 = t158 * qJ(3) + t155;
t225 = t159 * pkin(7) + t239;
t36 = -t158 * t212 + t225 - t289;
t279 = (t161 * t36 + t13) * t158 + t159 * (t161 * t35 - t12);
t164 = sin(qJ(1));
t264 = pkin(1) * qJD(1);
t231 = t164 * t264;
t10 = t12 - t231;
t166 = cos(qJ(1));
t230 = t166 * t264;
t11 = t13 - t230;
t268 = pkin(1) * t164;
t33 = t35 - t268;
t160 = t166 * pkin(1);
t34 = t160 + t36;
t278 = t158 * (t161 * t34 + t11) + t159 * (t161 * t33 - t10);
t277 = 2 * m(3);
t276 = 2 * m(4);
t275 = 2 * m(5);
t274 = 2 * m(6);
t209 = rSges(5,1) * t163 + rSges(5,2) * t165;
t110 = t209 * qJD(4);
t269 = m(5) * t110;
t267 = t301 * t248 + t289;
t245 = t159 * t165;
t265 = rSges(6,2) * t158;
t266 = t301 * t245 + t265 - t290;
t263 = t158 * rSges(5,3);
t250 = t158 * t161;
t141 = rSges(5,2) * t245;
t240 = rSges(5,1) * t246 + t141;
t238 = t158 ^ 2 + t159 ^ 2;
t233 = -rSges(5,3) + t273;
t232 = m(6) * t237;
t229 = rSges(5,2) * t237;
t226 = t296 * rSges(5,1) - t161 * t141;
t82 = rSges(5,1) * t249 + rSges(5,2) * t248 + t159 * rSges(5,3);
t97 = t159 * rSges(3,1) - rSges(3,2) * t158;
t211 = t110 * t238;
t86 = -rSges(3,1) * t247 + rSges(3,2) * t250;
t96 = -rSges(3,1) * t158 - rSges(3,2) * t159;
t203 = t163 * t69 + t165 * t77;
t70 = Icges(6,6) * t158 - t159 * t191;
t78 = Icges(6,4) * t158 - t159 * t195;
t202 = t163 * t78 - t165 * t70;
t201 = t163 * t70 + t165 * t78;
t199 = -t163 * t75 + t165 * t79;
t76 = Icges(5,6) * t158 - t159 * t194;
t80 = Icges(5,5) * t158 - t159 * t196;
t198 = t163 * t80 + t165 * t76;
t197 = t163 * t76 - t165 * t80;
t65 = t159 * rSges(4,3) + t302 * t158 + t147;
t66 = -rSges(4,2) * t159 + t158 * rSges(4,3) + t239;
t60 = t225 + t82;
t186 = t123 * t236 + ((-t195 - t196) * t165 - t308) * qJD(4);
t85 = t96 * t161;
t185 = t202 * t158;
t184 = t198 * t158;
t180 = t193 * t161;
t179 = t192 * t158;
t57 = rSges(4,3) * t247 + t302 * t250 + t243;
t59 = t158 * t233 + t147 + t240;
t58 = rSges(4,2) * t247 + t145 + (-t155 + (-rSges(4,3) - qJ(3)) * t158) * t161;
t25 = (t161 * t233 - t229) * t158 - t226 + t243;
t23 = t25 - t231;
t119 = rSges(5,1) * t221;
t26 = -rSges(5,2) * t222 + t119 + t145 + (t233 * t159 + (-qJ(3) - t209) * t158) * t161;
t24 = t26 - t230;
t55 = t59 - t268;
t56 = t160 + t60;
t170 = m(5) * ((t161 * t55 - t23) * t159 + (t161 * t56 + t24) * t158);
t169 = m(5) * ((t161 * t59 - t25) * t159 + (t161 * t60 + t26) * t158);
t168 = -t197 * t247 / 0.2e1 + (t300 * t158 - t298 * t159 + t201) * t247 / 0.2e1 + (t183 * t305 + (-t198 - t202) * qJD(4) + t304 * t158) * t158 / 0.2e1 + (t196 * t247 * t165 + (-t200 - t204) * qJD(4) + t304 * t159) * t159 / 0.2e1 - (t298 * t158 + t300 * t159 + t199 + t203) * t250 / 0.2e1;
t167 = t309 * qJD(4) + t186;
t135 = rSges(5,1) * t165 - rSges(5,2) * t163;
t90 = t160 + t97;
t89 = t96 - t268;
t84 = -t240 + t263;
t74 = Icges(6,2) * t158 - t159 * t193;
t73 = Icges(6,2) * t159 + t158 * t193;
t72 = Icges(5,3) * t158 - t159 * t192;
t71 = Icges(5,3) * t159 + t179;
t68 = t86 - t230;
t67 = t85 - t231;
t62 = t160 + t66;
t61 = t65 - t268;
t54 = t58 - t230;
t53 = t57 - t231;
t46 = t158 * t286 + t159 * t180;
t45 = t158 * t180 - t159 * t286;
t44 = t158 * t287 + t192 * t247;
t43 = -t159 * t287 + t161 * t179;
t22 = t158 * t72 - t198 * t159;
t21 = t158 * t71 - t294;
t20 = t158 * t74 - t202 * t159;
t19 = t158 * t73 - t293;
t18 = t159 * t72 + t184;
t17 = t200 * t158 + t159 * t71;
t16 = t159 * t74 + t185;
t15 = t204 * t158 + t159 * t73;
t14 = t158 * t267 + t159 * t266;
t1 = t172 * t159 + t171 * t158 + ((t152 + t267) * t159 + (t265 + (t212 + t292) * t159 - t266) * t158) * t161;
t2 = [-t126 * t236 + (t10 * t34 + t11 * t33) * t274 + (t23 * t56 + t24 * t55) * t275 + (t53 * t62 + t54 * t61) * t276 + (t67 * t90 + t68 * t89) * t277 + t186 - t307 * t237; m(6) * (t10 * t36 + t11 * t35 + t12 * t34 + t13 * t33) + m(5) * (t23 * t60 + t24 * t59 + t25 * t56 + t26 * t55) + m(4) * (t53 * t66 + t54 * t65 + t57 * t62 + t58 * t61) + m(3) * (t67 * t97 + t68 * t96 + t85 * t90 + t86 * t89) + t167; (t12 * t36 + t13 * t35) * t274 + (t25 * t60 + t26 * t59) * t275 + (t57 * t66 + t58 * t65) * t276 + (t85 * t97 + t86 * t96) * t277 + t167; m(6) * t278 + t170 + m(4) * ((t161 * t61 - t53) * t159 + (t161 * t62 + t54) * t158); m(6) * t279 + t169 + m(4) * ((t161 * t65 - t57) * t159 + (t161 * t66 + t58) * t158); 0; t168 + m(6) * (-t10 * t64 + t11 * t63 + t27 * t34 + t28 * t33) + t135 * t170 - (t158 * t55 - t159 * t56) * t269; t168 + m(6) * (-t12 * t64 + t13 * t63 + t27 * t36 + t28 * t35) + t135 * t169 - (t158 * t59 - t159 * t60) * t269; -m(5) * t211 + m(6) * t280; ((-t158 * t82 + t159 * t84) * ((-t161 * t82 - t119 + (rSges(5,3) * t161 + t229) * t159) * t159 + (rSges(5,2) * t224 + (t159 * t209 + t263 - t84) * t161 + t226) * t158) - t135 * t211) * t275 + t159 * ((t159 * t44 + (t18 + t294) * t161) * t159 + (-t17 * t161 + (t236 * t80 - t237 * t76) * t158 + (t199 * qJD(4) + t161 * t198 + t43) * t159) * t158) + t158 * ((t158 * t43 + (-t21 + t184) * t161) * t158 + (t22 * t161 + (-t236 * t79 + t237 * t75) * t159 + (t197 * qJD(4) + t200 * t161 + t44) * t158) * t159) + (t1 * t14 - t27 * t64 + t28 * t63) * t274 + t159 * ((t159 * t46 + (t16 + t293) * t161) * t159 + (-t15 * t161 + (t236 * t78 + t237 * t70) * t158 + (t203 * qJD(4) + t161 * t202 + t45) * t159) * t158) + t158 * ((t158 * t45 + (-t19 + t185) * t161) * t158 + (t20 * t161 + (-t236 * t77 - t237 * t69) * t159 + (-t201 * qJD(4) + t204 * t161 + t46) * t158) * t159) + ((-t15 - t17) * t159 + (-t16 - t18) * t158) * t250 + ((t19 + t21) * t159 + (t20 + t22) * t158) * t247; m(6) * ((t158 * t33 - t159 * t34) * t237 - t278 * t165); m(6) * ((t158 * t35 - t159 * t36) * t237 - t279 * t165); t238 * t232; m(6) * ((t1 + (t158 * t63 + t159 * t64) * qJD(4)) * t163 + (qJD(4) * t14 - t280) * t165); 0.2e1 * (0.1e1 - t238) * t165 * t232;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
