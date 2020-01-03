% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPPR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:38
% EndTime: 2019-12-31 17:47:43
% DurationCPUTime: 2.91s
% Computational Cost: add. (6810->242), mult. (17311->370), div. (0->0), fcn. (20222->8), ass. (0->168)
t181 = sin(pkin(7));
t186 = cos(qJ(1));
t240 = sin(pkin(8));
t199 = t186 * t240;
t182 = cos(pkin(8));
t184 = sin(qJ(1));
t224 = t184 * t182;
t166 = t181 * t224 + t199;
t200 = t184 * t240;
t223 = t186 * t182;
t167 = -t181 * t200 + t223;
t177 = t186 * qJ(2);
t214 = t186 * pkin(3) + t177;
t183 = cos(pkin(7));
t202 = qJ(3) * t181 + pkin(1);
t247 = pkin(2) + qJ(4);
t277 = t247 * t183 + t202;
t188 = t167 * pkin(4) - t184 * t277 + t214;
t185 = cos(qJ(5));
t249 = sin(qJ(5));
t206 = t183 * t249;
t144 = t167 * t185 - t184 * t206;
t226 = t183 * t184;
t145 = t167 * t249 + t185 * t226;
t198 = t144 * rSges(6,1) - t145 * rSges(6,2);
t67 = (rSges(6,3) + pkin(6)) * t166 + t188 + t198;
t109 = -t166 * rSges(6,3) - t198;
t69 = t166 * pkin(6) - t109 + t188;
t246 = t67 - t69;
t294 = m(6) * t246;
t201 = t183 * t240;
t162 = t181 * t185 + t249 * t201;
t163 = -t181 * t249 + t185 * t201;
t227 = t182 * t183;
t129 = -rSges(6,1) * t163 + rSges(6,2) * t162 + rSges(6,3) * t227;
t291 = -t109 * t227 - t129 * t166;
t238 = Icges(6,4) * t144;
t103 = -Icges(6,2) * t145 + Icges(6,6) * t166 + t238;
t138 = Icges(6,4) * t145;
t106 = Icges(6,1) * t144 + Icges(6,5) * t166 - t138;
t229 = t181 * t184;
t165 = t181 * t199 + t224;
t225 = t183 * t186;
t141 = -t165 * t249 + t185 * t225;
t142 = t165 * t185 + t186 * t206;
t122 = rSges(6,1) * t141 - rSges(6,2) * t142;
t123 = rSges(6,1) * t145 + rSges(6,2) * t144;
t237 = Icges(6,4) * t163;
t127 = Icges(6,2) * t162 + Icges(6,6) * t227 - t237;
t132 = Icges(6,5) * t162 + Icges(6,6) * t163;
t134 = Icges(6,1) * t162 + t237;
t205 = t227 / 0.2e1;
t164 = -t181 * t223 + t200;
t191 = t142 * rSges(6,1) + t141 * rSges(6,2) + t164 * rSges(6,3);
t203 = (pkin(3) + qJ(2)) * t184;
t68 = t165 * pkin(4) + t164 * pkin(6) + t186 * t277 + t191 + t203;
t288 = -(t134 / 0.2e1 - t127 / 0.2e1) * t163 + m(6) * (t122 * t68 - t123 * t67) + t132 * t205;
t287 = -t246 * t184 - t186 * t68;
t286 = m(5) + m(6);
t284 = t164 / 0.2e1;
t283 = -t166 / 0.2e1;
t282 = m(6) * (t122 * t184 - t123 * t186);
t159 = Icges(6,4) * t162;
t128 = -Icges(6,1) * t163 + Icges(6,5) * t227 + t159;
t133 = Icges(6,2) * t163 + t159;
t215 = t128 + t133;
t281 = t215 * t162;
t279 = -(rSges(4,3) + qJ(3)) * t181 - pkin(1);
t278 = (rSges(5,3) + t247) * t183 + t202;
t213 = -t184 ^ 2 - t186 ^ 2;
t168 = t213 * t181;
t276 = 0.2e1 * t168;
t169 = t213 * t183;
t275 = 0.2e1 * t169;
t274 = 2 * qJD(1);
t273 = 4 * qJD(1);
t272 = m(5) / 0.2e1;
t270 = m(6) / 0.2e1;
t130 = t186 * rSges(4,1) + t177 + ((rSges(4,2) - pkin(2)) * t183 + t279) * t184;
t131 = -rSges(4,2) * t225 + (rSges(4,1) + qJ(2)) * t184 + (pkin(2) * t183 - t279) * t186;
t228 = t181 * t186;
t269 = m(4) * (-t130 * t229 + t131 * t228);
t268 = m(4) * (t130 * t186 + t131 * t184);
t187 = t167 * rSges(5,1) - t166 * rSges(5,2) - t184 * t278 + t214;
t112 = t165 * rSges(5,1) - t164 * rSges(5,2) + t186 * t278 + t203;
t95 = t112 * t228;
t267 = m(5) * (-t187 * t229 + t95);
t96 = t112 * t225;
t266 = m(5) * (-t187 * t226 + t96);
t265 = m(5) * (t112 * t184 + t187 * t186);
t70 = t164 * t129 - t191 * t227;
t192 = t186 * t70;
t61 = t70 * t228;
t264 = m(6) * (t192 * t181 - t61);
t62 = t70 * t225;
t263 = m(6) * (t192 * t183 - t62);
t59 = t68 * t228;
t261 = m(6) * (-t67 * t229 + t59);
t60 = t68 * t225;
t260 = m(6) * (-t67 * t226 + t60);
t259 = m(6) * (-t229 * t291 - t61);
t258 = m(6) * (-t226 * t291 - t62);
t63 = t67 * t186;
t257 = m(6) * (t68 * t184 + t63);
t256 = m(6) * (-t70 * t184 + t291 * t186);
t255 = t181 * t282;
t254 = t183 * t282;
t253 = m(6) * (-t122 * t186 - t184 * t123);
t248 = m(3) * ((rSges(3,2) * t229 + t186 * rSges(3,3) + t177) * t186 + (-rSges(3,2) * t228 + (rSges(3,3) + qJ(2)) * t184) * t184);
t244 = m(6) * qJD(5);
t239 = Icges(6,4) * t142;
t232 = t112 * t186;
t211 = t270 + t272;
t88 = (m(4) / 0.2e1 + t211) * t276;
t222 = t88 * qJD(1);
t102 = Icges(6,2) * t141 + Icges(6,6) * t164 + t239;
t221 = -Icges(6,1) * t141 + t102 + t239;
t220 = -Icges(6,1) * t145 - t103 - t238;
t137 = Icges(6,4) * t141;
t105 = Icges(6,1) * t142 + Icges(6,5) * t164 + t137;
t219 = -Icges(6,2) * t142 + t105 + t137;
t218 = Icges(6,2) * t144 - t106 + t138;
t216 = t127 - t134;
t114 = t211 * t275;
t212 = t114 * qJD(1);
t210 = m(6) / 0.4e1 + m(5) / 0.4e1;
t208 = (t284 - t164 / 0.2e1) * (t141 * t102 + t142 * t105 + t164 * (Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t164));
t207 = (t166 / 0.2e1 + t283) * ((Icges(6,5) * t144 - Icges(6,6) * t145 + Icges(6,3) * t166) * t166 - t103 * t145 + t106 * t144);
t135 = rSges(6,1) * t162 + rSges(6,2) * t163;
t79 = t122 * t227 - t135 * t164;
t80 = -t123 * t227 - t135 * t166;
t197 = t184 * t79 + t186 * t80;
t116 = Icges(6,5) * t141 - Icges(6,6) * t142;
t117 = Icges(6,5) * t145 + Icges(6,6) * t144;
t196 = t164 * t116 - t166 * t117;
t190 = -t221 * t164 + t220 * t166;
t189 = t219 * t164 - t218 * t166;
t113 = t210 * t275 - t286 * t169 / 0.2e1;
t87 = (m(4) / 0.4e1 + t210) * t276 - (m(4) + t286) * t168 / 0.2e1;
t77 = t253 / 0.2e1;
t74 = t254 / 0.2e1;
t73 = t255 / 0.2e1;
t58 = t122 * t166 + t123 * t164;
t48 = (t132 * t227 + t216 * t163 + t281) * t227;
t42 = t256 / 0.2e1;
t37 = t258 / 0.2e1;
t36 = t259 / 0.2e1;
t31 = -t166 * t132 + t144 * t216 + t215 * t145;
t30 = t164 * t132 + t215 * t141 - t216 * t142;
t29 = t117 * t227 + t218 * t162 + t220 * t163;
t28 = t116 * t227 + t219 * t162 + t221 * t163;
t26 = t260 + t266;
t25 = t261 + t267 + t269;
t22 = t248 + t257 + t265 + t268;
t21 = (t128 / 0.2e1 + t133 / 0.2e1) * t162 + t288;
t18 = t263 / 0.2e1;
t17 = t264 / 0.2e1;
t13 = t42 - t253 / 0.2e1;
t12 = t77 + t42;
t11 = t77 - t256 / 0.2e1;
t9 = t37 + t18 - t254 / 0.2e1;
t8 = t36 + t17 - t255 / 0.2e1;
t7 = t74 + t37 - t263 / 0.2e1;
t6 = t74 + t18 - t258 / 0.2e1;
t5 = t73 + t36 - t264 / 0.2e1;
t4 = t73 + t17 - t259 / 0.2e1;
t1 = t207 * t164 + t208 * t166;
t2 = [-t68 * t273 * t294 / 0.4e1 + t22 * qJD(2) + t25 * qJD(3) + t26 * qJD(4) + t21 * qJD(5), qJD(1) * t22 + qJD(3) * t87 + qJD(4) * t113 + qJD(5) * t12, qJD(1) * t25 + qJD(2) * t87 + qJD(5) * t5, qJD(1) * t26 + qJD(2) * t113 + qJD(5) * t7, t21 * qJD(1) + t12 * qJD(2) + t5 * qJD(3) + t7 * qJD(4) + (t48 + m(6) * (-t122 * t70 - t123 * t291 + t67 * t80 + t68 * t79) + (-t29 / 0.2e1 - t31 / 0.2e1 - t208) * t166 + (t28 / 0.2e1 + t30 / 0.2e1 - t207) * t164) * qJD(5); t88 * qJD(3) + t114 * qJD(4) + t11 * qJD(5) + (-t257 / 0.4e1 - t268 / 0.4e1 - t265 / 0.4e1 - t248 / 0.4e1) * t273 + (-t186 * t69 + t63) * t270 * t274, 0, t222, t212, t11 * qJD(1) + (t80 * t184 - t186 * t79) * t244; -t88 * qJD(2) + t4 * qJD(5) + (-t261 / 0.4e1 - t269 / 0.4e1 - t267 / 0.4e1) * t273 + ((t287 * t181 + t59) * t270 + (-t232 * t181 + t95) * t272) * t274, -t222, 0, 0, t4 * qJD(1) + (t197 * t181 - t58 * t183) * t244; -t114 * qJD(2) + t6 * qJD(5) + (-t260 / 0.4e1 - t266 / 0.4e1) * t273 + ((-t232 * t183 + t96) * t272 + (t287 * t183 + t60) * t270) * t274, -t212, 0, 0, t6 * qJD(1) + (t58 * t181 + t197 * t183) * t244; t13 * qJD(2) + t8 * qJD(3) + t9 * qJD(4) + t1 * qJD(5) + (t70 * t294 - t281 / 0.2e1 - t288) * qJD(1), t13 * qJD(1), t8 * qJD(1), t9 * qJD(1), t1 * qJD(1) + (m(6) * (t291 * t80 + (t164 * t109 + t166 * t191) * t58 - t70 * t79) + (t189 * t141 + t190 * t142 + t196 * t164 + t30 * t227) * t284 + (-t144 * t190 + t189 * t145 - t196 * t166 + t31 * t227) * t283 + (t28 * t164 - t29 * t166 + t48) * t205) * qJD(5);];
Cq = t2;
