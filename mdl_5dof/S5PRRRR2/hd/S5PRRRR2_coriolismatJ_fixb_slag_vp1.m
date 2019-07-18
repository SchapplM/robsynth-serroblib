% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:30:12
% EndTime: 2019-07-18 13:30:15
% DurationCPUTime: 2.10s
% Computational Cost: add. (16470->210), mult. (11345->287), div. (0->0), fcn. (9588->8), ass. (0->143)
t208 = qJ(2) + qJ(3);
t207 = qJ(4) + t208;
t202 = sin(t207);
t203 = cos(t207);
t169 = -rSges(5,1) * t202 - rSges(5,2) * t203;
t204 = sin(t208);
t274 = pkin(3) * t204;
t157 = t169 - t274;
t170 = t203 * rSges(5,1) - t202 * rSges(5,2);
t205 = cos(t208);
t273 = pkin(3) * t205;
t158 = t170 + t273;
t241 = t170 * t157 - t158 * t169;
t211 = cos(qJ(5));
t252 = t202 * t211;
t209 = sin(qJ(5));
t253 = t202 * t209;
t220 = rSges(6,1) * t252 - rSges(6,2) * t253 - t203 * rSges(6,3);
t139 = t203 * pkin(6) - t220;
t129 = t139 - t274;
t193 = rSges(6,1) * t211 - rSges(6,2) * t209;
t229 = t193 * t203;
t138 = (rSges(6,3) + pkin(6)) * t202 + t229;
t130 = t138 + t273;
t269 = t138 * t129 - t130 * t139;
t305 = m(6) / 0.2e1;
t306 = m(5) / 0.2e1;
t275 = sin(qJ(2)) * pkin(2);
t123 = t129 - t275;
t272 = cos(qJ(2)) * pkin(2);
t124 = t130 + t272;
t55 = -t138 * t123 + t124 * t139;
t155 = t157 - t275;
t156 = t158 + t272;
t83 = -t170 * t155 + t156 * t169;
t270 = (t269 + t55) * t305 + (t241 + t83) * t306;
t271 = (-t269 + t55) * t305 + (-t241 + t83) * t306;
t4 = t271 - t270;
t320 = t4 * qJD(2);
t206 = Icges(6,4) * t211;
t188 = -Icges(6,2) * t209 + t206;
t189 = Icges(6,1) * t209 + t206;
t319 = t188 + t189;
t278 = m(4) * (t272 * (-rSges(4,1) * t204 - rSges(4,2) * t205) + (t205 * rSges(4,1) - t204 * rSges(4,2)) * t275);
t265 = Icges(6,4) * t209;
t187 = Icges(6,2) * t211 + t265;
t190 = Icges(6,1) * t211 - t265;
t230 = t319 * t211 / 0.2e1 + (-t187 / 0.2e1 + t190 / 0.2e1) * t209;
t191 = rSges(6,1) * t209 + rSges(6,2) * t211;
t165 = t191 * t202;
t115 = t123 * t165;
t166 = t191 * t203;
t63 = -t124 * t166 + t115;
t69 = -t138 * t166 + t139 * t165;
t295 = m(6) * (t69 + t63);
t316 = t230 + t295 / 0.2e1;
t64 = t129 * t165 - t130 * t166;
t294 = m(6) * (t69 + t64);
t315 = t230 + t294 / 0.2e1;
t79 = -t158 * t155 + t156 * t157;
t313 = t202 ^ 2;
t312 = t203 ^ 2;
t311 = 0.4e1 * qJD(2);
t308 = 2 * qJD(4);
t302 = m(5) * t79;
t300 = m(5) * t83;
t299 = m(5) * t241;
t296 = m(6) * (t64 + t63);
t293 = m(6) * (t115 + (-t129 * t202 + (-t124 + t130) * t203) * t191);
t292 = m(6) * (t115 + (-t139 * t202 + (-t124 + t138) * t203) * t191);
t291 = m(6) * ((-t130 + t138) * t203 + (t129 - t139) * t202) * t191;
t50 = -t130 * t123 + t124 * t129;
t290 = m(6) * t50;
t288 = m(6) * t55;
t287 = m(6) * t269;
t285 = m(6) * t63;
t284 = m(6) * t64;
t283 = m(6) * t69;
t282 = -t202 / 0.2e1;
t281 = t202 / 0.2e1;
t280 = -t203 / 0.2e1;
t146 = Icges(6,5) * t252 - Icges(6,6) * t253 - Icges(6,3) * t203;
t149 = Icges(6,6) * t202 + t188 * t203;
t248 = t209 * t149;
t231 = -t146 + t248;
t151 = Icges(6,5) * t202 + t190 * t203;
t244 = t211 * t151;
t135 = t202 * t244;
t186 = Icges(6,5) * t211 - Icges(6,6) * t209;
t250 = t203 * t186;
t147 = Icges(6,3) * t202 + t250;
t232 = t203 * t147 - t135;
t239 = t202 * t147 + t203 * t244;
t183 = Icges(6,4) * t253;
t150 = Icges(6,1) * t252 - Icges(6,5) * t203 - t183;
t245 = t211 * t150;
t240 = -t202 * t146 - t203 * t245;
t148 = Icges(6,4) * t252 - Icges(6,2) * t253 - Icges(6,6) * t203;
t249 = t209 * t148;
t72 = -t203 * t249 - t240;
t73 = -t203 * t248 + t239;
t11 = (t203 * t231 - t239 + t73) * t203 + (t202 * t231 + t232 + t72) * t202;
t71 = -t202 * t248 - t232;
t12 = (t71 - t135 + (t147 + t249) * t203 + t240) * t203 + t239 * t202;
t40 = t202 * t71 - t203 * (-t202 * (-t245 + t249) - t203 * t146);
t41 = t202 * t73 - t203 * t72;
t2 = (t41 / 0.2e1 - t12 / 0.2e1) * t203 + (t11 / 0.2e1 + t40 / 0.2e1) * t202;
t279 = t2 * qJD(5);
t238 = t189 * t202 + t148;
t237 = -t189 * t203 - t149;
t236 = -Icges(6,2) * t252 + t150 - t183;
t235 = -t187 * t203 + t151;
t228 = t296 / 0.2e1 + t230;
t224 = Icges(6,5) * t209 + Icges(6,6) * t211;
t219 = (-t165 * t203 + t166 * t202) * t191;
t214 = (-t187 + t190) * t211 - t319 * t209;
t218 = t203 * t12 / 0.2e1 + (t11 + t40) * t282 + (t202 * t186 + t203 * t214 + t209 * t237 + t211 * t235) * t281 + (t202 * t214 - t209 * t238 + t211 * t236 - t250 + t41) * t280;
t217 = -t230 + (t281 + t282) * (t211 * t148 + t209 * t150);
t216 = t209 * t236 + t211 * t238;
t215 = -t209 * t235 + t211 * t237;
t160 = t203 * t224;
t159 = t224 * t202;
t118 = -t165 * t202 - t166 * t203;
t51 = t230 + t283;
t49 = t230 + t284;
t48 = t230 + t285;
t44 = t291 / 0.2e1;
t42 = t292 / 0.2e1;
t37 = -t287 - t299;
t36 = t293 / 0.2e1;
t30 = t288 + t300;
t19 = t278 + t290 + t302;
t18 = -t291 / 0.2e1 + t315;
t17 = t44 + t315;
t16 = -t292 / 0.2e1 + t316;
t15 = t42 + t316;
t14 = -t293 / 0.2e1 + t228;
t13 = t36 + t228;
t8 = t44 - t294 / 0.2e1 + t217;
t7 = t42 - t295 / 0.2e1 + t217;
t6 = t36 - t296 / 0.2e1 + t217;
t3 = t270 + t271;
t1 = [0, 0, 0, 0, m(6) * t118 * qJD(5); 0, qJD(3) * t19 + qJD(4) * t30 + qJD(5) * t48, t19 * qJD(2) + t3 * qJD(4) + t13 * qJD(5) + 0.2e1 * (t278 / 0.2e1 + t50 * t305 + t79 * t306) * qJD(3), t30 * qJD(2) + t3 * qJD(3) + t15 * qJD(5) + (t55 * t305 + t83 * t306) * t308, t48 * qJD(2) + t13 * qJD(3) + t15 * qJD(4) + (((-t123 * t203 - t124 * t202) * t193 + t219) * m(6) + t218) * qJD(5); 0, t4 * qJD(4) + t14 * qJD(5) + (-t290 / 0.4e1 - t302 / 0.4e1 - t278 / 0.4e1) * t311, qJD(4) * t37 + qJD(5) * t49, t320 + t37 * qJD(3) + t17 * qJD(5) + (-t241 * t306 - t269 * t305) * t308, t14 * qJD(2) + t49 * qJD(3) + t17 * qJD(4) + (((-t129 * t203 - t130 * t202) * t193 + t219) * m(6) + t218) * qJD(5); 0, -t4 * qJD(3) + t16 * qJD(5) + (-t288 / 0.4e1 - t300 / 0.4e1) * t311, -t320 + t18 * qJD(5) + 0.4e1 * (t299 / 0.4e1 + t287 / 0.4e1) * qJD(3), qJD(5) * t51, t16 * qJD(2) + t18 * qJD(3) + t51 * qJD(4) + (((-t138 * t202 - t139 * t203) * t193 + t219) * m(6) + t218) * qJD(5); 0, (t217 - t285) * qJD(2) + t6 * qJD(3) + t7 * qJD(4) + t279, t6 * qJD(2) + (t217 - t284) * qJD(3) + t8 * qJD(4) + t279, t7 * qJD(2) + t8 * qJD(3) + (t217 - t283) * qJD(4) + t279, (m(6) * ((t202 * t220 + t203 * (t202 * rSges(6,3) + t229)) * t118 + (t312 + t313) * t193 * t191) + (-t313 * t160 + (t216 * t203 + (t159 + t215) * t202) * t203) * t281 + (-t312 * t159 + (t215 * t202 + (t160 + t216) * t203) * t202) * t280) * qJD(5) + (qJD(2) + qJD(3) + qJD(4)) * t2;];
Cq  = t1;
