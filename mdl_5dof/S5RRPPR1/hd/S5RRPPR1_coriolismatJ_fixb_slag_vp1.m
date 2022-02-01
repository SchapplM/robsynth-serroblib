% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:19
% EndTime: 2022-01-20 09:51:25
% DurationCPUTime: 1.90s
% Computational Cost: add. (14256->192), mult. (7992->260), div. (0->0), fcn. (6878->10), ass. (0->136)
t285 = m(6) / 0.2e1;
t299 = m(5) / 0.2e1;
t195 = qJ(1) + qJ(2);
t191 = pkin(8) + t195;
t186 = sin(t191);
t187 = cos(t191);
t194 = pkin(9) + qJ(5);
t189 = sin(t194);
t236 = t187 * t189;
t217 = -rSges(6,2) * t236 + t186 * rSges(6,3);
t196 = cos(pkin(9));
t190 = cos(t194);
t251 = t190 * rSges(6,1);
t218 = pkin(4) * t196 + pkin(3) + t251;
t257 = -pkin(7) - qJ(4);
t193 = cos(t195);
t259 = pkin(2) * t193;
t106 = -t186 * t257 + t218 * t187 + t217 + t259;
t261 = cos(qJ(1)) * pkin(1);
t100 = t106 + t261;
t238 = t186 * t189;
t224 = rSges(6,2) * t238 + t187 * rSges(6,3);
t192 = sin(t195);
t260 = pkin(2) * t192;
t105 = -t218 * t186 - t187 * t257 + t224 - t260;
t258 = sin(qJ(1)) * pkin(1);
t99 = t105 - t258;
t46 = t100 * t186 + t99 * t187;
t49 = t105 * t187 + t106 * t186;
t296 = rSges(5,2) * sin(pkin(9)) - rSges(5,1) * t196 - pkin(3);
t297 = rSges(5,3) + qJ(4);
t113 = t296 * t186 + t297 * t187 - t260;
t110 = t113 - t258;
t114 = t297 * t186 - t296 * t187 + t259;
t111 = t114 + t261;
t58 = t110 * t187 + t111 * t186;
t67 = t113 * t187 + t114 * t186;
t255 = (t49 + t46) * t285 + (t67 + t58) * t299;
t256 = (t46 - t49) * t285 + (t58 - t67) * t299;
t3 = t256 - t255;
t301 = t3 * qJD(1);
t184 = Icges(6,4) * t190;
t160 = -Icges(6,2) * t189 + t184;
t161 = Icges(6,1) * t189 + t184;
t300 = t160 + t161;
t284 = m(4) * (t261 * (-rSges(4,1) * t186 - rSges(4,2) * t187 - t260) + (t187 * rSges(4,1) - t186 * rSges(4,2) + t259) * t258);
t280 = m(5) * (-t114 * t110 + t111 * t113);
t263 = m(3) * (t261 * (-rSges(3,1) * t192 - rSges(3,2) * t193) + (t193 * rSges(3,1) - t192 * rSges(3,2)) * t258);
t223 = qJD(1) + qJD(2);
t163 = rSges(6,1) * t189 + rSges(6,2) * t190;
t144 = t163 * t186;
t145 = t163 * t187;
t208 = t186 * t144 + t187 * t145;
t202 = t208 * t285;
t292 = t187 ^ 2;
t293 = t186 ^ 2;
t294 = t163 * (t292 + t293);
t207 = m(6) * t294;
t69 = t202 + t207 / 0.2e1;
t295 = t223 * t69;
t249 = Icges(6,4) * t189;
t159 = Icges(6,2) * t190 + t249;
t162 = Icges(6,1) * t190 - t249;
t214 = t300 * t190 / 0.2e1 + (-t159 / 0.2e1 + t162 / 0.2e1) * t189;
t291 = 0.4e1 * qJD(1);
t278 = m(5) * t58;
t277 = m(5) * t67;
t40 = -t100 * t145 + t99 * t144;
t41 = t105 * t144 - t106 * t145;
t276 = m(6) * (t41 + t40);
t275 = m(6) * ((-t100 + t106) * t187 + (-t105 + t99) * t186) * t163;
t32 = t100 * t105 - t106 * t99;
t272 = m(6) * t32;
t270 = m(6) * t40;
t269 = m(6) * t41;
t268 = m(6) * t46;
t267 = m(6) * t49;
t266 = -t186 / 0.2e1;
t265 = t186 / 0.2e1;
t264 = -t187 / 0.2e1;
t68 = t202 - t207 / 0.2e1;
t254 = t68 * qJD(4);
t237 = t186 * t190;
t127 = Icges(6,4) * t237 - Icges(6,2) * t238 - Icges(6,6) * t187;
t242 = t127 * t189;
t158 = Icges(6,5) * t190 - Icges(6,6) * t189;
t240 = t158 * t187;
t235 = t187 * t190;
t125 = Icges(6,5) * t237 - Icges(6,6) * t238 - Icges(6,3) * t187;
t169 = Icges(6,4) * t238;
t129 = Icges(6,1) * t237 - Icges(6,5) * t187 - t169;
t230 = -t186 * t125 - t129 * t235;
t126 = Icges(6,3) * t186 + t240;
t130 = Icges(6,5) * t186 + t162 * t187;
t229 = t186 * t126 + t130 * t235;
t228 = t161 * t186 + t127;
t128 = Icges(6,6) * t186 + t160 * t187;
t227 = -t161 * t187 - t128;
t226 = -Icges(6,2) * t237 + t129 - t169;
t225 = -t159 * t187 + t130;
t215 = t128 * t189 - t125;
t118 = t130 * t237;
t216 = t126 * t187 - t118;
t52 = -t127 * t236 - t230;
t53 = -t128 * t236 + t229;
t11 = (t215 * t187 - t229 + t53) * t187 + (t215 * t186 + t216 + t52) * t186;
t51 = -t128 * t238 - t216;
t12 = (t51 - t118 + (t126 + t242) * t187 + t230) * t187 + t229 * t186;
t27 = t186 * t51 - t187 * (-(-t129 * t190 + t242) * t186 - t125 * t187);
t28 = t186 * t53 - t187 * t52;
t2 = (t28 / 0.2e1 - t12 / 0.2e1) * t187 + (t11 / 0.2e1 + t27 / 0.2e1) * t186;
t222 = -t69 * qJD(4) + t2 * qJD(5);
t213 = t276 / 0.2e1 + t214;
t211 = Icges(6,5) * t189 + Icges(6,6) * t190;
t205 = (-t144 * t187 + t145 * t186) * t163;
t199 = (-t159 + t162) * t190 - t300 * t189;
t204 = t187 * t12 / 0.2e1 + (t11 + t27) * t266 + (t158 * t186 + t199 * t187 + t227 * t189 + t225 * t190) * t265 + (t199 * t186 - t228 * t189 + t226 * t190 - t240 + t28) * t264;
t203 = -t214 + (t265 + t266) * (t190 * t127 + t189 * t129);
t201 = t226 * t189 + t228 * t190;
t200 = -t225 * t189 + t227 * t190;
t164 = -rSges(6,2) * t189 + t251;
t139 = t187 * t211;
t138 = t211 * t186;
t64 = t69 * qJD(5);
t62 = t68 * qJD(5);
t31 = t267 + t277;
t30 = t214 + t269;
t29 = t214 + t270;
t26 = t268 + t278;
t18 = t275 / 0.2e1;
t13 = t263 + t272 + t280 + t284;
t8 = -t275 / 0.2e1 + t213;
t7 = t18 + t213;
t6 = t18 - t276 / 0.2e1 + t203;
t4 = t255 + t256;
t1 = [t13 * qJD(2) + t26 * qJD(4) + t29 * qJD(5), t13 * qJD(1) + t4 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t280 / 0.2e1 + t263 / 0.2e1 + t32 * t285 + t284 / 0.2e1) * qJD(2), 0, qJD(1) * t26 + qJD(2) * t4 + t62, t29 * qJD(1) + t7 * qJD(2) + ((-t46 * t164 + t205) * m(6) + t204) * qJD(5) + t254; -t3 * qJD(4) + t8 * qJD(5) + (-t272 / 0.4e1 - t284 / 0.4e1 - t280 / 0.4e1 - t263 / 0.4e1) * t291, qJD(4) * t31 + qJD(5) * t30, 0, qJD(2) * t31 - t301 + t62, t8 * qJD(1) + t30 * qJD(2) + ((-t49 * t164 + t205) * m(6) + t204) * qJD(5) + t254; 0, 0, 0, 0, -m(6) * t208 * qJD(5); t3 * qJD(2) + t64 + (-t268 / 0.4e1 - t278 / 0.4e1) * t291, t301 + t64 + 0.4e1 * (-t267 / 0.4e1 - t277 / 0.4e1) * qJD(2), 0, 0, t295; (t203 - t270) * qJD(1) + t6 * qJD(2) + t222, t6 * qJD(1) + (t203 - t269) * qJD(2) + t222, 0, -t295, (m(6) * (t164 * t294 - t208 * (t186 * (rSges(6,1) * t237 - t224) + t187 * (rSges(6,1) * t235 + t217))) + (-t293 * t139 + (t201 * t187 + (t138 + t200) * t186) * t187) * t265 + (-t292 * t138 + (t200 * t186 + (t139 + t201) * t187) * t186) * t264) * qJD(5) + t223 * t2;];
Cq = t1;
