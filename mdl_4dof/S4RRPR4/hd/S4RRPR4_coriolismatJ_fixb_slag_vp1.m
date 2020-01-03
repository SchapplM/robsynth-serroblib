% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:26
% EndTime: 2019-12-31 17:02:29
% DurationCPUTime: 1.56s
% Computational Cost: add. (9292->178), mult. (7217->247), div. (0->0), fcn. (6236->8), ass. (0->129)
t257 = m(5) / 0.2e1;
t271 = m(4) / 0.2e1;
t174 = qJ(1) + qJ(2);
t171 = sin(t174);
t172 = cos(t174);
t233 = sin(qJ(1)) * pkin(1);
t175 = cos(pkin(7));
t173 = pkin(7) + qJ(4);
t170 = cos(t173);
t228 = t170 * rSges(5,1);
t197 = pkin(3) * t175 + pkin(2) + t228;
t169 = sin(t173);
t214 = t169 * t171;
t202 = rSges(5,2) * t214 + t172 * rSges(5,3);
t231 = -pkin(6) - qJ(3);
t95 = -t197 * t171 - t172 * t231 + t202;
t93 = t95 - t233;
t234 = cos(qJ(1)) * pkin(1);
t213 = t169 * t172;
t196 = -rSges(5,2) * t213 + t171 * rSges(5,3);
t96 = -t171 * t231 + t197 * t172 + t196;
t94 = t96 + t234;
t48 = t94 * t171 + t93 * t172;
t55 = t96 * t171 + t95 * t172;
t268 = rSges(4,2) * sin(pkin(7)) - rSges(4,1) * t175 - pkin(2);
t269 = rSges(4,3) + qJ(3);
t108 = t268 * t171 + t269 * t172;
t103 = t108 - t233;
t109 = t269 * t171 - t268 * t172;
t104 = t109 + t234;
t66 = t103 * t172 + t104 * t171;
t69 = t108 * t172 + t109 * t171;
t229 = (t55 + t48) * t257 + (t69 + t66) * t271;
t230 = (t48 - t55) * t257 + (t66 - t69) * t271;
t4 = t230 - t229;
t273 = t4 * qJD(1);
t160 = Icges(5,4) * t170;
t143 = -Icges(5,2) * t169 + t160;
t144 = Icges(5,1) * t169 + t160;
t272 = t143 + t144;
t256 = m(3) * (t234 * (-rSges(3,1) * t171 - rSges(3,2) * t172) + (t172 * rSges(3,1) - t171 * rSges(3,2)) * t233);
t252 = m(4) * (-t109 * t103 + t104 * t108);
t31 = -t96 * t93 + t94 * t95;
t201 = qJD(1) + qJD(2);
t146 = rSges(5,1) * t169 + rSges(5,2) * t170;
t129 = t146 * t171;
t130 = t146 * t172;
t185 = t171 * t129 + t172 * t130;
t179 = t185 * t257;
t263 = t172 ^ 2;
t264 = t171 ^ 2;
t266 = t146 * (t263 + t264);
t184 = m(5) * t266;
t61 = t179 + t184 / 0.2e1;
t267 = t201 * t61;
t224 = Icges(5,4) * t169;
t142 = Icges(5,2) * t170 + t224;
t145 = Icges(5,1) * t170 - t224;
t191 = t272 * t170 / 0.2e1 + (-t142 / 0.2e1 + t145 / 0.2e1) * t169;
t119 = Icges(5,5) * t171 + t145 * t172;
t203 = -t142 * t172 + t119;
t152 = Icges(5,4) * t214;
t210 = t170 * t171;
t118 = Icges(5,1) * t210 - Icges(5,5) * t172 - t152;
t204 = -Icges(5,2) * t210 + t118 - t152;
t117 = Icges(5,6) * t171 + t143 * t172;
t205 = -t144 * t172 - t117;
t116 = Icges(5,4) * t210 - Icges(5,2) * t214 - Icges(5,6) * t172;
t206 = t144 * t171 + t116;
t265 = (-t203 * t171 + t204 * t172) * t169 + (t205 * t171 + t206 * t172) * t170;
t262 = 0.4e1 * qJD(1);
t250 = m(4) * t66;
t249 = m(4) * t69;
t40 = t93 * t129 - t94 * t130;
t41 = t95 * t129 - t96 * t130;
t248 = m(5) * (t41 + t40);
t247 = m(5) * ((-t94 + t96) * t172 + (t93 - t95) * t171) * t146;
t244 = m(5) * t31;
t242 = m(5) * t40;
t241 = m(5) * t41;
t240 = m(5) * t48;
t239 = m(5) * t55;
t238 = -t171 / 0.2e1;
t237 = t171 / 0.2e1;
t236 = -t172 / 0.2e1;
t114 = Icges(5,5) * t210 - Icges(5,6) * t214 - Icges(5,3) * t172;
t192 = t117 * t169 - t114;
t105 = t119 * t210;
t141 = Icges(5,5) * t170 - Icges(5,6) * t169;
t218 = t141 * t172;
t115 = Icges(5,3) * t171 + t218;
t193 = t172 * t115 - t105;
t209 = t170 * t172;
t207 = t171 * t115 + t119 * t209;
t208 = -t171 * t114 - t118 * t209;
t51 = -t116 * t213 - t208;
t52 = -t117 * t213 + t207;
t11 = (t192 * t172 - t207 + t52) * t172 + (t192 * t171 + t193 + t51) * t171;
t219 = t116 * t169;
t50 = -t117 * t214 - t193;
t12 = (t50 - t105 + (t115 + t219) * t172 + t208) * t172 + t207 * t171;
t22 = t171 * t50 - t172 * (-(-t118 * t170 + t219) * t171 - t114 * t172);
t23 = t171 * t52 - t172 * t51;
t2 = (t23 / 0.2e1 - t12 / 0.2e1) * t172 + (t11 / 0.2e1 + t22 / 0.2e1) * t171;
t235 = -t61 * qJD(3) + t2 * qJD(4);
t190 = t248 / 0.2e1 + t191;
t188 = Icges(5,5) * t169 + Icges(5,6) * t170;
t182 = (-t129 * t172 + t130 * t171) * t146;
t178 = (-t142 + t145) * t170 - t272 * t169;
t181 = t172 * t12 / 0.2e1 + (t11 + t22) * t238 + (t141 * t171 + t205 * t169 + t203 * t170 + t178 * t172) * t237 + (-t206 * t169 + t204 * t170 + t178 * t171 - t218 + t23) * t236;
t180 = -t191 + (t237 + t238) * (t116 * t170 + t118 * t169);
t147 = -rSges(5,2) * t169 + t228;
t124 = t172 * t188;
t123 = t188 * t171;
t60 = t179 - t184 / 0.2e1;
t58 = t61 * qJD(4);
t57 = t60 * qJD(3);
t56 = t60 * qJD(4);
t32 = t239 + t249;
t30 = t240 + t250;
t25 = t191 + t241;
t24 = t191 + t242;
t18 = t247 / 0.2e1;
t13 = t244 + t252 + t256;
t8 = -t247 / 0.2e1 + t190;
t7 = t18 + t190;
t5 = t229 + t230;
t3 = t18 - t248 / 0.2e1 + t180;
t1 = [qJD(2) * t13 + qJD(3) * t30 + qJD(4) * t24, t13 * qJD(1) + t5 * qJD(3) + t7 * qJD(4) + 0.2e1 * (t252 / 0.2e1 + t31 * t257 + t256 / 0.2e1) * qJD(2), qJD(1) * t30 + qJD(2) * t5 + t56, t24 * qJD(1) + t7 * qJD(2) + t57 + ((-t48 * t147 + t182) * m(5) + t181) * qJD(4); -t4 * qJD(3) + t8 * qJD(4) + (-t256 / 0.4e1 - t252 / 0.4e1 - t244 / 0.4e1) * t262, qJD(3) * t32 + qJD(4) * t25, qJD(2) * t32 - t273 + t56, t8 * qJD(1) + t25 * qJD(2) + t57 + ((-t55 * t147 + t182) * m(5) + t181) * qJD(4); t4 * qJD(2) + t58 + (-t250 / 0.4e1 - t240 / 0.4e1) * t262, t273 + t58 + 0.4e1 * (-t249 / 0.4e1 - t239 / 0.4e1) * qJD(2), 0, t267; (t180 - t242) * qJD(1) + t3 * qJD(2) + t235, t3 * qJD(1) + (t180 - t241) * qJD(2) + t235, -t267, (m(5) * (t147 * t266 - (t171 * (rSges(5,1) * t210 - t202) + t172 * (rSges(5,1) * t209 + t196)) * t185) + (-t264 * t124 + (t171 * t123 + t265) * t172) * t237 + (-t263 * t123 + (t172 * t124 + t265) * t171) * t236) * qJD(4) + t201 * t2;];
Cq = t1;
