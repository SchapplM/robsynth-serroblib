% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:33:38
% EndTime: 2020-01-03 11:33:46
% DurationCPUTime: 1.98s
% Computational Cost: add. (13630->193), mult. (7600->252), div. (0->0), fcn. (6548->10), ass. (0->132)
t181 = qJ(1) + pkin(8);
t178 = qJ(3) + t181;
t171 = sin(t178);
t172 = cos(t178);
t182 = cos(pkin(9));
t273 = -rSges(5,2) * sin(pkin(9)) + pkin(3) + rSges(5,1) * t182;
t274 = -qJ(4) - rSges(5,3);
t110 = t273 * t171 + t274 * t172;
t196 = sin(qJ(1)) * pkin(1) + pkin(2) * sin(t181);
t100 = t110 + t196;
t173 = pkin(4) * t182 + pkin(3);
t183 = -pkin(7) - qJ(4);
t180 = pkin(9) + qJ(5);
t174 = sin(t180);
t218 = t171 * t174;
t204 = -rSges(6,2) * t218 - t172 * rSges(6,3);
t176 = cos(t180);
t233 = rSges(6,1) * t176;
t102 = t172 * t183 + (t173 + t233) * t171 + t204;
t111 = -t274 * t171 + t273 * t172;
t104 = t111 * t171;
t276 = m(6) / 0.2e1;
t277 = m(5) / 0.2e1;
t215 = t172 * t176;
t216 = t172 * t174;
t197 = rSges(6,1) * t215 - rSges(6,2) * t216;
t103 = t172 * t173 + (rSges(6,3) - t183) * t171 + t197;
t203 = pkin(2) * cos(t181) + cos(qJ(1)) * pkin(1);
t96 = t103 + t203;
t91 = t96 * t171;
t95 = t102 + t196;
t101 = t111 + t203;
t98 = t101 * t171;
t99 = t103 * t171;
t238 = (t91 + t99 + (-t102 - t95) * t172) * t276 + (t104 + t98 + (-t100 - t110) * t172) * t277;
t230 = t102 - t95;
t239 = (t230 * t172 + t91 - t99) * t276 + (-t104 + t98 + (-t100 + t110) * t172) * t277;
t4 = t239 - t238;
t279 = t4 * qJD(1);
t169 = Icges(6,4) * t176;
t147 = -Icges(6,2) * t174 + t169;
t268 = Icges(6,1) * t174 + t169;
t278 = t147 + t268;
t243 = m(6) * (-t102 * t172 + t99);
t244 = m(6) * (-t172 * t95 + t91);
t260 = m(4) * (-t203 * (rSges(4,1) * t171 + rSges(4,2) * t172) + t196 * (t172 * rSges(4,1) - rSges(4,2) * t171));
t256 = m(5) * (t100 * t111 - t110 * t101);
t248 = m(6) * (-t102 * t96 + t95 * t103);
t150 = rSges(6,1) * t174 + rSges(6,2) * t176;
t131 = t150 * t171;
t132 = t150 * t172;
t97 = -t171 * t131 - t172 * t132;
t189 = m(6) * t97;
t202 = qJD(1) + qJD(3);
t263 = t172 ^ 2;
t264 = t171 ^ 2;
t269 = t150 * (t263 + t264);
t188 = -m(6) * t269 / 0.2e1;
t63 = t189 / 0.2e1 + t188;
t272 = t202 * t63;
t228 = Icges(6,4) * t174;
t146 = Icges(6,2) * t176 + t228;
t149 = Icges(6,1) * t176 - t228;
t267 = (t146 - t149) * t176 + t278 * t174;
t154 = Icges(6,4) * t216;
t119 = Icges(6,1) * t215 + t171 * Icges(6,5) - t154;
t207 = -Icges(6,2) * t215 + t119 - t154;
t117 = Icges(6,4) * t215 - Icges(6,2) * t216 + t171 * Icges(6,6);
t209 = t172 * t268 + t117;
t266 = -t207 * t174 - t176 * t209;
t118 = -Icges(6,5) * t172 + t149 * t171;
t208 = -t146 * t171 + t118;
t116 = -Icges(6,6) * t172 + t147 * t171;
t210 = t171 * t268 + t116;
t265 = -t208 * t174 - t210 * t176;
t198 = t278 * t176 / 0.2e1 + (-t146 / 0.2e1 + t149 / 0.2e1) * t174;
t262 = 0.4e1 * qJD(1);
t254 = m(5) * (-t100 * t172 + t98);
t253 = m(5) * (-t110 * t172 + t104);
t39 = -t95 * t131 - t96 * t132;
t41 = -t102 * t131 - t103 * t132;
t252 = m(6) * (t41 + t39);
t251 = m(6) * ((t103 - t96) * t172 + t230 * t171) * t150;
t246 = m(6) * t39;
t245 = m(6) * t41;
t242 = -t171 / 0.2e1;
t241 = -t172 / 0.2e1;
t240 = t172 / 0.2e1;
t64 = -t189 / 0.2e1 + t188;
t236 = t64 * qJD(4);
t223 = t116 * t174;
t222 = t117 * t174;
t221 = t118 * t176;
t145 = Icges(6,5) * t176 - Icges(6,6) * t174;
t219 = t145 * t171;
t217 = t171 * t176;
t107 = t118 * t217;
t108 = t119 * t217;
t109 = t116 * t216;
t114 = -Icges(6,3) * t172 + t219;
t190 = t119 * t176 - t222;
t49 = -t114 * t171 - t118 * t215 + t109;
t115 = Icges(6,5) * t215 - Icges(6,6) * t216 + Icges(6,3) * t171;
t50 = t115 * t171 + t172 * t190;
t11 = (t49 + t108 - t109 + (t114 - t222) * t171) * t171 + (-t107 - t50 + (t114 + t190) * t172 + (t221 + t223) * t171) * t172;
t47 = -t114 * t172 - t116 * t218 + t107;
t48 = t115 * t172 + t117 * t218 - t108;
t12 = (t109 - t48 + (t115 - t221) * t172) * t172 + (-t107 + t47 + (t115 + t223) * t171) * t171;
t23 = -t171 * t48 - t172 * t47;
t24 = -t171 * t50 - t172 * t49;
t2 = (-t12 / 0.2e1 - t24 / 0.2e1) * t172 + (t23 / 0.2e1 - t11 / 0.2e1) * t171;
t201 = t63 * qJD(4) + t2 * qJD(5);
t194 = t252 / 0.2e1 + t198;
t192 = Icges(6,5) * t174 + Icges(6,6) * t176;
t187 = t171 * t11 / 0.2e1 + (-t145 * t172 - t267 * t171 - t210 * t174 + t208 * t176) * t241 + (t12 + t24) * t240 + (t267 * t172 + t174 * t209 - t176 * t207 - t219 + t23) * t242;
t186 = -t198 + (t240 + t241) * (t176 * t117 + t174 * t119);
t151 = -rSges(6,2) * t174 + t233;
t126 = t192 * t172;
t125 = t171 * t192;
t62 = t63 * qJD(5);
t60 = t64 * qJD(5);
t32 = t243 + t253;
t30 = t198 + t245;
t27 = t198 + t246;
t22 = t244 + t254;
t18 = t251 / 0.2e1;
t13 = t248 + t256 + t260;
t8 = -t251 / 0.2e1 + t194;
t7 = t18 + t194;
t6 = t238 + t239;
t3 = t18 - t252 / 0.2e1 + t186;
t1 = [t13 * qJD(3) + t22 * qJD(4) + t27 * qJD(5), 0, t13 * qJD(1) + t6 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t260 / 0.2e1 + t256 / 0.2e1 + t248 / 0.2e1) * qJD(3), qJD(1) * t22 + qJD(3) * t6 + t60, t27 * qJD(1) + t7 * qJD(3) + (-t151 * t244 + t187) * qJD(5) + t236; 0, 0, 0, 0, qJD(5) * t189; -t4 * qJD(4) + t8 * qJD(5) + (-t260 / 0.4e1 - t256 / 0.4e1 - t248 / 0.4e1) * t262, 0, qJD(4) * t32 + qJD(5) * t30, qJD(3) * t32 - t279 + t60, t8 * qJD(1) + t30 * qJD(3) + (-t151 * t243 + t187) * qJD(5) + t236; (-t254 / 0.4e1 - t244 / 0.4e1) * t262 + t4 * qJD(3) - t62, 0, t279 - t62 + 0.4e1 * (-t243 / 0.4e1 - t253 / 0.4e1) * qJD(3), 0, -t272; (t186 - t246) * qJD(1) + t3 * qJD(3) + t201, 0, t3 * qJD(1) + (t186 - t245) * qJD(3) + t201, t272, (m(6) * (t151 * t269 + (t172 * (rSges(6,3) * t171 + t197) + t171 * (rSges(6,1) * t217 + t204)) * t97) + (-t263 * t125 + (t266 * t171 + (t126 - t265) * t172) * t171) * t241 + (t264 * t126 + (t265 * t172 + (-t125 - t266) * t171) * t172) * t242) * qJD(5) + t202 * t2;];
Cq = t1;
