% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR6
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:45
% DurationCPUTime: 1.51s
% Computational Cost: add. (11394->175), mult. (7097->246), div. (0->0), fcn. (6110->8), ass. (0->121)
t269 = m(6) / 0.2e1;
t270 = m(5) / 0.2e1;
t179 = qJ(1) + pkin(8);
t178 = qJ(3) + t179;
t174 = sin(t178);
t175 = cos(t178);
t163 = t175 * qJ(4);
t180 = sin(qJ(5));
t182 = cos(qJ(5));
t195 = rSges(6,1) * t180 + rSges(6,2) * t182;
t184 = -t174 * rSges(6,3) + t195 * t175;
t253 = -pkin(3) - pkin(7);
t101 = t253 * t174 + t163 + t184;
t199 = -sin(qJ(1)) * pkin(1) - pkin(2) * sin(t179);
t93 = t101 + t199;
t102 = (rSges(6,3) - t253) * t175 + (qJ(4) + t195) * t174;
t198 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t179);
t94 = t102 + t198;
t42 = t94 * t174 + t93 * t175;
t49 = t101 * t175 + t102 * t174;
t267 = pkin(3) - rSges(5,2);
t122 = t175 * rSges(5,3) - t174 * t267 + t163;
t107 = t122 + t199;
t123 = (rSges(5,3) + qJ(4)) * t174 + t267 * t175;
t108 = t123 + t198;
t60 = t107 * t175 + t108 * t174;
t73 = t122 * t175 + t123 * t174;
t229 = (t49 + t42) * t269 + (t73 + t60) * t270;
t230 = (t42 - t49) * t269 + (t60 - t73) * t270;
t4 = t230 - t229;
t272 = t4 * qJD(1);
t224 = Icges(6,4) * t182;
t156 = -Icges(6,2) * t180 + t224;
t194 = Icges(6,1) * t180 + t224;
t271 = t156 + t194;
t248 = m(5) * (-t123 * t107 + t108 * t122);
t240 = m(6) * (t94 * t101 - t102 * t93);
t191 = Icges(6,5) * t180 + Icges(6,6) * t182;
t266 = t191 * t174;
t265 = t191 * t175;
t252 = m(4) * (t198 * (-rSges(4,1) * t174 - rSges(4,2) * t175) - (t175 * rSges(4,1) - t174 * rSges(4,2)) * t199);
t161 = rSges(6,1) * t182 - rSges(6,2) * t180;
t140 = t161 * t175;
t139 = t161 * t174;
t222 = t102 * t139;
t40 = t101 * t140 + t222;
t264 = m(6) * t40;
t263 = m(6) * t195;
t225 = Icges(6,4) * t180;
t193 = Icges(6,2) * t182 + t225;
t127 = -Icges(6,6) * t174 + t193 * t175;
t129 = -Icges(6,5) * t174 + t194 * t175;
t262 = (t127 * t182 + t129 * t180) * t175;
t158 = Icges(6,1) * t182 - t225;
t261 = t271 * t180 - (-t193 + t158) * t182;
t207 = t156 * t175 + t129;
t209 = -t158 * t175 + t127;
t260 = t209 * t180 - t207 * t182;
t216 = t174 * t182;
t153 = Icges(6,4) * t216;
t217 = t174 * t180;
t128 = Icges(6,1) * t217 + Icges(6,5) * t175 + t153;
t208 = -Icges(6,2) * t217 + t128 + t153;
t126 = Icges(6,6) * t175 + t193 * t174;
t210 = -t158 * t174 + t126;
t259 = t210 * t180 - t208 * t182;
t200 = -t271 * t182 / 0.2e1 + (t193 / 0.2e1 - t158 / 0.2e1) * t180;
t172 = t174 ^ 2;
t173 = t175 ^ 2;
t258 = 0.4e1 * qJD(1);
t246 = m(5) * t60;
t245 = m(5) * t73;
t226 = t94 * t139;
t38 = t93 * t140 + t226;
t244 = m(6) * (t40 + t38);
t243 = m(6) * (t226 - t222 + (-t101 + t93) * t140);
t238 = m(6) * t38;
t237 = m(6) * t42;
t235 = m(6) * t49;
t98 = -t139 * t175 + t140 * t174;
t234 = m(6) * t98;
t233 = t174 / 0.2e1;
t232 = -t175 / 0.2e1;
t231 = t175 / 0.2e1;
t203 = t234 / 0.2e1;
t228 = qJD(4) * t203;
t124 = Icges(6,3) * t175 + t266;
t117 = t174 * t124;
t125 = -Icges(6,3) * t174 + t265;
t190 = -t126 * t182 - t128 * t180;
t52 = t175 * t124 + t126 * t216 + t128 * t217;
t53 = -t175 * t125 - t127 * t216 - t129 * t217;
t54 = t190 * t175 + t117;
t55 = -t174 * t125 + t262;
t11 = (-t54 + t117 + t53) * t174 + (t55 - t262 + (t125 + t190) * t174 + t52) * t175;
t12 = t172 * t125 + (-t117 + t53 + (t125 - t190) * t175) * t175;
t27 = t174 * t53 + t175 * t52;
t28 = t174 * t55 + t175 * t54;
t2 = (t12 / 0.2e1 + t28 / 0.2e1) * t175 + (-t27 / 0.2e1 + t11 / 0.2e1) * t174;
t204 = -qJD(4) * t234 / 0.2e1 + t2 * qJD(5);
t202 = qJD(1) / 0.4e1 + qJD(3) / 0.4e1;
t201 = (t172 + t173) * t195;
t196 = t244 / 0.2e1 + t200;
t192 = Icges(6,5) * t182 - Icges(6,6) * t180;
t133 = t174 * t192;
t186 = -t174 * t11 / 0.2e1 + (t12 + t28) * t232 + (-t261 * t174 - t208 * t180 - t210 * t182 - t265) * t231 + (t261 * t175 + t207 * t180 + t209 * t182 - t266 + t27) * t233;
t185 = -t200 + (t231 + t232) * (t180 * t127 - t182 * t129);
t134 = t192 * t175;
t99 = -t139 * t174 - t140 * t175;
t88 = qJD(5) * t203;
t33 = t235 + t245;
t32 = t200 + t264;
t30 = t200 + t238;
t26 = t237 + t246;
t18 = t243 / 0.2e1;
t13 = t240 + t248 + t252;
t8 = -t243 / 0.2e1 + t196;
t7 = t18 + t196;
t5 = t229 + t230;
t3 = t18 - t244 / 0.2e1 + t185;
t1 = [qJD(3) * t13 + qJD(4) * t26 + qJD(5) * t30, 0, t13 * qJD(1) + t5 * qJD(4) + t7 * qJD(5) + 0.2e1 * (t248 / 0.2e1 + t240 / 0.2e1 + t252 / 0.2e1) * qJD(3), qJD(1) * t26 + qJD(3) * t5 + t88, t30 * qJD(1) + t7 * qJD(3) + (-(t174 * t93 - t175 * t94) * t263 + t186) * qJD(5) + t228; 0, 0, 0, 0, m(6) * t99 * qJD(5); -t4 * qJD(4) + t8 * qJD(5) + (-t248 / 0.4e1 - t240 / 0.4e1 - t252 / 0.4e1) * t258, 0, qJD(4) * t33 + qJD(5) * t32, qJD(3) * t33 - t272 + t88, t8 * qJD(1) + t32 * qJD(3) + (-(t101 * t174 - t102 * t175) * t263 + t186) * qJD(5) + t228; t4 * qJD(3) + t88 + (-t246 / 0.4e1 - t237 / 0.4e1) * t258, 0, t272 + t88 + 0.4e1 * (-t245 / 0.4e1 - t235 / 0.4e1) * qJD(3), 0, 0.2e1 * (t202 * t98 - qJD(5) * t201 / 0.2e1) * m(6); (t185 - t238) * qJD(1) + t3 * qJD(3) + t204, 0, t3 * qJD(1) + (t185 - t264) * qJD(3) + t204, -0.2e1 * t202 * t234, (m(6) * (-t161 * t201 + (-t175 * t184 + (-t175 * rSges(6,3) - t195 * t174) * t174) * t99) + (t173 * t133 + (t260 * t174 + (-t134 - t259) * t175) * t174) * t231 + (-t172 * t134 + (t259 * t175 + (t133 - t260) * t174) * t175) * t233) * qJD(5) + (qJD(1) + qJD(3)) * t2;];
Cq = t1;
