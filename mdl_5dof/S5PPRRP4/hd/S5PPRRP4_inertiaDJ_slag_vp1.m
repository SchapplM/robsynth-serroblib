% Calculate time derivative of joint inertia matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:33
% DurationCPUTime: 3.95s
% Computational Cost: add. (2094->187), mult. (5120->281), div. (0->0), fcn. (5228->6), ass. (0->114)
t98 = cos(qJ(4));
t140 = qJD(4) * t98;
t143 = sin(pkin(7));
t144 = cos(pkin(7));
t169 = sin(qJ(3));
t170 = cos(qJ(3));
t71 = -t143 * t169 - t144 * t170;
t67 = t71 * qJD(3);
t97 = sin(qJ(4));
t156 = t97 * t67;
t72 = -t143 * t170 + t144 * t169;
t107 = t140 * t72 - t156;
t236 = Icges(5,4) + Icges(6,4);
t232 = t236 * t97;
t234 = Icges(5,2) + Icges(6,2);
t193 = t234 * t98 + t232;
t238 = t107 * t193;
t233 = t236 * t98;
t235 = -Icges(5,1) - Icges(6,1);
t237 = t235 * t97 - t233;
t220 = Icges(5,5) + Icges(6,5);
t219 = Icges(5,6) + Icges(6,6);
t231 = t234 * t97 - t233;
t230 = t235 * t98 + t232;
t229 = -t219 * t71 + t231 * t72;
t228 = t219 * t72 + t231 * t71;
t196 = -t220 * t71 + t230 * t72;
t195 = t220 * t72 + t230 * t71;
t218 = Icges(5,3) + Icges(6,3);
t225 = t219 * t97 - t220 * t98;
t209 = -t195 * t98 + t228 * t97;
t208 = -t196 * t98 + t229 * t97;
t153 = rSges(6,3) + qJ(5) + pkin(6);
t212 = -t218 * t71 + t225 * t72;
t199 = t218 * t72 + t225 * t71;
t222 = t153 * t72;
t141 = qJD(4) * t97;
t138 = t71 * t141;
t66 = t72 * qJD(3);
t108 = -t66 * t98 + t138;
t137 = t71 * t140;
t109 = t66 * t97 + t137;
t183 = 2 * m(5);
t168 = rSges(5,2) * t97;
t128 = -rSges(5,1) * t98 + t168;
t81 = t128 * qJD(4);
t90 = -rSges(5,1) * t97 - rSges(5,2) * t98;
t217 = t81 * t90 * t183 + (t108 * t220 + t109 * t219 - t218 * t67) * t72;
t214 = t219 * t98 + t220 * t97;
t167 = rSges(6,2) * t97;
t127 = -rSges(6,1) * t98 + t167;
t95 = pkin(4) * t98 + pkin(3);
t110 = -t127 + t95;
t213 = pkin(3) - t110;
t211 = t209 * t72;
t210 = t208 * t71;
t139 = t72 * t141;
t163 = t67 * t98;
t106 = t139 + t163;
t205 = 0.2e1 * t72;
t166 = rSges(6,2) * t98;
t186 = (rSges(6,1) + pkin(4)) * t97 + t166;
t100 = qJD(4) * t186;
t70 = t72 * pkin(6);
t151 = t213 * t71 + t222 - t70;
t173 = pkin(6) * t71;
t160 = t72 * t98;
t161 = t72 * t97;
t29 = -rSges(6,1) * t160 + rSges(6,2) * t161 - t153 * t71 - t72 * t95;
t152 = pkin(3) * t72 + t173 + t29;
t172 = t66 * pkin(6);
t191 = t153 * t67;
t64 = t67 * pkin(6);
t99 = rSges(6,1) * t163 + t100 * t72 - t153 * t66 + t67 * t95;
t1 = -t152 * t67 + t151 * t66 + (t172 - (pkin(3) + t167) * t67 + t99) * t72 + (t71 * t100 + t213 * t66 - t191 + t64) * t71;
t182 = 2 * m(6);
t204 = t1 * t182;
t101 = rSges(5,1) * t106 - t66 * rSges(5,3);
t164 = t67 * rSges(5,3);
t150 = -rSges(5,1) * t160 - rSges(5,3) * t71;
t48 = rSges(5,2) * t161 + t150;
t162 = t72 * rSges(5,3);
t50 = t128 * t71 + t162;
t2 = -t67 * t48 + t72 * t101 + t66 * t50 + t71 * (rSges(5,1) * t108 - t164) + (t107 * t72 + t109 * t71) * rSges(5,2);
t203 = t183 * t2;
t201 = -t106 * t220 - t107 * t219 + t218 * t66;
t194 = t225 * qJD(4);
t190 = t90 ^ 2 * t183;
t188 = t237 * t98;
t187 = t199 * t71 - t212 * t72 - t210 - t211;
t185 = t193 * t97 + t188;
t184 = t185 * t66;
t177 = -t71 / 0.2e1;
t174 = t72 / 0.2e1;
t142 = qJD(4) * t71;
t89 = -rSges(6,1) * t97 - t166;
t135 = pkin(4) * t97 - t89;
t122 = t66 * t71 - t67 * t72;
t111 = pkin(3) - t128;
t80 = t127 * qJD(4);
t54 = rSges(4,1) * t67 + rSges(4,2) * t66;
t53 = rSges(4,1) * t66 - rSges(4,2) * t67;
t52 = t135 * t71;
t51 = t135 * t72;
t32 = t111 * t71 - t162 - t70;
t31 = -t173 + (-pkin(3) + t168) * t72 + t150;
t30 = t110 * t71 - t222;
t28 = pkin(4) * t107 + t67 * t89 - t72 * t80;
t27 = pkin(4) * t109 - t66 * t89 - t71 * t80;
t14 = rSges(5,2) * t107 + pkin(3) * t67 + t101 - t172;
t13 = t111 * t66 + t142 * t90 + t164 + t64;
t12 = -rSges(6,2) * t156 - qJD(5) * t71 + t99;
t11 = -qJD(5) * t72 + t110 * t66 - t142 * t186 + t191;
t3 = [0; 0; 0; 0; m(4) * (t143 * t54 - t144 * t53) + m(5) * (-t13 * t144 + t14 * t143) + m(6) * (-t11 * t144 + t12 * t143); 0.2e1 * m(4) * ((-rSges(4,1) * t72 + rSges(4,2) * t71) * t54 + (rSges(4,1) * t71 + rSges(4,2) * t72) * t53) + (t13 * t32 + t14 * t31) * t183 + (t11 * t30 + t12 * t29) * t182 + ((-t237 - t231) * t98 + (-t193 - t230) * t97) * qJD(4); m(5) * t2 + m(6) * t1; m(5) * ((-t143 * t66 - t144 * t67) * t90 + (-t143 * t71 + t144 * t72) * t81) + m(6) * (t143 * t27 - t144 * t28); m(6) * (t11 * t51 + t12 * t52 + t27 * t29 + t28 * t30) + (-t196 * t97 - t229 * t98) * t66 / 0.2e1 + (-t195 * t97 - t228 * t98) * t67 / 0.2e1 + (-t139 * t237 - t188 * t67 + t238) * t177 + (t208 * qJD(4) + t106 * t237 + t214 * t66 - t238) * t71 / 0.2e1 - (t209 * qJD(4) + t108 * t237 - t109 * t193 + t214 * t67) * t72 / 0.2e1 + (t137 * t193 - t138 * t237 + t184) * t174 + t214 * t122 + (-t184 / 0.2e1 - t194 * t174) * t72 + (-t185 * t67 / 0.2e1 + t194 * t177) * t71 + ((-t31 * t71 - t32 * t72) * t81 + (-t13 * t72 - t14 * t71 - t31 * t66 + t32 * t67) * t90) * m(5); (t52 * t27 + t51 * t28) * t182 + (t48 * t203 + t152 * t204 + t217 * t72 - (0.2e1 * t209 * t71 + t190 + (t205 + t72) * t199) * t67 + (t187 + t211) * t66) * t72 + (t50 * t203 + t151 * t204 + (t208 * t205 + t190) * t66 - (t187 + t210) * t67 + (t201 * t72 + t212 * t67 + t199 * t66 + (-t195 * t67 + t196 * t66) * t98 + (t228 * t67 - t229 * t66) * t97) * t72 + (t201 * t71 - 0.3e1 * t212 * t66 + t217) * t71) * t71; 0; m(6) * (-t143 * t67 + t144 * t66); m(6) * (-t11 * t71 + t12 * t72 - t29 * t67 - t30 * t66); m(6) * (t27 * t72 - t28 * t71 - t51 * t66 - t52 * t67); t122 * t182;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
