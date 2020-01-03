% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:36
% DurationCPUTime: 2.15s
% Computational Cost: add. (4248->145), mult. (4363->206), div. (0->0), fcn. (3960->6), ass. (0->100)
t141 = qJ(1) + pkin(6);
t137 = sin(t141);
t138 = cos(t141);
t142 = sin(qJ(3));
t144 = cos(qJ(3));
t229 = rSges(5,2) * t144 + (rSges(5,1) + pkin(3)) * t142;
t233 = t229 * t138;
t234 = t229 * t137;
t243 = m(5) * (t137 * t234 + t138 * t233);
t35 = -t243 / 0.2e1;
t30 = t243 / 0.2e1;
t245 = -Icges(4,4) - Icges(5,4);
t244 = Icges(4,5) + Icges(5,5);
t241 = Icges(4,2) + Icges(5,2);
t239 = Icges(4,6) + Icges(5,6);
t240 = Icges(4,3) + Icges(5,3);
t180 = Icges(5,4) * t142;
t104 = Icges(5,1) * t144 - t180;
t181 = Icges(4,4) * t142;
t106 = Icges(4,1) * t144 - t181;
t235 = (t104 + t106) * t138 + t244 * t137;
t175 = t137 * t142;
t238 = t245 * t175;
t221 = t245 * t144;
t97 = Icges(5,5) * t144 - Icges(5,6) * t142;
t98 = Icges(4,5) * t144 - Icges(4,6) * t142;
t236 = (t97 + t98) * t138 + t240 * t137;
t174 = t137 * t144;
t228 = t240 * t138 - t174 * t244 + t239 * t175;
t231 = t239 * t138 + t245 * t174 + t241 * t175;
t227 = t238 + (Icges(4,1) + Icges(5,1)) * t174 - t244 * t138;
t232 = t235 * t174;
t178 = Icges(5,2) * t142;
t179 = Icges(4,2) * t142;
t230 = (-t178 - t179 - t221) * t138 + t239 * t137;
t226 = t236 * t138 - t232;
t225 = t231 * t142;
t172 = t138 * t144;
t218 = t236 * t137 + t235 * t172;
t224 = t228 * t137 - t227 * t172;
t198 = pkin(3) * t144;
t134 = pkin(2) + t198;
t197 = -qJ(4) - pkin(5);
t146 = -rSges(5,1) * t174 + rSges(5,2) * t175 - t137 * t134 + (rSges(5,3) - t197) * t138;
t201 = sin(qJ(1)) * pkin(1);
t44 = t146 - t201;
t189 = rSges(5,1) * t144;
t157 = t134 + t189;
t173 = t138 * t142;
t170 = -rSges(5,2) * t173 - t137 * t197;
t200 = cos(qJ(1)) * pkin(1);
t45 = t137 * rSges(5,3) + t157 * t138 + t170 + t200;
t223 = m(5) * (t45 * t137 + t138 * t44);
t222 = -t230 * t173 + t218;
t220 = -t142 * t244 - t239 * t144;
t219 = t230 * t142 + t228;
t217 = t231 * t173 - t230 * t175 - t224 - t226;
t216 = -t137 / 0.2e1;
t204 = t137 / 0.2e1;
t203 = -t138 / 0.2e1;
t215 = t138 / 0.2e1;
t214 = qJD(1) * t223;
t101 = Icges(4,2) * t144 + t181;
t182 = Icges(5,1) * t142;
t183 = Icges(4,1) * t142;
t133 = t138 * pkin(5);
t190 = rSges(4,1) * t144;
t164 = pkin(2) + t190;
t171 = rSges(4,2) * t175 + t138 * rSges(4,3);
t46 = -t164 * t137 + t133 + t171 - t201;
t122 = rSges(4,2) * t173;
t47 = t200 - t122 + t164 * t138 + (rSges(4,3) + pkin(5)) * t137;
t108 = rSges(4,1) * t142 + rSges(4,2) * t144;
t92 = t108 * t137;
t93 = t108 * t138;
t99 = Icges(5,2) * t144 + t180;
t6 = (t106 / 0.2e1 - t101 / 0.2e1 + t104 / 0.2e1 - t99 / 0.2e1) * t142 + (t183 / 0.2e1 - t179 / 0.2e1 + t182 / 0.2e1 - t178 / 0.2e1 - t221) * t144 + m(5) * (-t233 * t45 + t234 * t44) + m(4) * (t46 * t92 - t47 * t93);
t213 = t6 * qJD(1);
t135 = t137 ^ 2;
t136 = t138 ^ 2;
t169 = t135 + t136;
t212 = -t182 - t183 + t221;
t191 = m(5) * qJD(3);
t168 = qJD(3) * t137;
t167 = t220 * t216;
t166 = t220 * t215;
t165 = t98 / 0.2e1 + t97 / 0.2e1;
t162 = rSges(5,2) * t142 - t189 - t198;
t150 = (t219 * t138 - t218 + t222) * t215 + (t228 * t138 + (t227 * t144 + t225) * t137) * t203 + (t219 * t137 + t217 + t226) * t204;
t149 = t218 * t216 + t222 * t204 + ((-t225 + t236) * t138 + t217 + t224 - t232) * t203;
t111 = -rSges(4,2) * t142 + t190;
t79 = t162 * t138;
t77 = t162 * t137;
t43 = -t137 * t92 - t138 * t93;
t31 = t229 * t169;
t16 = 0.2e1 * t35;
t15 = t35 + t30;
t14 = 0.2e1 * t30;
t1 = t150 * t137 + t149 * t138;
t2 = [t6 * qJD(3) + qJD(4) * t223, 0, t213 + (t44 * t79 + t45 * t77) * t191 + t15 * qJD(4) + (m(4) * (t108 * t93 - t111 * t47) + t165 * t137 - t150) * t168 + ((m(4) * (-t108 * t92 - t111 * t46) + t165 * t138 - t149) * t138 + (((-t101 - t99) * t138 + t235) * t204 + (-t241 * t174 + t227 + t238) * t203) * t144 + ((t212 * t138 - t230) * t204 + (t212 * t137 + t231) * t203) * t142) * qJD(3), t15 * qJD(3) + t214; 0, 0, 0.2e1 * (m(4) * t43 / 0.2e1 - m(5) * t31 / 0.2e1) * qJD(3), 0; t1 * qJD(3) + t16 * qJD(4) - t213, 0, t1 * qJD(1) + (t167 * t137 + t166 * t138) * t168 * t138 + (m(4) * (t169 * t111 * t108 + (t137 * (rSges(4,1) * t174 - t171) + t138 * (rSges(4,1) * t172 + t137 * rSges(4,3) - t122)) * t43) + m(5) * (-((-pkin(2) * t137 + t133 - t146) * t137 + ((-pkin(2) + t157) * t138 + (rSges(5,3) - pkin(5)) * t137 + t170) * t138) * t31 - t234 * t77 - t233 * t79) + t166 * t137 * t135 + t167 * t136 * t138) * qJD(3), t16 * qJD(1); t14 * qJD(3) - t214, 0, t14 * qJD(1) + (t137 * t79 - t138 * t77) * t191, 0;];
Cq = t2;
