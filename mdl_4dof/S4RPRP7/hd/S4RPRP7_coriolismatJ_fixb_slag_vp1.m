% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP7_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:03
% EndTime: 2019-12-31 16:47:07
% DurationCPUTime: 2.51s
% Computational Cost: add. (2486->214), mult. (5769->306), div. (0->0), fcn. (5199->4), ass. (0->141)
t163 = sin(qJ(1));
t164 = cos(qJ(3));
t165 = cos(qJ(1));
t202 = t164 * t165;
t162 = sin(qJ(3));
t204 = t162 * t165;
t86 = Icges(5,5) * t204 - Icges(5,6) * t163 - Icges(5,3) * t202;
t209 = Icges(4,4) * t162;
t175 = Icges(4,2) * t164 + t209;
t92 = -Icges(4,6) * t163 + t175 * t165;
t273 = t86 - t92;
t142 = Icges(5,5) * t202;
t94 = Icges(5,1) * t204 - Icges(5,4) * t163 - t142;
t208 = Icges(4,4) * t164;
t176 = Icges(4,1) * t162 + t208;
t96 = -Icges(4,5) * t163 + t176 * t165;
t272 = -t94 - t96;
t228 = rSges(5,1) + pkin(3);
t262 = rSges(5,3) + qJ(4);
t266 = t262 * t162 + t228 * t164;
t269 = t266 * t165;
t271 = t163 * t269;
t171 = Icges(4,5) * t162 + Icges(4,6) * t164;
t88 = -Icges(4,3) * t163 + t171 * t165;
t90 = Icges(5,4) * t204 - Icges(5,2) * t163 - Icges(5,6) * t202;
t270 = -t88 - t90;
t195 = t228 * t162;
t268 = (Icges(5,4) + Icges(4,5)) * t164 + (-Icges(4,6) + Icges(5,6)) * t162;
t267 = (t272 * t162 + t273 * t164) * t165;
t160 = t163 ^ 2;
t265 = m(5) / 0.2e1;
t203 = t163 * t164;
t205 = t162 * t163;
t264 = t270 * t165 + t273 * t203 + t272 * t205;
t223 = m(5) * qJD(4);
t161 = t165 ^ 2;
t198 = t160 + t161;
t77 = (0.1e1 - t198) * t164 * t162;
t263 = t77 * t223;
t173 = Icges(5,4) * t162 - Icges(5,6) * t164;
t154 = Icges(5,5) * t164;
t210 = Icges(5,1) * t162;
t93 = Icges(5,4) * t165 + (-t154 + t210) * t163;
t226 = t165 * (Icges(5,2) * t165 + t173 * t163) + t93 * t205;
t261 = t270 * t163 + t226 - t267;
t91 = Icges(4,6) * t165 + t175 * t163;
t143 = Icges(4,4) * t203;
t95 = Icges(4,1) * t205 + Icges(4,5) * t165 + t143;
t179 = -t162 * t95 - t164 * t91;
t259 = t179 * t165 - t93 * t204;
t257 = t163 / 0.2e1;
t255 = t165 / 0.2e1;
t149 = t165 * qJ(2);
t244 = pkin(1) + pkin(5);
t182 = -t244 * t163 + t149;
t196 = -t163 * rSges(5,2) - t262 * t202;
t50 = t165 * t195 + t182 + t196;
t199 = t262 * t203;
t51 = (rSges(5,2) + t244) * t165 + (qJ(2) + t195) * t163 - t199;
t253 = t51 * t163 + t165 * t50;
t67 = t266 * t163;
t252 = -t165 * t67 + t271;
t249 = t268 * t163;
t248 = t268 * t165;
t134 = rSges(4,1) * t164 - rSges(4,2) * t162;
t113 = t134 * t163;
t115 = t134 * t165;
t64 = t228 * t203 + t262 * t205;
t43 = -t165 * t64 + t271;
t227 = t43 * t265 + m(4) * (-t113 * t165 + t163 * t115) / 0.2e1;
t123 = -Icges(4,2) * t162 + t208;
t207 = Icges(5,5) * t162;
t125 = Icges(5,1) * t164 + t207;
t127 = Icges(4,1) * t164 - t209;
t206 = Icges(5,3) * t162;
t247 = (t127 / 0.2e1 - t175 / 0.2e1 + t125 / 0.2e1 - Icges(5,3) * t164 / 0.2e1 + t207 / 0.2e1) * t162 - (-t176 / 0.2e1 - t123 / 0.2e1 + t154 - t210 / 0.2e1 + t206 / 0.2e1) * t164;
t245 = 4 * qJD(1);
t243 = m(3) * (t165 * (rSges(3,3) * t165 + t149) + (rSges(3,3) + qJ(2)) * t160);
t181 = rSges(4,1) * t162 + rSges(4,2) * t164;
t166 = -t163 * rSges(4,3) + t181 * t165;
t61 = t166 + t182;
t62 = (rSges(4,3) + t244) * t165 + (qJ(2) + t181) * t163;
t242 = m(4) * (t113 * t62 + t115 * t61);
t241 = m(4) * (t62 * t163 + t165 * t61);
t183 = -t51 * t204 + t50 * t205;
t238 = m(5) * (-t43 * t164 + t183);
t237 = m(5) * (t252 * t164 + t183);
t236 = m(5) * (t269 * t50 + t51 * t64);
t234 = m(5) * t253;
t232 = m(5) * t252;
t224 = m(5) * qJD(3);
t222 = t162 * t93;
t218 = (t154 + t206) * t163 - t93;
t217 = -Icges(5,3) * t204 - t142 + t94;
t216 = -Icges(4,2) * t205 + t143 + t95;
t215 = -t123 * t165 - t96;
t141 = Icges(5,5) * t205;
t85 = Icges(5,6) * t165 - Icges(5,3) * t203 + t141;
t214 = Icges(5,1) * t203 + t141 + t85;
t213 = -t125 * t165 - t86;
t212 = t127 * t163 - t91;
t211 = -t127 * t165 + t92;
t201 = t262 * t164 - t195;
t34 = t165 * (Icges(4,3) * t165 + t171 * t163) + t91 * t203 + t95 * t205;
t27 = t253 * t164;
t197 = m(5) * t27 * qJD(1);
t194 = -t171 / 0.2e1 - t173 / 0.2e1;
t193 = -t164 * t85 + t90;
t192 = t218 * t165;
t191 = t217 * t163;
t190 = t216 * t165;
t189 = t215 * t163;
t188 = t214 * t165;
t187 = t213 * t163;
t186 = t212 * t165;
t185 = t211 * t163;
t184 = t198 * t181;
t66 = t201 * t163;
t68 = t201 * t165;
t177 = t66 * t163 + t165 * t68;
t32 = -t85 * t203 + t226;
t170 = -t264 * t163 / 0.2e1 + ((t193 - t222) * t165 - t259 + t264) * t257 - (t34 + t32) * t165 / 0.2e1 + ((t179 + t88) * t163 + t34 + t261 + t267) * t255;
t169 = t160 * t88 / 0.2e1 + (t193 * t163 + t261 - t32) * t257 + ((-t179 + t222 - t270) * t165 + t259 + t264) * t255;
t116 = t198 * t162;
t60 = t67 * t205;
t44 = -t232 / 0.2e1;
t40 = -t161 * t266 - t64 * t163;
t31 = (-t228 * t204 - t196) * t165 + (-rSges(5,2) * t165 - t228 * t205 + t199) * t163;
t19 = t164 * t198 * t31 + t204 * t269 + t60;
t17 = t237 / 0.2e1;
t15 = t238 / 0.2e1;
t14 = t234 + t241 + t243;
t13 = t232 / 0.2e1 + t227;
t12 = t44 + t227;
t11 = t44 - t227;
t9 = t236 + t242 - t247;
t4 = t17 - t238 / 0.2e1;
t3 = t17 + t15;
t2 = t15 - t237 / 0.2e1;
t1 = t170 * t163 + t169 * t165;
t5 = [t14 * qJD(2) + t9 * qJD(3) - t27 * t223, qJD(1) * t14 + qJD(3) * t12, t9 * qJD(1) + t12 * qJD(2) + t3 * qJD(4) + (m(5) * (t50 * t66 - t51 * t68 + (-t64 + t67) * t269) + (m(4) * (-t113 * t134 + t181 * t62) + t194 * t165 - t169) * t165 + (m(4) * (t115 * t134 - t181 * t61) + t194 * t163 - t170) * t163 + (-t190 / 0.2e1 + t192 / 0.2e1 - t189 / 0.2e1 + t191 / 0.2e1) * t162 + (t186 / 0.2e1 + t188 / 0.2e1 + t185 / 0.2e1 + t187 / 0.2e1) * t164) * qJD(3), t3 * qJD(3) - t197; t13 * qJD(3) + (-t243 / 0.4e1 - t241 / 0.4e1 - t234 / 0.4e1) * t245, 0, t13 * qJD(1) + 0.2e1 * (t177 * t265 - m(4) * t184 / 0.2e1) * qJD(3) + t116 * t223, t116 * t224; t11 * qJD(2) + t1 * qJD(3) + t4 * qJD(4) + (-t236 / 0.4e1 - t242 / 0.4e1) * t245 + t247 * qJD(1), t11 * qJD(1), t1 * qJD(1) + (m(4) * ((-t165 * t166 + (-t165 * rSges(4,3) - t181 * t163) * t163) * (-t163 * t113 - t115 * t165) - t134 * t184) + m(5) * (t269 * t68 + t31 * t40 + t66 * t67) + (((-t189 - t190 + t191 + t192) * t164 + t249 * t163 + ((-t212 - t214) * t165 + (-t211 - t213) * t163) * t162) * t165 - t248 * t160) * t257 + ((-t248 * t165 + (t186 + t185 + t188 + t187) * t162 + ((t215 - t217) * t163 + (t216 - t218) * t165) * t164) * t163 + t249 * t161) * t255) * qJD(3) + t19 * t223, t4 * qJD(1) + t19 * t224 - t263; t2 * qJD(3) + t197, 0, t2 * qJD(1) + (t60 + (t165 * t269 + t40) * t162 + (-t177 + t31) * t164 - t19) * t224 + t263, t77 * t224;];
Cq = t5;
