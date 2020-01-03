% Calculate vector of inverse dynamics joint torques for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:38
% EndTime: 2019-12-31 17:35:44
% DurationCPUTime: 3.95s
% Computational Cost: add. (8430->316), mult. (9097->428), div. (0->0), fcn. (9622->8), ass. (0->182)
t171 = sin(pkin(8));
t172 = cos(pkin(8));
t248 = qJ(3) + qJ(4);
t224 = sin(t248);
t225 = cos(t248);
t137 = -t171 * t225 + t172 * t224;
t170 = qJD(3) + qJD(4);
t133 = t137 * pkin(7);
t184 = t171 * t224 + t172 * t225;
t173 = sin(qJ(5));
t257 = t184 * t173;
t220 = rSges(6,2) * t257 - t137 * rSges(6,3);
t175 = cos(qJ(5));
t256 = t184 * t175;
t73 = -rSges(6,1) * t256 + t220;
t187 = -pkin(4) * t184 - t133 + t73;
t214 = rSges(6,1) * t173 + rSges(6,2) * t175;
t239 = qJD(5) * t214;
t299 = -t137 * t239 - t170 * t187;
t298 = qJD(3) ^ 2;
t264 = Icges(6,4) * t175;
t202 = -Icges(6,2) * t173 + t264;
t69 = Icges(6,6) * t137 + t184 * t202;
t265 = Icges(6,4) * t173;
t204 = Icges(6,1) * t175 - t265;
t72 = Icges(6,5) * t137 + t184 * t204;
t205 = t173 * t69 - t175 * t72;
t200 = Icges(6,5) * t175 - Icges(6,6) * t173;
t66 = Icges(6,3) * t137 + t184 * t200;
t27 = t137 * t66 - t184 * t205;
t115 = t184 * t170;
t112 = t115 * pkin(7);
t114 = t137 * t170;
t238 = qJD(5) * t173;
t191 = -t115 * t175 + t137 * t238;
t237 = qJD(5) * t175;
t192 = t115 * t173 + t137 * t237;
t46 = t191 * rSges(6,1) + rSges(6,2) * t192 - t114 * rSges(6,3);
t177 = -t115 * pkin(4) - t114 * pkin(7) + t46;
t230 = t184 * t238;
t290 = -t114 * t173 + t184 * t237;
t195 = rSges(6,1) * t230 + t290 * rSges(6,2) - t115 * rSges(6,3);
t253 = t137 * t175;
t254 = t137 * t173;
t74 = -rSges(6,1) * t253 + rSges(6,2) * t254 + rSges(6,3) * t184;
t270 = t137 * pkin(4) - pkin(7) * t184 - t74;
t196 = t270 * t170 + t184 * t239;
t268 = rSges(6,1) * t175;
t226 = pkin(4) + t268;
t174 = sin(qJ(3));
t249 = t172 * t174;
t234 = pkin(3) * t249;
t176 = cos(qJ(3));
t296 = pkin(3) * t176;
t288 = t296 * t171;
t120 = -t234 + t288;
t167 = qJD(2) * t171;
t247 = qJD(3) * t120 + t167;
t30 = -t196 + t247;
t251 = t171 * t174;
t121 = pkin(3) * t251 + t296 * t172;
t244 = qJD(2) * t172;
t197 = -qJD(3) * t121 - t244;
t31 = t197 - t299;
t297 = (t114 * t226 - t112 + t195 - t196) * t31 + (t177 + t299) * t30;
t65 = Icges(6,3) * t184 - t137 * t200;
t295 = t137 * t65;
t92 = -rSges(5,1) * t184 + t137 * rSges(5,2);
t293 = t170 * t92;
t201 = Icges(6,2) * t175 + t265;
t203 = Icges(6,1) * t173 + t264;
t198 = -t173 * t201 + t175 * t203;
t199 = Icges(6,5) * t173 + Icges(6,6) * t175;
t81 = t199 * t137;
t285 = t184 * t198 + t81;
t291 = t285 * t170;
t206 = t173 * t72 + t175 * t69;
t82 = t199 * t184;
t250 = t171 * t176;
t163 = pkin(3) * t250;
t141 = t163 - t234;
t93 = -t137 * rSges(5,1) - rSges(5,2) * t184;
t267 = t170 * t93;
t78 = -t115 * rSges(5,1) + t114 * rSges(5,2);
t284 = (t78 - t293) * (t247 + t267);
t71 = Icges(6,5) * t184 - t137 * t204;
t272 = t201 * t137 + t71;
t68 = Icges(6,6) * t184 - t137 * t202;
t274 = -t203 * t137 + t68;
t282 = t173 * t272 + t175 * t274;
t169 = -qJDD(3) - qJDD(4);
t279 = t169 / 0.2e1;
t278 = -t170 / 0.2e1;
t276 = -t184 * t65 + t71 * t253;
t275 = t71 * t256 + t295;
t273 = t184 * t203 + t69;
t271 = -t184 * t201 + t72;
t266 = t173 * t68;
t259 = t114 * t175;
t252 = t200 * t170;
t246 = -t201 + t204;
t245 = -t202 - t203;
t242 = qJD(5) * t184;
t241 = qJD(5) * t137;
t157 = rSges(6,2) * t173 - t268;
t147 = t157 * qJD(5);
t240 = qJD(5) * t147;
t236 = m(3) + m(4) + m(5);
t235 = qJDD(2) * t172;
t233 = -m(6) - t236;
t166 = qJDD(2) * t171;
t142 = -t172 * t176 - t251;
t190 = t142 * pkin(3);
t231 = qJDD(3) * t120 + t298 * t190 + t166;
t223 = -t242 / 0.2e1;
t222 = -t241 / 0.2e1;
t221 = t241 / 0.2e1;
t208 = t173 * t71 + t175 * t68;
t42 = Icges(6,4) * t191 + Icges(6,2) * t192 - Icges(6,6) * t114;
t44 = Icges(6,1) * t191 + Icges(6,4) * t192 - Icges(6,5) * t114;
t181 = qJD(5) * t208 + t173 * t42 - t175 * t44;
t193 = -t230 - t259;
t41 = Icges(6,4) * t193 - Icges(6,2) * t290 + Icges(6,6) * t115;
t43 = Icges(6,1) * t193 - Icges(6,4) * t290 + Icges(6,5) * t115;
t182 = qJD(5) * t206 + t173 * t41 - t175 * t43;
t207 = -t175 * t71 + t266;
t39 = Icges(6,5) * t193 - Icges(6,6) * t290 + Icges(6,3) * t115;
t40 = Icges(6,5) * t191 + Icges(6,6) * t192 - Icges(6,3) * t114;
t219 = t184 * (t114 * t207 + t115 * t65 + t137 * t40 - t181 * t184) + t137 * (t114 * t205 + t115 * t66 + t137 * t39 - t182 * t184);
t218 = t184 * (-t114 * t65 + t115 * t207 + t137 * t181 + t184 * t40) + t137 * (-t114 * t66 + t115 * t205 + t137 * t182 + t184 * t39);
t143 = t249 - t250;
t102 = rSges(4,1) * t142 + rSges(4,2) * t143;
t77 = -t114 * rSges(5,1) - t115 * rSges(5,2);
t24 = t254 * t68 - t276;
t25 = t184 * t66 - t72 * t253 + t254 * t69;
t213 = t137 * t25 + t184 * t24;
t26 = -t257 * t68 + t275;
t212 = t137 * t27 + t184 * t26;
t211 = -t137 * t31 + t184 * t30;
t45 = -rSges(6,1) * t259 - t195;
t210 = t137 * t46 - t184 * t45;
t209 = t137 * t74 + t184 * t73;
t108 = -t173 * t203 - t175 * t201;
t140 = t142 * qJD(3);
t189 = t173 * t271 + t175 * t273;
t188 = -qJDD(3) * t121 - t298 * t141 - t235;
t186 = -t184 * t226 - t133 + t220;
t185 = (t173 * t245 + t175 * t246) * t170;
t180 = -t270 - t234;
t145 = t202 * qJD(5);
t146 = t204 * qJD(5);
t179 = qJD(5) * t108 - t145 * t173 + t146 * t175;
t10 = qJD(5) * t212 + t291;
t17 = qJD(5) * t207 - t173 * t44 - t175 * t42;
t18 = qJD(5) * t205 - t173 * t43 - t175 * t41;
t144 = t200 * qJD(5);
t22 = t114 * t198 - t115 * t199 - t137 * t144 - t179 * t184;
t23 = t114 * t199 + t115 * t198 + t137 * t179 - t144 * t184;
t49 = t137 * t198 - t82;
t48 = t49 * t170;
t79 = qJD(5) * t115 + qJDD(5) * t137;
t80 = -qJD(5) * t114 + qJDD(5) * t184;
t9 = qJD(5) * t213 - t48;
t178 = -Icges(5,3) * t169 + 0.2e1 * t108 * t279 + 0.2e1 * (-qJD(5) * t198 - t145 * t175 - t146 * t173) * t278 + t9 * t221 - (-t206 - t285) * t79 / 0.2e1 - (-t208 + t49) * t80 / 0.2e1 + (-t48 + ((t25 - t26 + t275) * t137 - t276 * t184) * qJD(5) + t18 + t22) * t222 + (t17 + ((-t24 + (-t66 + t266) * t137 - t276) * t137 - (-t25 - (t207 - t66) * t184 + t295) * t184) * qJD(5) + t23 + t10 - t291) * t223;
t139 = t143 * qJD(3);
t103 = -rSges(4,1) * t143 + rSges(4,2) * t142;
t99 = rSges(4,1) * t140 + rSges(4,2) * t139;
t98 = -rSges(4,1) * t139 + rSges(4,2) * t140;
t91 = qJD(3) * t102 - t244;
t90 = qJD(3) * t103 + t167;
t88 = t214 * t184;
t87 = t214 * t137;
t60 = t197 + t293;
t58 = -qJD(3) * t98 + qJDD(3) * t102 - t235;
t57 = qJD(3) * t99 + qJDD(3) * t103 + t166;
t34 = -t169 * t92 - t170 * t77 + t188;
t33 = -t169 * t93 + t170 * t78 + t231;
t32 = qJD(5) * t209 + qJD(1);
t14 = -t137 * t240 + t79 * t214 + (t114 * pkin(4) - t112 - t45) * t170 - t187 * t169 + t188;
t13 = t169 * t270 + t170 * t177 + t184 * t240 - t214 * t80 + t231;
t12 = qJD(5) * t210 + t73 * t80 + t74 * t79 + qJDD(1);
t1 = [m(6) * t12 + (m(2) + t236) * qJDD(1) + (-m(2) + t233) * g(3); t233 * (g(1) * t171 - g(2) * t172) + m(4) * (t171 * t57 - t172 * t58) + m(5) * (t33 * t171 - t34 * t172) + m(6) * (t13 * t171 - t14 * t172) + m(3) * (t171 ^ 2 + t172 ^ 2) * qJDD(2); Icges(4,3) * qJDD(3) + t178 + (-g(1) * (t163 + t180) - g(2) * (t190 + t187) + t13 * (t180 + t288) + t14 * (-t121 + t186) + t297) * m(6) + (-g(1) * (t141 + t93) - g(2) * (t190 + t92) + t33 * (t120 + t93) + t34 * (-t121 + t92) + (-t77 + t267) * t60 + t284) * m(5) + (t90 * t99 - t91 * t98 + (qJD(3) * t91 - g(1) + t57) * t103 + (-qJD(3) * t90 - g(2) + t58) * t102) * m(4); t178 + (-g(2) * t187 + t14 * t186 - (-g(1) + t13) * t270 + t297) * m(6) + (-t60 * t77 + t284 + (t170 * t60 - g(1) + t33) * t93 + (-g(2) + t34) * t92) * m(5); t115 * t10 / 0.2e1 + t137 * (qJD(5) * t219 - t169 * t285 - t22 * t170 + t26 * t80 + t27 * t79) / 0.2e1 + t79 * t212 / 0.2e1 + (-t114 * t26 + t115 * t27 + t219) * t221 - t114 * t9 / 0.2e1 + t184 * (qJD(5) * t218 + t49 * t169 - t23 * t170 + t24 * t80 + t25 * t79) / 0.2e1 + t80 * t213 / 0.2e1 + (-t114 * t24 + t115 * t25 + t218) * t242 / 0.2e1 + (-t137 * t206 - t184 * t208) * t279 + (t114 * t208 - t115 * t206 + t18 * t137 + t17 * t184) * t278 + ((-t82 * t241 + t252) * t137 - (-t185 + (t282 * t184 + (-t81 + t189) * t137) * qJD(5)) * t184) * t222 + (-(-t81 * t242 - t252) * t184 + (-t185 + (t189 * t137 - (-t282 + t82) * t184) * qJD(5)) * t137) * t223 + t170 * (-(t173 * t246 - t175 * t245) * t170 + ((-t137 * t271 - t184 * t272) * t175 + (t137 * t273 + t184 * t274) * t173) * qJD(5)) / 0.2e1 + (t12 * t209 + t32 * (-t114 * t73 + t115 * t74 + t210) + t211 * t147 - (-t30 * t114 - t31 * t115 + t13 * t184 - t14 * t137) * t214 - (t30 * t87 + t31 * t88) * t170 - (t32 * (t137 * t87 + t184 * t88) + t211 * t157) * qJD(5) + g(1) * t88 - g(2) * t87 - g(3) * t157) * m(6);];
tau = t1;
