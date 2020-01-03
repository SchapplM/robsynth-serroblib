% Calculate kinetic energy for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:28
% EndTime: 2019-12-31 16:52:28
% DurationCPUTime: 0.53s
% Computational Cost: add. (510->122), mult. (528->208), div. (0->0), fcn. (440->8), ass. (0->75)
t239 = cos(pkin(7));
t271 = t239 * pkin(2);
t237 = pkin(7) + qJ(3);
t231 = sin(t237);
t270 = Icges(4,4) * t231;
t232 = cos(t237);
t269 = Icges(4,4) * t232;
t233 = qJ(4) + t237;
t228 = sin(t233);
t268 = Icges(5,4) * t228;
t229 = cos(t233);
t267 = Icges(5,4) * t229;
t241 = sin(qJ(1));
t242 = cos(qJ(1));
t225 = pkin(1) * t241 - qJ(2) * t242;
t265 = pkin(5) * t242 - t271 * t241 - t225;
t264 = pkin(3) * t232;
t262 = qJD(3) * t241;
t261 = qJD(3) * t242;
t260 = qJD(3) + qJD(4);
t259 = pkin(3) * qJD(3) * t231;
t222 = qJD(1) * (pkin(1) * t242 + qJ(2) * t241);
t258 = -qJD(2) * t242 + qJD(1) * (pkin(5) * t241 + t271 * t242) + t222;
t238 = sin(pkin(7));
t257 = rSges(3,1) * t239 - rSges(3,2) * t238;
t256 = rSges(4,1) * t232 - rSges(4,2) * t231;
t255 = rSges(5,1) * t229 - rSges(5,2) * t228;
t254 = Icges(4,1) * t232 - t270;
t253 = Icges(5,1) * t229 - t268;
t252 = -Icges(4,2) * t231 + t269;
t251 = -Icges(5,2) * t228 + t267;
t250 = Icges(4,5) * t232 - Icges(4,6) * t231;
t249 = Icges(5,5) * t229 - Icges(5,6) * t228;
t207 = -Icges(4,6) * t242 + t252 * t241;
t209 = -Icges(4,5) * t242 + t254 * t241;
t248 = t207 * t231 - t209 * t232;
t208 = Icges(4,6) * t241 + t252 * t242;
t210 = Icges(4,5) * t241 + t254 * t242;
t247 = -t208 * t231 + t210 * t232;
t218 = Icges(4,2) * t232 + t270;
t219 = Icges(4,1) * t231 + t269;
t246 = -t218 * t231 + t219 * t232;
t223 = t260 * t241;
t224 = t260 * t242;
t245 = (Icges(5,5) * t228 + Icges(5,6) * t229) * qJD(1) - (-Icges(5,3) * t242 + t249 * t241) * t224 + (Icges(5,3) * t241 + t249 * t242) * t223;
t198 = -Icges(5,6) * t242 + t251 * t241;
t199 = Icges(5,6) * t241 + t251 * t242;
t200 = -Icges(5,5) * t242 + t253 * t241;
t201 = Icges(5,5) * t241 + t253 * t242;
t214 = Icges(5,2) * t229 + t268;
t215 = Icges(5,1) * t228 + t267;
t244 = (-t199 * t228 + t201 * t229) * t223 - (-t198 * t228 + t200 * t229) * t224 + (-t214 * t228 + t215 * t229) * qJD(1);
t234 = qJD(2) * t241;
t227 = rSges(2,1) * t242 - rSges(2,2) * t241;
t226 = rSges(2,1) * t241 + rSges(2,2) * t242;
t220 = rSges(4,1) * t231 + rSges(4,2) * t232;
t217 = Icges(4,5) * t231 + Icges(4,6) * t232;
t216 = rSges(5,1) * t228 + rSges(5,2) * t229;
t212 = rSges(4,3) * t241 + t256 * t242;
t211 = -rSges(4,3) * t242 + t256 * t241;
t206 = Icges(4,3) * t241 + t250 * t242;
t205 = -Icges(4,3) * t242 + t250 * t241;
t203 = rSges(5,3) * t241 + t255 * t242;
t202 = -rSges(5,3) * t242 + t255 * t241;
t194 = qJD(1) * t241 * rSges(3,3) + t222 + (qJD(1) * t257 - qJD(2)) * t242;
t193 = t234 + (t242 * rSges(3,3) - t257 * t241 - t225) * qJD(1);
t192 = pkin(6) * t241 + t264 * t242;
t191 = -pkin(6) * t242 + t264 * t241;
t190 = (t211 * t241 + t212 * t242) * qJD(3);
t189 = qJD(1) * t212 - t220 * t262 + t258;
t188 = -t220 * t261 + t234 + (-t211 + t265) * qJD(1);
t187 = -t241 * t259 - t216 * t223 + (t192 + t203) * qJD(1) + t258;
t186 = -t242 * t259 - t216 * t224 + t234 + (-t191 - t202 + t265) * qJD(1);
t185 = t202 * t223 + t203 * t224 + (t191 * t241 + t192 * t242) * qJD(3);
t1 = m(3) * (t193 ^ 2 + t194 ^ 2) / 0.2e1 + m(4) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + ((t241 * t217 + t246 * t242) * qJD(1) + (t241 ^ 2 * t206 + (t248 * t242 + (-t205 + t247) * t241) * t242) * qJD(3)) * t262 / 0.2e1 - ((-t242 * t217 + t246 * t241) * qJD(1) + (t242 ^ 2 * t205 + (t247 * t241 + (-t206 + t248) * t242) * t241) * qJD(3)) * t261 / 0.2e1 + m(5) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + t223 * (t245 * t241 + t244 * t242) / 0.2e1 - t224 * (t244 * t241 - t245 * t242) / 0.2e1 + (((t208 * t232 + t210 * t231) * t241 - (t207 * t232 + t231 * t209) * t242) * qJD(3) + (t199 * t229 + t201 * t228) * t223 - (t198 * t229 + t200 * t228) * t224 + (t214 * t229 + t215 * t228 + t218 * t232 + t219 * t231) * qJD(1)) * qJD(1) / 0.2e1 + (m(2) * (t226 ^ 2 + t227 ^ 2) + Icges(2,3) + Icges(3,2) * t239 ^ 2 + (Icges(3,1) * t238 + 0.2e1 * Icges(3,4) * t239) * t238) * qJD(1) ^ 2 / 0.2e1;
T = t1;
