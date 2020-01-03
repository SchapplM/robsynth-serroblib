% Calculate kinetic energy for
% S5RRPPR2
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:20
% EndTime: 2020-01-03 11:57:20
% DurationCPUTime: 0.41s
% Computational Cost: add. (611->107), mult. (528->182), div. (0->0), fcn. (498->10), ass. (0->57)
t249 = sin(pkin(9));
t248 = qJ(1) + qJ(2);
t243 = pkin(8) + t248;
t239 = sin(t243);
t240 = cos(t243);
t253 = cos(qJ(5));
t250 = cos(pkin(9));
t251 = sin(qJ(5));
t266 = t250 * t251;
t224 = -t239 * t266 - t240 * t253;
t265 = t250 * t253;
t225 = t239 * t265 - t240 * t251;
t226 = -t239 * t253 + t240 * t266;
t227 = -t239 * t251 - t240 * t265;
t267 = t240 * t249;
t268 = t239 * t249;
t257 = (Icges(6,5) * t225 + Icges(6,6) * t224 + Icges(6,3) * t268) * t239 - (Icges(6,5) * t227 + Icges(6,6) * t226 - Icges(6,3) * t267) * t240;
t275 = t257 * t249;
t247 = qJD(1) + qJD(2);
t262 = qJD(5) * t249;
t274 = t247 * (pkin(4) * t250 + pkin(7) * t249) - (-rSges(6,3) * t250 + (rSges(6,1) * t253 - rSges(6,2) * t251) * t249) * t262 - qJD(4);
t273 = -qJD(4) + t247 * (rSges(5,1) * t250 - rSges(5,2) * t249);
t270 = pkin(2) * t247;
t269 = pkin(1) * qJD(1);
t252 = sin(qJ(1));
t241 = t252 * t269;
t244 = sin(t248);
t264 = t244 * t270 + t241;
t254 = cos(qJ(1));
t242 = t254 * t269;
t245 = cos(t248);
t263 = t245 * t270 + t242;
t261 = t247 * (pkin(3) * t239 - qJ(4) * t240) + t264;
t255 = qJD(3) ^ 2;
t238 = -t250 * qJD(5) + t247;
t237 = -rSges(2,1) * t254 + rSges(2,2) * t252;
t236 = rSges(2,1) * t252 + rSges(2,2) * t254;
t233 = -pkin(3) * t240 - qJ(4) * t239;
t230 = -Icges(6,5) * t250 + (Icges(6,1) * t253 - Icges(6,4) * t251) * t249;
t229 = -Icges(6,6) * t250 + (Icges(6,4) * t253 - Icges(6,2) * t251) * t249;
t228 = -Icges(6,3) * t250 + (Icges(6,5) * t253 - Icges(6,6) * t251) * t249;
t223 = t242 - t247 * (-rSges(3,1) * t245 + rSges(3,2) * t244);
t222 = t241 + t247 * (rSges(3,1) * t244 + rSges(3,2) * t245);
t221 = -t247 * (-rSges(4,1) * t240 + rSges(4,2) * t239) + t263;
t220 = t247 * (rSges(4,1) * t239 + rSges(4,2) * t240) + t264;
t219 = rSges(6,1) * t227 + rSges(6,2) * t226 - rSges(6,3) * t267;
t218 = rSges(6,1) * t225 + rSges(6,2) * t224 + rSges(6,3) * t268;
t217 = Icges(6,1) * t227 + Icges(6,4) * t226 - Icges(6,5) * t267;
t216 = Icges(6,1) * t225 + Icges(6,4) * t224 + Icges(6,5) * t268;
t215 = Icges(6,4) * t227 + Icges(6,2) * t226 - Icges(6,6) * t267;
t214 = Icges(6,4) * t225 + Icges(6,2) * t224 + Icges(6,6) * t268;
t211 = (rSges(5,3) * t239 - t233) * t247 + t273 * t240 + t263;
t210 = -t247 * t240 * rSges(5,3) + t273 * t239 + t261;
t209 = qJD(3) + (t218 * t240 + t219 * t239) * t262;
t208 = -t238 * t219 - t247 * t233 + t274 * t240 + t263;
t207 = t238 * t218 + t274 * t239 + t261;
t1 = m(3) * (t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(4) * (t220 ^ 2 + t221 ^ 2 + t255) / 0.2e1 + m(5) * (t210 ^ 2 + t211 ^ 2 + t255) / 0.2e1 + m(6) * (t207 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + t238 * ((-t250 * t228 + (-t229 * t251 + t230 * t253) * t249) * t238 + (((-t214 * t251 + t216 * t253) * t239 - (-t215 * t251 + t217 * t253) * t240) * t249 - t257 * t250) * t262) / 0.2e1 + t239 * ((t224 * t229 + t225 * t230 + t228 * t268) * t238 + (-(t215 * t224 + t217 * t225) * t240 + (t224 * t214 + t225 * t216 + t275) * t239) * t262) * t262 / 0.2e1 - t240 * ((t226 * t229 + t227 * t230 - t228 * t267) * t238 + ((t214 * t226 + t216 * t227) * t239 + (-t226 * t215 - t227 * t217 - t275) * t240) * t262) * t262 / 0.2e1 + (m(2) * (t236 ^ 2 + t237 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (Icges(3,3) + Icges(4,3) + Icges(5,2) * t250 ^ 2 + (Icges(5,1) * t249 + 0.2e1 * Icges(5,4) * t250) * t249) * t247 ^ 2 / 0.2e1;
T = t1;
