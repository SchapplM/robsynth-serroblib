% Calculate kinetic energy for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:28
% EndTime: 2019-12-31 17:56:28
% DurationCPUTime: 0.46s
% Computational Cost: add. (472->73), mult. (501->123), div. (0->0), fcn. (512->8), ass. (0->44)
t260 = qJ(1) + pkin(8);
t257 = sin(t260);
t258 = cos(t260);
t265 = sin(qJ(4));
t266 = cos(qJ(4));
t224 = -t257 * t265 - t258 * t266;
t273 = t224 ^ 2;
t225 = -t257 * t266 + t258 * t265;
t272 = t225 ^ 2;
t241 = sin(qJ(5));
t243 = cos(qJ(5));
t271 = -Icges(6,5) * t241 - Icges(6,6) * t243;
t269 = t224 * t225;
t240 = qJD(1) - qJD(4);
t268 = t271 * t240;
t242 = sin(qJ(1));
t264 = t242 * pkin(1);
t261 = qJD(5) * (-t241 * rSges(6,1) - t243 * rSges(6,2));
t259 = -pkin(2) * t257 + qJ(3) * t258 - t264;
t256 = -rSges(6,1) * t243 + rSges(6,2) * t241;
t253 = -Icges(6,5) * t243 + Icges(6,6) * t241;
t244 = cos(qJ(1));
t239 = qJD(1) * t244 * pkin(1);
t249 = -qJD(3) * t258 + qJD(1) * (pkin(2) * t258 + qJ(3) * t257) + t239;
t248 = qJD(1) * t258 * pkin(3) + t249;
t238 = qJD(3) * t257;
t247 = t238 + (-pkin(3) * t257 + t259) * qJD(1);
t245 = qJD(2) ^ 2;
t233 = t244 * rSges(2,1) - t242 * rSges(2,2);
t232 = t242 * rSges(2,1) + t244 * rSges(2,2);
t223 = t239 + qJD(1) * (rSges(3,1) * t258 - rSges(3,2) * t257);
t222 = (-rSges(3,1) * t257 - rSges(3,2) * t258 - t264) * qJD(1);
t221 = qJD(1) * (rSges(4,1) * t258 + rSges(4,3) * t257) + t249;
t220 = t238 + (-rSges(4,1) * t257 + rSges(4,3) * t258 + t259) * qJD(1);
t219 = t225 * rSges(6,3) + t224 * t256;
t218 = -t224 * rSges(6,3) + t225 * t256;
t213 = Icges(6,3) * t225 + t224 * t253;
t212 = -Icges(6,3) * t224 + t225 * t253;
t211 = t240 * (-t224 * rSges(5,1) - t225 * rSges(5,2)) + t248;
t210 = -t240 * (-t225 * rSges(5,1) + t224 * rSges(5,2)) + t247;
t209 = qJD(2) + (t218 * t225 + t219 * t224) * qJD(5);
t208 = -t225 * t261 + (-t224 * pkin(4) + t225 * pkin(7) + t219) * t240 + t248;
t207 = -t224 * t261 + (t225 * pkin(4) + t224 * pkin(7) - t218) * t240 + t247;
t1 = m(3) * (t222 ^ 2 + t223 ^ 2 + t245) / 0.2e1 + m(4) * (t220 ^ 2 + t221 ^ 2 + t245) / 0.2e1 + m(5) * (t210 ^ 2 + t211 ^ 2 + t245) / 0.2e1 + t240 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t207 ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + qJD(5) * t225 * (t225 * t268 + (-t212 * t269 + t272 * t213) * qJD(5)) / 0.2e1 - qJD(5) * t224 * (-t224 * t268 + (t273 * t212 - t213 * t269) * qJD(5)) / 0.2e1 + t240 * ((t243 ^ 2 * Icges(6,2) + (Icges(6,1) * t241 + 0.2e1 * Icges(6,4) * t243) * t241) * t240 + (t272 + t273) * t271 * qJD(5)) / 0.2e1 + (m(2) * (t232 ^ 2 + t233 ^ 2) + Icges(2,3) + Icges(3,3) + Icges(4,2)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
