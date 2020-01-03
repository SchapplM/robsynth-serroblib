% Calculate kinetic energy for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:36
% EndTime: 2019-12-31 17:39:37
% DurationCPUTime: 0.45s
% Computational Cost: add. (460->65), mult. (475->115), div. (0->0), fcn. (500->6), ass. (0->39)
t257 = pkin(8) + qJ(2);
t255 = sin(t257);
t256 = cos(t257);
t261 = sin(qJ(4));
t262 = cos(qJ(4));
t225 = -t255 * t261 - t256 * t262;
t269 = t225 ^ 2;
t226 = -t255 * t262 + t256 * t261;
t268 = t226 ^ 2;
t241 = sin(qJ(5));
t242 = cos(qJ(5));
t267 = -Icges(6,5) * t241 - Icges(6,6) * t242;
t265 = t225 * t226;
t240 = qJD(2) - qJD(4);
t264 = t267 * t240;
t258 = qJD(5) * (-t241 * rSges(6,1) - t242 * rSges(6,2));
t254 = -rSges(6,1) * t242 + rSges(6,2) * t241;
t251 = -Icges(6,5) * t242 + Icges(6,6) * t241;
t247 = -qJD(3) * t256 + qJD(2) * (t256 * pkin(2) + t255 * qJ(3));
t246 = qJD(2) * t256 * pkin(3) + t247;
t228 = t255 * pkin(2) - t256 * qJ(3);
t239 = qJD(3) * t255;
t245 = t239 + (-t255 * pkin(3) - t228) * qJD(2);
t244 = qJD(1) ^ 2;
t243 = qJD(2) ^ 2;
t230 = t256 * rSges(3,1) - t255 * rSges(3,2);
t229 = t255 * rSges(3,1) + t256 * rSges(3,2);
t224 = qJD(2) * (t256 * rSges(4,1) + t255 * rSges(4,3)) + t247;
t223 = t239 + (-t255 * rSges(4,1) + t256 * rSges(4,3) - t228) * qJD(2);
t222 = t226 * rSges(6,3) + t254 * t225;
t221 = -t225 * rSges(6,3) + t254 * t226;
t216 = Icges(6,3) * t226 + t251 * t225;
t215 = -Icges(6,3) * t225 + t251 * t226;
t214 = t240 * (-t225 * rSges(5,1) - t226 * rSges(5,2)) + t246;
t213 = -t240 * (-t226 * rSges(5,1) + t225 * rSges(5,2)) + t245;
t212 = qJD(1) + (t221 * t226 + t222 * t225) * qJD(5);
t211 = -t226 * t258 + (-t225 * pkin(4) + t226 * pkin(7) + t222) * t240 + t246;
t210 = -t225 * t258 + (t226 * pkin(4) + t225 * pkin(7) - t221) * t240 + t245;
t1 = m(2) * t244 / 0.2e1 + m(3) * (t244 + (t229 ^ 2 + t230 ^ 2) * t243) / 0.2e1 + m(4) * (t223 ^ 2 + t224 ^ 2 + t244) / 0.2e1 + m(5) * (t213 ^ 2 + t214 ^ 2 + t244) / 0.2e1 + t240 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + qJD(5) * t226 * (t226 * t264 + (-t215 * t265 + t268 * t216) * qJD(5)) / 0.2e1 - qJD(5) * t225 * (-t225 * t264 + (t269 * t215 - t216 * t265) * qJD(5)) / 0.2e1 + t240 * ((t242 ^ 2 * Icges(6,2) + (Icges(6,1) * t241 + 0.2e1 * Icges(6,4) * t242) * t241) * t240 + (t268 + t269) * t267 * qJD(5)) / 0.2e1 + (Icges(4,2) + Icges(3,3)) * t243 / 0.2e1;
T = t1;
