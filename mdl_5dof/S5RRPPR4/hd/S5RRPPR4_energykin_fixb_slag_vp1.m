% Calculate kinetic energy for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:39
% EndTime: 2019-12-31 19:27:40
% DurationCPUTime: 0.46s
% Computational Cost: add. (486->72), mult. (499->123), div. (0->0), fcn. (512->8), ass. (0->44)
t260 = qJ(1) + qJ(2);
t256 = sin(t260);
t257 = cos(t260);
t263 = sin(pkin(8));
t264 = cos(pkin(8));
t223 = -t256 * t263 - t257 * t264;
t273 = t223 ^ 2;
t224 = -t256 * t264 + t257 * t263;
t272 = t224 ^ 2;
t240 = sin(qJ(5));
t242 = cos(qJ(5));
t271 = -Icges(6,5) * t240 - Icges(6,6) * t242;
t269 = t224 * t223;
t239 = qJD(1) + qJD(2);
t268 = t271 * t239;
t265 = pkin(1) * qJD(1);
t259 = qJD(5) * (-rSges(6,1) * t240 - rSges(6,2) * t242);
t241 = sin(qJ(1));
t258 = t241 * t265;
t255 = qJD(3) * t256 - t258;
t254 = -rSges(6,1) * t242 + rSges(6,2) * t240;
t251 = -Icges(6,5) * t242 + Icges(6,6) * t240;
t225 = t256 * pkin(2) - t257 * qJ(3);
t247 = -t256 * pkin(3) - t225;
t243 = cos(qJ(1));
t237 = t243 * t265;
t246 = -qJD(3) * t257 + t239 * (t257 * pkin(2) + t256 * qJ(3)) + t237;
t245 = t239 * t257 * pkin(3) + t246;
t232 = rSges(2,1) * t243 - rSges(2,2) * t241;
t231 = rSges(2,1) * t241 + rSges(2,2) * t243;
t221 = t237 + t239 * (t257 * rSges(3,1) - t256 * rSges(3,2));
t220 = -t258 - t239 * (t256 * rSges(3,1) + t257 * rSges(3,2));
t219 = t239 * (t257 * rSges(4,1) + t256 * rSges(4,3)) + t246;
t218 = (-t256 * rSges(4,1) + t257 * rSges(4,3) - t225) * t239 + t255;
t217 = rSges(6,3) * t224 + t254 * t223;
t216 = -rSges(6,3) * t223 + t254 * t224;
t211 = Icges(6,3) * t224 + t251 * t223;
t210 = -Icges(6,3) * t223 + t251 * t224;
t209 = t239 * (-rSges(5,1) * t223 - rSges(5,2) * t224) + t245;
t208 = (t224 * rSges(5,1) - t223 * rSges(5,2) + t247) * t239 + t255;
t207 = -qJD(4) + (t216 * t224 + t217 * t223) * qJD(5);
t206 = -t224 * t259 + (-pkin(4) * t223 + pkin(7) * t224 + t217) * t239 + t245;
t205 = -t223 * t259 + (t224 * pkin(4) + t223 * pkin(7) - t216 + t247) * t239 + t255;
t1 = m(3) * (t220 ^ 2 + t221 ^ 2) / 0.2e1 + m(4) * (t218 ^ 2 + t219 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t208 ^ 2 + t209 ^ 2) / 0.2e1 + m(6) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + qJD(5) * t224 * (t224 * t268 + (-t210 * t269 + t272 * t211) * qJD(5)) / 0.2e1 - qJD(5) * t223 * (-t223 * t268 + (t273 * t210 - t211 * t269) * qJD(5)) / 0.2e1 + t239 * ((t242 ^ 2 * Icges(6,2) + (Icges(6,1) * t240 + 0.2e1 * Icges(6,4) * t242) * t240) * t239 + (t272 + t273) * t271 * qJD(5)) / 0.2e1 + (m(2) * (t231 ^ 2 + t232 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (Icges(3,3) + Icges(4,2) + Icges(5,3)) * t239 ^ 2 / 0.2e1;
T = t1;
