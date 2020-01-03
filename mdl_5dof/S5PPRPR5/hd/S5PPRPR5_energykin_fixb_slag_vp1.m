% Calculate kinetic energy for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energykin_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:24
% EndTime: 2019-12-31 17:33:24
% DurationCPUTime: 0.41s
% Computational Cost: add. (228->58), mult. (484->103), div. (0->0), fcn. (534->6), ass. (0->34)
t243 = sin(pkin(7));
t244 = cos(pkin(7));
t264 = sin(qJ(3));
t265 = cos(qJ(3));
t233 = -t243 * t264 - t244 * t265;
t271 = t233 ^ 2;
t234 = -t243 * t265 + t244 * t264;
t270 = t234 ^ 2;
t245 = sin(qJ(5));
t246 = cos(qJ(5));
t269 = Icges(6,5) * t246 - Icges(6,6) * t245;
t268 = t234 * t233;
t267 = t269 * qJD(3);
t260 = qJD(2) * t244;
t259 = qJD(5) * (-rSges(6,1) * t246 + rSges(6,2) * t245);
t242 = qJD(2) * t243;
t258 = qJD(3) * (-pkin(3) * t234 - qJ(4) * t233) + qJD(4) * t234 + t242;
t257 = -qJD(4) * t233 - t260;
t256 = rSges(6,1) * t245 + rSges(6,2) * t246;
t253 = Icges(6,5) * t245 + Icges(6,6) * t246;
t249 = qJD(1) ^ 2;
t230 = -pkin(3) * t233 + qJ(4) * t234;
t228 = -t260 - qJD(3) * (-rSges(4,1) * t233 - rSges(4,2) * t234);
t227 = t242 + qJD(3) * (-rSges(4,1) * t234 + rSges(4,2) * t233);
t226 = -t234 * rSges(6,3) - t256 * t233;
t225 = -t233 * rSges(6,3) + t256 * t234;
t220 = -Icges(6,3) * t234 - t253 * t233;
t219 = -Icges(6,3) * t233 + t253 * t234;
t218 = (-rSges(5,2) * t233 - rSges(5,3) * t234 - t230) * qJD(3) + t257;
t217 = qJD(3) * (rSges(5,2) * t234 - rSges(5,3) * t233) + t258;
t216 = qJD(1) + (t225 * t234 - t226 * t233) * qJD(5);
t215 = t233 * t259 + (pkin(6) * t233 - t225 - t230) * qJD(3) + t257;
t214 = -t234 * t259 + (-pkin(6) * t234 + t226) * qJD(3) + t258;
t1 = m(2) * t249 / 0.2e1 + m(3) * (t249 + (t243 ^ 2 + t244 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t227 ^ 2 + t228 ^ 2 + t249) / 0.2e1 + m(5) * (t217 ^ 2 + t218 ^ 2 + t249) / 0.2e1 + m(6) * (t214 ^ 2 + t215 ^ 2 + t216 ^ 2) / 0.2e1 - qJD(3) * (-(t246 ^ 2 * Icges(6,1) + (-0.2e1 * Icges(6,4) * t246 + Icges(6,2) * t245) * t245) * qJD(3) + (-t270 - t271) * t269 * qJD(5)) / 0.2e1 + (Icges(4,3) + Icges(5,1)) * qJD(3) ^ 2 / 0.2e1 - (t233 * (-t233 * t267 + (t271 * t219 + t220 * t268) * qJD(5)) + t234 * (-t234 * t267 + (t219 * t268 + t270 * t220) * qJD(5))) * qJD(5) / 0.2e1;
T = t1;
