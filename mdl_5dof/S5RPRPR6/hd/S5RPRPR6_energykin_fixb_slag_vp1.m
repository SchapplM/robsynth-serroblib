% Calculate kinetic energy for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_energykin_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR6_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR6_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:41
% EndTime: 2019-12-31 18:17:41
% DurationCPUTime: 0.28s
% Computational Cost: add. (383->89), mult. (299->144), div. (0->0), fcn. (212->8), ass. (0->53)
t236 = sin(qJ(1));
t255 = t236 * pkin(1);
t235 = sin(qJ(5));
t254 = Icges(6,4) * t235;
t237 = cos(qJ(5));
t253 = Icges(6,4) * t237;
t238 = cos(qJ(1));
t228 = qJD(1) * t238 * pkin(1);
t234 = qJ(1) + pkin(8);
t230 = cos(t234);
t252 = qJD(1) * pkin(2) * t230 + t228;
t231 = qJ(3) + t234;
t226 = sin(t231);
t251 = qJD(5) * t226;
t227 = cos(t231);
t233 = qJD(1) + qJD(3);
t250 = t233 * (t227 * pkin(3) + t226 * qJ(4)) + t252;
t249 = rSges(6,1) * t235 + rSges(6,2) * t237;
t248 = Icges(6,1) * t235 + t253;
t247 = Icges(6,2) * t237 + t254;
t246 = Icges(6,5) * t235 + Icges(6,6) * t237;
t208 = Icges(6,6) * t227 + t226 * t247;
t210 = Icges(6,5) * t227 + t226 * t248;
t245 = -t208 * t237 - t210 * t235;
t209 = Icges(6,6) * t226 - t227 * t247;
t211 = Icges(6,5) * t226 - t227 * t248;
t244 = t209 * t237 + t211 * t235;
t219 = -Icges(6,2) * t235 + t253;
t220 = Icges(6,1) * t237 - t254;
t243 = t219 * t237 + t220 * t235;
t229 = sin(t234);
t242 = (-pkin(2) * t229 - t255) * qJD(1);
t241 = qJD(4) * t226 + t242;
t239 = qJD(2) ^ 2;
t223 = t238 * rSges(2,1) - t236 * rSges(2,2);
t222 = t237 * rSges(6,1) - t235 * rSges(6,2);
t221 = t236 * rSges(2,1) + t238 * rSges(2,2);
t218 = Icges(6,5) * t237 - Icges(6,6) * t235;
t217 = t226 * pkin(3) - t227 * qJ(4);
t215 = t228 + qJD(1) * (t230 * rSges(3,1) - t229 * rSges(3,2));
t214 = (-t229 * rSges(3,1) - t230 * rSges(3,2) - t255) * qJD(1);
t213 = t226 * rSges(6,3) - t227 * t249;
t212 = t227 * rSges(6,3) + t226 * t249;
t207 = Icges(6,3) * t226 - t227 * t246;
t206 = Icges(6,3) * t227 + t226 * t246;
t205 = t233 * (t227 * rSges(4,1) - t226 * rSges(4,2)) + t252;
t204 = -t233 * (t226 * rSges(4,1) + t227 * rSges(4,2)) + t242;
t203 = -qJD(4) * t227 + t233 * (-t227 * rSges(5,2) + t226 * rSges(5,3)) + t250;
t202 = (t226 * rSges(5,2) + t227 * rSges(5,3) - t217) * t233 + t241;
t201 = qJD(2) + (-t212 * t226 + t213 * t227) * qJD(5);
t200 = t233 * t212 + (pkin(7) * t233 - qJD(5) * t222 - qJD(4)) * t227 + t250;
t199 = t222 * t251 + (-pkin(7) * t226 - t213 - t217) * t233 + t241;
t1 = m(3) * (t214 ^ 2 + t215 ^ 2 + t239) / 0.2e1 + m(4) * (t204 ^ 2 + t205 ^ 2 + t239) / 0.2e1 + m(5) * (t202 ^ 2 + t203 ^ 2 + t239) / 0.2e1 + m(6) * (t199 ^ 2 + t200 ^ 2 + t201 ^ 2) / 0.2e1 + qJD(5) * t227 * ((t227 * t218 + t226 * t243) * t233 + (t227 ^ 2 * t206 + (t244 * t226 + (t207 - t245) * t227) * t226) * qJD(5)) / 0.2e1 + ((t226 * t218 - t243 * t227) * t233 + (t226 ^ 2 * t207 + (t245 * t227 + (t206 - t244) * t226) * t227) * qJD(5)) * t251 / 0.2e1 + t233 * ((-t235 * t219 + t237 * t220) * t233 + ((-t235 * t208 + t237 * t210) * t227 + (-t235 * t209 + t237 * t211) * t226) * qJD(5)) / 0.2e1 + (Icges(4,3) + Icges(5,1)) * t233 ^ 2 / 0.2e1 + (m(2) * (t221 ^ 2 + t223 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
