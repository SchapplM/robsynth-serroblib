% Calculate kinetic energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:15
% EndTime: 2019-03-08 18:26:16
% DurationCPUTime: 0.20s
% Computational Cost: add. (119->72), mult. (285->100), div. (0->0), fcn. (259->6), ass. (0->40)
t228 = cos(pkin(4));
t249 = t228 ^ 2;
t246 = rSges(5,3) + qJ(4);
t226 = sin(pkin(4));
t229 = sin(qJ(1));
t245 = t226 * t229;
t230 = cos(qJ(1));
t244 = t226 * t230;
t225 = sin(pkin(6));
t243 = t229 * t225;
t227 = cos(pkin(6));
t242 = t229 * t227;
t241 = t230 * t225;
t240 = t230 * t227;
t212 = -t228 * t240 + t243;
t213 = t228 * t241 + t242;
t218 = t229 * pkin(1) - qJ(2) * t244;
t239 = -t213 * pkin(2) - t212 * qJ(3) - t218;
t214 = t228 * t242 + t241;
t237 = qJD(2) * t226;
t223 = t229 * t237;
t238 = qJD(3) * t214 + t223;
t236 = qJD(3) * t227;
t235 = 0.2e1 * t228;
t234 = t226 * (rSges(5,1) + pkin(3));
t233 = -t230 * t237 + qJD(1) * (t230 * pkin(1) + qJ(2) * t245);
t215 = -t228 * t243 + t240;
t232 = qJD(1) * (t215 * pkin(2) + t214 * qJ(3)) + qJD(3) * t212 + t233;
t224 = qJD(2) * t228;
t220 = t230 * rSges(2,1) - t229 * rSges(2,2);
t219 = t229 * rSges(2,1) + t230 * rSges(2,2);
t217 = -t226 * t236 + t224;
t209 = t224 + (qJD(4) * t225 - t236) * t226;
t206 = qJD(1) * (t215 * rSges(3,1) - t214 * rSges(3,2) + rSges(3,3) * t245) + t233;
t205 = t223 + (-t213 * rSges(3,1) + t212 * rSges(3,2) + rSges(3,3) * t244 - t218) * qJD(1);
t204 = qJD(1) * (rSges(4,1) * t245 - t215 * rSges(4,2) + t214 * rSges(4,3)) + t232;
t203 = (rSges(4,1) * t244 + t213 * rSges(4,2) - t212 * rSges(4,3) + t239) * qJD(1) + t238;
t202 = qJD(4) * t213 + (t214 * rSges(5,2) + t246 * t215 + t229 * t234) * qJD(1) + t232;
t201 = qJD(4) * t215 + (-t212 * rSges(5,2) - t246 * t213 + t230 * t234 + t239) * qJD(1) + t238;
t1 = m(3) * (qJD(2) ^ 2 * t249 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + m(4) * (t203 ^ 2 + t204 ^ 2 + t217 ^ 2) / 0.2e1 + m(5) * (t201 ^ 2 + t202 ^ 2 + t209 ^ 2) / 0.2e1 + (m(2) * (t219 ^ 2 + t220 ^ 2) + Icges(2,3) + (Icges(3,3) + Icges(5,1) + Icges(4,1)) * t249 + ((Icges(3,2) + Icges(5,2) + Icges(4,3)) * t226 * t227 ^ 2 + (Icges(3,6) - Icges(4,5) - Icges(5,4)) * t235 * t227 + ((-0.2e1 * Icges(5,6) * t227 + (Icges(3,1) + Icges(4,2) + Icges(5,3)) * t225 + 0.2e1 * (Icges(3,4) + Icges(4,6)) * t227) * t226 + (Icges(3,5) + Icges(5,5) - Icges(4,4)) * t235) * t225) * t226) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
