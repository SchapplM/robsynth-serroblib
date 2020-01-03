% Calculate kinetic energy for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR6_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:37
% DurationCPUTime: 0.31s
% Computational Cost: add. (218->97), mult. (531->160), div. (0->0), fcn. (540->6), ass. (0->49)
t233 = cos(pkin(6));
t252 = t233 ^ 2;
t250 = pkin(3) * t233;
t235 = sin(qJ(1));
t237 = cos(qJ(1));
t226 = t235 * pkin(1) - t237 * qJ(2);
t232 = sin(pkin(6));
t240 = pkin(2) * t233 + qJ(3) * t232;
t249 = -t240 * t235 - t226;
t231 = qJD(2) * t235;
t247 = qJD(3) * t232;
t248 = t237 * t247 + t231;
t246 = qJD(4) * t235;
t245 = qJD(4) * t237;
t244 = t235 * qJD(1);
t225 = qJD(1) * (t237 * pkin(1) + t235 * qJ(2));
t243 = qJD(1) * t240 * t237 + t235 * t247 + t225;
t242 = rSges(3,1) * t233 - rSges(3,2) * t232;
t241 = rSges(4,1) * t233 + rSges(4,3) * t232;
t234 = sin(qJ(4));
t236 = cos(qJ(4));
t224 = t232 * t236 - t233 * t234;
t239 = t232 * t234 + t233 * t236;
t228 = t237 * rSges(2,1) - t235 * rSges(2,2);
t227 = t235 * rSges(2,1) + t237 * rSges(2,2);
t220 = t239 * t237;
t219 = t224 * t237;
t218 = t239 * t235;
t217 = t224 * t235;
t216 = t224 * rSges(5,1) - rSges(5,2) * t239;
t215 = Icges(5,1) * t224 - Icges(5,4) * t239;
t214 = Icges(5,4) * t224 - Icges(5,2) * t239;
t213 = Icges(5,5) * t224 - Icges(5,6) * t239;
t212 = rSges(3,3) * t244 + t225 + (qJD(1) * t242 - qJD(2)) * t237;
t211 = t231 + (t237 * rSges(3,3) - t242 * t235 - t226) * qJD(1);
t210 = t220 * rSges(5,1) + t219 * rSges(5,2) - t235 * rSges(5,3);
t209 = t218 * rSges(5,1) + t217 * rSges(5,2) + t237 * rSges(5,3);
t208 = Icges(5,1) * t220 + Icges(5,4) * t219 - Icges(5,5) * t235;
t207 = Icges(5,1) * t218 + Icges(5,4) * t217 + Icges(5,5) * t237;
t206 = Icges(5,4) * t220 + Icges(5,2) * t219 - Icges(5,6) * t235;
t205 = Icges(5,4) * t218 + Icges(5,2) * t217 + Icges(5,6) * t237;
t204 = Icges(5,5) * t220 + Icges(5,6) * t219 - Icges(5,3) * t235;
t203 = Icges(5,5) * t218 + Icges(5,6) * t217 + Icges(5,3) * t237;
t202 = rSges(4,2) * t244 + (qJD(1) * t241 - qJD(2)) * t237 + t243;
t201 = (t237 * rSges(4,2) - t241 * t235 + t249) * qJD(1) + t248;
t200 = -qJD(3) * t233 + (-t209 * t235 - t210 * t237) * qJD(4);
t199 = t216 * t246 - qJD(2) * t237 + (-t235 * pkin(5) + t237 * t250 + t210) * qJD(1) + t243;
t198 = t216 * t245 + (-t237 * pkin(5) - t235 * t250 - t209 + t249) * qJD(1) + t248;
t1 = m(3) * (t211 ^ 2 + t212 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t252 + t201 ^ 2 + t202 ^ 2) / 0.2e1 + m(5) * (t198 ^ 2 + t199 ^ 2 + t200 ^ 2) / 0.2e1 - ((-t235 * t213 + t219 * t214 + t220 * t215) * qJD(1) + (-(-t235 * t204 + t219 * t206 + t220 * t208) * t235 + (-t235 * t203 + t219 * t205 + t220 * t207) * t237) * qJD(4)) * t246 / 0.2e1 + ((t237 * t213 + t217 * t214 + t218 * t215) * qJD(1) + (-(t237 * t204 + t217 * t206 + t218 * t208) * t235 + (t237 * t203 + t217 * t205 + t218 * t207) * t237) * qJD(4)) * t245 / 0.2e1 + qJD(1) * ((-t214 * t239 + t224 * t215) * qJD(1) + (-(-t206 * t239 + t224 * t208) * t235 + (-t205 * t239 + t224 * t207) * t237) * qJD(4)) / 0.2e1 + (m(2) * (t227 ^ 2 + t228 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t252 + ((Icges(3,1) + Icges(4,1)) * t232 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t233) * t232) * qJD(1) ^ 2 / 0.2e1;
T = t1;
