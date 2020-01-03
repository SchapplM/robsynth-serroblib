% Calculate kinetic energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:52:01
% EndTime: 2020-01-03 11:52:02
% DurationCPUTime: 0.32s
% Computational Cost: add. (474->87), mult. (293->144), div. (0->0), fcn. (206->10), ass. (0->60)
t229 = qJD(1) + qJD(3);
t254 = pkin(3) * t229;
t253 = pkin(1) * qJD(1);
t252 = pkin(2) * qJD(1);
t231 = sin(qJ(5));
t251 = Icges(6,4) * t231;
t233 = cos(qJ(5));
t250 = Icges(6,4) * t233;
t232 = sin(qJ(1));
t223 = t232 * t253;
t230 = qJ(1) + pkin(9);
t225 = sin(t230);
t249 = t225 * t252 + t223;
t234 = cos(qJ(1));
t224 = t234 * t253;
t226 = cos(t230);
t248 = t226 * t252 + t224;
t228 = qJ(3) + t230;
t222 = qJ(4) + t228;
t218 = sin(t222);
t247 = qJD(5) * t218;
t219 = cos(t222);
t246 = qJD(5) * t219;
t220 = sin(t228);
t245 = t220 * t254 + t249;
t221 = cos(t228);
t244 = t221 * t254 + t248;
t243 = rSges(6,1) * t233 - rSges(6,2) * t231;
t242 = Icges(6,1) * t233 - t251;
t241 = -Icges(6,2) * t231 + t250;
t240 = Icges(6,5) * t233 - Icges(6,6) * t231;
t200 = -Icges(6,6) * t219 + t241 * t218;
t202 = -Icges(6,5) * t219 + t242 * t218;
t239 = -t200 * t231 + t202 * t233;
t201 = -Icges(6,6) * t218 - t241 * t219;
t203 = -Icges(6,5) * t218 - t242 * t219;
t238 = t201 * t231 - t203 * t233;
t211 = Icges(6,2) * t233 + t251;
t212 = Icges(6,1) * t231 + t250;
t237 = t211 * t231 - t212 * t233;
t235 = qJD(2) ^ 2;
t227 = qJD(4) + t229;
t215 = -t234 * rSges(2,1) + t232 * rSges(2,2);
t214 = t232 * rSges(2,1) + t234 * rSges(2,2);
t213 = t231 * rSges(6,1) + t233 * rSges(6,2);
t210 = Icges(6,5) * t231 + Icges(6,6) * t233;
t207 = t224 - qJD(1) * (-t226 * rSges(3,1) + t225 * rSges(3,2));
t206 = t223 + qJD(1) * (t225 * rSges(3,1) + t226 * rSges(3,2));
t205 = -t218 * rSges(6,3) - t243 * t219;
t204 = -t219 * rSges(6,3) + t243 * t218;
t199 = -Icges(6,3) * t218 - t240 * t219;
t198 = -Icges(6,3) * t219 + t240 * t218;
t197 = -t229 * (-t221 * rSges(4,1) + t220 * rSges(4,2)) + t248;
t196 = t229 * (t220 * rSges(4,1) + t221 * rSges(4,2)) + t249;
t195 = -t227 * (-t219 * rSges(5,1) + t218 * rSges(5,2)) + t244;
t194 = t227 * (t218 * rSges(5,1) + t219 * rSges(5,2)) + t245;
t193 = qJD(2) + (t204 * t218 - t205 * t219) * qJD(5);
t192 = -t213 * t247 + (t219 * pkin(4) + t218 * pkin(8) - t205) * t227 + t244;
t191 = t213 * t246 + (t218 * pkin(4) - t219 * pkin(8) + t204) * t227 + t245;
t1 = m(3) * (t206 ^ 2 + t207 ^ 2 + t235) / 0.2e1 + m(4) * (t196 ^ 2 + t197 ^ 2 + t235) / 0.2e1 + t229 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t194 ^ 2 + t195 ^ 2 + t235) / 0.2e1 + t227 ^ 2 * Icges(5,3) / 0.2e1 + m(6) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + t227 * ((t233 * t211 + t231 * t212) * t227 + (-(t233 * t200 + t231 * t202) * t219 - (t233 * t201 + t231 * t203) * t218) * qJD(5)) / 0.2e1 - ((-t219 * t210 - t237 * t218) * t227 + (t219 ^ 2 * t198 + (t238 * t218 + (t199 - t239) * t219) * t218) * qJD(5)) * t246 / 0.2e1 - ((-t218 * t210 + t237 * t219) * t227 + (t218 ^ 2 * t199 + (t239 * t219 + (t198 - t238) * t218) * t219) * qJD(5)) * t247 / 0.2e1 + (m(2) * (t214 ^ 2 + t215 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
