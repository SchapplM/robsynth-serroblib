% Calculate kinetic energy for
% S4PPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:39
% EndTime: 2019-12-31 16:19:40
% DurationCPUTime: 0.52s
% Computational Cost: add. (230->89), mult. (575->161), div. (0->0), fcn. (574->6), ass. (0->50)
t213 = cos(pkin(6));
t235 = t213 ^ 2;
t212 = sin(pkin(6));
t236 = t212 ^ 2;
t237 = t235 + t236;
t239 = qJD(3) * t237;
t238 = t212 * t213;
t217 = cos(qJ(3));
t234 = t212 * t217;
t233 = t213 * t217;
t214 = sin(qJ(4));
t215 = sin(qJ(3));
t232 = t214 * t215;
t216 = cos(qJ(4));
t231 = t215 * t216;
t230 = qJD(3) * t212;
t229 = qJD(4) * t215;
t228 = qJD(4) * t217;
t222 = Icges(4,5) * t215 + Icges(4,6) * t217;
t219 = qJD(1) ^ 2;
t211 = qJD(2) * t212;
t210 = t217 * pkin(3) + t215 * pkin(5);
t209 = t217 * rSges(4,1) - t215 * rSges(4,2);
t208 = qJD(3) * t213 - t212 * t228;
t207 = t213 * t228 + t230;
t206 = t212 * t214 - t213 * t231;
t205 = t212 * t216 + t213 * t232;
t204 = t212 * t231 + t213 * t214;
t203 = -t212 * t232 + t213 * t216;
t202 = t215 * rSges(5,3) + (rSges(5,1) * t216 - rSges(5,2) * t214) * t217;
t201 = Icges(5,5) * t215 + (Icges(5,1) * t216 - Icges(5,4) * t214) * t217;
t200 = Icges(5,6) * t215 + (Icges(5,4) * t216 - Icges(5,2) * t214) * t217;
t199 = Icges(5,3) * t215 + (Icges(5,5) * t216 - Icges(5,6) * t214) * t217;
t198 = (-qJD(3) * t209 - qJD(2)) * t213;
t197 = t209 * t230 + t211;
t192 = Icges(4,3) * t212 - t222 * t213;
t191 = Icges(4,3) * t213 + t222 * t212;
t190 = t206 * rSges(5,1) + t205 * rSges(5,2) + rSges(5,3) * t233;
t189 = t204 * rSges(5,1) + t203 * rSges(5,2) - rSges(5,3) * t234;
t188 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t233;
t187 = Icges(5,1) * t204 + Icges(5,4) * t203 - Icges(5,5) * t234;
t186 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t233;
t185 = Icges(5,4) * t204 + Icges(5,2) * t203 - Icges(5,6) * t234;
t184 = Icges(5,5) * t206 + Icges(5,6) * t205 + Icges(5,3) * t233;
t183 = Icges(5,5) * t204 + Icges(5,6) * t203 - Icges(5,3) * t234;
t182 = qJD(1) - (rSges(4,1) * t215 + rSges(4,2) * t217) * t239;
t181 = t189 * t229 - t208 * t202 + (-qJD(3) * t210 - qJD(2)) * t213;
t180 = -t190 * t229 + t207 * t202 + t210 * t230 + t211;
t179 = -t207 * t189 + t208 * t190 + qJD(1) + (-pkin(3) * t215 + pkin(5) * t217) * t239;
t1 = m(2) * t219 / 0.2e1 + m(3) * (t237 * qJD(2) ^ 2 + t219) / 0.2e1 + m(4) * (t182 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(5) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + t208 * ((-t183 * t234 + t203 * t185 + t204 * t187) * t208 + (-t184 * t234 + t203 * t186 + t204 * t188) * t207 + (-t199 * t234 + t203 * t200 + t204 * t201) * t229) / 0.2e1 + t207 * ((t183 * t233 + t205 * t185 + t206 * t187) * t208 + (t184 * t233 + t205 * t186 + t206 * t188) * t207 + (t199 * t233 + t205 * t200 + t206 * t201) * t229) / 0.2e1 + ((t183 * t208 + t184 * t207 + t199 * t229) * t215 + ((-t185 * t214 + t187 * t216) * t208 + (-t186 * t214 + t188 * t216) * t207 + (-t200 * t214 + t201 * t216) * t229) * t217) * t229 / 0.2e1 + (t213 * (t235 * t191 + t192 * t238) + t212 * (t191 * t238 + t236 * t192)) * qJD(3) ^ 2 / 0.2e1;
T = t1;
