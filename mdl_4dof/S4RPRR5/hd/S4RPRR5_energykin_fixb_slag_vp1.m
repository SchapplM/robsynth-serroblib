% Calculate kinetic energy for
% S4RPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR5_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:51:31
% EndTime: 2019-12-31 16:51:32
% DurationCPUTime: 0.41s
% Computational Cost: add. (229->59), mult. (469->110), div. (0->0), fcn. (500->6), ass. (0->36)
t212 = sin(qJ(3));
t213 = sin(qJ(1));
t214 = cos(qJ(3));
t215 = cos(qJ(1));
t180 = -t213 * t212 - t215 * t214;
t222 = t180 ^ 2;
t181 = t215 * t212 - t213 * t214;
t221 = t181 ^ 2;
t196 = sin(qJ(4));
t197 = cos(qJ(4));
t220 = -Icges(5,5) * t196 - Icges(5,6) * t197;
t218 = t181 * t180;
t195 = qJD(1) - qJD(3);
t217 = t220 * t195;
t209 = qJD(4) * (-rSges(5,1) * t196 - rSges(5,2) * t197);
t208 = -qJD(2) * t215 + qJD(1) * (t215 * pkin(1) + t213 * qJ(2));
t207 = -rSges(5,1) * t197 + rSges(5,2) * t196;
t204 = -Icges(5,5) * t197 + Icges(5,6) * t196;
t200 = qJD(1) * t215 * pkin(2) + t208;
t187 = t213 * pkin(1) - t215 * qJ(2);
t194 = qJD(2) * t213;
t199 = t194 + (-t213 * pkin(2) - t187) * qJD(1);
t189 = t215 * rSges(2,1) - t213 * rSges(2,2);
t188 = t213 * rSges(2,1) + t215 * rSges(2,2);
t179 = qJD(1) * (t215 * rSges(3,1) + t213 * rSges(3,3)) + t208;
t178 = t194 + (-t213 * rSges(3,1) + t215 * rSges(3,3) - t187) * qJD(1);
t177 = rSges(5,3) * t181 + t207 * t180;
t176 = -rSges(5,3) * t180 + t207 * t181;
t171 = Icges(5,3) * t181 + t204 * t180;
t170 = -Icges(5,3) * t180 + t204 * t181;
t169 = t195 * (-rSges(4,1) * t180 - rSges(4,2) * t181) + t200;
t168 = -t195 * (-rSges(4,1) * t181 + rSges(4,2) * t180) + t199;
t167 = (t176 * t181 + t177 * t180) * qJD(4);
t166 = -t181 * t209 + (-pkin(3) * t180 + pkin(6) * t181 + t177) * t195 + t200;
t165 = -t180 * t209 + (pkin(3) * t181 + pkin(6) * t180 - t176) * t195 + t199;
t1 = m(3) * (t178 ^ 2 + t179 ^ 2) / 0.2e1 + m(4) * (t168 ^ 2 + t169 ^ 2) / 0.2e1 + t195 ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + qJD(4) * t181 * (t181 * t217 + (-t170 * t218 + t221 * t171) * qJD(4)) / 0.2e1 - qJD(4) * t180 * (-t180 * t217 + (t222 * t170 - t171 * t218) * qJD(4)) / 0.2e1 + t195 * ((t197 ^ 2 * Icges(5,2) + (Icges(5,1) * t196 + 0.2e1 * Icges(5,4) * t197) * t196) * t195 + (t221 + t222) * t220 * qJD(4)) / 0.2e1 + (m(2) * (t188 ^ 2 + t189 ^ 2) + Icges(2,3) + Icges(3,2)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
