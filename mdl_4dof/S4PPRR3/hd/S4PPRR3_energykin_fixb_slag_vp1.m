% Calculate kinetic energy for
% S4PPRR3
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:21
% EndTime: 2019-12-31 16:17:21
% DurationCPUTime: 0.36s
% Computational Cost: add. (196->45), mult. (426->90), div. (0->0), fcn. (476->6), ass. (0->29)
t188 = sin(pkin(6));
t189 = cos(pkin(6));
t205 = sin(qJ(3));
t206 = cos(qJ(3));
t178 = -t188 * t205 - t189 * t206;
t212 = t178 ^ 2;
t179 = -t188 * t206 + t189 * t205;
t211 = t179 ^ 2;
t190 = sin(qJ(4));
t191 = cos(qJ(4));
t210 = -Icges(5,5) * t190 - Icges(5,6) * t191;
t208 = t179 * t178;
t207 = t210 * qJD(3);
t202 = qJD(2) * t189;
t201 = qJD(4) * (-rSges(5,1) * t190 - rSges(5,2) * t191);
t200 = -rSges(5,1) * t191 + rSges(5,2) * t190;
t197 = -Icges(5,5) * t191 + Icges(5,6) * t190;
t193 = qJD(1) ^ 2;
t187 = qJD(2) * t188;
t177 = -t202 - qJD(3) * (-rSges(4,1) * t178 - rSges(4,2) * t179);
t176 = t187 + qJD(3) * (-rSges(4,1) * t179 + rSges(4,2) * t178);
t175 = rSges(5,3) * t179 + t178 * t200;
t174 = -rSges(5,3) * t178 + t179 * t200;
t169 = Icges(5,3) * t179 + t178 * t197;
t168 = -Icges(5,3) * t178 + t179 * t197;
t167 = -t179 * t201 - t202 + (pkin(3) * t178 - pkin(5) * t179 - t175) * qJD(3);
t166 = -t178 * t201 + t187 + (-pkin(3) * t179 - pkin(5) * t178 + t174) * qJD(3);
t165 = qJD(1) + (t174 * t179 + t175 * t178) * qJD(4);
t1 = m(2) * t193 / 0.2e1 + m(3) * (t193 + (t188 ^ 2 + t189 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t176 ^ 2 + t177 ^ 2 + t193) / 0.2e1 + qJD(3) ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + qJD(4) * t179 * (-t179 * t207 + (-t168 * t208 + t211 * t169) * qJD(4)) / 0.2e1 - qJD(4) * t178 * (t178 * t207 + (t212 * t168 - t169 * t208) * qJD(4)) / 0.2e1 - qJD(3) * (-(t191 ^ 2 * Icges(5,2) + (Icges(5,1) * t190 + 0.2e1 * Icges(5,4) * t191) * t190) * qJD(3) + (t211 + t212) * t210 * qJD(4)) / 0.2e1;
T = t1;
