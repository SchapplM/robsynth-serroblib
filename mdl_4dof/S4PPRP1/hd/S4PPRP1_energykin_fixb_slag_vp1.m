% Calculate kinetic energy for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:04
% EndTime: 2019-03-08 18:12:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (39->25), mult. (80->41), div. (0->0), fcn. (62->4), ass. (0->16)
t177 = rSges(5,1) + pkin(3);
t176 = cos(qJ(3));
t175 = sin(qJ(3));
t174 = -rSges(5,3) - qJ(4);
t169 = cos(pkin(5));
t173 = qJD(2) * t169;
t172 = qJD(1) ^ 2;
t168 = sin(pkin(5));
t167 = qJD(2) * t168;
t163 = -t168 * t176 + t169 * t175;
t162 = -t168 * t175 - t169 * t176;
t161 = -t173 - qJD(3) * (-t162 * rSges(4,1) - t163 * rSges(4,2));
t160 = t167 + qJD(3) * (-t163 * rSges(4,1) + t162 * rSges(4,2));
t159 = -t173 - qJD(4) * t162 + (t177 * t162 + t174 * t163) * qJD(3);
t158 = qJD(4) * t163 + t167 + (t174 * t162 - t177 * t163) * qJD(3);
t1 = m(2) * t172 / 0.2e1 + m(3) * (t172 + (t168 ^ 2 + t169 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t172) / 0.2e1 + m(5) * (t158 ^ 2 + t159 ^ 2 + t172) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * qJD(3) ^ 2;
T  = t1;
