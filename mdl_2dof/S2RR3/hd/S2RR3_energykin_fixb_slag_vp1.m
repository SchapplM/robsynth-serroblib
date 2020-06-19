% Calculate kinetic energy for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S2RR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_energykin_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_energykin_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR3_energykin_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR3_energykin_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:22
% EndTime: 2020-06-19 09:14:22
% DurationCPUTime: 0.05s
% Computational Cost: add. (18->13), mult. (32->27), div. (0->0), fcn. (10->4), ass. (0->12)
t66 = pkin(1) * qJD(1);
t64 = cos(qJ(1));
t63 = sin(qJ(1));
t62 = qJ(1) + qJ(2);
t61 = qJD(1) + qJD(2);
t60 = cos(t62);
t59 = sin(t62);
t58 = t64 * rSges(2,1) - t63 * rSges(2,2);
t57 = t63 * rSges(2,1) + t64 * rSges(2,2);
t56 = t64 * t66 + t61 * (t60 * rSges(3,1) - t59 * rSges(3,2));
t55 = -t63 * t66 - t61 * (t59 * rSges(3,1) + t60 * rSges(3,2));
t1 = m(3) * (t55 ^ 2 + t56 ^ 2) / 0.2e1 + t61 ^ 2 * Icges(3,3) / 0.2e1 + (m(2) * (t57 ^ 2 + t58 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2;
T = t1;
