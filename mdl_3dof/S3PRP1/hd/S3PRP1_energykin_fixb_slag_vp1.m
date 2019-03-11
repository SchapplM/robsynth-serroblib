% Calculate kinetic energy for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S3PRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_energykin_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_energykin_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_energykin_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRP1_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRP1_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:02:55
% EndTime: 2019-03-08 18:02:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (18->16), mult. (38->28), div. (0->0), fcn. (14->2), ass. (0->10)
t83 = rSges(4,1) + pkin(2);
t82 = rSges(4,3) + qJ(3);
t81 = qJD(2) ^ 2;
t80 = cos(qJ(2));
t79 = sin(qJ(2));
t78 = rSges(3,1) * t79 + rSges(3,2) * t80;
t77 = qJD(1) + qJD(2) * (rSges(3,1) * t80 - rSges(3,2) * t79);
t76 = qJD(3) * t79 + (-t79 * t83 + t80 * t82) * qJD(2);
t75 = -qJD(3) * t80 + qJD(1) + (t79 * t82 + t80 * t83) * qJD(2);
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t78 ^ 2 * t81 + t77 ^ 2) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2) / 0.2e1 + (Icges(3,3) + Icges(4,2)) * t81 / 0.2e1;
T  = t1;
