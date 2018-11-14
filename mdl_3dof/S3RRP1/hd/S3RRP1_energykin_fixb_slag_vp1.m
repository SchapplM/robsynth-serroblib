% Calculate kinetic energy for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3RRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_energykin_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_energykin_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_energykin_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRP1_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:07
% EndTime: 2018-11-14 10:15:07
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->24), mult. (57->40), div. (0->0), fcn. (22->4), ass. (0->18)
t121 = rSges(4,1) + pkin(2);
t120 = pkin(1) * qJD(1);
t119 = rSges(4,3) + qJ(3);
t115 = sin(qJ(1));
t118 = t115 * t120;
t116 = cos(qJ(1));
t114 = qJ(1) + qJ(2);
t113 = qJD(1) + qJD(2);
t111 = cos(t114);
t110 = sin(t114);
t109 = t116 * t120;
t108 = t116 * rSges(2,1) - t115 * rSges(2,2);
t107 = t115 * rSges(2,1) + t116 * rSges(2,2);
t106 = t109 + t113 * (t111 * rSges(3,1) - t110 * rSges(3,2));
t105 = -t118 - t113 * (t110 * rSges(3,1) + t111 * rSges(3,2));
t104 = -qJD(3) * t111 + t109 + (t119 * t110 + t121 * t111) * t113;
t103 = -t118 + qJD(3) * t110 + (-t121 * t110 + t119 * t111) * t113;
t1 = m(3) * (t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2) / 0.2e1 + (Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * t113 ^ 2 + (m(2) * (t107 ^ 2 + t108 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
