% Calculate kinetic energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
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
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3RRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energykin_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_energykin_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRR1_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:48
% EndTime: 2018-11-14 10:15:48
% DurationCPUTime: 0.04s
% Computational Cost: add. (45->24), mult. (53->43), div. (0->0), fcn. (18->6), ass. (0->21)
t110 = qJD(1) + qJD(2);
t117 = pkin(2) * t110;
t116 = pkin(1) * qJD(1);
t111 = qJ(1) + qJ(2);
t112 = sin(qJ(1));
t115 = t112 * t116;
t113 = cos(qJ(1));
t109 = qJ(3) + t111;
t108 = cos(t111);
t107 = sin(t111);
t106 = qJD(3) + t110;
t105 = t113 * t116;
t104 = cos(t109);
t103 = sin(t109);
t102 = t113 * rSges(2,1) - t112 * rSges(2,2);
t101 = t112 * rSges(2,1) + t113 * rSges(2,2);
t100 = t105 + t110 * (t108 * rSges(3,1) - t107 * rSges(3,2));
t99 = -t115 - t110 * (t107 * rSges(3,1) + t108 * rSges(3,2));
t98 = t105 + t108 * t117 + t106 * (t104 * rSges(4,1) - t103 * rSges(4,2));
t97 = -t115 - t107 * t117 - t106 * (t103 * rSges(4,1) + t104 * rSges(4,2));
t1 = m(3) * (t100 ^ 2 + t99 ^ 2) / 0.2e1 + t110 ^ 2 * Icges(3,3) / 0.2e1 + m(4) * (t97 ^ 2 + t98 ^ 2) / 0.2e1 + t106 ^ 2 * Icges(4,3) / 0.2e1 + (m(2) * (t101 ^ 2 + t102 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
