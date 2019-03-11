% Calculate kinetic energy for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S3RPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_energykin_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_energykin_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_energykin_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPR1_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:51
% EndTime: 2019-03-08 18:05:51
% DurationCPUTime: 0.03s
% Computational Cost: add. (35->29), mult. (69->50), div. (0->0), fcn. (38->4), ass. (0->17)
t115 = cos(qJ(1));
t114 = cos(qJ(3));
t113 = sin(qJ(1));
t112 = sin(qJ(3));
t111 = qJD(1) - qJD(3);
t110 = qJD(2) * t113;
t109 = t115 * rSges(2,1) - t113 * rSges(2,2);
t108 = t113 * rSges(2,1) + t115 * rSges(2,2);
t107 = t113 * pkin(1) - t115 * qJ(2);
t106 = qJD(1) * (t115 * pkin(1) + t113 * qJ(2));
t105 = -t115 * t112 + t113 * t114;
t104 = -t113 * t112 - t115 * t114;
t103 = t106 - qJD(2) * t115 + qJD(1) * (t115 * rSges(3,1) + t113 * rSges(3,3));
t102 = t110 + (-t113 * rSges(3,1) + t115 * rSges(3,3) - t107) * qJD(1);
t101 = t106 + t111 * (-t104 * rSges(4,1) + t105 * rSges(4,2)) + (qJD(1) * pkin(2) - qJD(2)) * t115;
t100 = t110 - t111 * (t105 * rSges(4,1) + t104 * rSges(4,2)) + (-t113 * pkin(2) - t107) * qJD(1);
t1 = m(3) * (t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2) / 0.2e1 + t111 ^ 2 * Icges(4,3) / 0.2e1 + (m(2) * (t108 ^ 2 + t109 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
