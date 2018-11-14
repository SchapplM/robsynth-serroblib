% Calculate kinetic energy for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3RPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_energykin_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_energykin_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_energykin_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPP1_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPP1_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:33
% EndTime: 2018-11-14 10:13:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (30->26), mult. (63->42), div. (0->0), fcn. (28->2), ass. (0->13)
t102 = rSges(4,3) + qJ(3);
t98 = sin(qJ(1));
t99 = cos(qJ(1));
t101 = -qJD(2) * t99 + qJD(1) * (pkin(1) * t99 + qJ(2) * t98);
t97 = qJD(2) * t98;
t96 = rSges(2,1) * t99 - rSges(2,2) * t98;
t95 = rSges(2,1) * t98 + rSges(2,2) * t99;
t94 = pkin(1) * t98 - qJ(2) * t99;
t92 = qJD(1) * (-rSges(3,2) * t99 + rSges(3,3) * t98) + t101;
t91 = t97 + (rSges(3,2) * t98 + rSges(3,3) * t99 - t94) * qJD(1);
t90 = qJD(3) * t98 + (t98 * rSges(4,2) + t102 * t99) * qJD(1) + t101;
t89 = qJD(3) * t99 + t97 + (t99 * rSges(4,2) - t102 * t98 - t94) * qJD(1);
t1 = m(3) * (t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2) / 0.2e1 + (m(2) * (t95 ^ 2 + t96 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(4,1) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
