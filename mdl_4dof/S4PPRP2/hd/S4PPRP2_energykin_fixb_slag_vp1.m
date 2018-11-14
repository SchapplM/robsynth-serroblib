% Calculate kinetic energy for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:57
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:57:21
% EndTime: 2018-11-14 13:57:21
% DurationCPUTime: 0.03s
% Computational Cost: add. (36->21), mult. (44->31), div. (0->0), fcn. (14->2), ass. (0->13)
t131 = rSges(5,1) + pkin(3);
t130 = rSges(5,3) + qJ(4);
t129 = qJD(1) ^ 2;
t128 = qJD(2) ^ 2;
t127 = qJD(3) ^ 2;
t126 = pkin(5) + qJ(3);
t125 = cos(t126);
t124 = sin(t126);
t123 = t124 * rSges(4,1) + t125 * rSges(4,2);
t122 = qJD(1) + qJD(3) * (t125 * rSges(4,1) - t124 * rSges(4,2));
t121 = qJD(4) * t124 + (-t131 * t124 + t130 * t125) * qJD(3);
t120 = -qJD(4) * t125 + qJD(1) + (t130 * t124 + t131 * t125) * qJD(3);
t1 = m(2) * t129 / 0.2e1 + m(3) * (t128 + t129) / 0.2e1 + m(4) * (t127 * t123 ^ 2 + t122 ^ 2 + t128) / 0.2e1 + m(5) * (t120 ^ 2 + t121 ^ 2 + t128) / 0.2e1 + (Icges(4,3) + Icges(5,2)) * t127 / 0.2e1;
T  = t1;
