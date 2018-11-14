% Calculate kinetic energy for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:14
% EndTime: 2018-11-14 14:11:14
% DurationCPUTime: 0.06s
% Computational Cost: add. (43->28), mult. (57->47), div. (0->0), fcn. (18->6), ass. (0->20)
t125 = sin(qJ(2));
t130 = t125 * pkin(2);
t126 = cos(qJ(2));
t129 = qJD(2) * t126 * pkin(2) + qJD(1);
t124 = qJ(2) + pkin(6);
t128 = qJD(2) ^ 2;
t127 = qJD(3) ^ 2;
t123 = -qJD(2) - qJD(4);
t122 = qJ(4) + t124;
t121 = cos(t124);
t120 = sin(t124);
t118 = cos(t122);
t117 = sin(t122);
t116 = t125 * rSges(3,1) + t126 * rSges(3,2);
t115 = qJD(1) + qJD(2) * (t126 * rSges(3,1) - t125 * rSges(3,2));
t114 = (-t120 * rSges(4,1) - t121 * rSges(4,2) - t130) * qJD(2);
t113 = qJD(2) * (t121 * rSges(4,1) - t120 * rSges(4,2)) + t129;
t112 = t123 * (t117 * rSges(5,1) + t118 * rSges(5,2)) + (-pkin(3) * t120 - t130) * qJD(2);
t111 = qJD(2) * pkin(3) * t121 - t123 * (t118 * rSges(5,1) - t117 * rSges(5,2)) + t129;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t128 * t116 ^ 2 + t115 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t127) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t127) / 0.2e1 + t123 ^ 2 * Icges(5,3) / 0.2e1 + (Icges(3,3) + Icges(4,3)) * t128 / 0.2e1;
T  = t1;
