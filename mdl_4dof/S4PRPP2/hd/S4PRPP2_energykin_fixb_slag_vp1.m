% Calculate kinetic energy for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2018-11-14 14:00
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:00:15
% EndTime: 2018-11-14 14:00:15
% DurationCPUTime: 0.04s
% Computational Cost: add. (44->28), mult. (61->44), div. (0->0), fcn. (22->4), ass. (0->18)
t137 = m(3) / 0.2e1;
t136 = rSges(5,1) + pkin(3);
t129 = sin(qJ(2));
t135 = t129 * pkin(2);
t134 = rSges(5,3) + qJ(4);
t130 = cos(qJ(2));
t133 = qJD(2) * t130 * pkin(2) + qJD(1);
t131 = qJD(3) ^ 2;
t128 = qJ(2) + pkin(5);
t127 = cos(t128);
t126 = sin(t128);
t124 = t129 * rSges(3,1) + t130 * rSges(3,2);
t123 = qJD(1) + qJD(2) * (t130 * rSges(3,1) - t129 * rSges(3,2));
t122 = (-t126 * rSges(4,1) - t127 * rSges(4,2) - t135) * qJD(2);
t121 = qJD(2) * (t127 * rSges(4,1) - t126 * rSges(4,2)) + t133;
t120 = qJD(4) * t126 + (-t136 * t126 + t134 * t127 - t135) * qJD(2);
t119 = -qJD(4) * t127 + (t134 * t126 + t136 * t127) * qJD(2) + t133;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + t123 ^ 2 * t137 + m(4) * (t121 ^ 2 + t122 ^ 2 + t131) / 0.2e1 + m(5) * (t119 ^ 2 + t120 ^ 2 + t131) / 0.2e1 + (t124 ^ 2 * t137 + Icges(3,3) / 0.2e1 + Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * qJD(2) ^ 2;
T  = t1;
