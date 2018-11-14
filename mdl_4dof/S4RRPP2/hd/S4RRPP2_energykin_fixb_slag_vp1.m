% Calculate kinetic energy for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RRPP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:24
% EndTime: 2018-11-14 13:52:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (81->37), mult. (87->58), div. (0->0), fcn. (36->4), ass. (0->22)
t155 = rSges(5,1) + pkin(3);
t154 = pkin(1) * qJD(1);
t149 = cos(qJ(1));
t142 = t149 * t154;
t147 = qJ(1) + qJ(2);
t143 = sin(t147);
t144 = cos(t147);
t146 = qJD(1) + qJD(2);
t153 = t146 * (t144 * pkin(2) + t143 * qJ(3)) + t142;
t148 = sin(qJ(1));
t152 = t148 * t154;
t151 = qJD(3) * t143 - t152;
t140 = t149 * rSges(2,1) - t148 * rSges(2,2);
t139 = t148 * rSges(2,1) + t149 * rSges(2,2);
t138 = t143 * pkin(2) - t144 * qJ(3);
t136 = t142 + t146 * (t144 * rSges(3,1) - t143 * rSges(3,2));
t135 = -t152 - t146 * (t143 * rSges(3,1) + t144 * rSges(3,2));
t134 = -qJD(3) * t144 + t146 * (t144 * rSges(4,1) + t143 * rSges(4,3)) + t153;
t133 = (-t143 * rSges(4,1) + t144 * rSges(4,3) - t138) * t146 + t151;
t132 = t146 * t143 * rSges(5,2) + (t155 * t146 - qJD(3)) * t144 + t153;
t131 = (t144 * rSges(5,2) - t155 * t143 - t138) * t146 + t151;
t1 = m(3) * (t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(4) * (t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + (m(2) * (t139 ^ 2 + t140 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2 + (Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1 + Icges(5,3) / 0.2e1) * t146 ^ 2;
T  = t1;
