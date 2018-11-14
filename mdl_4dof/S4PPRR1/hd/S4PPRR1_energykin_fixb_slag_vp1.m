% Calculate kinetic energy for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:08
% EndTime: 2018-11-14 13:40:08
% DurationCPUTime: 0.05s
% Computational Cost: add. (46->29), mult. (72->53), div. (0->0), fcn. (48->6), ass. (0->22)
t170 = cos(qJ(3));
t175 = t170 * pkin(3);
t167 = sin(pkin(6));
t169 = sin(qJ(3));
t174 = t167 * t169;
t168 = cos(pkin(6));
t173 = t168 * t169;
t172 = qJD(1) ^ 2;
t166 = qJ(3) + qJ(4);
t165 = -qJD(3) - qJD(4);
t164 = cos(t166);
t163 = sin(t166);
t162 = qJD(2) * t167;
t160 = t167 * t170 - t173;
t159 = -t168 * t170 - t174;
t158 = -t168 * t163 + t167 * t164;
t157 = -t167 * t163 - t168 * t164;
t156 = -qJD(2) * t168 - qJD(3) * (-t159 * rSges(4,1) + t160 * rSges(4,2));
t155 = t162 + qJD(3) * (t160 * rSges(4,1) + t159 * rSges(4,2));
t154 = -qJD(3) * pkin(3) * t174 + t165 * (-t157 * rSges(5,1) + t158 * rSges(5,2)) + (-qJD(3) * t175 - qJD(2)) * t168;
t153 = t162 + qJD(3) * (-pkin(3) * t173 + t175 * t167) - t165 * (t158 * rSges(5,1) + t157 * rSges(5,2));
t1 = m(2) * t172 / 0.2e1 + m(3) * (t172 + (t167 ^ 2 + t168 ^ 2) * qJD(2) ^ 2) / 0.2e1 + m(4) * (t155 ^ 2 + t156 ^ 2 + t172) / 0.2e1 + qJD(3) ^ 2 * Icges(4,3) / 0.2e1 + m(5) * (t153 ^ 2 + t154 ^ 2 + t172) / 0.2e1 + t165 ^ 2 * Icges(5,3) / 0.2e1;
T  = t1;
