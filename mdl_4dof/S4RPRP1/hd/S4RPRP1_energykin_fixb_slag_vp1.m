% Calculate kinetic energy for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:27
% EndTime: 2018-11-14 13:48:27
% DurationCPUTime: 0.05s
% Computational Cost: add. (82->37), mult. (85->55), div. (0->0), fcn. (32->6), ass. (0->25)
t171 = rSges(5,1) + pkin(3);
t163 = sin(qJ(1));
t170 = t163 * pkin(1);
t169 = rSges(5,3) + qJ(4);
t164 = cos(qJ(1));
t156 = qJD(1) * t164 * pkin(1);
t162 = qJ(1) + pkin(6);
t158 = cos(t162);
t168 = qJD(1) * pkin(2) * t158 + t156;
t157 = sin(t162);
t167 = (-pkin(2) * t157 - t170) * qJD(1);
t165 = qJD(2) ^ 2;
t161 = qJD(1) + qJD(3);
t159 = qJ(3) + t162;
t155 = cos(t159);
t154 = sin(t159);
t152 = t164 * rSges(2,1) - t163 * rSges(2,2);
t151 = t163 * rSges(2,1) + t164 * rSges(2,2);
t150 = t156 + qJD(1) * (t158 * rSges(3,1) - t157 * rSges(3,2));
t149 = (-t157 * rSges(3,1) - t158 * rSges(3,2) - t170) * qJD(1);
t148 = t161 * (t155 * rSges(4,1) - t154 * rSges(4,2)) + t168;
t147 = -t161 * (t154 * rSges(4,1) + t155 * rSges(4,2)) + t167;
t146 = -qJD(4) * t155 + (t169 * t154 + t171 * t155) * t161 + t168;
t145 = qJD(4) * t154 + (-t171 * t154 + t169 * t155) * t161 + t167;
t1 = m(3) * (t149 ^ 2 + t150 ^ 2 + t165) / 0.2e1 + m(4) * (t147 ^ 2 + t148 ^ 2 + t165) / 0.2e1 + m(5) * (t145 ^ 2 + t146 ^ 2 + t165) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t161 ^ 2 + (m(2) * (t151 ^ 2 + t152 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
