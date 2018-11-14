% Calculate kinetic energy for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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

function T = S4RRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:29
% EndTime: 2018-11-14 13:51:29
% DurationCPUTime: 0.05s
% Computational Cost: add. (88->37), mult. (84->55), div. (0->0), fcn. (32->6), ass. (0->26)
t174 = rSges(5,1) + pkin(3);
t164 = qJ(1) + qJ(2);
t160 = sin(t164);
t173 = pkin(2) * t160;
t172 = pkin(1) * qJD(1);
t171 = rSges(5,3) + qJ(4);
t166 = cos(qJ(1));
t158 = t166 * t172;
t161 = cos(t164);
t163 = qJD(1) + qJD(2);
t170 = t163 * pkin(2) * t161 + t158;
t165 = sin(qJ(1));
t169 = t165 * t172;
t167 = qJD(3) ^ 2;
t159 = pkin(6) + t164;
t157 = cos(t159);
t156 = sin(t159);
t155 = t166 * rSges(2,1) - t165 * rSges(2,2);
t154 = t165 * rSges(2,1) + t166 * rSges(2,2);
t152 = t158 + t163 * (t161 * rSges(3,1) - t160 * rSges(3,2));
t151 = -t169 - t163 * (t160 * rSges(3,1) + t161 * rSges(3,2));
t150 = t163 * (t157 * rSges(4,1) - t156 * rSges(4,2)) + t170;
t149 = -t169 + (-t156 * rSges(4,1) - t157 * rSges(4,2) - t173) * t163;
t148 = -qJD(4) * t157 + (t171 * t156 + t174 * t157) * t163 + t170;
t147 = -t169 + qJD(4) * t156 + (-t174 * t156 + t171 * t157 - t173) * t163;
t1 = m(3) * (t151 ^ 2 + t152 ^ 2) / 0.2e1 + m(4) * (t149 ^ 2 + t150 ^ 2 + t167) / 0.2e1 + m(5) * (t147 ^ 2 + t148 ^ 2 + t167) / 0.2e1 + (m(2) * (t154 ^ 2 + t155 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2 + (Icges(3,3) / 0.2e1 + Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t163 ^ 2;
T  = t1;
