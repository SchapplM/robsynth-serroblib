% Calculate kinetic energy for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2018-11-14 13:50
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPRP2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP2_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP2_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:49:22
% EndTime: 2018-11-14 13:49:22
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->40), mult. (115->67), div. (0->0), fcn. (72->4), ass. (0->24)
t171 = cos(qJ(3));
t179 = t171 * pkin(3);
t169 = sin(qJ(3));
t170 = sin(qJ(1));
t178 = t170 * t169;
t172 = cos(qJ(1));
t177 = t172 * t169;
t176 = -qJD(2) * t172 + qJD(1) * (t172 * pkin(1) + t170 * qJ(2));
t175 = qJD(1) * t172 * pkin(2) + t176;
t159 = t170 * pkin(1) - t172 * qJ(2);
t166 = qJD(2) * t170;
t174 = t166 + (-t170 * pkin(2) - t159) * qJD(1);
t168 = qJD(1) - qJD(3);
t161 = t172 * rSges(2,1) - t170 * rSges(2,2);
t160 = t170 * rSges(2,1) + t172 * rSges(2,2);
t157 = t170 * t171 - t177;
t156 = -t172 * t171 - t178;
t155 = qJD(1) * (t172 * rSges(3,1) + t170 * rSges(3,3)) + t176;
t154 = t166 + (-t170 * rSges(3,1) + t172 * rSges(3,3) - t159) * qJD(1);
t153 = t168 * (-t156 * rSges(4,1) + t157 * rSges(4,2)) + t175;
t152 = -t168 * (t157 * rSges(4,1) + t156 * rSges(4,2)) + t174;
t151 = (-t156 * rSges(5,1) + t157 * rSges(5,2) + pkin(3) * t178 + t179 * t172) * t168 + t175;
t150 = (-t157 * rSges(5,1) - t156 * rSges(5,2) + pkin(3) * t177 - t179 * t170) * t168 + t174;
t1 = m(3) * (t154 ^ 2 + t155 ^ 2) / 0.2e1 + m(4) * (t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + (Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t168 ^ 2 + (m(2) * (t160 ^ 2 + t161 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * qJD(1) ^ 2;
T  = t1;
