% Calculate kinetic energy for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2018-11-14 13:41
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP1_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP1_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:59
% EndTime: 2018-11-14 13:40:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (62->31), mult. (69->45), div. (0->0), fcn. (28->2), ass. (0->15)
t144 = rSges(5,3) + qJ(4);
t140 = pkin(5) + qJ(2);
t138 = sin(t140);
t139 = cos(t140);
t143 = -qJD(3) * t139 + qJD(2) * (t139 * pkin(2) + t138 * qJ(3));
t142 = qJD(1) ^ 2;
t137 = qJD(3) * t138;
t136 = t139 * rSges(3,1) - t138 * rSges(3,2);
t135 = t138 * rSges(3,1) + t139 * rSges(3,2);
t134 = t138 * pkin(2) - t139 * qJ(3);
t132 = qJD(2) * (-t139 * rSges(4,2) + t138 * rSges(4,3)) + t143;
t131 = t137 + (t138 * rSges(4,2) + t139 * rSges(4,3) - t134) * qJD(2);
t130 = qJD(4) * t138 + (t138 * rSges(5,2) + t144 * t139) * qJD(2) + t143;
t129 = qJD(4) * t139 + t137 + (t139 * rSges(5,2) - t144 * t138 - t134) * qJD(2);
t1 = m(4) * (t131 ^ 2 + t132 ^ 2 + t142) / 0.2e1 + m(5) * (t129 ^ 2 + t130 ^ 2 + t142) / 0.2e1 + (m(3) * (t135 ^ 2 + t136 ^ 2) / 0.2e1 + Icges(3,3) / 0.2e1 + Icges(4,1) / 0.2e1 + Icges(5,1) / 0.2e1) * qJD(2) ^ 2 + (m(2) + m(3)) * t142 / 0.2e1;
T  = t1;
