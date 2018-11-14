% Calculate kinetic energy for
% S4PPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
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
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR5_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR5_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR5_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR5_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:00
% EndTime: 2018-11-14 14:05:00
% DurationCPUTime: 0.02s
% Computational Cost: add. (19->16), mult. (39->31), div. (0->0), fcn. (20->4), ass. (0->12)
t112 = cos(pkin(5));
t110 = -qJD(3) * t112 + qJD(1);
t116 = qJD(1) ^ 2;
t115 = qJD(2) ^ 2;
t114 = cos(qJ(4));
t113 = sin(qJ(4));
t111 = sin(pkin(5));
t109 = t111 * t114 - t112 * t113;
t108 = -t111 * t113 - t112 * t114;
t107 = qJD(3) * t111 + qJD(4) * (t109 * rSges(5,1) + t108 * rSges(5,2));
t106 = -qJD(4) * (-t108 * rSges(5,1) + t109 * rSges(5,2)) + t110;
t1 = m(2) * t116 / 0.2e1 + m(3) * (t115 + t116) / 0.2e1 + m(4) * (qJD(3) ^ 2 * t111 ^ 2 + t110 ^ 2 + t115) / 0.2e1 + m(5) * (t106 ^ 2 + t107 ^ 2 + t115) / 0.2e1 + qJD(4) ^ 2 * Icges(5,3) / 0.2e1;
T  = t1;
