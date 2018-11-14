% Calculate kinetic energy for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3PRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_energykin_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_energykin_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_energykin_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRR1_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRR1_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:56
% EndTime: 2018-11-14 10:11:56
% DurationCPUTime: 0.03s
% Computational Cost: add. (21->16), mult. (34->32), div. (0->0), fcn. (10->4), ass. (0->13)
t78 = pkin(2) * qJD(2);
t77 = qJD(2) ^ 2;
t76 = cos(qJ(2));
t75 = sin(qJ(2));
t74 = qJ(2) + qJ(3);
t73 = qJD(2) + qJD(3);
t72 = cos(t74);
t71 = sin(t74);
t70 = t75 * rSges(3,1) + t76 * rSges(3,2);
t69 = qJD(1) + qJD(2) * (t76 * rSges(3,1) - t75 * rSges(3,2));
t68 = -t75 * t78 - t73 * (t71 * rSges(4,1) + t72 * rSges(4,2));
t67 = qJD(1) + t76 * t78 + t73 * (t72 * rSges(4,1) - t71 * rSges(4,2));
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t77 * t70 ^ 2 + t69 ^ 2) / 0.2e1 + t77 * Icges(3,3) / 0.2e1 + m(4) * (t67 ^ 2 + t68 ^ 2) / 0.2e1 + t73 ^ 2 * Icges(4,3) / 0.2e1;
T  = t1;
