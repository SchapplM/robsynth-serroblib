% Calculate kinetic energy for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S2RR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_energykin_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_energykin_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_energykin_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_energykin_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:06
% EndTime: 2020-01-03 11:19:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (86->49), mult. (213->96), div. (0->0), fcn. (174->4), ass. (0->33)
t79 = sin(qJ(2));
t94 = Icges(3,4) * t79;
t81 = cos(qJ(2));
t93 = Icges(3,4) * t81;
t80 = sin(qJ(1));
t92 = qJD(2) * t80;
t82 = cos(qJ(1));
t91 = qJD(2) * t82;
t90 = rSges(3,1) * t81 - rSges(3,2) * t79;
t85 = -Icges(3,2) * t79 + t93;
t67 = Icges(3,6) * t82 + t85 * t80;
t86 = Icges(3,1) * t81 - t94;
t69 = Icges(3,5) * t82 + t86 * t80;
t89 = -t67 * t79 + t69 * t81;
t68 = -Icges(3,6) * t80 + t85 * t82;
t70 = -Icges(3,5) * t80 + t86 * t82;
t88 = t68 * t79 - t70 * t81;
t74 = -Icges(3,2) * t81 - t94;
t75 = -Icges(3,1) * t79 - t93;
t87 = -t74 * t79 + t75 * t81;
t84 = Icges(3,5) * t81 - Icges(3,6) * t79;
t78 = t82 * rSges(2,1) - t80 * rSges(2,2);
t77 = t80 * rSges(2,1) + t82 * rSges(2,2);
t76 = -t79 * rSges(3,1) - t81 * rSges(3,2);
t73 = -Icges(3,5) * t79 - Icges(3,6) * t81;
t72 = -t80 * rSges(3,3) + t90 * t82;
t71 = t82 * rSges(3,3) + t90 * t80;
t66 = -Icges(3,3) * t80 + t84 * t82;
t65 = Icges(3,3) * t82 + t84 * t80;
t64 = t76 * t91 + (-pkin(1) * t82 - t71) * qJD(1);
t63 = t76 * t92 + (-pkin(1) * t80 + t72) * qJD(1);
t62 = (-t71 * t80 - t72 * t82) * qJD(2);
t1 = m(3) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + ((t82 * t73 + t87 * t80) * qJD(1) + (t82 ^ 2 * t65 + (t88 * t80 + (-t66 + t89) * t82) * t80) * qJD(2)) * t91 / 0.2e1 + qJD(1) * ((-t81 * t74 - t79 * t75) * qJD(1) + ((-t81 * t67 - t79 * t69) * t82 - (-t81 * t68 - t79 * t70) * t80) * qJD(2)) / 0.2e1 - ((-t80 * t73 + t87 * t82) * qJD(1) + (t80 ^ 2 * t66 + (t89 * t82 + (-t65 + t88) * t80) * t82) * qJD(2)) * t92 / 0.2e1 + (m(2) * (t77 ^ 2 + t78 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
