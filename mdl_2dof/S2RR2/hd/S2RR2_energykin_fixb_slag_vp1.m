% Calculate kinetic energy for
% S2RR2
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
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S2RR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_energykin_fixb_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_energykin_fixb_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_energykin_fixb_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_energykin_fixb_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR2_energykin_fixb_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR2_energykin_fixb_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:42
% EndTime: 2019-03-08 18:00:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (86->49), mult. (213->96), div. (0->0), fcn. (174->4), ass. (0->33)
t82 = sin(qJ(2));
t97 = Icges(3,4) * t82;
t84 = cos(qJ(2));
t96 = Icges(3,4) * t84;
t83 = sin(qJ(1));
t95 = qJD(2) * t83;
t85 = cos(qJ(1));
t94 = qJD(2) * t85;
t93 = rSges(3,1) * t84 - rSges(3,2) * t82;
t88 = -Icges(3,2) * t82 + t96;
t70 = Icges(3,6) * t85 - t83 * t88;
t89 = Icges(3,1) * t84 - t97;
t72 = Icges(3,5) * t85 - t83 * t89;
t92 = -t70 * t82 + t72 * t84;
t71 = Icges(3,6) * t83 + t85 * t88;
t73 = Icges(3,5) * t83 + t85 * t89;
t91 = t71 * t82 - t73 * t84;
t77 = Icges(3,2) * t84 + t97;
t78 = Icges(3,1) * t82 + t96;
t90 = t77 * t82 - t78 * t84;
t87 = Icges(3,5) * t84 - Icges(3,6) * t82;
t81 = t85 * rSges(2,1) - t83 * rSges(2,2);
t80 = -t83 * rSges(2,1) - t85 * rSges(2,2);
t79 = t82 * rSges(3,1) + t84 * rSges(3,2);
t76 = Icges(3,5) * t82 + Icges(3,6) * t84;
t75 = t83 * rSges(3,3) + t85 * t93;
t74 = t85 * rSges(3,3) - t83 * t93;
t69 = Icges(3,3) * t83 + t85 * t87;
t68 = Icges(3,3) * t85 - t83 * t87;
t67 = -t79 * t94 + (pkin(1) * t85 + t74) * qJD(1);
t66 = t79 * t95 + (-pkin(1) * t83 - t75) * qJD(1);
t65 = (-t74 * t83 + t75 * t85) * qJD(2);
t1 = m(3) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + ((t83 * t76 - t85 * t90) * qJD(1) + (t83 ^ 2 * t69 + (t92 * t85 + (t68 - t91) * t83) * t85) * qJD(2)) * t95 / 0.2e1 + qJD(1) * ((t84 * t77 + t82 * t78) * qJD(1) + ((t84 * t71 + t82 * t73) * t83 + (t84 * t70 + t82 * t72) * t85) * qJD(2)) / 0.2e1 + ((t85 * t76 + t83 * t90) * qJD(1) + (t85 ^ 2 * t68 + (t91 * t83 + (t69 - t92) * t85) * t83) * qJD(2)) * t94 / 0.2e1 + (m(2) * (t80 ^ 2 + t81 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
