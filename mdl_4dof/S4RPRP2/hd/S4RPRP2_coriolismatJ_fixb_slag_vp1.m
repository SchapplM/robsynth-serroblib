% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:52
% EndTime: 2019-03-08 18:30:52
% DurationCPUTime: 0.35s
% Computational Cost: add. (1344->60), mult. (2492->90), div. (0->0), fcn. (2530->4), ass. (0->48)
t113 = m(4) / 0.2e1;
t68 = sin(qJ(1));
t91 = cos(qJ(3));
t77 = t68 * t91;
t67 = sin(qJ(3));
t69 = cos(qJ(1));
t87 = t69 * t67;
t57 = -t77 + t87;
t88 = t68 * t67;
t72 = t69 * t91 + t88;
t45 = -rSges(4,1) * t72 + t57 * rSges(4,2);
t76 = -t57 * rSges(4,1) - rSges(4,2) * t72;
t120 = t45 * t68 + t69 * t76;
t126 = t120 * t113;
t111 = m(5) / 0.2e1;
t65 = t91 * pkin(3) + pkin(2);
t74 = t57 * rSges(5,1) + rSges(5,2) * t72 + pkin(3) * t87;
t119 = t68 * t65 - t74;
t124 = -t68 * pkin(2) + t119;
t66 = t69 * qJ(2);
t122 = -t68 * pkin(1) - t119 + t66;
t107 = pkin(1) + pkin(2);
t42 = t68 * qJ(2) + t107 * t69 - t45;
t70 = -t107 * t68 + t66 - t76;
t121 = t42 * t76 - t45 * t70;
t75 = rSges(5,1) * t72 - t57 * rSges(5,2);
t34 = -pkin(3) * t88 + (pkin(2) - t65) * t69 - t75;
t86 = (t124 * t69 + t34 * t68) * t111 + t126;
t116 = 4 * qJD(1);
t115 = 2 * qJD(3);
t106 = m(3) * ((t69 * rSges(3,3) + t66) * t69 + (rSges(3,3) + qJ(2)) * t68 ^ 2);
t105 = m(4) * t121;
t102 = m(4) * (t42 * t68 + t69 * t70);
t101 = m(4) * t120;
t29 = (pkin(1) + t65) * t69 + (pkin(3) * t67 + qJ(2)) * t68 + t75;
t38 = t72 * pkin(3) + t75;
t39 = -pkin(3) * t77 + t74;
t99 = m(5) * (t122 * t38 - t39 * t29);
t98 = m(5) * (t122 * t69 + t29 * t68);
t96 = m(5) * (t68 * t38 + t39 * t69);
t85 = t96 / 0.2e1 - t101 / 0.2e1;
t78 = -t96 / 0.2e1;
t6 = t102 + t98 + t106;
t4 = t99 + t105;
t3 = t101 / 0.2e1 + t78 + t86;
t2 = t85 + t86;
t1 = t85 - t86;
t5 = [t6 * qJD(2) + t4 * qJD(3), t6 * qJD(1) + t2 * qJD(3), t4 * qJD(1) + t2 * qJD(2) + (-t121 * t113 + ((t29 - t34) * t39 + (-t122 + t124) * t38) * t111) * t115, 0; t1 * qJD(3) + (-t106 / 0.4e1 - t102 / 0.4e1 - t98 / 0.4e1) * t116, 0, t1 * qJD(1) + (t78 + t126) * t115, 0; t3 * qJD(2) + (-t105 / 0.4e1 - t99 / 0.4e1) * t116, t3 * qJD(1) (-t124 * t38 + t34 * t39) * m(5) * qJD(3), 0; 0, 0, 0, 0;];
Cq  = t5;
