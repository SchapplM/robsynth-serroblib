% Calculate vector of inverse dynamics joint torques for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:49
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:48:37
% EndTime: 2018-11-14 13:48:38
% DurationCPUTime: 0.81s
% Computational Cost: add. (1202->99), mult. (907->114), div. (0->0), fcn. (462->6), ass. (0->59)
t65 = cos(qJ(1));
t60 = t65 * pkin(1);
t108 = qJDD(1) * t60 - g(2);
t64 = sin(qJ(1));
t66 = qJD(1) ^ 2;
t107 = (-qJDD(1) * t64 - t65 * t66) * pkin(1) - g(1);
t63 = qJ(1) + pkin(6);
t59 = qJ(3) + t63;
t53 = cos(t59);
t44 = t53 * qJ(4);
t52 = sin(t59);
t98 = rSges(5,1) + pkin(3);
t106 = t53 * rSges(5,3) - t52 * t98 + t44;
t58 = cos(t63);
t51 = pkin(2) * t58;
t57 = sin(t63);
t90 = pkin(1) * t64;
t77 = -pkin(2) * t57 - t90;
t105 = -qJDD(1) * t51 - t77 * t66 - t108;
t104 = (-qJDD(1) * t57 - t58 * t66) * pkin(2) + t107;
t40 = qJD(4) * t52;
t62 = qJD(1) + qJD(3);
t103 = t106 * t62 + t40;
t48 = t53 * rSges(4,1);
t30 = -rSges(4,2) * t52 + t48;
t61 = qJDD(1) + qJDD(3);
t27 = rSges(4,1) * t52 + rSges(4,2) * t53;
t84 = t62 * t27;
t102 = t30 * t61 - t62 * t84 - t105;
t100 = t98 * t53;
t97 = rSges(5,3) + qJ(4);
t14 = t97 * t52 + t100;
t74 = -t62 * t98 + qJD(4);
t85 = t53 * t62;
t78 = rSges(5,3) * t85 + t62 * t44 + t40;
t96 = qJDD(4) * t53 - t14 * t61 - (t74 * t52 + t78) * t62 + t105;
t41 = qJD(4) * t53;
t86 = t52 * t62;
t95 = qJDD(4) * t52 + t106 * t61 + (t74 * t53 - t86 * t97 + t41) * t62 + t104;
t36 = rSges(4,2) * t86;
t18 = rSges(4,1) * t85 - t36;
t101 = t18 * t62 + t27 * t61 - t104;
t73 = t77 * qJD(1);
t7 = t73 + t103;
t94 = t7 * t41;
t32 = t58 * rSges(3,1) - rSges(3,2) * t57;
t24 = t32 + t60;
t92 = t51 + t60;
t72 = t92 * qJD(1);
t12 = t30 * t62 + t72;
t87 = t12 * t27;
t81 = (Icges(5,2) + Icges(4,3)) * t61;
t39 = rSges(2,1) * t65 - rSges(2,2) * t64;
t38 = rSges(2,1) * t64 + rSges(2,2) * t65;
t31 = rSges(3,1) * t57 + rSges(3,2) * t58;
t11 = t73 - t84;
t8 = t14 * t62 - t41 + t72;
t67 = t94 + t8 * t78 + (-t7 * t100 + (-t7 * t97 - t8 * t98) * t52) * t62;
t1 = [t81 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t11 * t36 + (-t11 * t48 - t87) * t62 + (-t11 * t92 + t12 * t77) * qJD(1) + t102 * (t30 + t92) - t101 * (-t27 + t77)) * m(4) + ((qJDD(1) * t32 + t108) * t24 + (-qJDD(1) * t31 + (-0.2e1 * t32 + 0.2e1 * t24 - t60) * t66 + t107) * (-t31 - t90)) * m(3) + (g(1) * t38 - g(2) * t39 + (t38 ^ 2 + t39 ^ 2) * qJDD(1)) * m(2) + ((-t7 * t92 + t8 * t77) * qJD(1) + t67 + t95 * (t106 + t77) - t96 * (t14 + t92)) * m(5); (-g(3) + qJDD(2)) * (m(3) + m(4) + m(5)); t81 + (t67 - t94 - t8 * t103 + (t7 * t62 - t96) * t14 + t95 * t106) * m(5) + (-t11 * t18 - t12 * t84 + t87 * t62 + (t11 * t62 + t102) * t30 + t101 * t27) * m(4); (t95 * t52 + t96 * t53) * m(5);];
tau  = t1;
