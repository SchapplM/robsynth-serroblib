% Calculate vector of inverse dynamics joint torques for
% S4RPRR1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:29
% EndTime: 2018-11-14 13:50:30
% DurationCPUTime: 0.85s
% Computational Cost: add. (1127->95), mult. (832->112), div. (0->0), fcn. (392->8), ass. (0->61)
t61 = cos(qJ(1));
t55 = t61 * pkin(1);
t101 = qJDD(1) * t55 - g(2);
t60 = sin(qJ(1));
t62 = qJD(1) ^ 2;
t100 = (-qJDD(1) * t60 - t61 * t62) * pkin(1) - g(1);
t59 = qJ(1) + pkin(7);
t51 = sin(t59);
t52 = cos(t59);
t99 = (-qJDD(1) * t51 - t52 * t62) * pkin(2) + t100;
t44 = pkin(2) * t52;
t86 = pkin(1) * t60;
t71 = -pkin(2) * t51 - t86;
t98 = qJDD(1) * t44 + t71 * t62 + t101;
t54 = qJ(3) + t59;
t47 = qJ(4) + t54;
t40 = sin(t47);
t41 = cos(t47);
t22 = rSges(5,1) * t40 + rSges(5,2) * t41;
t45 = sin(t54);
t46 = cos(t54);
t57 = qJDD(1) + qJDD(3);
t50 = qJDD(4) + t57;
t58 = qJD(1) + qJD(3);
t53 = qJD(4) + t58;
t56 = t58 ^ 2;
t23 = t41 * rSges(5,1) - rSges(5,2) * t40;
t93 = t23 * t53;
t91 = -t93 * t53 - t22 * t50 + (-t45 * t57 - t46 * t56) * pkin(3) + t99;
t79 = t45 * t58;
t33 = rSges(4,2) * t79;
t78 = t46 * t58;
t19 = rSges(4,1) * t78 - t33;
t26 = rSges(4,1) * t45 + rSges(4,2) * t46;
t97 = t19 * t58 + t26 * t57 - t99;
t17 = t53 * t22;
t90 = -t17 * t53 + t23 * t50 + (-t45 * t56 + t46 * t57) * pkin(3) + t98;
t38 = t46 * rSges(4,1);
t27 = -rSges(4,2) * t45 + t38;
t77 = t58 * t26;
t96 = t27 * t57 - t58 * t77 + t98;
t95 = -pkin(3) * t78 - t93;
t92 = -pkin(3) * t79 - t17;
t29 = t52 * rSges(3,1) - rSges(3,2) * t51;
t25 = t29 + t55;
t89 = t44 + t55;
t68 = t71 * qJD(1);
t5 = t68 + t92;
t67 = t89 * qJD(1);
t10 = t27 * t58 + t67;
t81 = t10 * t26;
t42 = Icges(5,3) * t50;
t76 = Icges(4,3) * t57 + t42;
t16 = pkin(3) * t46 + t23;
t35 = rSges(2,1) * t61 - rSges(2,2) * t60;
t34 = rSges(2,1) * t60 + rSges(2,2) * t61;
t28 = rSges(3,1) * t51 + rSges(3,2) * t52;
t15 = -pkin(3) * t45 - t22;
t9 = t68 - t77;
t6 = t67 - t95;
t1 = [t76 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t5 * t95 + t6 * t92 + (-t5 * t89 + t6 * t71) * qJD(1) + t90 * (t16 + t89) + t91 * (t15 + t71)) * m(5) + (t9 * t33 + (-t9 * t38 - t81) * t58 + (t10 * t71 - t89 * t9) * qJD(1) + t96 * (t27 + t89) - t97 * (-t26 + t71)) * m(4) + ((qJDD(1) * t29 + t101) * t25 + (-qJDD(1) * t28 + (-0.2e1 * t29 + 0.2e1 * t25 - t55) * t62 + t100) * (-t28 - t86)) * m(3) + (g(1) * t34 - g(2) * t35 + (t34 ^ 2 + t35 ^ 2) * qJDD(1)) * m(2); (-g(3) + qJDD(2)) * (m(3) + m(4) + m(5)); t76 + (t91 * t15 + t90 * t16) * m(5) + (-t10 * t77 - t19 * t9 + t81 * t58 + (t9 * t58 + t96) * t27 + t97 * t26) * m(4); t42 + ((t5 * t53 + t90) * t23 + (t6 * t53 - t91) * t22 - t5 * t93 - t6 * t17) * m(5);];
tau  = t1;
