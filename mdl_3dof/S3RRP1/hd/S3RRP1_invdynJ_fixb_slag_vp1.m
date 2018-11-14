% Calculate vector of inverse dynamics joint torques for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S3RRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp1: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRP1_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:07
% EndTime: 2018-11-14 10:15:08
% DurationCPUTime: 0.49s
% Computational Cost: add. (740->76), mult. (687->87), div. (0->0), fcn. (358->4), ass. (0->48)
t51 = qJ(1) + qJ(2);
t46 = sin(t51);
t34 = qJD(3) * t46;
t50 = qJD(1) + qJD(2);
t47 = cos(t51);
t37 = t47 * qJ(3);
t82 = rSges(4,1) + pkin(2);
t83 = -t47 * rSges(4,3) + t46 * t82 - t37;
t74 = t83 * t50;
t86 = t34 - t74;
t85 = t82 * t47;
t53 = cos(qJ(1));
t48 = t53 * pkin(1);
t54 = qJD(1) ^ 2;
t52 = sin(qJ(1));
t72 = t52 * pkin(1);
t81 = -qJDD(1) * t48 + t54 * t72 + g(2);
t80 = (-qJDD(1) * t52 - t54 * t53) * pkin(1) - g(1);
t65 = pkin(1) * qJD(1);
t62 = t52 * t65;
t23 = t46 * rSges(3,1) + t47 * rSges(3,2);
t68 = t50 * t23;
t13 = -t62 - t68;
t79 = rSges(4,3) + qJ(3);
t69 = t47 * t50;
t70 = t46 * t50;
t16 = rSges(3,1) * t69 - rSges(3,2) * t70;
t49 = qJDD(1) + qJDD(2);
t78 = -t50 * t16 - t49 * t23 + t80;
t26 = t47 * rSges(3,1) - t46 * rSges(3,2);
t77 = t49 * t26 - t50 * t68 - t81;
t12 = t79 * t46 + t85;
t57 = -t50 * t82 + qJD(3);
t73 = rSges(4,3) * t69 + t50 * t37;
t63 = t34 + t73;
t76 = qJDD(3) * t47 - t12 * t49 - (t46 * t57 + t63) * t50 + t81;
t35 = qJD(3) * t47;
t75 = qJDD(3) * t46 - t83 * t49 + (t47 * t57 - t70 * t79 + t35) * t50 + t80;
t66 = (Icges(4,2) + Icges(3,3)) * t49;
t61 = t53 * t65;
t58 = t35 - t61;
t33 = t53 * rSges(2,1) - t52 * rSges(2,2);
t32 = t52 * rSges(2,1) + t53 * rSges(2,2);
t7 = -t62 + t86;
t8 = t12 * t50 - t58;
t55 = (-t7 * t85 + (-t7 * t79 - t8 * t82) * t46) * t50;
t14 = t50 * t26 + t61;
t1 = [Icges(2,3) * qJDD(1) + t66 + (t77 * (t26 + t48) + t78 * (-t23 - t72) + (-t16 - t61 + t14) * t13) * m(3) + (g(1) * t32 - g(2) * t33 + (t32 ^ 2 + t33 ^ 2) * qJDD(1)) * m(2) + (t7 * t58 + t55 + t75 * (-t83 - t72) - t76 * (t48 + t12) + (t7 + t73 + t74) * t8) * m(4); t66 + (t55 + (t63 - t86) * t8 - t75 * t83 + (t50 * t7 - t76) * t12) * m(4) + (-t13 * t16 - t14 * t68 + (t13 * t50 + t77) * t26 + (t14 * t50 - t78) * t23) * m(3); (t75 * t46 + t76 * t47) * m(4);];
tau  = t1;
