% Calculate vector of inverse dynamics joint torques for
% S4PRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:23
% DurationCPUTime: 0.41s
% Computational Cost: add. (748->84), mult. (620->93), div. (0->0), fcn. (288->6), ass. (0->54)
t48 = sin(qJ(2));
t49 = cos(qJ(2));
t50 = qJD(2) ^ 2;
t78 = (-qJDD(2) * t48 - t49 * t50) * pkin(1) - g(3);
t70 = t49 * pkin(1);
t71 = t48 * pkin(1);
t77 = -qJDD(2) * t70 + t50 * t71 - g(1);
t46 = qJD(2) + qJD(3);
t47 = qJ(2) + qJ(3);
t41 = sin(t47);
t42 = cos(t47);
t54 = t41 * rSges(4,1) + t42 * rSges(4,2);
t15 = t54 * t46;
t64 = pkin(1) * qJD(2);
t61 = t48 * t64;
t11 = t15 + t61;
t43 = qJ(4) + t47;
t36 = sin(t43);
t37 = cos(t43);
t20 = t37 * rSges(5,1) - t36 * rSges(5,2);
t45 = -qJDD(2) - qJDD(3);
t39 = -qJDD(4) + t45;
t40 = qJD(4) + t46;
t72 = t46 ^ 2;
t68 = t37 * t40;
t69 = t36 * t40;
t9 = rSges(5,1) * t69 + rSges(5,2) * t68;
t74 = t39 * t20 + t40 * t9 + (t41 * t72 + t42 * t45) * pkin(2) + t77;
t22 = t42 * rSges(4,1) - t41 * rSges(4,2);
t76 = t46 * t15 + t45 * t22 + t77;
t10 = -rSges(5,1) * t68 + rSges(5,2) * t69;
t19 = -t36 * rSges(5,1) - t37 * rSges(5,2);
t73 = t40 * t10 - t39 * t19 + (t41 * t45 - t42 * t72) * pkin(2) + t78;
t66 = t42 * t46;
t67 = t41 * t46;
t16 = -rSges(4,1) * t66 + rSges(4,2) * t67;
t75 = t46 * t16 + t45 * t54 + t78;
t65 = Icges(5,3) * t39;
t63 = pkin(2) * t67;
t30 = pkin(2) * t66;
t62 = t63 + t9;
t60 = t49 * t64;
t58 = -t40 * t19 + t63;
t57 = -t40 * t20 - t30;
t33 = t49 * rSges(3,1) - t48 * rSges(3,2);
t55 = t48 * rSges(3,1) + t49 * rSges(3,2);
t53 = -Icges(4,3) * t45 - t65;
t14 = -pkin(2) * t42 - t20;
t12 = -t46 * t22 - t60;
t6 = t57 - t60;
t13 = -pkin(2) * t41 + t19;
t51 = t10 - t30;
t5 = t58 + t61;
t1 = [(g(2) + qJDD(1)) * (m(2) + m(3) + m(4) + m(5)); Icges(3,3) * qJDD(2) + t53 + (t6 * (t61 + t62) - t5 * (t51 - t60) + t74 * (t14 - t70) + t73 * (t13 - t71)) * m(5) + (t76 * (-t22 - t70) + t75 * (-t54 - t71) + (t12 - t16 + t60) * t11) * m(4) + ((qJDD(2) * t55 + g(3)) * t55 + (qJDD(2) * t33 + g(1)) * t33) * m(3); t53 + ((-t58 + t62) * t6 + (-t51 + t57) * t5 + t74 * t14 + t73 * t13) * m(5) + (-t11 * t16 + t12 * t15 + (-t11 * t46 - t76) * t22 - (t12 * t46 + t75) * t54) * m(4); -t65 + (-t5 * t10 + t6 * t9 + (-t40 * t5 - t74) * t20 + (t40 * t6 + t73) * t19) * m(5);];
tau  = t1;
