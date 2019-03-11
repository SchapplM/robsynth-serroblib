% Calculate vector of inverse dynamics joint torques for
% S4PRPP2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:59
% EndTime: 2019-03-08 18:18:59
% DurationCPUTime: 0.42s
% Computational Cost: add. (429->65), mult. (523->67), div. (0->0), fcn. (255->4), ass. (0->36)
t42 = qJ(2) + pkin(5);
t39 = sin(t42);
t43 = sin(qJ(2));
t62 = t43 * pkin(2);
t70 = rSges(5,1) + pkin(3);
t76 = -t39 * t70 - t62;
t40 = cos(t42);
t67 = rSges(5,3) + qJ(4);
t75 = t67 * t40;
t49 = t75 + t76;
t44 = cos(qJ(2));
t41 = t44 * pkin(2);
t69 = t67 * t39;
t58 = t40 * t70 + t69;
t74 = t41 + t58;
t45 = qJD(2) ^ 2;
t73 = -t45 * t41 - g(1);
t29 = qJD(4) * t40;
t47 = -qJD(2) * t70 + qJD(4);
t72 = qJDD(4) * t39 + t49 * qJDD(2) + (-qJD(2) * t69 + t47 * t40 + t29) * qJD(2) + t73;
t55 = qJDD(2) * pkin(2);
t46 = t44 * t55 - t45 * t62 + qJDD(1);
t28 = qJD(4) * t39;
t51 = qJD(2) * t75 + t28;
t1 = -qJDD(4) * t40 + t58 * qJDD(2) + (t47 * t39 + t51) * qJD(2) + t46;
t71 = t1 - g(2);
t66 = m(4) + m(5);
t20 = t40 * rSges(4,1) - t39 * rSges(4,2);
t25 = t44 * rSges(3,1) - t43 * rSges(3,2);
t24 = t43 * rSges(3,1) + t44 * rSges(3,2);
t17 = t39 * rSges(4,1) + t40 * rSges(4,2);
t22 = qJD(2) * t24;
t21 = qJD(2) * t25 + qJD(1);
t10 = -qJD(2) * t22 + qJDD(2) * t25 + qJDD(1);
t5 = qJDD(2) * t20 - t17 * t45 + t46;
t2 = [m(2) * qJDD(1) + m(3) * t10 + m(4) * t5 + m(5) * t1 + (-m(2) - m(3) - t66) * g(2); (Icges(3,3) + Icges(4,3) + Icges(5,2)) * qJDD(2) + (-t21 * t22 + (t10 - g(2)) * t25 + (qJD(2) * t21 + qJDD(2) * t24 + t25 * t45 + g(1)) * t24) * m(3) + (t71 * t74 + t72 * t49 + (-t28 + t51 + (-t49 + t76) * qJD(2)) * (t74 * qJD(2) + qJD(1) - t29)) * m(5) + ((t5 - g(2)) * (t20 + t41) + (-qJDD(2) * t17 - t45 * t20 - t43 * t55 + t73) * (-t17 - t62)) * m(4); t66 * (-g(3) + qJDD(3)); (t39 * t72 - t71 * t40) * m(5);];
tau  = t2;
