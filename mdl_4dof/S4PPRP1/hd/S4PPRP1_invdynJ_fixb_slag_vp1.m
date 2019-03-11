% Calculate vector of inverse dynamics joint torques for
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:14
% EndTime: 2019-03-08 18:12:14
% DurationCPUTime: 0.32s
% Computational Cost: add. (398->64), mult. (844->78), div. (0->0), fcn. (794->4), ass. (0->38)
t41 = sin(pkin(5));
t42 = cos(pkin(5));
t56 = sin(qJ(3));
t57 = cos(qJ(3));
t33 = -t41 * t57 + t42 * t56;
t22 = qJD(3) * t33;
t45 = t41 * t56 + t42 * t57;
t23 = t45 * qJD(3);
t24 = qJD(4) * t33;
t58 = rSges(5,1) + pkin(3);
t63 = rSges(5,3) + qJ(4);
t43 = t58 * t22 - t23 * t63 - t24;
t51 = qJDD(2) * t42;
t59 = -t63 * t33 - t45 * t58;
t2 = t43 * qJD(3) - qJD(4) * t22 + qJDD(3) * t59 + qJDD(4) * t45 - t51;
t64 = -g(2) + t2;
t25 = t45 * qJD(4);
t62 = -qJD(3) * t59 - t25;
t39 = qJDD(2) * t41;
t44 = -t63 * t22 - t58 * t23 + t25;
t7 = -t58 * t33 + t63 * t45;
t1 = t44 * qJD(3) + qJD(4) * t23 + t7 * qJDD(3) + qJDD(4) * t33 + t39;
t61 = -g(1) + t1;
t60 = t7 * qJD(3) + t24;
t52 = qJD(2) * t42;
t50 = -m(3) - m(4) - m(5);
t12 = -rSges(4,1) * t45 + t33 * rSges(4,2);
t40 = qJD(2) * t41;
t15 = -t33 * rSges(4,1) - rSges(4,2) * t45;
t11 = -t23 * rSges(4,1) + t22 * rSges(4,2);
t10 = -t22 * rSges(4,1) - t23 * rSges(4,2);
t9 = qJD(3) * t12 - t52;
t8 = qJD(3) * t15 + t40;
t6 = -qJD(3) * t10 + qJDD(3) * t12 - t51;
t5 = qJD(3) * t11 + qJDD(3) * t15 + t39;
t4 = -t52 - t62;
t3 = t40 + t60;
t13 = [(-g(3) + qJDD(1)) * (m(2) - t50); t50 * (g(1) * t41 - g(2) * t42) + m(4) * (t5 * t41 - t6 * t42) + m(5) * (t1 * t41 - t2 * t42) + m(3) * (t41 ^ 2 + t42 ^ 2) * qJDD(2); (Icges(4,3) + Icges(5,2)) * qJDD(3) + (t61 * t7 + (t43 + t60) * t4 + (t44 + t62) * t3 + t64 * t59) * m(5) + (-t9 * t10 + t8 * t11 + (qJD(3) * t9 - g(1) + t5) * t15 + (-qJD(3) * t8 - g(2) + t6) * t12) * m(4); (-t4 * t22 + t3 * t23 + (qJD(3) * t4 + t61) * t33 - (qJD(3) * t3 - t64) * t45) * m(5);];
tau  = t13;
