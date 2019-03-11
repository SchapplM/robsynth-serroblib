% Calculate vector of inverse dynamics joint torques for
% S4PPRP2
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
%   pkin=[a2,a3,a4,d3,theta2]';
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
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:15
% EndTime: 2019-03-08 18:13:15
% DurationCPUTime: 0.24s
% Computational Cost: add. (353->46), mult. (351->46), div. (0->0), fcn. (175->2), ass. (0->25)
t31 = pkin(5) + qJ(3);
t30 = cos(t31);
t46 = rSges(5,3) + qJ(4);
t41 = t46 * t30;
t29 = sin(t31);
t47 = rSges(5,1) + pkin(3);
t42 = t47 * t29;
t38 = t41 - t42;
t43 = t46 * t29;
t9 = t47 * t30 + t43;
t32 = -qJD(3) * t47 + qJD(4);
t21 = qJD(4) * t29;
t33 = qJD(3) * t41 + t21;
t1 = -qJDD(4) * t30 + qJDD(1) + t9 * qJDD(3) + (t32 * t29 + t33) * qJD(3);
t45 = g(2) - t1;
t22 = qJD(4) * t30;
t44 = qJDD(4) * t29 + t38 * qJDD(3) + (-qJD(3) * t43 + t32 * t30 + t22) * qJD(3) - g(1);
t18 = t30 * rSges(4,1) - t29 * rSges(4,2);
t36 = qJD(3) * t18;
t34 = m(3) + m(4) + m(5);
t15 = t29 * rSges(4,1) + t30 * rSges(4,2);
t11 = qJD(3) * t15;
t10 = qJD(1) + t36;
t7 = -qJD(3) * t11 + qJDD(3) * t18 + qJDD(1);
t2 = [(m(2) + m(3)) * qJDD(1) + m(4) * t7 + m(5) * t1 + (-m(2) - t34) * g(2); (-g(3) + qJDD(2)) * t34; (Icges(4,3) + Icges(5,2)) * qJDD(3) + (-t10 * t11 + (-g(2) + t7) * t18 + (qJDD(3) * t15 + g(1) + (t10 + t36) * qJD(3)) * t15) * m(4) + (-t45 * t9 + t44 * t38 + (-t21 + t33 + (-t38 - t42) * qJD(3)) * (t9 * qJD(3) + qJD(1) - t22)) * m(5); (t44 * t29 + t45 * t30) * m(5);];
tau  = t2;
