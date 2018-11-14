% Calculate vector of inverse dynamics joint torques for
% S4PRRR1
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRRR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:18
% EndTime: 2018-11-14 13:44:19
% DurationCPUTime: 0.46s
% Computational Cost: add. (1030->73), mult. (620->84), div. (0->0), fcn. (288->6), ass. (0->47)
t49 = pkin(7) + qJ(2);
t44 = cos(t49);
t59 = pkin(2) * qJD(2);
t57 = t44 * t59;
t46 = qJ(3) + t49;
t39 = cos(t46);
t50 = qJD(2) + qJD(3);
t63 = t39 * t50;
t40 = qJ(4) + t46;
t34 = sin(t40);
t35 = cos(t40);
t21 = t35 * rSges(5,1) - t34 * rSges(5,2);
t45 = qJD(4) + t50;
t74 = t45 * t21;
t6 = pkin(3) * t63 + t57 + t74;
t43 = sin(t49);
t51 = qJD(2) ^ 2;
t76 = (-qJDD(2) * t43 - t44 * t51) * pkin(2) - g(1);
t37 = pkin(2) * t44;
t68 = pkin(2) * t43;
t75 = qJDD(2) * t37 - t51 * t68 - g(2);
t58 = t43 * t59;
t38 = sin(t46);
t24 = t38 * rSges(4,1) + t39 * rSges(4,2);
t61 = t50 * t24;
t9 = -t58 - t61;
t20 = t34 * rSges(5,1) + t35 * rSges(5,2);
t48 = qJDD(2) + qJDD(3);
t42 = qJDD(4) + t48;
t47 = t50 ^ 2;
t73 = -t45 * t74 - t42 * t20 + (-t38 * t48 - t39 * t47) * pkin(3) + t76;
t64 = t38 * t50;
t17 = rSges(4,1) * t63 - rSges(4,2) * t64;
t72 = -t50 * t17 - t48 * t24 + t76;
t15 = t45 * t20;
t71 = -t45 * t15 + t42 * t21 + (-t38 * t47 + t39 * t48) * pkin(3) + t75;
t25 = t39 * rSges(4,1) - t38 * rSges(4,2);
t70 = t48 * t25 - t50 * t61 + t75;
t36 = Icges(5,3) * t42;
t60 = Icges(4,3) * t48 + t36;
t14 = pkin(3) * t39 + t21;
t27 = t44 * rSges(3,1) - t43 * rSges(3,2);
t26 = t43 * rSges(3,1) + t44 * rSges(3,2);
t13 = -pkin(3) * t38 - t20;
t10 = t50 * t25 + t57;
t5 = -pkin(3) * t64 - t15 - t58;
t1 = [(-g(3) + qJDD(1)) * (m(2) + m(3) + m(4) + m(5)); Icges(3,3) * qJDD(2) + t60 + (t70 * (t25 + t37) + t72 * (-t24 - t68) + (-t17 - t57 + t10) * t9) * m(4) + (g(1) * t26 - g(2) * t27 + (t26 ^ 2 + t27 ^ 2) * qJDD(2)) * m(3) + (t71 * (t14 + t37) + t73 * (t13 - t68)) * m(5); t60 + (t73 * t13 + t71 * t14) * m(5) + (-t10 * t61 - t9 * t17 + (t50 * t9 + t70) * t25 + (t10 * t50 - t72) * t24) * m(4); t36 + ((t45 * t5 + t71) * t21 + (t45 * t6 - t73) * t20 - t5 * t74 - t6 * t15) * m(5);];
tau  = t1;
