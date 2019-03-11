% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:20:52
% EndTime: 2019-03-08 18:20:53
% DurationCPUTime: 0.20s
% Computational Cost: add. (1164->25), mult. (1070->41), div. (0->0), fcn. (1066->4), ass. (0->27)
t39 = pkin(6) + qJ(2);
t37 = sin(t39);
t38 = cos(t39);
t51 = sin(qJ(4));
t52 = cos(qJ(4));
t29 = -t37 * t51 - t38 * t52;
t30 = -t37 * t52 + t38 * t51;
t22 = t29 * rSges(5,1) + t30 * rSges(5,2);
t41 = -t30 * rSges(5,1) + t29 * rSges(5,2);
t61 = t22 * t37 + t38 * t41;
t65 = m(5) * t61;
t10 = -t65 / 0.2e1;
t11 = t65 / 0.2e1;
t42 = m(5) * qJD(4);
t57 = pkin(2) + pkin(3);
t19 = t37 * qJ(3) + t57 * t38 - t22;
t36 = t38 * qJ(3);
t40 = -t57 * t37 + t36 - t41;
t5 = t19 * t41 - t22 * t40;
t64 = t5 * t42;
t4 = m(4) * ((t38 * rSges(4,3) + t36) * t38 + (rSges(4,3) + qJ(3)) * t37 ^ 2) + m(5) * (t19 * t37 + t38 * t40);
t63 = t4 * qJD(2);
t62 = m(5) * t5 * qJD(2);
t3 = 0.2e1 * t11;
t2 = t10 + t11;
t1 = 0.2e1 * t10;
t6 = [0, 0, 0, 0; 0, t4 * qJD(3) + t64, t2 * qJD(4) + t63, t2 * qJD(3) + t62 - t64; 0, t1 * qJD(4) - t63, 0, t1 * qJD(2) + t61 * t42; 0, t3 * qJD(3) - t62, t3 * qJD(2), 0;];
Cq  = t6;
