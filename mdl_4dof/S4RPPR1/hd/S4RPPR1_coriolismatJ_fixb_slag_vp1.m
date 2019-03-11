% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RPPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:30
% EndTime: 2019-03-08 18:27:31
% DurationCPUTime: 0.20s
% Computational Cost: add. (1212->28), mult. (1126->43), div. (0->0), fcn. (1116->6), ass. (0->28)
t41 = qJ(1) + pkin(6);
t39 = sin(t41);
t40 = cos(t41);
t56 = sin(qJ(4));
t57 = cos(qJ(4));
t31 = -t39 * t56 - t40 * t57;
t32 = -t39 * t57 + t40 * t56;
t22 = t31 * rSges(5,1) + t32 * rSges(5,2);
t45 = -t32 * rSges(5,1) + t31 * rSges(5,2);
t67 = t22 * t39 + t40 * t45;
t71 = m(5) * t67;
t10 = -t71 / 0.2e1;
t11 = t71 / 0.2e1;
t47 = m(5) * qJD(4);
t59 = cos(qJ(1)) * pkin(1);
t63 = pkin(2) + pkin(3);
t19 = t39 * qJ(3) + t63 * t40 - t22 + t59;
t46 = -sin(qJ(1)) * pkin(1) + t40 * qJ(3);
t44 = -t63 * t39 - t45 + t46;
t5 = t19 * t45 - t22 * t44;
t70 = t5 * t47;
t4 = m(4) * ((t40 * rSges(4,3) + t46) * t40 + (t59 + (rSges(4,3) + qJ(3)) * t39) * t39) + m(5) * (t19 * t39 + t40 * t44);
t69 = t4 * qJD(1);
t68 = m(5) * t5 * qJD(1);
t3 = 0.2e1 * t11;
t2 = t10 + t11;
t1 = 0.2e1 * t10;
t6 = [t4 * qJD(3) + t70, 0, t2 * qJD(4) + t69, t2 * qJD(3) + t68 - t70; 0, 0, 0, 0; t1 * qJD(4) - t69, 0, 0, t1 * qJD(1) + t67 * t47; t3 * qJD(3) - t68, 0, t3 * qJD(1), 0;];
Cq  = t6;
