% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RPPR2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:26
% EndTime: 2018-11-14 13:47:27
% DurationCPUTime: 0.26s
% Computational Cost: add. (1110->37), mult. (1308->61), div. (0->0), fcn. (1308->6), ass. (0->34)
t51 = sin(qJ(1));
t52 = cos(qJ(1));
t58 = pkin(6) + qJ(4);
t56 = sin(t58);
t57 = cos(t58);
t35 = -t51 * t56 - t52 * t57;
t36 = -t51 * t57 + t52 * t56;
t28 = t35 * rSges(5,1) + t36 * rSges(5,2);
t55 = -t36 * rSges(5,1) + t35 * rSges(5,2);
t78 = t28 * t51 + t52 * t55;
t81 = m(5) * t78;
t11 = -t81 / 0.2e1;
t12 = t81 / 0.2e1;
t50 = sin(pkin(6));
t59 = cos(pkin(6));
t68 = t59 * pkin(3) + pkin(1) + pkin(2);
t19 = t68 * t52 + (pkin(3) * t50 + qJ(2)) * t51 - t28;
t49 = t52 * qJ(2);
t65 = t52 * t50;
t53 = pkin(3) * t65 - t68 * t51 + t49 - t55;
t5 = t19 * t55 - t28 * t53;
t60 = m(5) * qJD(4);
t80 = t5 * t60;
t79 = m(5) * t5 * qJD(1);
t75 = m(5) * (t19 * t51 + t52 * t53);
t73 = m(3) * ((t52 * rSges(3,3) + t49) * t52 + (rSges(3,3) + qJ(2)) * t51 ^ 2);
t37 = -t51 * t50 - t52 * t59;
t38 = -t51 * t59 + t65;
t72 = m(4) * (t52 * (t38 * rSges(4,1) - t37 * rSges(4,2) + t49) + (-t37 * rSges(4,1) - t38 * rSges(4,2) + t51 * qJ(2)) * t51);
t4 = t72 + t73 + t75;
t3 = 0.2e1 * t12;
t2 = t11 + t12;
t1 = 0.2e1 * t11;
t6 = [t4 * qJD(2) + t80, t4 * qJD(1) + t2 * qJD(4), 0, t2 * qJD(2) + t79 - t80; t1 * qJD(4) + 0.4e1 * (-t73 / 0.4e1 - t72 / 0.4e1 - t75 / 0.4e1) * qJD(1), 0, 0, t1 * qJD(1) + t78 * t60; 0, 0, 0, 0; t3 * qJD(2) - t79, t3 * qJD(1), 0, 0;];
Cq  = t6;
