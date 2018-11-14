% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2018-11-14 13:54
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RRPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:53:18
% EndTime: 2018-11-14 13:53:19
% DurationCPUTime: 0.30s
% Computational Cost: add. (3399->41), mult. (1707->55), div. (0->0), fcn. (1172->8), ass. (0->39)
t69 = qJ(1) + qJ(2);
t66 = pkin(7) + t69;
t73 = qJ(4) + t66;
t62 = cos(t73);
t72 = sin(t73);
t52 = -t72 * rSges(5,1) - t62 * rSges(5,2);
t64 = sin(t66);
t67 = sin(t69);
t87 = pkin(2) * t67;
t42 = -pkin(3) * t64 + t52 - t87;
t53 = t62 * rSges(5,1) - t72 * rSges(5,2);
t65 = cos(t66);
t68 = cos(t69);
t86 = pkin(2) * t68;
t43 = pkin(3) * t65 + t53 + t86;
t77 = t53 * t42 - t43 * t52;
t109 = m(5) * qJD(2) * t77;
t85 = sin(qJ(1)) * pkin(1);
t40 = t42 - t85;
t84 = cos(qJ(1)) * pkin(1);
t41 = t43 + t84;
t101 = t53 * t40 - t41 * t52;
t108 = t101 * m(5) * qJD(1);
t107 = m(3) * (t84 * (-t67 * rSges(3,1) - t68 * rSges(3,2)) + (t68 * rSges(3,1) - t67 * rSges(3,2)) * t85);
t106 = m(4) * (t84 * (-t64 * rSges(4,1) - t65 * rSges(4,2) - t87) + (t65 * rSges(4,1) - t64 * rSges(4,2) + t86) * t85);
t105 = m(5) * (-t43 * t40 + t41 * t42);
t74 = m(5) * qJD(4);
t104 = t101 * t74;
t103 = t77 * t74;
t4 = t105 + t106 + t107;
t102 = t4 * qJD(1);
t97 = m(5) * (-t77 - t101);
t92 = m(5) * (t77 - t101);
t8 = t92 / 0.2e1;
t7 = t97 / 0.2e1;
t3 = t8 - t97 / 0.2e1;
t2 = t7 - t92 / 0.2e1;
t1 = t7 + t8;
t5 = [t4 * qJD(2) - t104, t102 + t1 * qJD(4) + 0.2e1 * (t107 / 0.2e1 + t105 / 0.2e1 + t106 / 0.2e1) * qJD(2), 0, t1 * qJD(2) - t104 - t108; t2 * qJD(4) - t102, -t103, 0, t2 * qJD(1) - t103 - t109; 0, 0, 0, 0; t3 * qJD(2) + t108, t3 * qJD(1) + t109, 0, 0;];
Cq  = t5;
