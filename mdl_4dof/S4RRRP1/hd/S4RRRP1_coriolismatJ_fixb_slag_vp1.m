% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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

function Cq = S4RRRP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:20
% EndTime: 2018-11-14 13:54:21
% DurationCPUTime: 0.43s
% Computational Cost: add. (4725->53), mult. (2725->78), div. (0->0), fcn. (1892->6), ass. (0->49)
t129 = m(5) / 0.2e1;
t130 = m(4) / 0.2e1;
t110 = sin(qJ(1)) * pkin(1);
t89 = qJ(1) + qJ(2);
t87 = sin(t89);
t112 = pkin(2) * t87;
t137 = rSges(5,1) + pkin(3);
t95 = qJ(3) + t89;
t86 = cos(t95);
t93 = sin(t95);
t71 = -t86 * rSges(5,2) - t137 * t93;
t65 = t71 - t112;
t59 = t65 - t110;
t109 = cos(qJ(1)) * pkin(1);
t88 = cos(t89);
t111 = pkin(2) * t88;
t72 = -t93 * rSges(5,2) + t137 * t86;
t66 = t72 + t111;
t60 = t66 + t109;
t25 = -t72 * t59 + t60 * t71;
t77 = -t93 * rSges(4,1) - t86 * rSges(4,2);
t73 = t77 - t112;
t69 = t73 - t110;
t78 = t86 * rSges(4,1) - t93 * rSges(4,2);
t74 = t78 + t111;
t70 = t74 + t109;
t31 = -t78 * t69 + t70 * t77;
t96 = t78 * t73 - t74 * t77;
t97 = t72 * t65 - t66 * t71;
t107 = (t96 + t31) * t130 + (t97 + t25) * t129;
t108 = (-t96 + t31) * t130 + (-t97 + t25) * t129;
t138 = t107 - t108;
t139 = t138 * qJD(1);
t127 = m(3) * (t109 * (-t87 * rSges(3,1) - t88 * rSges(3,2)) + (t88 * rSges(3,1) - t87 * rSges(3,2)) * t110);
t23 = -t66 * t59 + t60 * t65;
t28 = -t74 * t69 + t70 * t73;
t135 = 0.4e1 * qJD(1);
t132 = 2 * qJD(3);
t123 = m(4) * t28;
t121 = m(4) * t31;
t120 = m(4) * t96;
t117 = m(5) * t23;
t115 = m(5) * t25;
t114 = m(5) * t97;
t14 = -t114 - t120;
t13 = t115 + t121;
t4 = t117 + t123 + t127;
t1 = t107 + t108;
t2 = [t4 * qJD(2) + t13 * qJD(3), t4 * qJD(1) + t1 * qJD(3) + 0.2e1 * (t23 * t129 + t28 * t130 + t127 / 0.2e1) * qJD(2), t13 * qJD(1) + t1 * qJD(2) + (t25 * t129 + t31 * t130) * t132, 0; -t138 * qJD(3) + (-t127 / 0.4e1 - t123 / 0.4e1 - t117 / 0.4e1) * t135, t14 * qJD(3), -t139 + t14 * qJD(2) + (-t129 * t97 - t130 * t96) * t132, 0; t138 * qJD(2) + (-t121 / 0.4e1 - t115 / 0.4e1) * t135, t139 + 0.4e1 * (t120 / 0.4e1 + t114 / 0.4e1) * qJD(2), 0, 0; 0, 0, 0, 0;];
Cq  = t2;
