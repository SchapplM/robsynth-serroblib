% Calculate vector of inverse dynamics joint torques for
% S4PRPR1
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4PRPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:12
% EndTime: 2018-11-14 13:42:13
% DurationCPUTime: 0.57s
% Computational Cost: add. (998->108), mult. (968->128), div. (0->0), fcn. (698->4), ass. (0->57)
t61 = pkin(6) + qJ(2);
t58 = sin(t61);
t59 = cos(t61);
t36 = t59 * pkin(2) + t58 * qJ(3);
t51 = qJD(3) * t59;
t24 = t36 * qJD(2) - t51;
t53 = t59 * qJ(3);
t33 = t58 * pkin(2) - t53;
t62 = qJD(2) - qJD(4);
t63 = sin(qJ(4));
t80 = cos(qJ(4));
t65 = t58 * t63 + t59 * t80;
t16 = t62 * t65;
t68 = t58 * t80;
t71 = qJD(2) * t59;
t87 = -t59 * t63 + t68;
t17 = -qJD(2) * t68 + t87 * qJD(4) + t63 * t71;
t5 = t16 * rSges(5,1) - t17 * rSges(5,2);
t60 = qJDD(2) - qJDD(4);
t64 = qJD(2) ^ 2;
t70 = qJD(2) * qJD(3);
t76 = qJDD(3) * t58 + t59 * t70;
t78 = -rSges(5,1) * t87 + rSges(5,2) * t65;
t1 = -qJD(2) * t24 - qJDD(2) * t33 + t60 * t78 - t62 * t5 + (-qJDD(2) * t58 - t59 * t64) * pkin(3) + t76;
t92 = -g(1) + t1;
t18 = -rSges(5,1) * t65 - rSges(5,2) * t87;
t6 = t17 * rSges(5,1) + t16 * rSges(5,2);
t72 = qJD(2) * t58;
t50 = qJD(3) * t58;
t75 = qJ(3) * t71 + t50;
t66 = -qJDD(3) * t59 + qJD(2) * (-pkin(2) * t72 + t75) + qJDD(2) * t36 + t58 * t70;
t2 = -t60 * t18 + t62 * t6 + (qJDD(2) * t59 - t58 * t64) * pkin(3) + t66;
t91 = -g(2) + t2;
t90 = t62 * t78;
t88 = t59 * pkin(3) + t36;
t37 = t59 * rSges(4,1) + t58 * rSges(4,3);
t86 = t58 / 0.2e1;
t85 = -t59 / 0.2e1;
t84 = -m(4) - m(5);
t83 = -pkin(2) - pkin(3);
t81 = -rSges(4,1) - pkin(2);
t55 = t59 * rSges(4,3);
t34 = t58 * rSges(4,1) - t55;
t77 = -t33 - t34;
t22 = t36 + t37;
t74 = -qJD(2) * t33 + t50;
t73 = Icges(5,3) * t60;
t38 = t59 * rSges(3,1) - t58 * rSges(3,2);
t35 = t58 * rSges(3,1) + t59 * rSges(3,2);
t48 = rSges(4,3) * t71;
t15 = t22 * qJD(2) - t51;
t14 = t77 * qJD(2) + t50;
t10 = t88 * qJD(2) - t18 * t62 - t51;
t9 = t90 + t50 + (-t58 * pkin(3) - t33) * qJD(2);
t4 = qJDD(2) * t37 + qJD(2) * (-rSges(4,1) * t72 + t48) + t66;
t3 = t77 * qJDD(2) + (-t37 * qJD(2) - t24) * qJD(2) + t76;
t7 = [(-g(3) + qJDD(1)) * (m(2) + m(3) - t84); -m(3) * (-g(1) * t35 + g(2) * t38) + t73 + (Icges(3,3) + Icges(4,2) + m(3) * (t35 ^ 2 + t38 ^ 2)) * qJDD(2) + (-(-pkin(3) * t72 + t74 - t9 + t90) * t10 + t9 * (-t5 + t51) + t10 * (t6 + t75) + (t9 * t83 * t59 + (-t9 * qJ(3) + t10 * t83) * t58) * qJD(2) + t91 * (-t18 + t88) + t92 * (t83 * t58 + t53 + t78)) * m(5) + (-(-qJD(2) * t34 - t14 + t74) * t15 + t14 * t51 + t15 * (t48 + t75) + (t14 * t81 * t59 + (t14 * (-rSges(4,3) - qJ(3)) + t15 * t81) * t58) * qJD(2) + (-g(2) + t4) * t22 + (-g(1) + t3) * (t81 * t58 + t53 + t55)) * m(4); t84 * (g(1) * t58 - g(2) * t59) + 0.2e1 * (t1 * t86 + t2 * t85) * m(5) + 0.2e1 * (t3 * t86 + t4 * t85) * m(4); -t73 + (-t10 * t6 + t9 * t5 - (-t10 * t62 + t92) * t78 + (t9 * t62 + t91) * t18) * m(5);];
tau  = t7;
