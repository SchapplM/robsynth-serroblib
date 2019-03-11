% Calculate vector of inverse dynamics joint torques for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% qJDD [3x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [3x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S3RPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp1: qJDD has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPR1_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:05:51
% EndTime: 2019-03-08 18:05:52
% DurationCPUTime: 0.52s
% Computational Cost: add. (493->104), mult. (960->127), div. (0->0), fcn. (698->4), ass. (0->55)
t61 = sin(qJ(1));
t62 = cos(qJ(1));
t36 = t62 * pkin(1) + t61 * qJ(2);
t51 = qJD(2) * t62;
t24 = t36 * qJD(1) - t51;
t53 = t62 * qJ(2);
t33 = t61 * pkin(1) - t53;
t59 = qJD(1) - qJD(3);
t60 = sin(qJ(3));
t78 = cos(qJ(3));
t64 = t61 * t60 + t62 * t78;
t14 = t59 * t64;
t67 = t61 * t78;
t69 = qJD(1) * t62;
t84 = -t62 * t60 + t67;
t15 = -qJD(1) * t67 + t84 * qJD(3) + t60 * t69;
t5 = t14 * rSges(4,1) - t15 * rSges(4,2);
t58 = qJDD(1) - qJDD(3);
t63 = qJD(1) ^ 2;
t68 = qJD(1) * qJD(2);
t74 = qJDD(2) * t61 + t62 * t68;
t76 = -rSges(4,1) * t84 + rSges(4,2) * t64;
t1 = -qJD(1) * t24 - qJDD(1) * t33 + t58 * t76 - t59 * t5 + (-qJDD(1) * t61 - t62 * t63) * pkin(2) + t74;
t89 = -g(1) + t1;
t16 = -rSges(4,1) * t64 - rSges(4,2) * t84;
t6 = t15 * rSges(4,1) + t14 * rSges(4,2);
t70 = qJD(1) * t61;
t50 = qJD(2) * t61;
t73 = qJ(2) * t69 + t50;
t65 = -qJDD(2) * t62 + qJD(1) * (-pkin(1) * t70 + t73) + qJDD(1) * t36 + t61 * t68;
t2 = -t58 * t16 + t59 * t6 + (qJDD(1) * t62 - t61 * t63) * pkin(2) + t65;
t88 = -g(2) + t2;
t87 = t59 * t76;
t85 = t62 * pkin(2) + t36;
t37 = t62 * rSges(3,1) + t61 * rSges(3,3);
t83 = t61 / 0.2e1;
t82 = -t62 / 0.2e1;
t81 = -pkin(1) - pkin(2);
t79 = -rSges(3,1) - pkin(1);
t55 = t62 * rSges(3,3);
t34 = t61 * rSges(3,1) - t55;
t75 = -t33 - t34;
t23 = t36 + t37;
t72 = -qJD(1) * t33 + t50;
t71 = Icges(4,3) * t58;
t38 = t62 * rSges(2,1) - t61 * rSges(2,2);
t35 = t61 * rSges(2,1) + t62 * rSges(2,2);
t48 = rSges(3,3) * t69;
t20 = qJD(1) * t23 - t51;
t19 = qJD(1) * t75 + t50;
t10 = t85 * qJD(1) - t59 * t16 - t51;
t9 = t87 + t50 + (-t61 * pkin(2) - t33) * qJD(1);
t4 = qJDD(1) * t37 + qJD(1) * (-rSges(3,1) * t70 + t48) + t65;
t3 = t75 * qJDD(1) + (-t37 * qJD(1) - t24) * qJD(1) + t74;
t7 = [-m(2) * (-g(1) * t35 + g(2) * t38) + t71 + (Icges(2,3) + Icges(3,2) + m(2) * (t35 ^ 2 + t38 ^ 2)) * qJDD(1) + (-(-pkin(2) * t70 + t72 + t87 - t9) * t10 + t9 * (-t5 + t51) + t10 * (t6 + t73) + (t9 * t81 * t62 + (-t9 * qJ(2) + t10 * t81) * t61) * qJD(1) + t88 * (-t16 + t85) + t89 * (t81 * t61 + t53 + t76)) * m(4) + (-(-qJD(1) * t34 - t19 + t72) * t20 + t19 * t51 + t20 * (t48 + t73) + (t19 * t79 * t62 + (t19 * (-rSges(3,3) - qJ(2)) + t20 * t79) * t61) * qJD(1) + (-g(2) + t4) * t23 + (-g(1) + t3) * (t79 * t61 + t53 + t55)) * m(3); (-m(3) - m(4)) * (g(1) * t61 - g(2) * t62) + 0.2e1 * (t1 * t83 + t2 * t82) * m(4) + 0.2e1 * (t3 * t83 + t4 * t82) * m(3); -t71 + (-t10 * t6 + t9 * t5 - (-t10 * t59 + t89) * t76 + (t59 * t9 + t88) * t16) * m(4);];
tau  = t7;
