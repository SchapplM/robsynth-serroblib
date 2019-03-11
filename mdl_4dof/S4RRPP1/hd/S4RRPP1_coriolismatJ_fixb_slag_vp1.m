% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:32:59
% EndTime: 2019-03-08 18:32:59
% DurationCPUTime: 0.28s
% Computational Cost: add. (2470->37), mult. (1382->53), div. (0->0), fcn. (1024->6), ass. (0->34)
t61 = qJ(1) + qJ(2);
t59 = sin(t61);
t60 = cos(t61);
t72 = cos(qJ(1)) * pkin(1);
t73 = sin(qJ(1)) * pkin(1);
t94 = m(3) * (t72 * (-t59 * rSges(3,1) - t60 * rSges(3,2)) + (t60 * rSges(3,1) - t59 * rSges(3,2)) * t73);
t58 = pkin(6) + t61;
t56 = sin(t58);
t57 = cos(t58);
t74 = pkin(2) * t60;
t75 = pkin(2) * t59;
t93 = m(4) * (t72 * (-t56 * rSges(4,1) - t57 * rSges(4,2) - t75) + (t57 * rSges(4,1) - t56 * rSges(4,2) + t74) * t73);
t71 = rSges(5,1) + pkin(3);
t89 = rSges(5,3) + qJ(4);
t36 = -t71 * t56 + t89 * t57 - t75;
t34 = t36 - t73;
t37 = t89 * t56 + t71 * t57 + t74;
t35 = t37 + t72;
t92 = m(5) * (-t37 * t34 + t35 * t36);
t18 = t36 * t57 + t37 * t56;
t91 = m(5) * t18 * qJD(2);
t4 = t92 + t93 + t94;
t90 = t4 * qJD(1);
t14 = t34 * t57 + t35 * t56;
t88 = m(5) * t14 * qJD(1);
t84 = m(5) * (t14 - t18);
t79 = m(5) * (t18 + t14);
t64 = m(5) * qJD(4);
t8 = t79 / 0.2e1;
t7 = t84 / 0.2e1;
t3 = t8 - t84 / 0.2e1;
t2 = t8 + t7;
t1 = t7 - t79 / 0.2e1;
t5 = [t4 * qJD(2) + t14 * t64, t90 + t2 * qJD(4) + 0.2e1 * (t92 / 0.2e1 + t93 / 0.2e1 + t94 / 0.2e1) * qJD(2), 0, t2 * qJD(2) + t88; t3 * qJD(4) - t90, t18 * t64, 0, t3 * qJD(1) + t91; 0, 0, 0, 0; t1 * qJD(2) - t88, t1 * qJD(1) - t91, 0, 0;];
Cq  = t5;
