% Calculate vector of centrifugal and Coriolis load on the joints for
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:34
% EndTime: 2019-01-31 13:16:36
% DurationCPUTime: 0.49s
% Computational Cost: add. (988->90), mult. (739->107), div. (0->0), fcn. (332->8), ass. (0->68)
t48 = qJ(1) + qJ(2);
t42 = pkin(7) + t48;
t40 = qJ(4) + t42;
t34 = sin(t40);
t35 = cos(t40);
t19 = t34 * rSges(5,1) + t35 * rSges(5,2);
t47 = qJD(1) + qJD(2);
t41 = qJD(4) + t47;
t13 = t41 * t19;
t38 = cos(t42);
t44 = cos(t48);
t39 = pkin(2) * t44;
t60 = pkin(3) * t38 + t39;
t50 = cos(qJ(1));
t68 = pkin(1) * qJD(1);
t64 = t50 * t68;
t31 = t35 * rSges(5,1);
t20 = -t34 * rSges(5,2) + t31;
t71 = t41 * t20;
t4 = t60 * t47 + t64 + t71;
t81 = t4 * t13;
t32 = t38 * rSges(4,1);
t80 = -t32 - t39;
t49 = sin(qJ(1));
t65 = t49 * t68;
t43 = sin(t48);
t24 = t43 * rSges(3,1) + t44 * rSges(3,2);
t69 = t47 * t24;
t14 = -t65 - t69;
t37 = sin(t42);
t72 = t37 * rSges(4,2);
t79 = -t72 - t80;
t78 = pkin(2) * t43;
t61 = -pkin(3) * t37 - t78;
t3 = t61 * t47 - t13 - t65;
t46 = t47 ^ 2;
t77 = pkin(2) * t46;
t76 = t19 * t4;
t75 = t49 * pkin(1);
t45 = t50 * pkin(1);
t74 = t34 * t41;
t73 = t35 * t41;
t70 = t43 * t47;
t36 = t44 * rSges(3,1);
t51 = qJD(1) ^ 2;
t67 = t51 * t75;
t66 = t51 * t45;
t21 = t37 * rSges(4,1) + t38 * rSges(4,2);
t56 = -t21 - t78;
t25 = -t43 * rSges(3,2) + t36;
t18 = -rSges(3,2) * t70 + t47 * t36;
t23 = rSges(5,2) * t74;
t12 = rSges(5,1) * t73 - t23;
t58 = t20 + t60;
t55 = -t19 + t61;
t7 = t56 * t47 - t65;
t8 = t47 * t79 + t64;
t52 = (t8 * t56 + t7 * t80) * t47;
t27 = t47 * t72;
t16 = t47 * t21;
t15 = t47 * t25 + t64;
t10 = -t47 * t18 - t66;
t9 = -t47 * t69 - t67;
t6 = -t66 - t44 * t77 - t47 * (t47 * t32 - t27);
t5 = -t21 * t46 - t43 * t77 - t67;
t2 = -t41 * t12 - t60 * t46 - t66;
t1 = -t13 * t41 + t61 * t46 - t67;
t11 = [m(3) * (t10 * (-t24 - t75) + t9 * (t25 + t45) + (-t18 - t64 + t15) * t14) + (t2 * (t55 - t75) + t3 * (-t12 - t64) + t1 * (t45 + t58) + (-t3 * t60 + t4 * t61) * t47 + t4 * (-rSges(5,1) * t74 - rSges(5,2) * t73 - t65)) * m(5) + (t6 * (t56 - t75) + t7 * (t27 - t64) + t5 * (t45 + t79) + t52 - (-pkin(2) * t70 - t16 - t7) * t8) * m(4); (t8 * t16 - (-t7 * t79 - t8 * t78) * t47 + t7 * t27 + t5 * t79 + t6 * t56 + t52) * m(4) + (-(-t14 * t25 - t15 * t24) * t47 - t10 * t24 - t14 * t18 - t15 * t69 + t9 * t25) * m(3) + (t1 * t58 + t81 + t2 * t55 - t76 * t41 + (-t31 * t41 + t23 + t71) * t3) * m(5); 0; (t1 * t20 - t81 - t3 * t12 - t2 * t19 - (-t20 * t3 - t76) * t41) * m(5);];
tauc  = t11(:);
