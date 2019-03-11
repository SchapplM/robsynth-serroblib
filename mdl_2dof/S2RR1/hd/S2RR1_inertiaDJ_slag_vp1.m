% Calculate time derivative of joint inertia matrix for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [2x2]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S2RR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_inertiaDJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_inertiaDJ_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_inertiaDJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_inertiaDJ_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_inertiaDJ_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_inertiaDJ_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 17:59:55
% EndTime: 2019-03-08 17:59:56
% DurationCPUTime: 0.77s
% Computational Cost: add. (370->89), mult. (1022->151), div. (0->0), fcn. (786->4), ass. (0->59)
t37 = sin(qJ(2));
t92 = qJD(1) * t37;
t39 = cos(qJ(2));
t91 = qJD(1) * t39;
t52 = rSges(3,1) * t37 + rSges(3,2) * t39;
t90 = qJD(2) * t52;
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t67 = Icges(3,4) * t39;
t46 = -Icges(3,2) * t37 + t67;
t18 = Icges(3,6) * t38 - t46 * t40;
t68 = Icges(3,4) * t37;
t48 = Icges(3,1) * t39 - t68;
t20 = Icges(3,5) * t38 - t48 * t40;
t49 = t18 * t37 - t20 * t39;
t89 = t49 * t38;
t82 = Icges(3,5) * t40 + t48 * t38;
t83 = Icges(3,6) * t40 + t46 * t38;
t50 = -t37 * t83 + t39 * t82;
t88 = t50 * t40;
t87 = t38 * t90;
t74 = rSges(3,2) * t37;
t75 = rSges(3,1) * t39;
t53 = -t74 + t75;
t86 = t53 * t38;
t85 = t38 * rSges(3,3) - t40 * t75;
t44 = Icges(3,5) * t39 - Icges(3,6) * t37;
t84 = Icges(3,3) * t40 + t44 * t38;
t81 = 2 * m(3);
t80 = t38 ^ 2;
t79 = t40 ^ 2;
t78 = m(3) * t52;
t77 = rSges(3,3) + pkin(1);
t72 = t37 * t82;
t71 = t37 * t20;
t70 = t39 * t83;
t69 = t39 * t18;
t16 = Icges(3,3) * t38 - t44 * t40;
t63 = qJD(1) * t16;
t62 = qJD(1) * t38;
t61 = qJD(1) * t40;
t60 = qJD(2) * t37;
t58 = qJD(2) * t39;
t55 = rSges(3,3) * t61 + t40 * t90 + t62 * t75;
t22 = t40 * t74 + t85;
t41 = qJD(2) * (Icges(3,5) * t37 + Icges(3,6) * t39);
t26 = t53 * qJD(2);
t21 = -t40 * rSges(3,3) - t86;
t14 = t77 * t40 + t86;
t13 = t38 * pkin(1) + t22;
t8 = t38 * t41 + t63;
t7 = t84 * qJD(1) + t40 * t41;
t6 = (pkin(1) * t40 - t38 * t74) * qJD(1) + t55;
t5 = -t87 + (-t77 * t38 + t53 * t40) * qJD(1);
t4 = t38 * t16 + t49 * t40;
t3 = -t38 * t84 + t88;
t2 = -t16 * t40 + t89;
t1 = t50 * t38 + t40 * t84;
t9 = [(t13 * t6 + t14 * t5) * t81 + (-Icges(3,2) * t39 + t48 - t68) * t60 + (Icges(3,1) * t37 + t46 + t67) * t58; m(3) * (-(-t38 * t6 - t40 * t5) * t52 - (-t13 * t38 - t14 * t40) * t26) - (t50 * qJD(2) - t18 * t91 - t20 * t92) * t40 / 0.2e1 + (t49 * qJD(2) - t82 * t92 - t83 * t91) * t38 / 0.2e1 - (t79 / 0.2e1 + t80 / 0.2e1) * t44 * qJD(2) + ((t13 * t78 - t69 / 0.2e1 - t71 / 0.2e1) * t40 + (-t14 * t78 + t70 / 0.2e1 + t72 / 0.2e1) * t38) * qJD(1); ((t38 * t21 + t22 * t40) * ((qJD(1) * t21 + t55) * t40 + (t87 + (-t22 + t85) * qJD(1)) * t38) + (t79 + t80) * t52 * t26) * t81 + (-t1 * t40 + t2 * t38) * t62 - t40 * ((t40 * t8 + (t2 - t88) * qJD(1)) * t40 + (t1 * qJD(1) + (t18 * t58 + t20 * t60 + t63) * t38 + (-t7 + (t70 + t72) * qJD(2) + t49 * qJD(1)) * t40) * t38) + (-t3 * t40 + t4 * t38) * t61 + t38 * ((t38 * t7 + (t3 - t89) * qJD(1)) * t38 + (t4 * qJD(1) + (t58 * t83 + t60 * t82) * t40 + (-t8 + (t69 + t71) * qJD(2) + (t16 + t50) * qJD(1)) * t38) * t40);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t9(1) t9(2); t9(2) t9(3);];
Mq  = res;
