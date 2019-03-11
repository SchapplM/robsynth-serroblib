% Calculate time derivative of joint inertia matrix for
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:32:58
% EndTime: 2019-03-08 18:32:59
% DurationCPUTime: 0.19s
% Computational Cost: add. (689->56), mult. (462->74), div. (0->0), fcn. (256->6), ass. (0->45)
t47 = qJ(1) + qJ(2);
t42 = pkin(6) + t47;
t39 = sin(t42);
t40 = cos(t42);
t44 = cos(t47);
t41 = pkin(2) * t44;
t62 = rSges(5,1) + pkin(3);
t65 = rSges(5,3) + qJ(4);
t10 = t65 * t39 + t62 * t40 + t41;
t46 = qJD(1) + qJD(2);
t66 = t10 * t46;
t43 = sin(t47);
t58 = pkin(2) * t43;
t9 = -t62 * t39 + t65 * t40 - t58;
t63 = t46 * t9;
t24 = t44 * rSges(3,1) - t43 * rSges(3,2);
t61 = -t40 * rSges(4,1) + t39 * rSges(4,2) - t41;
t48 = sin(qJ(1));
t57 = t48 * pkin(1);
t53 = pkin(1) * qJD(1);
t52 = t48 * t53;
t49 = cos(qJ(1));
t51 = t49 * t53;
t20 = t24 * t46;
t23 = -t43 * rSges(3,1) - t44 * rSges(3,2);
t19 = t23 * t46;
t17 = -t39 * rSges(4,1) - t40 * rSges(4,2) - t58;
t12 = t61 * t46;
t11 = t17 * t46;
t3 = qJD(4) * t39 + t63;
t4 = qJD(4) * t40 - t66;
t45 = t49 * pkin(1);
t22 = t24 + t45;
t21 = t23 - t57;
t16 = -t20 - t51;
t15 = t19 - t52;
t14 = -t61 + t45;
t13 = t17 - t57;
t8 = t12 - t51;
t7 = t11 - t52;
t6 = t45 + t10;
t5 = t9 - t57;
t2 = t4 - t51;
t1 = t3 - t52;
t18 = [0.2e1 * m(3) * (t22 * t15 + t21 * t16) + 0.2e1 * m(4) * (t13 * t8 + t14 * t7) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); m(3) * (t24 * t15 + t23 * t16 + t19 * t22 - t20 * t21) + m(4) * (t11 * t14 + t12 * t13 + t17 * t8 - t61 * t7) + m(5) * (t10 * t1 + t9 * t2 + t3 * t6 + t4 * t5); 0.2e1 * m(3) * (t24 * t19 - t23 * t20) + 0.2e1 * m(4) * (-t11 * t61 + t17 * t12) + 0.2e1 * m(5) * (t10 * t3 + t9 * t4); 0; 0; 0; m(5) * ((t46 * t5 - t1) * t40 + (t46 * t6 + t2) * t39); m(5) * ((-t3 + t63) * t40 + (t4 + t66) * t39); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t18(1) t18(2) t18(4) t18(7); t18(2) t18(3) t18(5) t18(8); t18(4) t18(5) t18(6) t18(9); t18(7) t18(8) t18(9) t18(10);];
Mq  = res;
