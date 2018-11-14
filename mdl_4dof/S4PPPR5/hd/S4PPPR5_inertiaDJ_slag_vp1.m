% Calculate time derivative of joint inertia matrix for
% S4PPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:01
% EndTime: 2018-11-14 14:05:01
% DurationCPUTime: 0.07s
% Computational Cost: add. (23->8), mult. (78->22), div. (0->0), fcn. (58->4), ass. (0->11)
t10 = cos(qJ(4));
t7 = sin(pkin(5));
t8 = cos(pkin(5));
t9 = sin(qJ(4));
t5 = -t10 * t8 - t7 * t9;
t6 = t10 * t7 - t8 * t9;
t4 = t6 * qJD(4);
t3 = t5 * qJD(4);
t2 = -rSges(5,1) * t4 - rSges(5,2) * t3;
t1 = rSges(5,1) * t3 - rSges(5,2) * t4;
t11 = [0; 0; 0; 0; 0; 0; m(5) * t2; 0; m(5) * (t1 * t7 - t2 * t8); 0.2e1 * m(5) * ((rSges(5,1) * t5 - rSges(5,2) * t6) * t2 + (rSges(5,1) * t6 + rSges(5,2) * t5) * t1);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1) t11(2) t11(4) t11(7); t11(2) t11(3) t11(5) t11(8); t11(4) t11(5) t11(6) t11(9); t11(7) t11(8) t11(9) t11(10);];
Mq  = res;
