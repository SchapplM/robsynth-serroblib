% Calculate joint inertia matrix for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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
% Mq [3x3]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_inertiaJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_inertiaJ_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_inertiaJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPP1_inertiaJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPP1_inertiaJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:04:51
% EndTime: 2019-03-08 18:04:51
% DurationCPUTime: 0.08s
% Computational Cost: add. (54->28), mult. (104->37), div. (0->0), fcn. (62->2), ass. (0->14)
t12 = sin(qJ(1));
t13 = cos(qJ(1));
t15 = t13 * pkin(1) + t12 * qJ(2);
t14 = rSges(4,3) + qJ(3);
t10 = t13 * qJ(2);
t8 = t12 ^ 2 + t13 ^ 2;
t7 = rSges(2,1) * t13 - t12 * rSges(2,2);
t6 = -t12 * rSges(2,1) - rSges(2,2) * t13;
t5 = m(4) * t8;
t4 = -rSges(3,2) * t13 + t12 * rSges(3,3) + t15;
t3 = rSges(3,3) * t13 + t10 + (rSges(3,2) - pkin(1)) * t12;
t2 = t12 * rSges(4,2) + t14 * t13 + t15;
t1 = rSges(4,2) * t13 + t10 + (-pkin(1) - t14) * t12;
t9 = [Icges(2,3) + Icges(3,1) + Icges(4,1) + m(2) * (t6 ^ 2 + t7 ^ 2) + m(3) * (t3 ^ 2 + t4 ^ 2) + m(4) * (t1 ^ 2 + t2 ^ 2); m(3) * (t12 * t3 - t13 * t4) + m(4) * (t12 * t1 - t13 * t2); m(3) * t8 + t5; m(4) * (t1 * t13 + t12 * t2); 0; t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t9(1) t9(2) t9(4); t9(2) t9(3) t9(5); t9(4) t9(5) t9(6);];
Mq  = res;
