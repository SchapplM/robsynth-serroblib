% Calculate minimal parameter regressor of gravitation load for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% taug_reg [3x9]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S3RRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_gravloadJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_gravloadJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:31:25
% EndTime: 2019-05-04 18:31:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (45->14), mult. (36->18), div. (0->0), fcn. (32->4), ass. (0->10)
t8 = qJ(1) + qJ(2);
t6 = sin(t8);
t7 = cos(t8);
t12 = t7 * pkin(2) + t6 * qJ(3);
t11 = -t6 * pkin(2) + t7 * qJ(3);
t10 = cos(qJ(1));
t9 = sin(qJ(1));
t2 = g(1) * t7 + g(2) * t6;
t1 = g(1) * t6 - g(2) * t7;
t3 = [0, g(1) * t9 - g(2) * t10, g(1) * t10 + g(2) * t9, 0, t1, t2, t1, -t2, -g(1) * (-t9 * pkin(1) + t11) - g(2) * (t10 * pkin(1) + t12); 0, 0, 0, 0, t1, t2, t1, -t2, -g(1) * t11 - g(2) * t12; 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t3;
