% Calculate minimal parameter regressor of gravitation load for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% 
% Output:
% taug_reg [3x9]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S3RPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_gravloadJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3RPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_gravloadJ_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:28:55
% EndTime: 2019-05-04 18:28:55
% DurationCPUTime: 0.07s
% Computational Cost: add. (20->13), mult. (40->18), div. (0->0), fcn. (46->4), ass. (0->11)
t13 = cos(qJ(3));
t12 = sin(qJ(3));
t8 = sin(qJ(1));
t9 = cos(qJ(1));
t1 = -t8 * t12 - t9 * t13;
t2 = t9 * t12 - t8 * t13;
t11 = g(1) * t2 - g(2) * t1;
t10 = g(1) * t1 + g(2) * t2;
t4 = g(1) * t9 + g(2) * t8;
t3 = g(1) * t8 - g(2) * t9;
t5 = [0, t3, t4, t3, -t4, -g(1) * (-t8 * pkin(1) + t9 * qJ(2)) - g(2) * (t9 * pkin(1) + t8 * qJ(2)) 0, -t11, t10; 0, 0, 0, 0, 0, -t3, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t11, -t10;];
taug_reg  = t5;
