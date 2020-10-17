% Calculate minimal parameter regressor of gravitation load for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% 
% Output:
% taug_reg [3x7]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S3PRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_gravloadJ_regmin_slag_vp: qJ has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S3PRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_gravloadJ_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 18:23:54
% EndTime: 2019-05-04 18:23:54
% DurationCPUTime: 0.05s
% Computational Cost: add. (10->9), mult. (16->10), div. (0->0), fcn. (14->2), ass. (0->5)
t4 = cos(qJ(2));
t3 = sin(qJ(2));
t2 = g(1) * t4 + g(2) * t3;
t1 = g(1) * t3 - g(2) * t4;
t5 = [-g(2), 0, 0, 0, 0, 0, -g(2); 0, 0, t1, t2, t1, -t2, -g(1) * (-pkin(2) * t3 + qJ(3) * t4) - g(2) * (pkin(2) * t4 + qJ(3) * t3); 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
