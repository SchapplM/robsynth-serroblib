% Calculate minimal parameter regressor of gravitation load for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x12]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:26
% EndTime: 2019-12-31 16:17:26
% DurationCPUTime: 0.08s
% Computational Cost: add. (23->10), mult. (48->16), div. (0->0), fcn. (60->6), ass. (0->11)
t13 = cos(qJ(3));
t12 = sin(qJ(3));
t11 = sin(pkin(6));
t6 = cos(pkin(6));
t1 = -t11 * t12 - t6 * t13;
t2 = -t11 * t13 + t6 * t12;
t10 = g(1) * t2 - g(2) * t1;
t9 = -g(1) * t1 - g(2) * t2;
t8 = cos(qJ(4));
t7 = sin(qJ(4));
t3 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(1) * t11 + g(2) * t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, t10 * t8, -t10 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t8 + t9 * t7, -g(3) * t7 + t9 * t8;];
taug_reg = t3;
