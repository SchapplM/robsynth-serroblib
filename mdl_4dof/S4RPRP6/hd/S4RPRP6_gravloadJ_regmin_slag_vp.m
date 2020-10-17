% Calculate minimal parameter regressor of gravitation load for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_gravloadJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:14
% EndTime: 2019-12-31 16:46:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (28->21), mult. (60->24), div. (0->0), fcn. (53->4), ass. (0->12)
t8 = sin(qJ(3));
t14 = pkin(3) * t8;
t11 = cos(qJ(1));
t9 = sin(qJ(1));
t13 = t11 * pkin(1) + t9 * qJ(2);
t1 = g(1) * t9 - g(2) * t11;
t2 = g(1) * t11 + g(2) * t9;
t10 = cos(qJ(3));
t12 = g(3) * t8 - t1 * t10;
t7 = -qJ(4) - pkin(5);
t4 = t11 * qJ(2);
t3 = [0, t1, t2, -t1, -t2, -g(1) * (-t9 * pkin(1) + t4) - g(2) * t13, 0, 0, 0, 0, 0, -t2 * t8, -t2 * t10, t1, -g(1) * (t11 * t14 + t4 + (-pkin(1) + t7) * t9) - g(2) * (-t11 * t7 + t9 * t14 + t13); 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, g(3) * t10 + t1 * t8, 0, t12 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg = t3;
