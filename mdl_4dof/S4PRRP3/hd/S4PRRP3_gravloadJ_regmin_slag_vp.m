% Calculate minimal parameter regressor of gravitation load for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% taug_reg [4x15]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:05
% EndTime: 2021-01-14 22:27:06
% DurationCPUTime: 0.08s
% Computational Cost: add. (57->16), mult. (60->18), div. (0->0), fcn. (55->4), ass. (0->14)
t10 = pkin(6) + qJ(2);
t8 = sin(t10);
t9 = cos(t10);
t6 = g(1) * t9 + g(2) * t8;
t5 = g(1) * t8 - g(2) * t9;
t12 = sin(qJ(3));
t13 = cos(qJ(3));
t1 = -g(3) * t13 + t6 * t12;
t11 = -qJ(4) - pkin(5);
t7 = t13 * pkin(3) + pkin(2);
t4 = t5 * t13;
t3 = t5 * t12;
t2 = g(3) * t12 + t6 * t13;
t14 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(1) * (-t9 * t11 - t8 * t7) - g(2) * (-t8 * t11 + t9 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t14;
