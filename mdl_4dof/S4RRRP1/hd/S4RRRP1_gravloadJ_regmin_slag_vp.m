% Calculate minimal parameter regressor of gravitation load for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From joint_gravload_fixb_regressor_minpar_matlab.m
t12 = qJ(1) + qJ(2);
t10 = cos(t12);
t11 = qJ(3) + t12;
t8 = cos(t11);
t16 = pkin(2) * t10 + pkin(3) * t8;
t7 = sin(t11);
t9 = sin(t12);
t15 = -pkin(2) * t9 - pkin(3) * t7;
t1 = g(1) * t7 - g(2) * t8;
t14 = cos(qJ(1));
t13 = sin(qJ(1));
t4 = g(1) * t10 + g(2) * t9;
t3 = g(1) * t9 - g(2) * t10;
t2 = g(1) * t8 + g(2) * t7;
t5 = [0, g(1) * t13 - g(2) * t14, g(1) * t14 + g(2) * t13, 0, t3, t4, 0, t1, t2, -g(1) * (-t13 * pkin(1) + t15) - g(2) * (pkin(1) * t14 + t16); 0, 0, 0, 0, t3, t4, 0, t1, t2, -g(1) * t15 - g(2) * t16; 0, 0, 0, 0, 0, 0, 0, t1, t2, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3);];
taug_reg  = t5;
