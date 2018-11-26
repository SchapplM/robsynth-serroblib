% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:29
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRRP9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:28:53
% EndTime: 2018-11-23 16:28:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (130->58), mult. (129->56), div. (0->0), fcn. (184->8), ass. (0->40)
t27 = -pkin(9) - pkin(8);
t24 = cos(qJ(4));
t9 = t24 * pkin(4) + pkin(3);
t21 = sin(qJ(4));
t23 = sin(qJ(1));
t44 = t23 * t21;
t22 = sin(qJ(3));
t43 = t23 * t22;
t42 = t23 * t24;
t25 = cos(qJ(3));
t41 = t23 * t25;
t20 = qJ(4) + qJ(5);
t10 = sin(t20);
t40 = t25 * t10;
t26 = cos(qJ(1));
t39 = t26 * t21;
t38 = t26 * t22;
t37 = t26 * t24;
t19 = pkin(6) + 0;
t36 = t23 * pkin(1) + 0;
t35 = pkin(2) + t19;
t13 = t23 * pkin(7);
t34 = t13 + t36;
t33 = t26 * pkin(1) + t23 * qJ(2) + 0;
t32 = t26 * pkin(7) + t33;
t31 = pkin(3) * t22 - pkin(8) * t25;
t18 = -qJ(6) + t27;
t11 = cos(t20);
t5 = pkin(5) * t11 + t9;
t30 = t18 * t25 + t22 * t5;
t29 = t22 * t9 + t25 * t27;
t28 = -t26 * qJ(2) + t36;
t8 = t26 * t25;
t7 = t25 * t11;
t6 = t21 * pkin(4) + pkin(5) * t10;
t4 = t23 * t10 - t11 * t38;
t3 = t10 * t38 + t23 * t11;
t2 = t26 * t10 + t11 * t43;
t1 = -t10 * t43 + t26 * t11;
t12 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t23, 0, 0; t23, t26, 0, 0; 0, 0, 1, t19; 0, 0, 0, 1; 0, -t26, t23, t33; 0, -t23, -t26, t28; 1, 0, 0, t19; 0, 0, 0, 1; t43, t41, t26, t32; -t38, -t8, t23, t13 + t28; t25, -t22, 0, t35; 0, 0, 0, 1; t22 * t42 + t39, -t21 * t43 + t37, -t41, t31 * t23 + t32; -t22 * t37 + t44, t21 * t38 + t42, t8 (-qJ(2) - t31) * t26 + t34; t25 * t24, -t25 * t21, t22, t25 * pkin(3) + t22 * pkin(8) + t35; 0, 0, 0, 1; t2, t1, -t41, pkin(4) * t39 + t29 * t23 + t32; t4, t3, t8, pkin(4) * t44 + (-qJ(2) - t29) * t26 + t34; t7, -t40, t22, -t22 * t27 + t25 * t9 + t35; 0, 0, 0, 1; t2, t1, -t41, t30 * t23 + t26 * t6 + t32; t4, t3, t8, t23 * t6 + (-qJ(2) - t30) * t26 + t34; t7, -t40, t22, -t22 * t18 + t25 * t5 + t35; 0, 0, 0, 1;];
T_ges = t12;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
