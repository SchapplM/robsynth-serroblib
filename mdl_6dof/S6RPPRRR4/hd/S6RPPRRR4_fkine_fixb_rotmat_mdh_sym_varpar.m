% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 15:49
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPPRRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:49:02
% EndTime: 2018-11-23 15:49:02
% DurationCPUTime: 0.14s
% Computational Cost: add. (148->54), mult. (207->58), div. (0->0), fcn. (305->10), ass. (0->35)
t44 = cos(qJ(1));
t43 = sin(qJ(1));
t20 = qJ(5) + qJ(6);
t12 = sin(t20);
t24 = cos(qJ(4));
t42 = t12 * t24;
t13 = cos(t20);
t41 = t13 * t24;
t21 = sin(qJ(5));
t40 = t21 * t24;
t23 = cos(qJ(5));
t39 = t23 * t24;
t38 = cos(pkin(10));
t37 = sin(pkin(10));
t19 = pkin(6) + 0;
t36 = pkin(5) * t21 + pkin(7);
t11 = -qJ(3) + t19;
t35 = t44 * pkin(1) + t43 * qJ(2) + 0;
t34 = pkin(2) * t44 + t35;
t22 = sin(qJ(4));
t33 = -pkin(4) * t24 - pkin(8) * t22;
t10 = t23 * pkin(5) + pkin(4);
t25 = -pkin(9) - pkin(8);
t32 = -t10 * t24 + t22 * t25;
t5 = -t43 * t37 - t38 * t44;
t31 = -t5 * pkin(3) + t34;
t30 = t43 * pkin(1) - qJ(2) * t44 + 0;
t6 = t37 * t44 - t43 * t38;
t29 = t6 * pkin(7) + t31;
t28 = t43 * pkin(2) + t30;
t27 = -t6 * pkin(3) + t28;
t26 = -t5 * pkin(7) + t27;
t2 = t5 * t22;
t1 = t6 * t22;
t3 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t44, -t43, 0, 0; t43, t44, 0, 0; 0, 0, 1, t19; 0, 0, 0, 1; t44, 0, t43, t35; t43, 0, -t44, t30; 0, 1, 0, t19; 0, 0, 0, 1; -t5, -t6, 0, t34; -t6, t5, 0, t28; 0, 0, -1, t11; 0, 0, 0, 1; -t5 * t24, t2, t6, t29; -t6 * t24, t1, -t5, t26; -t22, -t24, 0, t11; 0, 0, 0, 1; t6 * t21 - t39 * t5, t6 * t23 + t40 * t5, -t2, t33 * t5 + t29; -t5 * t21 - t39 * t6, -t5 * t23 + t40 * t6, -t1, t33 * t6 + t26; -t22 * t23, t22 * t21, t24, -t22 * pkin(4) + t24 * pkin(8) + t11; 0, 0, 0, 1; t6 * t12 - t41 * t5, t6 * t13 + t42 * t5, -t2, t32 * t5 + t36 * t6 + t31; -t5 * t12 - t41 * t6, -t5 * t13 + t42 * t6, -t1, t32 * t6 - t36 * t5 + t27; -t22 * t13, t22 * t12, t24, -t22 * t10 - t24 * t25 + t11; 0, 0, 0, 1;];
T_ges = t3;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
