% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:16
% EndTime: 2019-12-31 17:29:16
% DurationCPUTime: 0.10s
% Computational Cost: add. (99->42), mult. (221->53), div. (0->0), fcn. (311->10), ass. (0->36)
t21 = sin(pkin(4));
t25 = sin(qJ(2));
t45 = t21 * t25;
t29 = cos(qJ(2));
t44 = t21 * t29;
t26 = sin(qJ(1));
t43 = t26 * t21;
t42 = t26 * t25;
t41 = t26 * t29;
t30 = cos(qJ(1));
t40 = t30 * t21;
t39 = t30 * t25;
t38 = t30 * t29;
t37 = pkin(5) + 0;
t22 = cos(pkin(4));
t36 = t22 * pkin(6) + t37;
t35 = t30 * pkin(1) + pkin(6) * t43 + 0;
t34 = t26 * pkin(1) - pkin(6) * t40 + 0;
t11 = t22 * t41 + t39;
t12 = -t22 * t42 + t38;
t33 = t12 * pkin(2) + t11 * pkin(7) + t35;
t32 = pkin(2) * t45 - pkin(7) * t44 + t36;
t10 = t22 * t39 + t41;
t9 = -t22 * t38 + t42;
t31 = t10 * pkin(2) + t9 * pkin(7) + t34;
t28 = cos(qJ(3));
t27 = cos(qJ(4));
t24 = sin(qJ(3));
t23 = sin(qJ(4));
t8 = t22 * t24 + t28 * t45;
t7 = -t22 * t28 + t24 * t45;
t4 = t12 * t28 + t24 * t43;
t3 = t12 * t24 - t28 * t43;
t2 = t10 * t28 - t24 * t40;
t1 = t10 * t24 + t28 * t40;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t30, -t26, 0, 0; t26, t30, 0, 0; 0, 0, 1, t37; 0, 0, 0, 1; t12, -t11, t43, t35; t10, -t9, -t40, t34; t45, t44, t22, t36; 0, 0, 0, 1; t4, -t3, t11, t33; t2, -t1, t9, t31; t8, -t7, -t44, t32; 0, 0, 0, 1; t11 * t23 + t4 * t27, t11 * t27 - t4 * t23, t3, t4 * pkin(3) + t3 * pkin(8) + t33; t2 * t27 + t9 * t23, -t2 * t23 + t9 * t27, t1, t2 * pkin(3) + t1 * pkin(8) + t31; -t23 * t44 + t8 * t27, -t8 * t23 - t27 * t44, t7, t8 * pkin(3) + t7 * pkin(8) + t32; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
