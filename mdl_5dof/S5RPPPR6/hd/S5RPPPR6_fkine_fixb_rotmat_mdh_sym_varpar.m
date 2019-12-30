% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 15:52:53
% EndTime: 2019-12-29 15:52:53
% DurationCPUTime: 0.24s
% Computational Cost: add. (86->49), mult. (141->50), div. (0->0), fcn. (199->8), ass. (0->34)
t24 = cos(pkin(7));
t26 = sin(qJ(1));
t12 = t26 * t24;
t22 = sin(pkin(7));
t38 = qJ(3) * t22;
t45 = pkin(2) * t12 + t26 * t38;
t21 = sin(pkin(8));
t44 = t24 * t21;
t23 = cos(pkin(8));
t43 = t24 * t23;
t42 = t26 * t22;
t41 = t26 * t23;
t28 = cos(qJ(1));
t40 = t28 * t22;
t39 = t28 * t23;
t13 = t28 * t24;
t37 = qJ(4) * t24;
t20 = pkin(5) + 0;
t36 = t26 * pkin(1) + 0;
t35 = t22 * pkin(2) + t20;
t34 = t28 * pkin(1) + t26 * qJ(2) + 0;
t33 = pkin(2) * t13 + t28 * t38 + t34;
t32 = -t28 * qJ(2) + t36;
t31 = -t24 * qJ(3) + t35;
t30 = t26 * pkin(3) + t28 * t37 + t33;
t29 = t26 * t37 + (-pkin(3) - qJ(2)) * t28 + t36 + t45;
t27 = cos(qJ(5));
t25 = sin(qJ(5));
t14 = t22 * qJ(4);
t4 = t21 * t42 - t39;
t3 = t28 * t21 + t22 * t41;
t2 = t21 * t40 + t41;
t1 = t26 * t21 - t22 * t39;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t26, 0, 0; t26, t28, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t13, -t40, t26, t34; t12, -t42, -t28, t32; t22, t24, 0, t20; 0, 0, 0, 1; t26, -t13, t40, t33; -t28, -t12, t42, t32 + t45; 0, -t22, -t24, t31; 0, 0, 0, 1; t2, -t1, t13, t30; t4, t3, t12, t29; -t44, -t43, t22, t14 + t31; 0, 0, 0, 1; t25 * t13 + t2 * t27, t27 * t13 - t2 * t25, t1, t2 * pkin(4) + t1 * pkin(6) + t30; t25 * t12 + t4 * t27, t27 * t12 - t4 * t25, -t3, t4 * pkin(4) - t3 * pkin(6) + t29; t22 * t25 - t27 * t44, t22 * t27 + t25 * t44, t43, t14 + (-pkin(4) * t21 + pkin(6) * t23 - qJ(3)) * t24 + t35; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
