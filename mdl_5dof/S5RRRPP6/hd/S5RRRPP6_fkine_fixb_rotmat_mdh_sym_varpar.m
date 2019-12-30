% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-29 19:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRPP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 19:47:02
% EndTime: 2019-12-29 19:47:03
% DurationCPUTime: 0.30s
% Computational Cost: add. (116->47), mult. (132->51), div. (0->0), fcn. (187->8), ass. (0->31)
t20 = qJ(3) + pkin(8);
t15 = sin(t20);
t24 = sin(qJ(2));
t39 = t24 * t15;
t23 = sin(qJ(3));
t25 = sin(qJ(1));
t38 = t25 * t23;
t12 = t25 * t24;
t27 = cos(qJ(2));
t37 = t25 * t27;
t28 = cos(qJ(1));
t13 = t28 * t24;
t36 = t28 * t27;
t21 = pkin(5) + 0;
t35 = t25 * pkin(1) + 0;
t34 = t28 * pkin(1) + t25 * pkin(6) + 0;
t26 = cos(qJ(3));
t14 = t26 * pkin(3) + pkin(2);
t22 = -qJ(4) - pkin(7);
t33 = t24 * t14 + t27 * t22 + t21;
t32 = pkin(2) * t27 + pkin(7) * t24;
t31 = -t28 * pkin(6) + t35;
t30 = pkin(3) * t38 - t22 * t13 + t14 * t36 + t34;
t29 = -t22 * t12 + t14 * t37 + (-pkin(3) * t23 - pkin(6)) * t28 + t35;
t16 = cos(t20);
t8 = t24 * t16;
t4 = t25 * t15 + t16 * t36;
t3 = t15 * t36 - t25 * t16;
t2 = -t28 * t15 + t16 * t37;
t1 = t15 * t37 + t28 * t16;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; t36, -t13, t25, t34; t37, -t12, -t28, t31; t24, t27, 0, t21; 0, 0, 0, 1; t26 * t36 + t38, -t23 * t36 + t25 * t26, t13, t32 * t28 + t34; -t28 * t23 + t26 * t37, -t23 * t37 - t28 * t26, t12, t32 * t25 + t31; t24 * t26, -t24 * t23, -t27, t24 * pkin(2) - t27 * pkin(7) + t21; 0, 0, 0, 1; t4, -t3, t13, t30; t2, -t1, t12, t29; t8, -t39, -t27, t33; 0, 0, 0, 1; t4, t13, t3, t4 * pkin(4) + t3 * qJ(5) + t30; t2, t12, t1, t2 * pkin(4) + t1 * qJ(5) + t29; t8, -t27, t39, (pkin(4) * t16 + qJ(5) * t15) * t24 + t33; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
