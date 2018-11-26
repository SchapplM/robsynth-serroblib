% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2018-11-23 17:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRRPRP9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:47:05
% EndTime: 2018-11-23 17:47:05
% DurationCPUTime: 0.14s
% Computational Cost: add. (163->55), mult. (316->55), div. (0->0), fcn. (430->8), ass. (0->37)
t36 = sin(qJ(3));
t37 = sin(qJ(2));
t56 = t37 * t36;
t40 = cos(qJ(3));
t25 = t37 * t40;
t38 = sin(qJ(1));
t26 = t38 * t37;
t41 = cos(qJ(2));
t55 = t38 * t41;
t42 = cos(qJ(1));
t28 = t42 * t37;
t54 = t42 * t41;
t34 = pkin(6) + 0;
t53 = t42 * pkin(1) + t38 * pkin(7) + 0;
t52 = t38 * pkin(1) - t42 * pkin(7) + 0;
t51 = pkin(2) * t54 + pkin(8) * t28 + t53;
t50 = t37 * pkin(2) - t41 * pkin(8) + t34;
t49 = pkin(2) * t55 + pkin(8) * t26 + t52;
t48 = pkin(3) * t25 + qJ(4) * t56 + t50;
t15 = t36 * t54 - t38 * t40;
t16 = t38 * t36 + t40 * t54;
t47 = t16 * pkin(3) + t15 * qJ(4) + t51;
t46 = pkin(4) * t25 + t41 * pkin(9) + t48;
t13 = t36 * t55 + t42 * t40;
t14 = -t42 * t36 + t40 * t55;
t45 = t14 * pkin(3) + t13 * qJ(4) + t49;
t44 = t16 * pkin(4) - pkin(9) * t28 + t47;
t43 = t14 * pkin(4) - pkin(9) * t26 + t45;
t39 = cos(qJ(5));
t35 = sin(qJ(5));
t8 = (t35 * t36 + t39 * t40) * t37;
t7 = t35 * t25 - t39 * t56;
t4 = t15 * t35 + t16 * t39;
t3 = -t15 * t39 + t16 * t35;
t2 = t13 * t35 + t14 * t39;
t1 = -t13 * t39 + t14 * t35;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t42, -t38, 0, 0; t38, t42, 0, 0; 0, 0, 1, t34; 0, 0, 0, 1; t54, -t28, t38, t53; t55, -t26, -t42, t52; t37, t41, 0, t34; 0, 0, 0, 1; t16, -t15, t28, t51; t14, -t13, t26, t49; t25, -t56, -t41, t50; 0, 0, 0, 1; t16, t28, t15, t47; t14, t26, t13, t45; t25, -t41, t56, t48; 0, 0, 0, 1; t4, -t3, -t28, t44; t2, -t1, -t26, t43; t8, -t7, t41, t46; 0, 0, 0, 1; t4, -t28, t3, t4 * pkin(5) + t3 * qJ(6) + t44; t2, -t26, t1, t2 * pkin(5) + t1 * qJ(6) + t43; t8, t41, t7, t8 * pkin(5) + t7 * qJ(6) + t46; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
