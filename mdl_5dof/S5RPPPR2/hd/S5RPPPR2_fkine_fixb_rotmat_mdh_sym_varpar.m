% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:50
% EndTime: 2019-12-05 17:30:50
% DurationCPUTime: 0.15s
% Computational Cost: add. (124->53), mult. (250->61), div. (0->0), fcn. (348->10), ass. (0->37)
t25 = sin(pkin(8));
t26 = sin(pkin(7));
t48 = t26 * t25;
t28 = cos(pkin(8));
t47 = t26 * t28;
t31 = sin(qJ(1));
t46 = t31 * t26;
t29 = cos(pkin(7));
t45 = t31 * t29;
t33 = cos(qJ(1));
t44 = t33 * t26;
t43 = t33 * t29;
t42 = qJ(3) * t26;
t23 = pkin(5) + 0;
t41 = t33 * qJ(2) + 0;
t40 = t33 * pkin(1) + t31 * qJ(2) + 0;
t39 = pkin(2) * t43 + t33 * t42 + t40;
t38 = t26 * pkin(2) - t29 * qJ(3) + t23;
t37 = pkin(3) * t47 + qJ(4) * t48 + t38;
t11 = t25 * t43 - t31 * t28;
t12 = t31 * t25 + t28 * t43;
t36 = t12 * pkin(3) + t11 * qJ(4) + t39;
t35 = (-pkin(2) * t29 - pkin(1) - t42) * t31 + t41;
t10 = t33 * t25 - t28 * t45;
t9 = t25 * t45 + t33 * t28;
t34 = t10 * pkin(3) - t9 * qJ(4) + t35;
t32 = cos(qJ(5));
t30 = sin(qJ(5));
t27 = cos(pkin(9));
t24 = sin(pkin(9));
t8 = -t29 * t24 + t27 * t47;
t7 = t24 * t47 + t29 * t27;
t4 = t12 * t27 + t24 * t44;
t3 = t12 * t24 - t27 * t44;
t2 = t10 * t27 - t24 * t46;
t1 = t10 * t24 + t27 * t46;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t23; -t31, -t33, 0, 0; t33, -t31, 0, 0; 0, 0, 0, 1; t26, t29, 0, t23; -t45, t46, t33, -t31 * pkin(1) + t41; t43, -t44, t31, t40; 0, 0, 0, 1; t47, -t48, -t29, t38; t10, t9, -t46, t35; t12, -t11, t44, t39; 0, 0, 0, 1; t8, -t7, t48, t37; t2, -t1, -t9, t34; t4, -t3, t11, t36; 0, 0, 0, 1; t30 * t48 + t8 * t32, -t8 * t30 + t32 * t48, t7, t8 * pkin(4) + t7 * pkin(6) + t37; t2 * t32 - t9 * t30, -t2 * t30 - t9 * t32, t1, t2 * pkin(4) + t1 * pkin(6) + t34; t11 * t30 + t4 * t32, t11 * t32 - t4 * t30, t3, t4 * pkin(4) + t3 * pkin(6) + t36; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
