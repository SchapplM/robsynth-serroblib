% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-29 18:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPPP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 18:07:25
% EndTime: 2019-12-29 18:07:26
% DurationCPUTime: 0.28s
% Computational Cost: add. (133->52), mult. (304->55), div. (0->0), fcn. (406->8), ass. (0->39)
t27 = sin(pkin(8));
t29 = cos(pkin(8));
t31 = sin(qJ(2));
t30 = cos(pkin(5));
t33 = cos(qJ(2));
t54 = t30 * t33;
t8 = t31 * t27 - t29 * t54;
t53 = t31 * t29;
t9 = t27 * t54 + t53;
t58 = t9 * pkin(3) + t8 * qJ(4);
t28 = sin(pkin(5));
t34 = cos(qJ(1));
t32 = sin(qJ(1));
t52 = t32 * t31;
t57 = t28 * t34 + t30 * t52;
t56 = t28 * t32;
t51 = t32 * t33;
t50 = t33 * t28;
t49 = t34 * t31;
t48 = t34 * t33;
t46 = qJ(3) * t28;
t45 = qJ(3) * t30;
t26 = pkin(6) + 0;
t44 = t32 * pkin(1) + 0;
t42 = t31 * t46;
t41 = t31 * pkin(2) + t26;
t40 = t34 * pkin(1) + t32 * pkin(7) + 0;
t39 = pkin(2) * t48 + t32 * t45 + t34 * t42 + t40;
t38 = -t33 * t46 + t41;
t5 = -t29 * t56 + (t27 * t33 + t30 * t53) * t34;
t6 = t29 * t48 + (-t30 * t49 + t56) * t27;
t37 = t6 * pkin(3) + t5 * qJ(4) + t39;
t36 = t32 * t42 + pkin(2) * t51 + (-pkin(7) - t45) * t34 + t44;
t3 = t27 * t51 + t57 * t29;
t4 = -t57 * t27 + t29 * t51;
t35 = t4 * pkin(3) + t3 * qJ(4) + t36;
t11 = t28 * t49 + t32 * t30;
t10 = t28 * t52 - t34 * t30;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t34, -t32, 0, 0; t32, t34, 0, 0; 0, 0, 1, t26; 0, 0, 0, 1; t48, -t49, t32, t40; t51, -t52, -t34, -t34 * pkin(7) + t44; t31, t33, 0, t26; 0, 0, 0, 1; t6, -t5, t11, t39; t4, -t3, t10, t36; t9, -t8, -t50, t38; 0, 0, 0, 1; t11, -t6, t5, t37; t10, -t4, t3, t35; -t50, -t9, t8, t38 + t58; 0, 0, 0, 1; t11, t5, t6, t11 * pkin(4) + t6 * qJ(5) + t37; t10, t3, t4, t10 * pkin(4) + t4 * qJ(5) + t35; -t50, t8, t9, t9 * qJ(5) + (-pkin(4) - qJ(3)) * t50 + t41 + t58; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
