% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2018-11-23 16:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRRPP2_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:11:57
% EndTime: 2018-11-23 16:11:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (195->45), mult. (163->41), div. (0->0), fcn. (228->8), ass. (0->36)
t28 = sin(qJ(4));
t29 = sin(qJ(3));
t19 = t29 * t28;
t31 = cos(qJ(4));
t20 = t29 * t31;
t50 = pkin(4) * t20 + qJ(5) * t19;
t27 = qJ(1) + pkin(9);
t21 = sin(t27);
t11 = t21 * t29;
t32 = cos(qJ(3));
t49 = t21 * t32;
t22 = cos(t27);
t13 = t22 * t29;
t48 = t22 * t32;
t47 = t28 * t32;
t46 = t31 * t32;
t45 = qJ(6) * t29;
t44 = pkin(6) + 0;
t30 = sin(qJ(1));
t43 = t30 * pkin(1) + 0;
t33 = cos(qJ(1));
t42 = t33 * pkin(1) + 0;
t23 = qJ(2) + t44;
t41 = t29 * pkin(3) + t23;
t40 = t22 * pkin(2) + t21 * pkin(7) + t42;
t39 = t21 * pkin(2) - t22 * pkin(7) + t43;
t38 = pkin(3) * t48 + pkin(8) * t13 + t40;
t37 = -t32 * pkin(8) + t41;
t36 = pkin(3) * t49 + pkin(8) * t11 + t39;
t5 = -t21 * t31 + t22 * t47;
t6 = t21 * t28 + t22 * t46;
t35 = t6 * pkin(4) + t5 * qJ(5) + t38;
t3 = t21 * t47 + t22 * t31;
t4 = t21 * t46 - t22 * t28;
t34 = t4 * pkin(4) + t3 * qJ(5) + t36;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t33, -t30, 0, 0; t30, t33, 0, 0; 0, 0, 1, t44; 0, 0, 0, 1; t22, -t21, 0, t42; t21, t22, 0, t43; 0, 0, 1, t23; 0, 0, 0, 1; t48, -t13, t21, t40; t49, -t11, -t22, t39; t29, t32, 0, t23; 0, 0, 0, 1; t6, -t5, t13, t38; t4, -t3, t11, t36; t20, -t19, -t32, t37; 0, 0, 0, 1; t6, t13, t5, t35; t4, t11, t3, t34; t20, -t32, t19, t37 + t50; 0, 0, 0, 1; t6, t5, -t13, t6 * pkin(5) - t22 * t45 + t35; t4, t3, -t11, t4 * pkin(5) - t21 * t45 + t34; t20, t19, t32, pkin(5) * t20 + (-pkin(8) + qJ(6)) * t32 + t41 + t50; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
