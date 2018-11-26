% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
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
% Datum: 2018-11-23 16:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RPRPRP10_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:02:26
% EndTime: 2018-11-23 16:02:27
% DurationCPUTime: 0.12s
% Computational Cost: add. (97->48), mult. (132->41), div. (0->0), fcn. (186->6), ass. (0->31)
t23 = sin(qJ(3));
t25 = cos(qJ(5));
t43 = t23 * t25;
t24 = sin(qJ(1));
t10 = t24 * t23;
t26 = cos(qJ(3));
t42 = t24 * t26;
t27 = cos(qJ(1));
t41 = t27 * t23;
t40 = t27 * t26;
t39 = qJ(4) * t26;
t21 = pkin(6) + 0;
t38 = t24 * pkin(1) + 0;
t37 = pkin(2) + t21;
t36 = t27 * pkin(1) + t24 * qJ(2) + 0;
t13 = t24 * pkin(7);
t35 = t27 * t39 + t13 + t38;
t34 = t27 * pkin(7) + t36;
t33 = t26 * pkin(3) + t23 * qJ(4) + t37;
t32 = -t27 * qJ(2) + t38;
t31 = t26 * pkin(8) + t33;
t30 = pkin(3) * t10 - t24 * t39 + t34;
t29 = t27 * pkin(4) + pkin(8) * t10 + t30;
t28 = t24 * pkin(4) + (-qJ(2) + (-pkin(3) - pkin(8)) * t23) * t27 + t35;
t22 = sin(qJ(5));
t9 = t23 * t22;
t4 = -t22 * t42 + t27 * t25;
t3 = t27 * t22 + t25 * t42;
t2 = t22 * t40 + t24 * t25;
t1 = t24 * t22 - t25 * t40;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t24, 0, 0; t24, t27, 0, 0; 0, 0, 1, t21; 0, 0, 0, 1; 0, -t27, t24, t36; 0, -t24, -t27, t32; 1, 0, 0, t21; 0, 0, 0, 1; t10, t42, t27, t34; -t41, -t40, t24, t13 + t32; t26, -t23, 0, t37; 0, 0, 0, 1; t27, -t10, -t42, t30; t24, t41, t40 (-pkin(3) * t23 - qJ(2)) * t27 + t35; 0, -t26, t23, t33; 0, 0, 0, 1; t4, -t3, t10, t29; t2, -t1, -t41, t28; t9, t43, t26, t31; 0, 0, 0, 1; t4, t10, t3, t4 * pkin(5) + t3 * qJ(6) + t29; t2, -t41, t1, t2 * pkin(5) + t1 * qJ(6) + t28; t9, t26, -t43 (pkin(5) * t22 - qJ(6) * t25) * t23 + t31; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
