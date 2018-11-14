% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RPPP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:21
% EndTime: 2018-11-14 13:45:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (162->41), mult. (157->34), div. (0->0), fcn. (177->10), ass. (0->32)
t21 = sin(pkin(6));
t25 = sin(qJ(1));
t26 = cos(qJ(1));
t39 = pkin(4) - pkin(6);
t33 = cos(t39) / 0.2e1;
t38 = pkin(4) + pkin(6);
t35 = cos(t38);
t31 = t33 + t35 / 0.2e1;
t3 = t25 * t21 - t26 * t31;
t23 = cos(pkin(6));
t32 = sin(t38) / 0.2e1;
t34 = sin(t39);
t30 = t32 - t34 / 0.2e1;
t4 = t25 * t23 + t26 * t30;
t45 = t4 * pkin(2) + t3 * qJ(3);
t22 = sin(pkin(4));
t17 = t25 * t22;
t44 = t26 * t22;
t42 = qJ(2) * t22;
t41 = pkin(5) + 0;
t40 = t25 * pkin(1) + 0;
t24 = cos(pkin(4));
t37 = t24 * qJ(2) + t41;
t36 = t26 * pkin(1) + t25 * t42 + 0;
t29 = -t26 * t42 + t40;
t10 = t32 + t34 / 0.2e1;
t11 = t33 - t35 / 0.2e1;
t28 = t11 * pkin(2) - t10 * qJ(3) + t37;
t5 = t26 * t21 + t25 * t31;
t6 = t26 * t23 - t25 * t30;
t27 = t6 * pkin(2) + t5 * qJ(3) + t36;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t25, 0, 0; t25, t26, 0, 0; 0, 0, 1, t41; 0, 0, 0, 1; t6, -t5, t17, t36; t4, -t3, -t44, t29; t11, t10, t24, t37; 0, 0, 0, 1; t17, -t6, t5, t27; -t44, -t4, t3, t29 + t45; t24, -t11, -t10, t28; 0, 0, 0, 1; t17, t5, t6, pkin(3) * t17 + t6 * qJ(4) + t27; -t44, t3, t4, t4 * qJ(4) + (-pkin(3) - qJ(2)) * t44 + t40 + t45; t24, -t10, t11, t24 * pkin(3) + t11 * qJ(4) + t28; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
