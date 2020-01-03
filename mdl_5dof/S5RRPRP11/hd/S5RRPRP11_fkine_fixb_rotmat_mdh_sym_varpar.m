% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRP11_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:14
% EndTime: 2019-12-31 20:12:15
% DurationCPUTime: 0.10s
% Computational Cost: add. (79->44), mult. (120->37), div. (0->0), fcn. (170->6), ass. (0->31)
t23 = sin(qJ(1));
t25 = cos(qJ(2));
t12 = t23 * t25;
t22 = sin(qJ(2));
t35 = qJ(3) * t22;
t42 = pkin(2) * t12 + t23 * t35;
t41 = t23 * t22;
t24 = cos(qJ(4));
t40 = t23 * t24;
t21 = sin(qJ(4));
t39 = t25 * t21;
t38 = t25 * t24;
t26 = cos(qJ(1));
t37 = t26 * t22;
t36 = t26 * t24;
t13 = t26 * t25;
t20 = pkin(5) + 0;
t34 = t23 * pkin(1) + 0;
t33 = t22 * pkin(2) + t20;
t32 = t26 * pkin(1) + t23 * pkin(6) + 0;
t31 = -t26 * pkin(6) + t34;
t30 = pkin(2) * t13 + t26 * t35 + t32;
t29 = -t25 * qJ(3) + t33;
t28 = t23 * pkin(3) + pkin(7) * t13 + t30;
t27 = pkin(7) * t12 + (-pkin(3) - pkin(6)) * t26 + t34 + t42;
t14 = t22 * pkin(7);
t4 = t21 * t41 - t36;
t3 = t26 * t21 + t22 * t40;
t2 = t21 * t37 + t40;
t1 = t23 * t21 - t22 * t36;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t26, -t23, 0, 0; t23, t26, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t13, -t37, t23, t32; t12, -t41, -t26, t31; t22, t25, 0, t20; 0, 0, 0, 1; t23, -t13, t37, t30; -t26, -t12, t41, t31 + t42; 0, -t22, -t25, t29; 0, 0, 0, 1; t2, -t1, t13, t28; t4, t3, t12, t27; -t39, -t38, t22, t14 + t29; 0, 0, 0, 1; t2, t13, t1, t2 * pkin(4) + t1 * qJ(5) + t28; t4, t12, -t3, t4 * pkin(4) - t3 * qJ(5) + t27; -t39, t22, t38, t14 + (-pkin(4) * t21 + qJ(5) * t24 - qJ(3)) * t25 + t33; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
