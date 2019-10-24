% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-10-24 10:25
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PRPRP6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-10-24 10:25:08
% EndTime: 2019-10-24 10:25:09
% DurationCPUTime: 0.11s
% Computational Cost: add. (79->44), mult. (120->39), div. (0->0), fcn. (170->6), ass. (0->31)
t21 = sin(pkin(7));
t26 = cos(qJ(2));
t12 = t21 * t26;
t24 = sin(qJ(2));
t35 = qJ(3) * t24;
t42 = pkin(2) * t12 + t21 * t35;
t41 = t21 * t24;
t22 = cos(pkin(7));
t40 = t22 * t24;
t13 = t22 * t26;
t23 = sin(qJ(4));
t39 = t23 * t24;
t25 = cos(qJ(4));
t38 = t24 * t25;
t37 = t26 * t23;
t36 = t26 * t25;
t34 = t21 * pkin(1) + 0;
t20 = qJ(1) + 0;
t33 = t22 * pkin(1) + t21 * pkin(5) + 0;
t32 = t24 * pkin(2) + t20;
t31 = -t22 * pkin(5) + t34;
t30 = pkin(2) * t13 + t22 * t35 + t33;
t29 = -t26 * qJ(3) + t32;
t28 = t21 * pkin(3) + pkin(6) * t13 + t30;
t27 = pkin(6) * t12 + (-pkin(3) - pkin(5)) * t22 + t34 + t42;
t18 = t24 * pkin(6);
t4 = t21 * t39 - t22 * t25;
t3 = t21 * t38 + t22 * t23;
t2 = t21 * t25 + t22 * t39;
t1 = t21 * t23 - t22 * t38;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t22, -t21, 0, 0; t21, t22, 0, 0; 0, 0, 1, t20; 0, 0, 0, 1; t13, -t40, t21, t33; t12, -t41, -t22, t31; t24, t26, 0, t20; 0, 0, 0, 1; t21, -t13, t40, t30; -t22, -t12, t41, t31 + t42; 0, -t24, -t26, t29; 0, 0, 0, 1; t2, -t1, t13, t28; t4, t3, t12, t27; -t37, -t36, t24, t18 + t29; 0, 0, 0, 1; t2, t13, t1, t2 * pkin(4) + t1 * qJ(5) + t28; t4, t12, -t3, t4 * pkin(4) - t3 * qJ(5) + t27; -t37, t24, t36, t18 + (-pkin(4) * t23 + qJ(5) * t25 - qJ(3)) * t26 + t32; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
