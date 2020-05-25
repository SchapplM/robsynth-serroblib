% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRPPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:37:55
% EndTime: 2019-12-31 17:37:55
% DurationCPUTime: 0.14s
% Computational Cost: add. (88->40), mult. (154->44), div. (0->0), fcn. (216->8), ass. (0->29)
t24 = sin(pkin(8));
t26 = cos(pkin(8));
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t45 = t31 * t24 - t29 * t26;
t25 = sin(pkin(7));
t44 = t25 * t29;
t15 = t25 * t31;
t27 = cos(pkin(7));
t43 = t27 * t29;
t16 = t27 * t31;
t40 = qJ(3) * t29;
t23 = qJ(1) + 0;
t39 = t27 * pkin(1) + t25 * pkin(5) + 0;
t5 = t29 * t24 + t31 * t26;
t38 = t25 * pkin(1) - t27 * pkin(5) + 0;
t37 = pkin(2) * t16 + t27 * t40 + t39;
t36 = t29 * pkin(2) - t31 * qJ(3) + t23;
t35 = pkin(2) * t15 + t25 * t40 + t38;
t34 = t29 * pkin(3) + t36;
t33 = pkin(3) * t15 + t27 * qJ(4) + t35;
t32 = pkin(3) * t16 - t25 * qJ(4) + t37;
t30 = cos(qJ(5));
t28 = sin(qJ(5));
t4 = t5 * t27;
t3 = t45 * t27;
t2 = t5 * t25;
t1 = t45 * t25;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t27, -t25, 0, 0; t25, t27, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t16, -t43, t25, t39; t15, -t44, -t27, t38; t29, t31, 0, t23; 0, 0, 0, 1; t16, t25, t43, t37; t15, -t27, t44, t35; t29, 0, -t31, t36; 0, 0, 0, 1; t4, -t3, -t25, t32; t2, -t1, t27, t33; -t45, -t5, 0, t34; 0, 0, 0, 1; -t25 * t28 + t4 * t30, -t25 * t30 - t4 * t28, t3, t4 * pkin(4) + t3 * pkin(6) + t32; t2 * t30 + t27 * t28, -t2 * t28 + t27 * t30, t1, t2 * pkin(4) + t1 * pkin(6) + t33; -t45 * t30, t45 * t28, t5, -pkin(4) * t45 + t5 * pkin(6) + t34; 0, 0, 0, 1;];
T_ges = t6;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
