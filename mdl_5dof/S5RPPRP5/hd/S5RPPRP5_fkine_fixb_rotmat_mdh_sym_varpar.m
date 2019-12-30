% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-29 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RPPRP5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 16:03:37
% EndTime: 2019-12-29 16:03:38
% DurationCPUTime: 0.24s
% Computational Cost: add. (81->37), mult. (134->36), div. (0->0), fcn. (188->6), ass. (0->27)
t24 = sin(pkin(7));
t27 = sin(qJ(1));
t40 = t27 * t24;
t25 = cos(pkin(7));
t15 = t27 * t25;
t29 = cos(qJ(1));
t39 = t29 * t24;
t16 = t29 * t25;
t38 = qJ(3) * t24;
t23 = pkin(5) + 0;
t37 = t29 * pkin(1) + t27 * qJ(2) + 0;
t26 = sin(qJ(4));
t28 = cos(qJ(4));
t5 = t24 * t26 + t25 * t28;
t36 = t27 * pkin(1) - t29 * qJ(2) + 0;
t35 = pkin(2) * t16 + t29 * t38 + t37;
t34 = t24 * pkin(2) - t25 * qJ(3) + t23;
t33 = t24 * pkin(3) + t34;
t32 = pkin(2) * t15 + t27 * t38 + t36;
t31 = pkin(3) * t16 - t27 * pkin(6) + t35;
t30 = pkin(3) * t15 + t29 * pkin(6) + t32;
t6 = t24 * t28 - t25 * t26;
t4 = t5 * t29;
t3 = t26 * t16 - t28 * t39;
t2 = t5 * t27;
t1 = t26 * t15 - t28 * t40;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t29, -t27, 0, 0; t27, t29, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t16, -t39, t27, t37; t15, -t40, -t29, t36; t24, t25, 0, t23; 0, 0, 0, 1; t16, t27, t39, t35; t15, -t29, t40, t32; t24, 0, -t25, t34; 0, 0, 0, 1; t4, -t3, -t27, t31; t2, -t1, t29, t30; t6, -t5, 0, t33; 0, 0, 0, 1; t4, -t27, t3, t4 * pkin(4) + t3 * qJ(5) + t31; t2, t29, t1, t2 * pkin(4) + t1 * qJ(5) + t30; t6, 0, t5, t6 * pkin(4) + t5 * qJ(5) + t33; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
