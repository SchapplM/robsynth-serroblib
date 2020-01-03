% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRPRR12_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:28:54
% EndTime: 2019-12-31 20:28:54
% DurationCPUTime: 0.12s
% Computational Cost: add. (88->41), mult. (154->46), div. (0->0), fcn. (216->8), ass. (0->30)
t26 = sin(qJ(2));
t27 = sin(qJ(1));
t43 = t27 * t26;
t30 = cos(qJ(2));
t15 = t27 * t30;
t25 = sin(qJ(4));
t42 = t30 * t25;
t31 = cos(qJ(1));
t41 = t31 * t26;
t16 = t31 * t30;
t40 = qJ(3) * t26;
t23 = pkin(5) + 0;
t39 = t31 * pkin(1) + t27 * pkin(6) + 0;
t29 = cos(qJ(4));
t5 = t26 * t25 + t30 * t29;
t38 = t27 * pkin(1) - t31 * pkin(6) + 0;
t37 = pkin(2) * t16 + t31 * t40 + t39;
t36 = t26 * pkin(2) - t30 * qJ(3) + t23;
t35 = pkin(2) * t15 + t27 * t40 + t38;
t34 = t26 * pkin(3) + t36;
t33 = pkin(3) * t15 + t31 * pkin(7) + t35;
t32 = pkin(3) * t16 - t27 * pkin(7) + t37;
t28 = cos(qJ(5));
t24 = sin(qJ(5));
t6 = t26 * t29 - t42;
t4 = t5 * t31;
t3 = t25 * t16 - t29 * t41;
t2 = t5 * t27;
t1 = t27 * t42 - t29 * t43;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t31, -t27, 0, 0; t27, t31, 0, 0; 0, 0, 1, t23; 0, 0, 0, 1; t16, -t41, t27, t39; t15, -t43, -t31, t38; t26, t30, 0, t23; 0, 0, 0, 1; t16, t27, t41, t37; t15, -t31, t43, t35; t26, 0, -t30, t36; 0, 0, 0, 1; t4, -t3, -t27, t32; t2, -t1, t31, t33; t6, -t5, 0, t34; 0, 0, 0, 1; -t27 * t24 + t4 * t28, -t4 * t24 - t27 * t28, t3, t4 * pkin(4) + t3 * pkin(8) + t32; t2 * t28 + t31 * t24, -t2 * t24 + t31 * t28, t1, t2 * pkin(4) + t1 * pkin(8) + t33; t6 * t28, -t6 * t24, t5, t6 * pkin(4) + t5 * pkin(8) + t34; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
