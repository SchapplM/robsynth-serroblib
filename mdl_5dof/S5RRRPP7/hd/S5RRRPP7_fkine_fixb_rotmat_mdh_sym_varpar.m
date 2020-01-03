% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5RRRPP7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:50
% EndTime: 2019-12-31 21:03:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (88->40), mult. (153->37), div. (0->0), fcn. (214->6), ass. (0->28)
t23 = sin(qJ(3));
t24 = sin(qJ(2));
t13 = t24 * t23;
t26 = cos(qJ(3));
t14 = t24 * t26;
t40 = pkin(3) * t14 + qJ(4) * t13;
t25 = sin(qJ(1));
t15 = t25 * t24;
t27 = cos(qJ(2));
t39 = t25 * t27;
t28 = cos(qJ(1));
t17 = t28 * t24;
t38 = t28 * t27;
t37 = qJ(5) * t24;
t22 = pkin(5) + 0;
t36 = t24 * pkin(2) + t22;
t35 = t28 * pkin(1) + t25 * pkin(6) + 0;
t34 = t25 * pkin(1) - t28 * pkin(6) + 0;
t33 = pkin(2) * t38 + pkin(7) * t17 + t35;
t32 = -t27 * pkin(7) + t36;
t31 = pkin(2) * t39 + pkin(7) * t15 + t34;
t5 = t23 * t38 - t25 * t26;
t6 = t25 * t23 + t26 * t38;
t30 = t6 * pkin(3) + t5 * qJ(4) + t33;
t3 = t23 * t39 + t28 * t26;
t4 = -t28 * t23 + t26 * t39;
t29 = t4 * pkin(3) + t3 * qJ(4) + t31;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t28, -t25, 0, 0; t25, t28, 0, 0; 0, 0, 1, t22; 0, 0, 0, 1; t38, -t17, t25, t35; t39, -t15, -t28, t34; t24, t27, 0, t22; 0, 0, 0, 1; t6, -t5, t17, t33; t4, -t3, t15, t31; t14, -t13, -t27, t32; 0, 0, 0, 1; t6, t17, t5, t30; t4, t15, t3, t29; t14, -t27, t13, t32 + t40; 0, 0, 0, 1; t6, t5, -t17, t6 * pkin(4) - t28 * t37 + t30; t4, t3, -t15, t4 * pkin(4) - t25 * t37 + t29; t14, t13, t27, pkin(4) * t14 + (-pkin(7) + qJ(5)) * t27 + t36 + t40; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
