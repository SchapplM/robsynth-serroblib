% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PPRPR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:20
% EndTime: 2019-12-31 17:33:20
% DurationCPUTime: 0.08s
% Computational Cost: add. (68->28), mult. (84->20), div. (0->0), fcn. (132->6), ass. (0->17)
t26 = cos(qJ(3));
t25 = sin(qJ(3));
t24 = sin(pkin(7));
t14 = qJ(1) + 0;
t15 = cos(pkin(7));
t23 = t15 * pkin(1) + t24 * qJ(2) + 0;
t9 = -pkin(5) + t14;
t22 = t15 * pkin(2) + t23;
t21 = t24 * pkin(1) - t15 * qJ(2) + 0;
t20 = t24 * pkin(2) + t21;
t3 = -t15 * t26 - t24 * t25;
t4 = t15 * t25 - t24 * t26;
t19 = -t3 * pkin(3) + t4 * qJ(4) + t22;
t18 = -t4 * pkin(3) - t3 * qJ(4) + t20;
t17 = cos(qJ(5));
t16 = sin(qJ(5));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t15, -t24, 0, 0; t24, t15, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t15, 0, t24, t23; t24, 0, -t15, t21; 0, 1, 0, t14; 0, 0, 0, 1; -t3, -t4, 0, t22; -t4, t3, 0, t20; 0, 0, -1, t9; 0, 0, 0, 1; 0, t3, t4, t19; 0, t4, -t3, t18; -1, 0, 0, t9; 0, 0, 0, 1; t4 * t16, t4 * t17, -t3, -t3 * pkin(6) + t19; -t3 * t16, -t3 * t17, -t4, -t4 * pkin(6) + t18; -t17, t16, 0, -pkin(4) + t9; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
