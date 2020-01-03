% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S5PPRRP4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:21
% EndTime: 2019-12-31 17:34:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (77->36), mult. (89->24), div. (0->0), fcn. (141->6), ass. (0->21)
t18 = sin(qJ(4));
t16 = cos(pkin(7));
t24 = sin(pkin(7));
t25 = sin(qJ(3));
t26 = cos(qJ(3));
t3 = -t16 * t26 - t24 * t25;
t28 = t3 * t18;
t4 = t16 * t25 - t24 * t26;
t27 = t4 * t18;
t15 = qJ(1) + 0;
t23 = t16 * pkin(1) + t24 * qJ(2) + 0;
t10 = -pkin(5) + t15;
t22 = t16 * pkin(2) + t23;
t21 = t24 * pkin(1) - t16 * qJ(2) + 0;
t20 = t24 * pkin(2) + t21;
t19 = cos(qJ(4));
t17 = -qJ(5) - pkin(6);
t8 = t19 * pkin(4) + pkin(3);
t2 = t3 * t19;
t1 = t4 * t19;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t24, 0, 0; t24, t16, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t16, 0, t24, t23; t24, 0, -t16, t21; 0, 1, 0, t15; 0, 0, 0, 1; -t3, -t4, 0, t22; -t4, t3, 0, t20; 0, 0, -1, t10; 0, 0, 0, 1; -t2, t28, t4, -t3 * pkin(3) + t4 * pkin(6) + t22; -t1, t27, -t3, -t4 * pkin(3) - t3 * pkin(6) + t20; -t18, -t19, 0, t10; 0, 0, 0, 1; -t2, t28, t4, -t4 * t17 - t3 * t8 + t22; -t1, t27, -t3, t3 * t17 - t4 * t8 + t20; -t18, -t19, 0, -t18 * pkin(4) + t10; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
