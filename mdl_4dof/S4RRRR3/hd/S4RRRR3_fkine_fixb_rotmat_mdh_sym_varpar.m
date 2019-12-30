% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% T_c_mdh [4x4x(4+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   5:  mdh base (link 0) -> mdh frame (5-1), link (5-1)
%   ...
%   4+1:  mdh base (link 0) -> mdh frame (4)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_c_mdh = S4RRRR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 15:05:52
% EndTime: 2019-12-29 15:05:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (62->31), mult. (33->28), div. (0->0), fcn. (61->8), ass. (0->17)
t17 = -pkin(6) - pkin(5);
t15 = cos(qJ(2));
t4 = t15 * pkin(2) + pkin(1);
t12 = qJ(2) + qJ(3);
t10 = pkin(4) + 0;
t13 = sin(qJ(2));
t18 = t13 * pkin(2) + t10;
t16 = cos(qJ(1));
t14 = sin(qJ(1));
t11 = -pkin(7) + t17;
t7 = qJ(4) + t12;
t6 = cos(t12);
t5 = sin(t12);
t3 = cos(t7);
t2 = sin(t7);
t1 = pkin(3) * t6 + t4;
t8 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t16, -t14, 0, 0; t14, t16, 0, 0; 0, 0, 1, t10; 0, 0, 0, 1; t16 * t15, -t16 * t13, t14, t16 * pkin(1) + t14 * pkin(5) + 0; t14 * t15, -t14 * t13, -t16, t14 * pkin(1) - t16 * pkin(5) + 0; t13, t15, 0, t10; 0, 0, 0, 1; t16 * t6, -t16 * t5, t14, -t14 * t17 + t16 * t4 + 0; t14 * t6, -t14 * t5, -t16, t14 * t4 + t16 * t17 + 0; t5, t6, 0, t18; 0, 0, 0, 1; t16 * t3, -t16 * t2, t14, t16 * t1 - t14 * t11 + 0; t14 * t3, -t14 * t2, -t16, t14 * t1 + t16 * t11 + 0; t2, t3, 0, pkin(3) * t5 + t18; 0, 0, 0, 1;];
T_ges = t8;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,4+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,4+1]); end % symbolisch
for i = 1:4+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
