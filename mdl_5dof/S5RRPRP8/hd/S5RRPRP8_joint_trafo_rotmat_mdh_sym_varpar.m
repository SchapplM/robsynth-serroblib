% Calculate homogenous joint transformation matrices for
% S5RRPRP8
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
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5RRPRP8_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 18:51:57
% EndTime: 2019-12-29 18:51:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (6->6), mult. (0->0), div. (0->0), fcn. (12->6), ass. (0->7)
t45 = cos(qJ(1));
t44 = cos(qJ(2));
t43 = cos(qJ(4));
t42 = sin(qJ(1));
t41 = sin(qJ(2));
t40 = sin(qJ(4));
t1 = [t45, -t42, 0, 0; t42, t45, 0, 0; 0, 0, 1, pkin(5); 0, 0, 0, 1; t44, -t41, 0, pkin(1); 0, 0, -1, -pkin(6); t41, t44, 0, 0; 0, 0, 0, 1; 1, 0, 0, pkin(2); 0, 0, -1, -qJ(3); 0, 1, 0, 0; 0, 0, 0, 1; t43, -t40, 0, pkin(3); 0, 0, -1, -pkin(7); t40, t43, 0, 0; 0, 0, 0, 1; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, qJ(5); 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
