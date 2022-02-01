% Calculate homogenous joint transformation matrices for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5RPRPR3_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:30
% EndTime: 2022-01-23 09:20:30
% DurationCPUTime: 0.02s
% Computational Cost: add. (7->7), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t43 = cos(qJ(1));
t42 = cos(qJ(3));
t41 = cos(qJ(5));
t40 = sin(qJ(1));
t39 = sin(qJ(3));
t38 = sin(qJ(5));
t37 = cos(pkin(8));
t36 = cos(pkin(9));
t35 = sin(pkin(8));
t34 = sin(pkin(9));
t1 = [t43, -t40, 0, 0; t40, t43, 0, 0; 0, 0, 1, pkin(5); t37, -t35, 0, pkin(1); t35, t37, 0, 0; 0, 0, 1, qJ(2); t42, -t39, 0, pkin(2); t39, t42, 0, 0; 0, 0, 1, pkin(6); t36, -t34, 0, pkin(3); 0, 0, -1, -qJ(4); t34, t36, 0, 0; t41, -t38, 0, pkin(4); 0, 0, -1, -pkin(7); t38, t41, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
