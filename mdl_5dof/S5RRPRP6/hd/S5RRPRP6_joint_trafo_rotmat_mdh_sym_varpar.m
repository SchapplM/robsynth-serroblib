% Calculate homogenous joint transformation matrices for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-03 17:50
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5RRPRP6_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-03 17:50:30
% EndTime: 2022-02-03 17:50:30
% DurationCPUTime: 0.03s
% Computational Cost: add. (6->6), mult. (0->0), div. (0->0), fcn. (16->8), ass. (0->9)
t48 = cos(qJ(1));
t47 = cos(qJ(2));
t46 = cos(qJ(4));
t45 = sin(qJ(1));
t44 = sin(qJ(2));
t43 = sin(qJ(4));
t42 = cos(pkin(8));
t41 = sin(pkin(8));
t1 = [t48, -t45, 0, 0; t45, t48, 0, 0; 0, 0, 1, pkin(5); t47, -t44, 0, pkin(1); 0, 0, -1, -pkin(6); t44, t47, 0, 0; t42, -t41, 0, pkin(2); t41, t42, 0, 0; 0, 0, 1, qJ(3); t46, -t43, 0, pkin(3); 0, 0, -1, -pkin(7); t43, t46, 0, 0; 1, 0, 0, pkin(4); 0, 1, 0, 0; 0, 0, 1, qJ(5);];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
