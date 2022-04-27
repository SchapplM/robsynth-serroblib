% Calculate homogenous joint transformation matrices for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-02-01 04:10
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5PRRRP8_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-02-01 04:10:23
% EndTime: 2022-02-01 04:10:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (9->9), mult. (6->6), div. (0->0), fcn. (26->10), ass. (0->11)
t101 = cos(qJ(2));
t100 = cos(qJ(3));
t99 = cos(qJ(4));
t98 = sin(qJ(2));
t97 = sin(qJ(3));
t96 = sin(qJ(4));
t95 = cos(pkin(5));
t94 = cos(pkin(9));
t93 = sin(pkin(5));
t92 = sin(pkin(9));
t1 = [t94, -t92, 0, 0; t92, t94, 0, 0; 0, 0, 1, qJ(1); t101, -t98, 0, pkin(1); t95 * t98, t95 * t101, -t93, -t93 * pkin(6); t93 * t98, t93 * t101, t95, t95 * pkin(6); t100, -t97, 0, pkin(2); 0, 0, -1, -pkin(7); t97, t100, 0, 0; t99, -t96, 0, pkin(3); 0, 0, -1, -pkin(8); t96, t99, 0, 0; 1, 0, 0, pkin(4); 0, 0, -1, -qJ(5); 0, 1, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
