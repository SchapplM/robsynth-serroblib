% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPRP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(7);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->3), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t12 = cos(pkin(8));
	t11 = sin(pkin(8));
	t10 = qJ(1) + pkin(7);
	t9 = cos(t10);
	t8 = sin(t10);
	t1 = [-t8 * t12, 0, 0, 0, 0; t9 * t12, 0, 0, 0, 0; 0, 0, 0, 0, 0; t8 * t11, 0, 0, 0, 0; -t9 * t11, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t8, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (28->14), div. (0->0), fcn. (48->6), ass. (0->14)
	t55 = cos(pkin(8));
	t56 = sin(qJ(4));
	t59 = t55 * t56;
	t57 = cos(qJ(4));
	t58 = t55 * t57;
	t54 = sin(pkin(8));
	t53 = qJ(1) + pkin(7);
	t52 = cos(t53);
	t51 = sin(t53);
	t50 = t51 * t56 + t52 * t58;
	t49 = t51 * t57 - t52 * t59;
	t48 = -t51 * t58 + t52 * t56;
	t47 = t51 * t59 + t52 * t57;
	t1 = [t48, 0, 0, t49, 0; t50, 0, 0, -t47, 0; 0, 0, 0, -t54 * t56, 0; t47, 0, 0, -t50, 0; t49, 0, 0, t48, 0; 0, 0, 0, -t54 * t57, 0; -t51 * t54, 0, 0, 0, 0; t52 * t54, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (28->14), div. (0->0), fcn. (48->6), ass. (0->14)
	t57 = cos(pkin(8));
	t58 = sin(qJ(4));
	t61 = t57 * t58;
	t59 = cos(qJ(4));
	t60 = t57 * t59;
	t56 = sin(pkin(8));
	t55 = qJ(1) + pkin(7);
	t54 = cos(t55);
	t53 = sin(t55);
	t52 = t53 * t58 + t54 * t60;
	t51 = t53 * t59 - t54 * t61;
	t50 = -t53 * t60 + t54 * t58;
	t49 = t53 * t61 + t54 * t59;
	t1 = [t50, 0, 0, t51, 0; t52, 0, 0, -t49, 0; 0, 0, 0, -t56 * t58, 0; t49, 0, 0, -t52, 0; t51, 0, 0, t50, 0; 0, 0, 0, -t56 * t59, 0; -t53 * t56, 0, 0, 0, 0; t54 * t56, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end