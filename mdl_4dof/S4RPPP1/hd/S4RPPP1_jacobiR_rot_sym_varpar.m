% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% JR_rot [9x4]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S4RPPP1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiR_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPP1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiR_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0; t9, 0, 0, 0; 0, 0, 0, 0; -t9, 0, 0, 0; -t8, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t38 = sin(pkin(6));
	t42 = sin(qJ(1));
	t47 = t42 * t38;
	t40 = cos(pkin(6));
	t46 = t42 * t40;
	t43 = cos(qJ(1));
	t45 = t43 * t38;
	t44 = t43 * t40;
	t41 = cos(pkin(4));
	t39 = sin(pkin(4));
	t1 = [-t41 * t45 - t46, 0, 0, 0; -t41 * t47 + t44, 0, 0, 0; 0, 0, 0, 0; -t41 * t44 + t47, 0, 0, 0; -t41 * t46 - t45, 0, 0, 0; 0, 0, 0, 0; t43 * t39, 0, 0, 0; t42 * t39, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t45 = sin(pkin(6));
	t49 = sin(qJ(1));
	t54 = t49 * t45;
	t47 = cos(pkin(6));
	t53 = t49 * t47;
	t50 = cos(qJ(1));
	t52 = t50 * t45;
	t51 = t50 * t47;
	t48 = cos(pkin(4));
	t46 = sin(pkin(4));
	t1 = [t50 * t46, 0, 0, 0; t49 * t46, 0, 0, 0; 0, 0, 0, 0; t48 * t52 + t53, 0, 0, 0; t48 * t54 - t51, 0, 0, 0; 0, 0, 0, 0; t48 * t51 - t54, 0, 0, 0; t48 * t53 + t52, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (4->4), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t45 = sin(pkin(6));
	t49 = sin(qJ(1));
	t54 = t49 * t45;
	t47 = cos(pkin(6));
	t53 = t49 * t47;
	t50 = cos(qJ(1));
	t52 = t50 * t45;
	t51 = t50 * t47;
	t48 = cos(pkin(4));
	t46 = sin(pkin(4));
	t1 = [t50 * t46, 0, 0, 0; t49 * t46, 0, 0, 0; 0, 0, 0, 0; t48 * t51 - t54, 0, 0, 0; t48 * t53 + t52, 0, 0, 0; 0, 0, 0, 0; -t48 * t52 - t53, 0, 0, 0; -t48 * t54 + t51, 0, 0, 0; 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,4);
end