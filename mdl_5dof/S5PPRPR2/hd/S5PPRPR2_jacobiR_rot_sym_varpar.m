% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRPR2
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
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:18
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPRPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:27
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:27
	% EndTime: 2019-10-24 10:18:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->6), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t11 = cos(pkin(7));
	t10 = sin(pkin(7));
	t9 = pkin(8) + qJ(3);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, -t11 * t7, 0, 0; 0, 0, -t10 * t7, 0, 0; 0, 0, t8, 0, 0; 0, 0, -t11 * t8, 0, 0; 0, 0, -t10 * t8, 0, 0; 0, 0, -t7, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->6)
	t34 = cos(pkin(7));
	t33 = sin(pkin(7));
	t32 = pkin(8) + qJ(3);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, t34 * t30, 0, 0; 0, 0, t33 * t30, 0, 0; 0, 0, -t31, 0, 0; 0, 0, t34 * t31, 0, 0; 0, 0, t33 * t31, 0, 0; 0, 0, t30, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:28
	% EndTime: 2019-10-24 10:18:28
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->8), mult. (26->18), div. (0->0), fcn. (45->6), ass. (0->12)
	t48 = sin(pkin(7));
	t50 = sin(qJ(5));
	t55 = t48 * t50;
	t51 = cos(qJ(5));
	t54 = t48 * t51;
	t49 = cos(pkin(7));
	t53 = t49 * t50;
	t52 = t49 * t51;
	t47 = pkin(8) + qJ(3);
	t46 = cos(t47);
	t45 = sin(t47);
	t1 = [0, 0, t46 * t53, 0, t45 * t52 - t55; 0, 0, t46 * t55, 0, t45 * t54 + t53; 0, 0, t45 * t50, 0, -t46 * t51; 0, 0, t46 * t52, 0, -t45 * t53 - t54; 0, 0, t46 * t54, 0, -t45 * t55 + t52; 0, 0, t45 * t51, 0, t46 * t50; 0, 0, -t49 * t45, 0, 0; 0, 0, -t48 * t45, 0, 0; 0, 0, t46, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end