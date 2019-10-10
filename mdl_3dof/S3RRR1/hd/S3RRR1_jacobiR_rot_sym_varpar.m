% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% JR_rot [9x3]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S3RRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_jacobiR_rot_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3RRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_jacobiR_rot_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:13:34
	% EndTime: 2019-10-09 20:13:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:13:34
	% EndTime: 2019-10-09 20:13:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0; t9, 0, 0; 0, 0, 0; -t9, 0, 0; -t8, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:13:34
	% EndTime: 2019-10-09 20:13:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(1) + qJ(2);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [-t14, -t14, 0; t15, t15, 0; 0, 0, 0; -t15, -t15, 0; -t14, -t14, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:13:34
	% EndTime: 2019-10-09 20:13:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t21 = qJ(1) + qJ(2) + qJ(3);
	t20 = cos(t21);
	t19 = sin(t21);
	t1 = [-t19, -t19, -t19; t20, t20, t20; 0, 0, 0; -t20, -t20, -t20; -t19, -t19, -t19; 0, 0, 0; 0, 0, 0; 0, 0, 0; 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,3);
end