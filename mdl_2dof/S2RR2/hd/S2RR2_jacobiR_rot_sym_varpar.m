% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% 
% Output:
% JR_rot [9x2]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S2RR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_jacobiR_rot_sym_varpar: qJ has to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S2RR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_jacobiR_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t3, 0; 0, 0; -t4, 0; -t4, 0; 0, 0; t3, 0; 0, 0; 0, 0; 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:02:34
	% EndTime: 2019-10-09 20:02:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t23 = sin(qJ(1));
	t24 = cos(qJ(2));
	t28 = t23 * t24;
	t22 = sin(qJ(2));
	t25 = cos(qJ(1));
	t27 = t25 * t22;
	t26 = t25 * t24;
	t21 = t23 * t22;
	t1 = [-t28, -t27; 0, t24; -t26, t21; t21, -t26; 0, -t22; t27, t28; t25, 0; 0, 0; -t23, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,2);
end