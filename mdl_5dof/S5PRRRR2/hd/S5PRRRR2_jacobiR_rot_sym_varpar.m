% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_jacobiR_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(2));
	t8 = sin(qJ(2));
	t1 = [0, -t8, 0, 0, 0; 0, t9, 0, 0, 0; 0, 0, 0, 0, 0; 0, -t9, 0, 0, 0; 0, -t8, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (14->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t16 = qJ(2) + qJ(3);
	t15 = cos(t16);
	t14 = sin(t16);
	t1 = [0, -t14, -t14, 0, 0; 0, t15, t15, 0, 0; 0, 0, 0, 0, 0; 0, -t15, -t15, 0, 0; 0, -t14, -t14, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (33->10), mult. (0->0), div. (0->0), fcn. (12->2), ass. (0->4)
	t21 = qJ(2) + qJ(3) + qJ(4);
	t20 = cos(t21);
	t19 = sin(t21);
	t1 = [0, -t19, -t19, -t19, 0; 0, t20, t20, t20, 0; 0, 0, 0, 0, 0; 0, -t20, -t20, -t20, 0; 0, -t19, -t19, -t19, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (55->12), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
	t29 = qJ(2) + qJ(3) + qJ(4);
	t27 = sin(t29);
	t31 = cos(qJ(5));
	t33 = t27 * t31;
	t28 = cos(t29);
	t30 = sin(qJ(5));
	t32 = t28 * t30;
	t26 = t28 * t31;
	t25 = t27 * t30;
	t1 = [0, -t33, -t33, -t33, -t32; 0, t26, t26, t26, -t25; 0, 0, 0, 0, t31; 0, t25, t25, t25, -t26; 0, -t32, -t32, -t32, -t33; 0, 0, 0, 0, -t30; 0, t28, t28, t28, 0; 0, t27, t27, t27, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end