% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRP5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPRP5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_jacobiR_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t10 = cos(qJ(1));
	t9 = sin(qJ(1));
	t8 = cos(pkin(7));
	t7 = sin(pkin(7));
	t1 = [-t9 * t8, 0, 0, 0, 0; t10 * t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9 * t7, 0, 0, 0, 0; -t10 * t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (10->4), ass. (0->5)
	t34 = cos(qJ(1));
	t33 = sin(qJ(1));
	t32 = cos(pkin(7));
	t31 = sin(pkin(7));
	t1 = [-t33 * t32, 0, 0, 0, 0; t34 * t32, 0, 0, 0, 0; 0, 0, 0, 0, 0; t34, 0, 0, 0, 0; t33, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t33 * t31, 0, 0, 0, 0; t34 * t31, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:09
	% EndTime: 2019-12-31 17:54:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (12->10), mult. (36->8), div. (0->0), fcn. (58->6), ass. (0->13)
	t21 = sin(pkin(7));
	t22 = cos(pkin(7));
	t23 = sin(qJ(4));
	t25 = cos(qJ(4));
	t28 = t21 * t25 - t22 * t23;
	t27 = t21 * t23 + t22 * t25;
	t26 = cos(qJ(1));
	t24 = sin(qJ(1));
	t20 = t27 * t26;
	t19 = t28 * t26;
	t18 = t27 * t24;
	t17 = t28 * t24;
	t1 = [-t18, 0, 0, t19, 0; t20, 0, 0, t17, 0; 0, 0, 0, -t27, 0; -t17, 0, 0, -t20, 0; t19, 0, 0, -t18, 0; 0, 0, 0, -t28, 0; -t26, 0, 0, 0, 0; -t24, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:54:10
	% EndTime: 2019-12-31 17:54:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (12->7), mult. (36->8), div. (0->0), fcn. (58->6), ass. (0->13)
	t62 = sin(pkin(7));
	t63 = cos(pkin(7));
	t64 = sin(qJ(4));
	t66 = cos(qJ(4));
	t69 = t62 * t66 - t63 * t64;
	t68 = t62 * t64 + t63 * t66;
	t67 = cos(qJ(1));
	t65 = sin(qJ(1));
	t61 = t68 * t67;
	t60 = t69 * t67;
	t59 = t68 * t65;
	t58 = t69 * t65;
	t1 = [-t59, 0, 0, t60, 0; t61, 0, 0, t58, 0; 0, 0, 0, -t68, 0; -t67, 0, 0, 0, 0; -t65, 0, 0, 0, 0; 0, 0, 0, 0, 0; t58, 0, 0, t61, 0; -t60, 0, 0, t59, 0; 0, 0, 0, t69, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end