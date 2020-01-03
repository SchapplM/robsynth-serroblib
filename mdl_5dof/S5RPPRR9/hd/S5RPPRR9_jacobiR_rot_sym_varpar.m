% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPPRR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_jacobiR_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:03:08
	% EndTime: 2019-12-31 18:03:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:03:08
	% EndTime: 2019-12-31 18:03:08
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
	% StartTime: 2019-12-31 18:03:08
	% EndTime: 2019-12-31 18:03:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [-t5, 0, 0, 0, 0; t6, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; t6, 0, 0, 0, 0; t5, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:03:08
	% EndTime: 2019-12-31 18:03:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->3), mult. (8->4), div. (0->0), fcn. (16->4), ass. (0->7)
	t26 = cos(qJ(1));
	t25 = sin(qJ(1));
	t24 = cos(pkin(8));
	t23 = sin(pkin(8));
	t22 = t25 * t23 + t26 * t24;
	t21 = t26 * t23 - t25 * t24;
	t1 = [t21, 0, 0, 0, 0; t22, 0, 0, 0, 0; 0, 0, 0, 0, 0; t22, 0, 0, 0, 0; -t21, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:03:08
	% EndTime: 2019-12-31 18:03:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->5), mult. (28->8), div. (0->0), fcn. (50->6), ass. (0->13)
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t28 = sin(pkin(8));
	t29 = cos(pkin(8));
	t21 = -t34 * t28 - t35 * t29;
	t26 = sin(qJ(4));
	t33 = t21 * t26;
	t27 = cos(qJ(4));
	t32 = t21 * t27;
	t22 = t35 * t28 - t34 * t29;
	t31 = t22 * t26;
	t30 = t22 * t27;
	t1 = [t30, 0, 0, t33, 0; -t32, 0, 0, t31, 0; 0, 0, 0, -t27, 0; -t31, 0, 0, t32, 0; t33, 0, 0, t30, 0; 0, 0, 0, t26, 0; t21, 0, 0, 0, 0; t22, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:03:08
	% EndTime: 2019-12-31 18:03:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (36->15), mult. (88->24), div. (0->0), fcn. (141->8), ass. (0->19)
	t96 = cos(qJ(1));
	t95 = sin(qJ(1));
	t83 = sin(qJ(5));
	t84 = sin(qJ(4));
	t94 = t84 * t83;
	t85 = cos(qJ(5));
	t93 = t84 * t85;
	t86 = cos(qJ(4));
	t92 = t86 * t83;
	t91 = t86 * t85;
	t90 = cos(pkin(8));
	t89 = sin(pkin(8));
	t78 = -t95 * t89 - t96 * t90;
	t79 = t96 * t89 - t95 * t90;
	t88 = t78 * t83 + t79 * t91;
	t87 = -t78 * t85 + t79 * t92;
	t77 = -t78 * t91 + t79 * t83;
	t76 = t78 * t92 + t79 * t85;
	t1 = [t88, 0, 0, t78 * t93, t76; t77, 0, 0, t79 * t93, t87; 0, 0, 0, -t91, t94; -t87, 0, 0, -t78 * t94, -t77; t76, 0, 0, -t79 * t94, t88; 0, 0, 0, t92, t93; t79 * t84, 0, 0, -t78 * t86, 0; -t78 * t84, 0, 0, -t79 * t86, 0; 0, 0, 0, -t84, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end