% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PPRPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (6->6), mult. (14->10), div. (0->0), fcn. (24->6), ass. (0->11)
	t22 = sin(pkin(7));
	t25 = sin(qJ(3));
	t30 = t22 * t25;
	t26 = cos(qJ(3));
	t29 = t22 * t26;
	t24 = cos(pkin(7));
	t28 = t24 * t25;
	t27 = t24 * t26;
	t23 = cos(pkin(8));
	t21 = sin(pkin(8));
	t1 = [0, 0, -t23 * t28 + t29, 0, 0; 0, 0, -t23 * t30 - t27, 0, 0; 0, 0, -t21 * t25, 0, 0; 0, 0, -t23 * t27 - t30, 0, 0; 0, 0, -t23 * t29 + t28, 0, 0; 0, 0, -t21 * t26, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (16->7), mult. (14->11), div. (0->0), fcn. (24->6), ass. (0->11)
	t29 = sin(pkin(7));
	t30 = cos(pkin(8));
	t34 = t29 * t30;
	t27 = qJ(3) + pkin(9);
	t25 = sin(t27);
	t31 = cos(pkin(7));
	t33 = t31 * t25;
	t26 = cos(t27);
	t32 = t31 * t26;
	t28 = sin(pkin(8));
	t1 = [0, 0, t29 * t26 - t30 * t33, 0, 0; 0, 0, -t25 * t34 - t32, 0, 0; 0, 0, -t28 * t25, 0, 0; 0, 0, -t29 * t25 - t30 * t32, 0, 0; 0, 0, -t26 * t34 + t33, 0, 0; 0, 0, -t28 * t26, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:05:41
	% EndTime: 2019-12-05 15:05:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (44->14), mult. (57->30), div. (0->0), fcn. (88->8), ass. (0->19)
	t73 = sin(pkin(8));
	t77 = sin(qJ(5));
	t83 = t73 * t77;
	t78 = cos(qJ(5));
	t82 = t73 * t78;
	t74 = sin(pkin(7));
	t75 = cos(pkin(8));
	t81 = t74 * t75;
	t72 = qJ(3) + pkin(9);
	t70 = sin(t72);
	t76 = cos(pkin(7));
	t80 = t76 * t70;
	t71 = cos(t72);
	t79 = t76 * t71;
	t69 = t74 * t70 + t75 * t79;
	t68 = t74 * t71 - t75 * t80;
	t67 = t71 * t81 - t80;
	t66 = -t70 * t81 - t79;
	t1 = [0, 0, t68 * t78, 0, -t69 * t77 + t76 * t82; 0, 0, t66 * t78, 0, -t67 * t77 + t74 * t82; 0, 0, -t70 * t82, 0, -t71 * t83 - t75 * t78; 0, 0, -t68 * t77, 0, -t69 * t78 - t76 * t83; 0, 0, -t66 * t77, 0, -t67 * t78 - t74 * t83; 0, 0, t70 * t83, 0, -t71 * t82 + t75 * t77; 0, 0, t69, 0, 0; 0, 0, t67, 0, 0; 0, 0, t73 * t71, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end