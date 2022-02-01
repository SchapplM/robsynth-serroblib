% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRPR3
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
%   Siehe auch: S5RPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RPRPR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_jacobiR_rot_sym_varpar: pkin has to be [9x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
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
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (7->4), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->4)
	t12 = qJ(1) + pkin(8);
	t11 = cos(t12);
	t10 = sin(t12);
	t1 = [-t10, 0, 0, 0, 0; t11, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t11, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->7), mult. (0->0), div. (0->0), fcn. (8->2), ass. (0->4)
	t17 = qJ(1) + pkin(8) + qJ(3);
	t16 = cos(t17);
	t15 = sin(t17);
	t1 = [-t15, 0, -t15, 0, 0; t16, 0, t16, 0, 0; 0, 0, 0, 0, 0; -t16, 0, -t16, 0, 0; -t15, 0, -t15, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->5), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->10)
	t24 = qJ(1) + pkin(8) + qJ(3);
	t22 = sin(t24);
	t26 = cos(pkin(9));
	t28 = t22 * t26;
	t23 = cos(t24);
	t25 = sin(pkin(9));
	t27 = t23 * t25;
	t21 = t23 * t26;
	t20 = t22 * t25;
	t1 = [-t28, 0, -t28, 0, 0; t21, 0, t21, 0, 0; 0, 0, 0, 0, 0; t20, 0, t20, 0, 0; -t27, 0, -t27, 0, 0; 0, 0, 0, 0, 0; t23, 0, t23, 0, 0; t22, 0, t22, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:21:20
	% EndTime: 2022-01-23 09:21:20
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (72->11), mult. (42->14), div. (0->0), fcn. (72->6), ass. (0->16)
	t70 = qJ(1) + pkin(8) + qJ(3);
	t68 = sin(t70);
	t71 = sin(pkin(9));
	t77 = t68 * t71;
	t72 = cos(pkin(9));
	t73 = sin(qJ(5));
	t76 = t72 * t73;
	t74 = cos(qJ(5));
	t75 = t72 * t74;
	t69 = cos(t70);
	t67 = t69 * t71;
	t66 = t68 * t73 + t69 * t75;
	t65 = t68 * t74 - t69 * t76;
	t64 = -t68 * t75 + t69 * t73;
	t63 = t68 * t76 + t69 * t74;
	t1 = [t64, 0, t64, 0, t65; t66, 0, t66, 0, -t63; 0, 0, 0, 0, -t71 * t73; t63, 0, t63, 0, -t66; t65, 0, t65, 0, t64; 0, 0, 0, 0, -t71 * t74; -t77, 0, -t77, 0, 0; t67, 0, t67, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
end