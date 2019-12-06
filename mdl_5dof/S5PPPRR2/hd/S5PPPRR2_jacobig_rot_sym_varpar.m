% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5PPPRR2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_jacobig_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (7->6), div. (0->0), fcn. (12->6), ass. (0->6)
	t24 = sin(pkin(9));
	t29 = t24 * cos(pkin(8));
	t28 = cos(pkin(7));
	t26 = cos(pkin(9));
	t25 = sin(pkin(7));
	t1 = [0, 0, 0, -t25 * t26 + t28 * t29, 0; 0, 0, 0, t25 * t29 + t28 * t26, 0; 0, 0, 0, sin(pkin(8)) * t24, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:59:52
	% EndTime: 2019-12-05 14:59:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (22->18), div. (0->0), fcn. (35->8), ass. (0->13)
	t70 = sin(pkin(8));
	t76 = cos(qJ(4));
	t80 = t70 * t76;
	t71 = sin(pkin(7));
	t73 = cos(pkin(8));
	t79 = t71 * t73;
	t69 = sin(pkin(9));
	t74 = cos(pkin(7));
	t78 = t74 * t69;
	t72 = cos(pkin(9));
	t77 = t74 * t72;
	t75 = sin(qJ(4));
	t1 = [0, 0, 0, -t71 * t72 + t73 * t78, (t71 * t69 + t73 * t77) * t75 - t74 * t80; 0, 0, 0, t69 * t79 + t77, (t72 * t79 - t78) * t75 - t71 * t80; 0, 0, 0, t70 * t69, t70 * t72 * t75 + t73 * t76;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,5);
end