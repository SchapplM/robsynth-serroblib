% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:31
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5PRRPR5_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR5_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_jacobig_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:23
	% EndTime: 2019-10-24 10:31:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:23
	% EndTime: 2019-10-24 10:31:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:23
	% EndTime: 2019-10-24 10:31:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(5));
	t1 = [0, sin(pkin(9)) * t18, 0, 0, 0; 0, -cos(pkin(9)) * t18, 0, 0, 0; 0, cos(pkin(5)), 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:24
	% EndTime: 2019-10-24 10:31:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t57 = cos(pkin(5));
	t59 = cos(qJ(2));
	t60 = t57 * t59;
	t58 = sin(qJ(2));
	t56 = cos(pkin(9));
	t55 = sin(pkin(5));
	t54 = sin(pkin(9));
	t1 = [0, t54 * t55, t54 * t60 + t56 * t58, 0, 0; 0, -t56 * t55, t54 * t58 - t56 * t60, 0, 0; 0, t57, -t55 * t59, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:24
	% EndTime: 2019-10-24 10:31:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t64 = cos(pkin(5));
	t66 = cos(qJ(2));
	t67 = t64 * t66;
	t65 = sin(qJ(2));
	t63 = cos(pkin(9));
	t62 = sin(pkin(5));
	t61 = sin(pkin(9));
	t1 = [0, t61 * t62, t61 * t67 + t63 * t65, 0, 0; 0, -t63 * t62, t61 * t65 - t63 * t67, 0, 0; 0, t64, -t62 * t66, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:24
	% EndTime: 2019-10-24 10:31:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (15->10), mult. (24->19), div. (0->0), fcn. (40->8), ass. (0->15)
	t95 = sin(pkin(9));
	t96 = sin(pkin(5));
	t105 = t95 * t96;
	t97 = cos(pkin(9));
	t104 = t97 * t96;
	t98 = cos(pkin(5));
	t99 = sin(qJ(2));
	t103 = t98 * t99;
	t100 = cos(qJ(2));
	t102 = t95 * t100;
	t101 = t97 * t100;
	t94 = qJ(3) + pkin(10);
	t93 = cos(t94);
	t92 = sin(t94);
	t1 = [0, t105, t98 * t102 + t97 * t99, 0, (-t95 * t103 + t101) * t92 - t93 * t105; 0, -t104, -t98 * t101 + t95 * t99, 0, (t97 * t103 + t102) * t92 + t93 * t104; 0, t98, -t96 * t100, 0, t96 * t99 * t92 - t98 * t93;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,5);
end