% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S4PRRR7
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% Jg_rot [3x4]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:35
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S4PRRR7_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR7_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_jacobig_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(4));
	t1 = [0, sin(pkin(8)) * t18, 0, 0; 0, -cos(pkin(8)) * t18, 0, 0; 0, cos(pkin(4)), 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t56 = cos(pkin(4));
	t58 = cos(qJ(2));
	t59 = t56 * t58;
	t57 = sin(qJ(2));
	t55 = cos(pkin(8));
	t54 = sin(pkin(4));
	t53 = sin(pkin(8));
	t1 = [0, t53 * t54, t53 * t59 + t55 * t57, 0; 0, -t55 * t54, t53 * t57 - t55 * t59, 0; 0, t56, -t54 * t58, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:59
	% EndTime: 2019-12-29 12:34:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t83 = sin(pkin(8));
	t84 = sin(pkin(4));
	t94 = t83 * t84;
	t85 = cos(pkin(8));
	t93 = t85 * t84;
	t86 = cos(pkin(4));
	t88 = sin(qJ(2));
	t92 = t86 * t88;
	t90 = cos(qJ(2));
	t91 = t86 * t90;
	t89 = cos(qJ(3));
	t87 = sin(qJ(3));
	t1 = [0, t94, t83 * t91 + t85 * t88, (-t83 * t92 + t85 * t90) * t87 - t89 * t94; 0, -t93, t83 * t88 - t85 * t91, (t83 * t90 + t85 * t92) * t87 + t89 * t93; 0, t86, -t84 * t90, t84 * t88 * t87 - t86 * t89;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,4);
end