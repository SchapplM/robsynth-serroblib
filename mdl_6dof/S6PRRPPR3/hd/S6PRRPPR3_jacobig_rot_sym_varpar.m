% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPPR3_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR3_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobig_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t57 = cos(pkin(6));
	t59 = cos(qJ(2));
	t60 = t57 * t59;
	t58 = sin(qJ(2));
	t56 = cos(pkin(10));
	t55 = sin(pkin(6));
	t54 = sin(pkin(10));
	t1 = [0, t54 * t55, t54 * t60 + t56 * t58, 0, 0, 0; 0, -t56 * t55, t54 * t58 - t56 * t60, 0, 0, 0; 0, t57, -t55 * t59, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t75 = cos(pkin(6));
	t77 = cos(qJ(2));
	t78 = t75 * t77;
	t76 = sin(qJ(2));
	t74 = cos(pkin(10));
	t73 = sin(pkin(6));
	t72 = sin(pkin(10));
	t1 = [0, t72 * t73, t72 * t78 + t74 * t76, 0, 0, 0; 0, -t74 * t73, t72 * t76 - t74 * t78, 0, 0, 0; 0, t75, -t73 * t77, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t67 = cos(pkin(6));
	t69 = cos(qJ(2));
	t70 = t67 * t69;
	t68 = sin(qJ(2));
	t66 = cos(pkin(10));
	t65 = sin(pkin(6));
	t64 = sin(pkin(10));
	t1 = [0, t64 * t65, t64 * t70 + t66 * t68, 0, 0, 0; 0, -t66 * t65, t64 * t68 - t66 * t70, 0, 0, 0; 0, t67, -t65 * t69, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:10:54
	% EndTime: 2019-10-09 22:10:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t85 = sin(pkin(10));
	t86 = sin(pkin(6));
	t96 = t85 * t86;
	t87 = cos(pkin(10));
	t95 = t87 * t86;
	t88 = cos(pkin(6));
	t90 = sin(qJ(2));
	t94 = t88 * t90;
	t92 = cos(qJ(2));
	t93 = t88 * t92;
	t91 = cos(qJ(3));
	t89 = sin(qJ(3));
	t1 = [0, t96, t85 * t93 + t87 * t90, 0, 0, (-t85 * t94 + t87 * t92) * t91 + t89 * t96; 0, -t95, t85 * t90 - t87 * t93, 0, 0, (t85 * t92 + t87 * t94) * t91 - t89 * t95; 0, t88, -t86 * t92, 0, 0, t86 * t90 * t91 + t88 * t89;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end