% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRPR6_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR6_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t42 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t42, 0, 0, 0, 0; 0, -cos(pkin(10)) * t42, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t57 = cos(pkin(6));
	t58 = sin(qJ(2));
	t60 = t57 * t58;
	t59 = cos(qJ(2));
	t56 = cos(pkin(10));
	t55 = sin(pkin(6));
	t54 = sin(pkin(10));
	t1 = [0, t54 * t55, 0, -t54 * t60 + t56 * t59, 0, 0; 0, -t56 * t55, 0, t54 * t59 + t56 * t60, 0, 0; 0, t57, 0, t55 * t58, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t84 = cos(pkin(6));
	t85 = sin(qJ(2));
	t87 = t84 * t85;
	t86 = cos(qJ(2));
	t83 = cos(pkin(10));
	t82 = sin(pkin(6));
	t81 = sin(pkin(10));
	t1 = [0, t81 * t82, 0, -t81 * t87 + t83 * t86, 0, 0; 0, -t83 * t82, 0, t81 * t86 + t83 * t87, 0, 0; 0, t84, 0, t82 * t85, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t89 = sin(pkin(10));
	t90 = sin(pkin(6));
	t100 = t89 * t90;
	t91 = cos(pkin(10));
	t99 = t91 * t90;
	t92 = cos(pkin(6));
	t94 = sin(qJ(2));
	t98 = t92 * t94;
	t96 = cos(qJ(2));
	t97 = t92 * t96;
	t95 = cos(qJ(4));
	t93 = sin(qJ(4));
	t1 = [0, t100, 0, -t89 * t98 + t91 * t96, 0, t93 * t100 - (t89 * t97 + t91 * t94) * t95; 0, -t99, 0, t89 * t96 + t91 * t98, 0, -t93 * t99 - (t89 * t94 - t91 * t97) * t95; 0, t92, 0, t90 * t94, 0, t90 * t95 * t96 + t92 * t93;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end