% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRP3_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP3_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
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
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t83 = cos(pkin(6));
	t85 = cos(qJ(2));
	t86 = t83 * t85;
	t84 = sin(qJ(2));
	t82 = cos(pkin(10));
	t81 = sin(pkin(6));
	t80 = sin(pkin(10));
	t1 = [0, t80 * t81, t80 * t86 + t82 * t84, 0, 0, 0; 0, -t82 * t81, t80 * t84 - t82 * t86, 0, 0, 0; 0, t83, -t81 * t85, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t90 = sin(pkin(10));
	t91 = sin(pkin(6));
	t101 = t90 * t91;
	t92 = cos(pkin(10));
	t100 = t92 * t91;
	t93 = cos(pkin(6));
	t95 = sin(qJ(2));
	t99 = t93 * t95;
	t97 = cos(qJ(2));
	t98 = t93 * t97;
	t96 = cos(qJ(3));
	t94 = sin(qJ(3));
	t1 = [0, t101, t90 * t98 + t92 * t95, 0, (-t90 * t99 + t92 * t97) * t94 - t96 * t101, 0; 0, -t100, t90 * t95 - t92 * t98, 0, (t90 * t97 + t92 * t99) * t94 + t96 * t100, 0; 0, t93, -t91 * t97, 0, t91 * t95 * t94 - t93 * t96, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:20:05
	% EndTime: 2019-10-09 22:20:05
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t119 = sin(pkin(10));
	t120 = sin(pkin(6));
	t130 = t119 * t120;
	t121 = cos(pkin(10));
	t129 = t121 * t120;
	t122 = cos(pkin(6));
	t124 = sin(qJ(2));
	t128 = t122 * t124;
	t126 = cos(qJ(2));
	t127 = t122 * t126;
	t125 = cos(qJ(3));
	t123 = sin(qJ(3));
	t1 = [0, t130, t119 * t127 + t121 * t124, 0, (-t119 * t128 + t121 * t126) * t123 - t125 * t130, 0; 0, -t129, t119 * t124 - t121 * t127, 0, (t119 * t126 + t121 * t128) * t123 + t125 * t129, 0; 0, t122, -t120 * t126, 0, t120 * t124 * t123 - t122 * t125, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end