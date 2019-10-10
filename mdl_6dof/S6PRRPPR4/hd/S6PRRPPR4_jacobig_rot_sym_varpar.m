% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPPR4_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR4_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
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
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t100 = cos(pkin(6));
	t102 = cos(qJ(2));
	t103 = t100 * t102;
	t101 = sin(qJ(2));
	t99 = cos(pkin(10));
	t98 = sin(pkin(6));
	t97 = sin(pkin(10));
	t1 = [0, t97 * t98, t99 * t101 + t97 * t103, 0, 0, 0; 0, -t99 * t98, t97 * t101 - t99 * t103, 0, 0, 0; 0, t100, -t98 * t102, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t109 = sin(pkin(10));
	t110 = sin(pkin(6));
	t120 = t109 * t110;
	t111 = cos(pkin(10));
	t119 = t111 * t110;
	t112 = cos(pkin(6));
	t114 = sin(qJ(2));
	t118 = t112 * t114;
	t116 = cos(qJ(2));
	t117 = t112 * t116;
	t115 = cos(qJ(3));
	t113 = sin(qJ(3));
	t1 = [0, t120, t109 * t117 + t111 * t114, 0, 0, -(-t109 * t118 + t111 * t116) * t113 + t115 * t120; 0, -t119, t109 * t114 - t111 * t117, 0, 0, -(t109 * t116 + t111 * t118) * t113 - t115 * t119; 0, t112, -t110 * t116, 0, 0, -t110 * t114 * t113 + t112 * t115;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end