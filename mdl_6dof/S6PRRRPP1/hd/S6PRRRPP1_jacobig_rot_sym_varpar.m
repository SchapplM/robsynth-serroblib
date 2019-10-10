% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPP1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(10)) * t18, 0, 0, 0, 0; 0, -cos(pkin(10)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
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
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t83 = sin(pkin(10));
	t84 = sin(pkin(6));
	t94 = t83 * t84;
	t85 = cos(pkin(10));
	t93 = t85 * t84;
	t86 = cos(pkin(6));
	t88 = sin(qJ(2));
	t92 = t86 * t88;
	t90 = cos(qJ(2));
	t91 = t86 * t90;
	t89 = cos(qJ(3));
	t87 = sin(qJ(3));
	t1 = [0, t94, t83 * t91 + t85 * t88, (-t83 * t92 + t85 * t90) * t87 - t89 * t94, 0, 0; 0, -t93, t83 * t88 - t85 * t91, (t83 * t90 + t85 * t92) * t87 + t89 * t93, 0, 0; 0, t86, -t84 * t90, t84 * t88 * t87 - t86 * t89, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t91 = sin(pkin(10));
	t92 = sin(pkin(6));
	t102 = t91 * t92;
	t93 = cos(pkin(10));
	t101 = t93 * t92;
	t94 = cos(pkin(6));
	t96 = sin(qJ(2));
	t100 = t94 * t96;
	t98 = cos(qJ(2));
	t99 = t94 * t98;
	t97 = cos(qJ(3));
	t95 = sin(qJ(3));
	t1 = [0, t102, t91 * t99 + t93 * t96, (-t91 * t100 + t93 * t98) * t95 - t97 * t102, 0, 0; 0, -t101, t91 * t96 - t93 * t99, (t93 * t100 + t91 * t98) * t95 + t97 * t101, 0, 0; 0, t94, -t92 * t98, t92 * t96 * t95 - t94 * t97, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t120 = sin(pkin(10));
	t121 = sin(pkin(6));
	t131 = t120 * t121;
	t122 = cos(pkin(10));
	t130 = t122 * t121;
	t123 = cos(pkin(6));
	t125 = sin(qJ(2));
	t129 = t123 * t125;
	t127 = cos(qJ(2));
	t128 = t123 * t127;
	t126 = cos(qJ(3));
	t124 = sin(qJ(3));
	t1 = [0, t131, t120 * t128 + t122 * t125, (-t120 * t129 + t122 * t127) * t124 - t126 * t131, 0, 0; 0, -t130, t120 * t125 - t122 * t128, (t120 * t127 + t122 * t129) * t124 + t126 * t130, 0, 0; 0, t123, -t121 * t127, t121 * t125 * t124 - t123 * t126, 0, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end