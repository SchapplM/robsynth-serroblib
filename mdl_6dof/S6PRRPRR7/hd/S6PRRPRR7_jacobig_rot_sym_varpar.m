% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR7_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR7_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t18, 0, 0, 0, 0; 0, -cos(pkin(11)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t57 = cos(pkin(6));
	t59 = cos(qJ(2));
	t60 = t57 * t59;
	t58 = sin(qJ(2));
	t56 = cos(pkin(11));
	t55 = sin(pkin(6));
	t54 = sin(pkin(11));
	t1 = [0, t54 * t55, t54 * t60 + t56 * t58, 0, 0, 0; 0, -t56 * t55, t54 * t58 - t56 * t60, 0, 0, 0; 0, t57, -t55 * t59, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:18
	% EndTime: 2019-10-09 22:37:18
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t72 = cos(pkin(6));
	t74 = cos(qJ(2));
	t75 = t72 * t74;
	t73 = sin(qJ(2));
	t71 = cos(pkin(11));
	t70 = sin(pkin(6));
	t69 = sin(pkin(11));
	t1 = [0, t69 * t70, t69 * t75 + t71 * t73, 0, 0, 0; 0, -t71 * t70, t69 * t73 - t71 * t75, 0, 0, 0; 0, t72, -t70 * t74, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->9), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->13)
	t86 = sin(pkin(11));
	t87 = sin(pkin(6));
	t97 = t86 * t87;
	t88 = cos(pkin(11));
	t96 = t88 * t87;
	t89 = cos(pkin(6));
	t91 = sin(qJ(2));
	t95 = t89 * t91;
	t93 = cos(qJ(2));
	t94 = t89 * t93;
	t92 = cos(qJ(3));
	t90 = sin(qJ(3));
	t1 = [0, t97, t86 * t94 + t88 * t91, 0, (-t86 * t95 + t88 * t93) * t92 + t90 * t97, 0; 0, -t96, t86 * t91 - t88 * t94, 0, (t86 * t93 + t88 * t95) * t92 - t90 * t96, 0; 0, t89, -t87 * t93, 0, t87 * t91 * t92 + t89 * t90, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:37:19
	% EndTime: 2019-10-09 22:37:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->9), mult. (39->20), div. (0->0), fcn. (63->8), ass. (0->16)
	t104 = sin(pkin(11));
	t105 = sin(pkin(6));
	t115 = t104 * t105;
	t106 = cos(pkin(11));
	t114 = t106 * t105;
	t107 = cos(pkin(6));
	t109 = sin(qJ(2));
	t113 = t107 * t109;
	t111 = cos(qJ(2));
	t112 = t107 * t111;
	t110 = cos(qJ(3));
	t108 = sin(qJ(3));
	t103 = t105 * t109 * t110 + t107 * t108;
	t102 = (-t104 * t113 + t106 * t111) * t110 + t108 * t115;
	t101 = (t104 * t111 + t106 * t113) * t110 - t108 * t114;
	t1 = [0, t115, t104 * t112 + t106 * t109, 0, t102, t102; 0, -t114, t104 * t109 - t106 * t112, 0, t101, t101; 0, t107, -t105 * t111, 0, t103, t103;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end